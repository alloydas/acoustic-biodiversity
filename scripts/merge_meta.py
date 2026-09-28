"""Left-join Xeno-canto metadata onto each continent's scores by 'id'.

Loads the flattened metadata CSV once, then for each continent writes
score_<continent>_meta.csv = score columns + metadata columns (minus the id key).
Score rows are the base score_<continent>.csv PLUS the gap-filled
score_<continent>_gap.csv PLUS the 10-year backfill score_<continent>_hist.csv
(all identical schema); duplicate ids are dropped, across continents too.

A fourth source, EXTRA (merged_metadata_all.csv: an externally computed
scores+metadata table, 1886-2023), adds every recording whose id is not already
in base/gap/hist. Overlapping ids keep the locally computed values. EXTRA rows
carry their own metadata (remapped to this file's column names) and are routed
to the continent the local files use for the same country and place (see
extra_continent).

A trailing 'index_source' column marks each row 'local' or 'extra'. The EXTRA
indices were computed with different settings (on the 461k shared recordings its
Bio_acoustic_Index median is ~324 vs ~9 here, rank rho 0.24), so they are not
comparable with local values; build_cells.py keeps them out of the index means.
"""
import os
import sys
import csv
import ast
import math
import collections

csv.field_size_limit(10 ** 7)

META = 'metadata_2015-2025.csv'
CONTINENTS = ['africa', 'america', 'asia', 'australia', 'europe']
EXTRA = 'merged_metadata_all.csv'
# EXTRA column name -> our column name (everything else matches by name)
EXTRA_RENAME = {'group': 'grp', 'lng': 'lon'}
# countries absent from the local data (continent can't be learned from it)
EXTRA_COUNTRIES = {'yemen': 'asia', 'syria': 'asia', 'afghanistan': 'asia', 'sudan': 'africa',
                   'kosovo': 'europe', 'vanuatu': 'australia', 'kiribati': 'australia'}


def load_meta(path):
    with open(path, newline='') as f:
        r = csv.reader(f)
        hdr = next(r)
        idx = hdr.index('id')
        cols = [c for i, c in enumerate(hdr) if i != idx]
        meta = {row[idx]: [v for i, v in enumerate(row) if i != idx] for row in r}
    return meta, cols


def fnum(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def geo_key(cnt, lat, lon):
    """(country, 5-degree cell) - the unit the continent vote is learned on."""
    return cnt, math.floor(lat / 5), math.floor(lon / 5)


def local_ids_and_country_map():
    """All locally scored ids (any continent); the owner continent of each id (the
    first in CONTINENTS order that has it); and the majority continent per
    (normalised) country and per (country, 5-degree cell), learned from them."""
    local, owner = set(), {}
    votes = collections.defaultdict(collections.Counter)
    gvotes = collections.defaultdict(collections.Counter)
    for c in CONTINENTS:
        for src in (f'score_{c}.csv', f'score_{c}_gap.csv', f'score_{c}_hist.csv'):
            if not os.path.exists(src):
                continue
            with open(src, newline='') as f:
                ids = {row[1] for row in csv.reader(f)}
            local |= ids
            for rid in ids:
                owner.setdefault(rid, c)
    # vote on what merge() writes: each id once, under its owner continent (the 1,569
    # ids in both asia and australia, nearly all Papua New Guinea, are written under asia)
    for rid, c in owner.items():
        m = META_ROWS.get(rid)
        if m is not None:
            k = m[META_CNT].strip().lower()
            votes[k][c] += 1
            lat, lon = fnum(m[META_LAT]), fnum(m[META_LON])
            if lat is not None and lon is not None:
                gvotes[geo_key(k, lat, lon)][c] += 1
    cmap = {k: v.most_common(1)[0][0] for k, v in votes.items() if k}
    cmap.update(EXTRA_COUNTRIES)
    gmap = {k: v.most_common(1)[0][0] for k, v in gvotes.items() if k[0]}
    return local, owner, cmap, gmap


def extra_continent(cnt, lat, lon, cmap, gmap):
    k = cnt.strip().lower()
    # where local recordings exist for the same country and 5-degree cell, follow
    # them: Xeno-canto files Siberia under asia, French Polynesia under australia,
    # the Caribbean Netherlands and UK territories under america
    if lat is not None and lon is not None:
        c = gmap.get(geo_key(k, lat, lon))
        if c is not None:
            return c
    c = cmap.get(k)
    if c == 'europe' and lat is not None and lon is not None:
        if k == 'russian federation':
            return 'asia' if lon > 60 else 'europe'
        # Dutch Caribbean (Aruba, Curacao, Bonaire, Sint Maarten) and Pitcairn are
        # filed under america / australia locally
        if k == 'netherlands' and lon < -30:
            return 'america'
        if k == 'united kingdom' and lat < 0 and lon < -100:
            return 'australia'
        # French overseas territories are filed under their geographic continent
        if k == 'france':
            if lat < 0 and lon < -100:
                return 'australia'
            if lon < -30:
                return 'america'
            if lon > 100:
                return 'australia'
            if lat < 0 and 30 < lon < 70:
                return 'africa'
    return c


def load_extra(scols, mcols, local, cmap, gmap):
    """EXTRA rows reshaped to scols + mcols, grouped by continent (first id wins).
    Ids already scored locally are dropped here, globally: EXTRA's country routing
    can disagree with the continent a local id was filed under."""
    out = collections.defaultdict(list)
    seen, skipped = set(), collections.Counter()
    with open(EXTRA, newline='') as f:
        r = csv.reader(f)
        hdr = [EXTRA_RENAME.get(h, h) for h in next(r)]
        ei = {h: i for i, h in enumerate(hdr)}
        missing = [c for c in scols + mcols if c not in ei]
        if missing:
            sys.exit(f'{EXTRA} lacks columns: {missing}')
        for row in r:
            rid = row[ei['id']]
            if rid in seen:
                skipped['dup id'] += 1
                continue
            seen.add(rid)
            if rid in local:
                skipped['already local'] += 1
                continue
            rec = {c: row[ei[c]] for c in scols + mcols}
            # 'also' is a python list repr there; ';'-joined here
            a = rec['also'].strip()
            if a.startswith('['):
                try:
                    rec['also'] = ';'.join(str(x).strip() for x in ast.literal_eval(a) if str(x).strip())
                except (ValueError, SyntaxError):
                    rec['also'] = ''
            if rec['smp'].endswith('.0'):
                rec['smp'] = rec['smp'][:-2]
            c = extra_continent(rec['cnt'], fnum(rec['lat']), fnum(rec['lon']), cmap, gmap)
            if c is None:
                skipped[f'no continent ({rec["cnt"].strip() or "blank"})'] += 1
                continue
            out[c].append([rec[k] for k in scols] + [rec[k] for k in mcols])
    print(f'{EXTRA}: {len(seen)} unique ids, skipped {dict(skipped)}')
    return out


def merge(continent, meta, mcols, extra, owner):
    base = f'score_{continent}.csv'
    gap = f'score_{continent}_gap.csv'
    hist = f'score_{continent}_hist.csv'
    out = f'score_{continent}_meta.csv'
    blank = [''] * len(mcols)
    matched = total = gap_rows = hist_rows = dupes = extra_rows = 0
    # 1,569 recordings sit in both the asia and australia score files (Xeno-canto's
    # areas overlap in Wallacea/New Guinea); keeping both double-counts them in every
    # cell they fall in. Each id is written only under its owner continent (the first
    # in CONTINENTS order that has it), so a partial rerun (argv) matches a full one.
    seen = set()
    with open(out, 'w', newline='') as fo:
        w = csv.writer(fo)
        shdr = sid = None
        # base, then gap, then the 10-year backfill (all same schema); skip each
        # file's header and any id already written (sources are designed to be
        # disjoint, but enforce it — dedup by id keeps the first occurrence).
        for src in (base, gap, hist):
            if not os.path.exists(src):
                continue
            with open(src, newline='') as f:
                r = csv.reader(f)
                hdr = next(r)
                if shdr is None:
                    shdr = hdr
                    sid = shdr.index('id')
                    w.writerow(shdr + mcols + ['index_source'])
                for row in r:
                    rid = row[sid]
                    if rid in seen or owner.get(rid, continent) != continent:
                        dupes += 1
                        continue
                    seen.add(rid)
                    total += 1
                    if src == gap:
                        gap_rows += 1
                    elif src == hist:
                        hist_rows += 1
                    m = meta.get(rid)
                    if m is not None:
                        matched += 1
                    w.writerow(row + (m if m is not None else blank) + ['local'])
        # EXTRA rows are already in shdr + mcols order, carry their own metadata and
        # were deduped against every local id in load_extra
        for row in extra.get(continent, ()):
            total += 1
            extra_rows += 1
            matched += 1
            w.writerow(row + ['extra'])
    pct = 100 * matched / total if total else 0
    print(f'{continent:10} {total:>7} rows (+{gap_rows} gap, +{hist_rows} hist, +{extra_rows} extra, {dupes} dup dropped)  '
          f'{matched:>7} matched ({pct:5.1f}%)  unmatched={total - matched}  -> {out}')


if __name__ == '__main__':
    meta, mcols = load_meta(META)
    print(f'loaded {len(meta)} metadata records, {len(mcols)} metadata columns')
    META_ROWS, META_CNT, META_LAT, META_LON = meta, mcols.index('cnt'), mcols.index('lat'), mcols.index('lon')
    with open(f'score_{CONTINENTS[0]}.csv', newline='') as f:
        scols = next(csv.reader(f))
    local, owner, cmap, gmap = local_ids_and_country_map()
    extra = load_extra(scols, mcols, local, cmap, gmap)
    print()
    for c in (sys.argv[1:] or CONTINENTS):
        merge(c, meta, mcols, extra, owner)
