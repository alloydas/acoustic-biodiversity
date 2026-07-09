"""Left-join Xeno-canto metadata onto each continent's scores by 'id'.

Loads the flattened metadata CSV once, then for each continent writes
score_<continent>_meta.csv = score columns + metadata columns (minus the id key).
Score rows are the base score_<continent>.csv PLUS the gap-filled
score_<continent>_gap.csv (identical schema); duplicate ids are dropped.
"""
import os
import sys
import csv

csv.field_size_limit(10 ** 7)

META = 'metadata_since_2023-01-01.csv'
CONTINENTS = ['africa', 'america', 'asia', 'australia', 'europe']


def load_meta(path):
    with open(path, newline='') as f:
        r = csv.reader(f)
        hdr = next(r)
        idx = hdr.index('id')
        cols = [c for i, c in enumerate(hdr) if i != idx]
        meta = {row[idx]: [v for i, v in enumerate(row) if i != idx] for row in r}
    return meta, cols


def merge(continent, meta, mcols):
    base = f'score_{continent}.csv'
    gap = f'score_{continent}_gap.csv'
    out = f'score_{continent}_meta.csv'
    blank = [''] * len(mcols)
    matched = total = gap_rows = dupes = 0
    seen = set()
    with open(out, 'w', newline='') as fo:
        w = csv.writer(fo)
        shdr = sid = None
        # base first, then gap (same schema); skip the gap file's header and
        # any id already written (gap is designed to be disjoint, but enforce it).
        for src in (base, gap):
            if not os.path.exists(src):
                continue
            with open(src, newline='') as f:
                r = csv.reader(f)
                hdr = next(r)
                if shdr is None:
                    shdr = hdr
                    sid = shdr.index('id')
                    w.writerow(shdr + mcols)
                for row in r:
                    rid = row[sid]
                    if rid in seen:
                        dupes += 1
                        continue
                    seen.add(rid)
                    total += 1
                    if src == gap:
                        gap_rows += 1
                    m = meta.get(rid)
                    if m is not None:
                        matched += 1
                    w.writerow(row + (m if m is not None else blank))
    pct = 100 * matched / total if total else 0
    print(f'{continent:10} {total:>7} rows (+{gap_rows} gap, {dupes} dup dropped)  '
          f'{matched:>7} matched ({pct:5.1f}%)  unmatched={total - matched}  -> {out}')


if __name__ == '__main__':
    meta, mcols = load_meta(META)
    print(f'loaded {len(meta)} metadata records, {len(mcols)} metadata columns\n')
    for c in (sys.argv[1:] or CONTINENTS):
        merge(c, meta, mcols)
