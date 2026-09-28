"""Build the JSON data blocks embedded in the Silent Signal artifact page.

Reads the pipeline outputs (grid_cells*.csv, biome/urban tables, the merged
score_<c>_meta.csv files) and writes four JSON files into OUT:

  payload.json   correlations, biome + urban series, completeness, biome map
  trend.json     turnover-vs-change tests (same method as make_report.py page 12)
  play.json      per-year cell layers for the playable time-lapse map
  csv.json       first 100 rows of each published table, for the data appendix

Year pages use the report's convention: years with >= MIN_YEAR_CELLS scored
cell-years (1990-2025 on the merged 1886-2025 data). Index statistics use only
rows with index_source == 'local'; the merged merged_metadata_all.csv rows carry
indices from different settings (see merge_meta.py / build_cells.py).

Run from the processing directory, after the rest of the pipeline: it reads
grid_cells.csv, grid_cells_yearly.csv, cell_biome.csv, the biome/urban yearly
and change tables, score_<c>_meta.csv, the raw score_<c>[_gap|_hist].csv (for
provenance), recordings_ecoregion.csv, recordings_urban_class.csv and
merged_metadata_all.csv (for the index-offset comparison); biodiversity_map.csv
and cell_change.csv go into the data appendix. Optionally it also reads
build_cells_20260924_allsources.log: build_cells.py's stdout with
INDEX_SOURCES = {'local', 'extra'}, used for the "mixed sources" rho.
"""
import collections
import csv
import json
import math
import os
import random
import re
import sys

csv.field_size_limit(10 ** 8)

OUT = sys.argv[1] if len(sys.argv) > 1 else 'silent_signal_data'
CONTINENTS = ['africa', 'america', 'asia', 'australia', 'europe']
MIN_YEAR_CELLS = 20
MINC = 5          # biome/urban-year needs this many scored cells for a median
PAIR_MIN = 100

INDICES = [
    'Acoustic_Complexity_Index__mean', 'Acoustic_Diversity_Index__main_value',
    'Acoustic_Evenness_Index__main_value', 'Bio_acoustic_Index__main_value',
    'Normalized_Difference_Sound_Index__main_value', 'Spectral_Entropy__main_value',
    'Temporal_Entropy__main_value', 'NB_peaks__main_value',
    'Acoustic_Diversity_Index_NR__main_value', 'Bio_acoustic_Index_NR__main_value',
    'Spectral_Entropy_NR__main_value',
]
SHORT = {
    'Acoustic_Complexity_Index__mean': 'ACI', 'Acoustic_Diversity_Index__main_value': 'ADI',
    'Acoustic_Evenness_Index__main_value': 'AEI', 'Bio_acoustic_Index__main_value': 'Bioacoustic',
    'Normalized_Difference_Sound_Index__main_value': 'NDSI', 'Spectral_Entropy__main_value': 'Spectral Entropy',
    'Temporal_Entropy__main_value': 'Temporal Entropy', 'NB_peaks__main_value': 'NB peaks',
    'Acoustic_Diversity_Index_NR__main_value': 'ADI (NR)', 'Bio_acoustic_Index_NR__main_value': 'Bioacoustic (NR)',
    'Spectral_Entropy_NR__main_value': 'Spectral Entropy (NR)',
}


def fnum(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def med(v):
    s = sorted(v); n = len(s)
    if not n:
        return None
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def ranks(vals):
    n = len(vals)
    order = sorted(range(n), key=lambda i: vals[i])
    rk = [0.0] * n; i = 0
    while i < n:
        j = i
        while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        for k in range(i, j + 1):
            rk[order[k]] = (i + j) / 2.0 + 1
        i = j + 1
    return rk


def spearman(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    n = len(pairs)
    if n < 10:
        return None, n
    xr = ranks([p[0] for p in pairs]); yr = ranks([p[1] for p in pairs])
    mx = sum(xr) / n; my = sum(yr) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xr, yr))
    den = math.sqrt(sum((a - mx) ** 2 for a in xr) * sum((b - my) ** 2 for b in yr))
    return (num / den if den else None), n


def theil_sen(xs, ys, z=1.959963985):
    """Same estimator + CI rule as make_report.py (verified there against scipy)."""
    n = len(xs); sl = []
    for i in range(n):
        for j in range(i + 1, n):
            if xs[j] != xs[i]:
                sl.append((ys[j] - ys[i]) / (xs[j] - xs[i]))
    sl.sort(); nt = len(sl)
    m = sl[nt // 2] if nt % 2 else 0.5 * (sl[nt // 2 - 1] + sl[nt // 2])
    sig = math.sqrt(n * (n - 1) * (2 * n + 5) / 18.0)
    ru = min(int(round((nt + z * sig) / 2.0)), nt - 1)
    rl = max(int(round((nt - z * sig) / 2.0)) - 1, 0)
    return m, sl[rl], sl[ru]


def kendall(xs, ys):
    n = len(xs); c = d = 0
    for i in range(n):
        for j in range(i + 1, n):
            s = (xs[j] - xs[i]) * (ys[j] - ys[i])
            if s > 0: c += 1
            elif s < 0: d += 1
    S = c - d; n0 = n * (n - 1) / 2.0
    tx = list(collections.Counter(xs).values()); ty = list(collections.Counter(ys).values())
    n1 = sum(t * (t - 1) / 2.0 for t in tx); n2 = sum(u * (u - 1) / 2.0 for u in ty)
    tau = S / math.sqrt((n0 - n1) * (n0 - n2))
    v0 = n * (n - 1) * (2 * n + 5)
    vt = sum(t * (t - 1) * (2 * t + 5) for t in tx); vu = sum(u * (u - 1) * (2 * u + 5) for u in ty)
    v1 = (sum(t * (t - 1) for t in tx) * sum(u * (u - 1) for u in ty)) / (2.0 * n * (n - 1))
    v2 = (sum(t * (t - 1) * (t - 2) for t in tx) * sum(u * (u - 1) * (u - 2) for u in ty)) \
        / (9.0 * n * (n - 1) * (n - 2))
    var = (v0 - vt - vu) / 18.0 + v1 + v2
    return tau, math.erfc(abs(S / math.sqrt(var)) / math.sqrt(2))


def log10_two_sided_p(zabs):
    """log10 of erfc(z/sqrt2), stable far into the tail."""
    x = zabs / math.sqrt(2)
    p = math.erfc(x)
    if p > 1e-300:
        return math.log10(p)
    # asymptotic erfc(x) ~ exp(-x^2) / (x sqrt(pi)) * (1 - 1/(2x^2))
    return (-x * x - math.log(x * math.sqrt(math.pi)) + math.log(1 - 1 / (2 * x * x))) / math.log(10)


def mann_whitney(a, b):
    """Two-sided normal-approximation Mann-Whitney U with tie correction."""
    n1, n2 = len(a), len(b)
    rk = ranks(a + b)
    r1 = sum(rk[:n1])
    u1 = r1 - n1 * (n1 + 1) / 2.0
    N = n1 + n2
    ties = collections.Counter(a + b).values()
    tc = sum(t ** 3 - t for t in ties) / (N * (N - 1))
    sd = math.sqrt(n1 * n2 / 12.0 * ((N + 1) - tc))
    z = (u1 - n1 * n2 / 2.0) / sd
    return z, log10_two_sided_p(abs(z))


def rd(v, k=3):
    return None if v is None else round(v, k)


# ------------------------------------------------------------------ grid cells
with open('grid_cells.csv', newline='') as f:
    cells = list(csv.DictReader(f))
usable = [r for r in cells if r['S_rare10'] != '']
rich = [float(r['S_rare10']) for r in usable]
corrs = []
for idx in INDICES:
    rho, n = spearman([fnum(r[idx]) for r in usable], rich)
    corrs.append({'name': SHORT[idx], 'rho': round(rho, 4), 'n': n})
corrs.sort(key=lambda c: -abs(c['rho']))
n_idx_cells = max(c['n'] for c in corrs)

# rho when the merged rows' indices were (wrongly) averaged in, from that run's log
# (optional: offsets[].rho_mixed is null without it)
MIXED_LOG = 'build_cells_20260924_allsources.log'
mixed = {}
if os.path.exists(MIXED_LOG):
    with open(MIXED_LOG) as f:
        for line in f:
            m = re.match(r'(\S+)\s+([+-]\d\.\d+)\s+(\d+)$', line.strip())
            if m and m.group(1) in SHORT:
                mixed[SHORT[m.group(1)]] = float(m.group(2))
else:
    print(f'note: {MIXED_LOG} not found; rho_mixed left empty')

# ------------------------------------------------------------------ years
yc = collections.Counter()
cy_rows = []
with open('grid_cells_yearly.csv', newline='') as f:
    for r in csv.DictReader(f):
        cy_rows.append(r)
        if r['S_rare10'] != '':
            yc[r['year']] += 1
RY = sorted(y for y, n in yc.items() if n >= MIN_YEAR_CELLS)
RYS = set(RY)

# ------------------------------------------------------------------ biomes
with open('biome_change.csv', newline='') as f:
    bchg = list(csv.DictReader(f))
with open('biome_yearly.csv', newline='') as f:
    byl = {(r['biome_name'], r['year']): r for r in csv.DictReader(f)}
biomes = []
raw_series = {}   # full-precision medians for the Kendall test (as make_report page 8)
for b in bchg:
    name = b['biome_name']
    series, effort, raw = [], [], []
    for y in RY:
        r = byl.get((name, y))
        ok = r and r['median_S_rare10'] != '' and int(r['n_scored_cells']) >= MINC
        raw.append(float(r['median_S_rare10']) if ok else None)
        series.append(round(raw[-1], 2) if ok else None)
        effort.append(int(r['n_recordings']) if r else 0)
    raw_series[name] = raw
    biomes.append({'name': name, 'series': series, 'effort': effort,
                   'delta': round(float(b['delta']), 3), 'first_year': b['first_year'],
                   'last_year': b['last_year'], 'first_median': round(float(b['first_median']), 2),
                   'last_median': round(float(b['last_median']), 2), 'rec': int(b['total_recordings']),
                   'up': int(b['cells_up']), 'down': int(b['cells_down']), 'flat': int(b['cells_flat'])})
biome_order = [b['biome_name'] for b in bchg]
btests = []
for b in biomes:
    pts = [(int(y), v) for y, v in zip(RY, raw_series[b['name']]) if v is not None]
    if len(pts) >= 5:
        tau, p = kendall([x for x, _ in pts], [v for _, v in pts])
        btests.append({'name': b['name'], 'n_years': len(pts), 'tau': round(tau, 3), 'p': round(p, 4)})
        b['tau'], b['p'] = round(tau, 3), round(p, 4)
bonf = 0.05 / len(btests)

# ------------------------------------------------------------------ urban
with open('urban_change.csv', newline='') as f:
    uchg = {r['urban_class']: r for r in csv.DictReader(f)}
with open('urban_yearly.csv', newline='') as f:
    uyl = list(csv.DictReader(f))
ulook = {(r['urban_class'], r['year']): r for r in uyl}
lab_years = sorted({r['year'] for r in uyl if r['class_carried_forward'] == '0'})
urban = []
for u in ['city', 'town', 'rural']:
    series, effort = [], []
    for y in RY:
        r = ulook.get((u, y))
        ok = r and r['median_S_rare10'] != '' and int(r['n_scored_cells']) >= MINC
        series.append(round(float(r['median_S_rare10']), 2) if ok else None)
        effort.append(int(r['n_recordings']) if r else 0)
    c = uchg[u]
    urban.append({'name': u, 'series': series, 'effort': effort,
                  'delta_lab': round(float(c['delta_labelled_window']), 3),
                  'up': int(c['cells_up']), 'down': int(c['cells_down']), 'flat': int(c['cells_flat'])})

# ------------------------------------------------------------------ per-recording pass
# index_source + 40-metric completeness + two urban index means, from the merged meta files
eco = {}
with open('recordings_ecoregion.csv', newline='') as f:
    for r in csv.DictReader(f):
        eco[r['id']] = (r['biome_name'], r['realm'], fnum(r['lat']), fnum(r['lon']))
uclass = {}
with open('recordings_urban_class.csv', newline='') as f:
    for r in csv.DictReader(f):
        uclass[r['id']] = r['urban_class']

hist_ids = set()
for c in CONTINENTS:
    if os.path.exists(f'score_{c}_hist.csv'):
        with open(f'score_{c}_hist.csv', newline='') as f:
            hist_ids |= {row[1] for row in csv.reader(f)}
prov = collections.defaultdict(collections.Counter)   # era -> source -> n (geolocated, dated)


def era_of(y):
    return None if y is None or y < 1990 or y > 2025 else ('1990-2014' if y < 2015 else '2015-2022' if y < 2023 else '2023-2025')


comp_b = collections.defaultdict(lambda: [0, 0])     # biome -> [n_local, complete]
comp_r = collections.defaultdict(lambda: [0, 0])
urb_ix = collections.defaultdict(lambda: {'aci': [], 'ndsi': []})
n_src = collections.Counter()
n_local_all = n_local_complete = 0
for c in CONTINENTS:
    with open(f'score_{c}_meta.csv', newline='') as f:
        rd_ = csv.reader(f)
        h = next(rd_); ci = {k: i for i, k in enumerate(h)}
        icols = [ci[k] for k in h[h.index('id') + 1: h.index('gen')]]
        assert len(icols) == 40, len(icols)
        isrc, iid = ci['index_source'], ci['id']
        iaci, indsi = ci['Acoustic_Complexity_Index__mean'], ci['Normalized_Difference_Sound_Index__main_value']
        idate, ilat, ilon = ci['date'], ci['lat'], ci['lon']
        for row in rd_:
            src = row[isrc]; n_src[src] += 1
            d_ = row[idate]
            era = era_of(int(d_[:4]) if d_[:4].isdigit() else None)
            if era and fnum(row[ilat]) is not None and fnum(row[ilon]) is not None:
                prov[era]['merged file' if src != 'local' else
                          'historical backfill' if row[iid] in hist_ids else 'base + gap runs'] += 1
            if src != 'local':
                continue
            complete = all(fnum(row[k]) is not None for k in icols)
            n_local_all += 1; n_local_complete += complete
            e = eco.get(row[iid])
            if e and e[0] and e[0] != 'N/A':
                comp_b[e[0]][0] += 1; comp_b[e[0]][1] += complete
                comp_r[e[1]][0] += 1; comp_r[e[1]][1] += complete
            u = uclass.get(row[iid])
            if u:
                a, nd = fnum(row[iaci]), fnum(row[indsi])
                if a is not None: urb_ix[u]['aci'].append(a)
                if nd is not None: urb_ix[u]['ndsi'].append(nd)
completeness = [{'name': b, 'n': comp_b[b][0], 'complete': comp_b[b][1],
                 'pct': round(100.0 * comp_b[b][1] / comp_b[b][0], 1)} for b in biome_order if comp_b[b][0]]
realms = sorted(({'name': k, 'n': v[0], 'complete': v[1], 'pct': round(100.0 * v[1] / v[0], 1)}
                 for k, v in comp_r.items() if k and k != 'N/A'), key=lambda r: -r['n'])
urban_ix = {u: {'n': len(v['aci']), 'aci_mean': round(sum(v['aci']) / len(v['aci']), 1),
                'ndsi_mean': round(sum(v['ndsi']) / len(v['ndsi']), 3)} for u, v in urb_ix.items()}

# ------------------------------------------------------------------ biome map (0.5 deg)
votes = collections.defaultdict(collections.Counter)
n_biome_labelled = 0
n_coord = sum(1 for v in eco.values() if v[2] is not None and v[3] is not None)
biome_rec = collections.Counter()
for bname, realm, la, lo in eco.values():
    if not bname or bname == 'N/A' or la is None or lo is None:
        continue
    n_biome_labelled += 1
    biome_rec[bname] += 1
    votes[(math.floor(la * 2) / 2, math.floor(lo * 2) / 2)][bname] += 1
bidx = {b: i for i, b in enumerate(biome_order)}
bmap = []
for (la, lo), cnt in votes.items():
    top = cnt.most_common(1)[0][0]
    if top in bidx:
        bmap.append([la, lo, bidx[top], sum(cnt.values())])
bmap.sort(key=lambda p: -p[3])
cellb = collections.Counter()
mixed_cells = n_cb = 0
with open('cell_biome.csv', newline='') as f:
    for r in csv.DictReader(f):
        cellb[r['biome_name']] += 1; n_cb += 1
        mixed_cells += int(r['n_biomes_in_cell']) > 1
map_counts = collections.Counter(p[2] for p in bmap)
map_rec = collections.Counter()
for p in bmap:
    map_rec[p[2]] += p[3]
tot_map_rec = sum(map_rec.values())

# ------------------------------------------------------------------ trend (page-12 method)
recs = []
for r in cy_rows:
    if r['S_rare10'] == '' or r['year'] not in RYS:
        continue
    recs.append((round(float(r['lat_cell']), 1), round(float(r['lon_cell']), 1),
                 int(r['year']), float(r['S_rare10']), int(r['n_rec'])))
TY = sorted({r[2] for r in recs})
byy = collections.defaultdict(list)
for r in recs:
    byy[r[2]].append(r[3])
raw = [med(byy[y]) for y in TY]
r_sl, r_lo, r_hi = theil_sen(TY, raw)
r_tau, r_p = kendall(TY, raw)
seen = collections.defaultdict(set)
for r in recs:
    seen[(r[0], r[1])].add(r[2])


def panel(min_y):
    keep = {c for c, ys in seen.items() if len(ys) >= min_y}
    b = collections.defaultdict(list)
    for r in recs:
        if (r[0], r[1]) in keep:
            b[r[2]].append(r[3])
    yy = [y for y in TY if len(b[y]) >= 5]
    return keep, yy, [med(b[y]) for y in yy]


forest = [{'lab': 'all scored cells', 'n': len(seen), 'est': rd(r_sl * 10), 'lo': rd(r_lo * 10),
           'hi': rd(r_hi * 10), 'raw': True, 'kind': 'slope'}]
for m in (2, 4, 6, 8):
    kk, yy, ss = panel(m)
    sl, lo, hi = theil_sen(yy, ss)
    forest.append({'lab': f'same cells, ≥{m} yrs', 'n': len(kk), 'est': rd(sl * 10),
                   'lo': rd(lo * 10), 'hi': rd(hi * 10), 'raw': False, 'kind': 'slope'})
k8, y8, s8 = panel(8)
p8 = theil_sen(y8, s8)

cyc = collections.Counter(r[2] for r in recs)
e_lo = next((y for y in TY if all(cyc[y + k] >= PAIR_MIN for k in range(3))), TY[0])
e_hi, l_lo, l_hi = e_lo + 2, TY[-3], TY[-1]
ear = collections.defaultdict(list); lat_ = collections.defaultdict(list)
for r in recs:
    if e_lo <= r[2] <= e_hi: ear[(r[0], r[1])].append(r[3])
    elif l_lo <= r[2] <= l_hi: lat_[(r[0], r[1])].append(r[3])
both = sorted(set(ear) & set(lat_))
dif = [sum(lat_[c]) / len(lat_[c]) - sum(ear[c]) / len(ear[c]) for c in both]
np_ = len(dif); p_mean = sum(dif) / np_
p_sd = math.sqrt(sum((x - p_mean) ** 2 for x in dif) / (np_ - 1))
p_ci = 1.959963985 * p_sd / math.sqrt(np_)
gap_dec = ((l_lo + l_hi) - (e_lo + e_hi)) / 2 / 10.0      # decades between window midpoints
forest.append({'lab': f'paired {e_lo}–{e_hi} vs {l_lo}–{l_hi}', 'n': np_,
               'est': rd(p_mean / gap_dec), 'lo': rd((p_mean - p_ci) / gap_dec), 'hi': rd((p_mean + p_ci) / gap_dec),
               'total': rd(p_mean), 'tlo': rd(p_mean - p_ci), 'thi': rd(p_mean + p_ci),
               'decades': round(gap_dec, 2), 'raw': False, 'kind': 'window'})

first = {}   # first scored year per cell over EVERY year: a cell scored in 1985 is not new in 1990
for r in cy_rows:
    if r['S_rare10'] != '':
        k = (round(float(r['lat_cell']), 1), round(float(r['lon_cell']), 1))
        first[k] = min(first.get(k, 9999), int(r['year']))
nw, rt, pnew, all_new, all_ret = [], [], [], [], []
for y in TY:
    a = [r[3] for r in recs if r[2] == y and first[(r[0], r[1])] == y]
    b = [r[3] for r in recs if r[2] == y and first[(r[0], r[1])] < y]
    nw.append(rd(med(a))); rt.append(rd(med(b)))
    pnew.append(round(100.0 * len(a) / (len(a) + len(b)), 1))
    all_new += a; all_ret += b
n_lower = sum(1 for a, b in zip(nw, rt) if a is not None and b is not None and a < b)
n_cmp = sum(1 for a, b in zip(nw, rt) if a is not None and b is not None)
mw_z, mw_log10p = mann_whitney(all_new, all_ret)
# sign test across years (each year one vote): immune to cells recurring across years
_k, _n = n_lower, n_cmp
sign_p = min(1.0, 2 * sum(math.comb(_n, i) for i in range(max(_k, _n - _k), _n + 1)) / 2 ** _n)

trend = {
    'years': [str(y) for y in TY], 'raw': [rd(v) for v in raw],
    'panel8': {'years': [str(y) for y in y8], 'series': [rd(v) for v in s8], 'n': len(k8),
               'dec': rd(p8[0] * 10), 'lo': rd(p8[1] * 10), 'hi': rd(p8[2] * 10)},
    'raw_dec': rd(r_sl * 10), 'raw_lo': rd(r_lo * 10), 'raw_hi': rd(r_hi * 10),
    'raw_tau': rd(r_tau), 'raw_p': round(r_p, 4), 'forest': forest,
    'new': nw, 'ret': rt, 'pct_new': pnew,
    'pct_new_range': [min(pnew), max(pnew)],
    'paired': {'n': np_, 'mean': rd(p_mean), 'lo': rd(p_mean - p_ci), 'hi': rd(p_mean + p_ci),
               'up': sum(1 for x in dif if x > 0), 'down': sum(1 for x in dif if x < 0),
               'early': [e_lo, e_hi], 'late': [l_lo, l_hi]},
    'n_lower': n_lower, 'n_cmp': n_cmp, 'sign_p': sign_p,
    'mw': {'z': round(mw_z, 2), 'log10p': round(mw_log10p, 1), 'new_med': rd(med(all_new)),
           'ret_med': rd(med(all_ret)), 'n_new': len(all_new), 'n_ret': len(all_ret)},
    'min_year_cells': MIN_YEAR_CELLS, 'pair_min': PAIR_MIN,
}

# ------------------------------------------------------------------ playable layers
# scored cell-years (0.1 deg) with a new/returning flag, and the effort footprint
# of every recorded cell-year aggregated to 0.5 deg.
yi = {y: i for i, y in enumerate(RY)}
play_scored = []   # [lat*10, lon*10, yearIdx, n_rec, S_rare10*100, isNew]
for r in recs:
    play_scored.append([int(round(r[0] * 10)), int(round(r[1] * 10)), yi[str(r[2])], r[4],
                        int(round(r[3] * 100)), 1 if first[(r[0], r[1])] == r[2] else 0])
eff = collections.Counter()
rec_year = collections.Counter(); cells_year = collections.Counter()
for r in cy_rows:
    if r['year'] not in RYS:
        continue
    la, lo = float(r['lat_cell']), float(r['lon_cell'])
    n = int(r['n_rec'])
    rec_year[r['year']] += n; cells_year[r['year']] += 1
    eff[(int(math.floor(la * 2)), int(math.floor(lo * 2)), yi[r['year']])] += n
play_eff = [[k[0], k[1], k[2], v] for k, v in eff.items()]   # lat*2, lon*2 (0.5-deg index)
_nn = sorted(r[3] for r in play_scored)
play = {
    'years': RY, 'n_cap': _nn[int(0.99 * (len(_nn) - 1))],
    'scored': play_scored,
    'effort': play_eff,
    'per_year': [{'year': y, 'recordings': rec_year[y], 'cells': cells_year[y], 'scored': yc[y],
                  'median': rd(raw[TY.index(int(y))]), 'new': nw[TY.index(int(y))],
                  'ret': rt[TY.index(int(y))], 'pct_new': pnew[TY.index(int(y))]} for y in RY],
}

# ------------------------------------------------------------------ index offsets between sources
loc = {}
for c in CONTINENTS:
    for src in (f'score_{c}.csv', f'score_{c}_gap.csv', f'score_{c}_hist.csv'):
        if not os.path.exists(src):
            continue
        with open(src, newline='') as f:
            for r in csv.DictReader(f):
                loc[r['id']] = [fnum(r[k]) for k in INDICES]
pairs = {k: [] for k in INDICES}; extra_only = {k: [] for k in INDICES}
with open('merged_metadata_all.csv', newline='') as f:
    for r in csv.DictReader(f):
        o = loc.get(r['id'])
        for i, k in enumerate(INDICES):
            b = fnum(r[k])
            if b is None:
                continue
            if o is None:
                extra_only[k].append(b)
            elif o[i] is not None:
                pairs[k].append((o[i], b))
offsets = []
for k in INDICES:
    p = pairs[k]
    smp = random.Random(0).sample(p, min(len(p), 100000))
    rho, _ = spearman([a for a, _ in smp], [b for _, b in smp])
    offsets.append({'name': SHORT[k], 'local_med': rd(med([a for a, _ in p]), 4),
                    'extra_med': rd(med([b for _, b in p]), 4), 'rank_rho': rd(rho),
                    'extra_only_med': rd(med(extra_only[k]), 4), 'n_pairs': len(p),
                    'rho_mixed': mixed.get(SHORT[k]),
                    'rho_local': next(c['rho'] for c in corrs if c['name'] == SHORT[k])})
offsets.sort(key=lambda o: -abs(o['rho_mixed'] or 0))

# ------------------------------------------------------------------ headline counts
n_rec_total = sum(n_src.values())
payload = {
    'corrs': corrs, 'n_cells': len(cells), 'n_usable': len(usable), 'n_idx_cells': n_idx_cells,
    'years': RY, 'biomes': biomes, 'biome_order': biome_order,
    'urban': urban, 'urban_first_labelled': lab_years[0], 'urban_last_labelled': lab_years[-1],
    'urban_ix': urban_ix,
    'completeness': completeness, 'realms': realms,
    'map': bmap, 'map_cells': len(bmap), 'map_counts': {biome_order[k]: v for k, v in map_counts.items()},
    'map_rec_share_top': round(100.0 * map_rec[0] / tot_map_rec, 1),
    'cell_biome': {'cells': n_cb, 'mixed': mixed_cells, 'top': cellb[biome_order[0]]},
    'offsets': offsets,
    'biome_tests': {'tests': btests, 'bonferroni': round(bonf, 5),
                    'n_nominal': sum(1 for t in btests if t['p'] < 0.05),
                    'n_bonf': sum(1 for t in btests if t['p'] < bonf),
                    'min_p': min(t['p'] for t in btests)},
    'provenance': {era: {k: v for k, v in cnt.items()} for era, cnt in sorted(prov.items())},
    'coverage': {'with_coords': n_coord, 'with_biome': n_biome_labelled,
                 'top_biome_share': round(100.0 * biome_rec[biome_order[0]] / n_biome_labelled, 1)},
    'counts': {'recordings': n_rec_total, 'local': n_src['local'], 'extra': n_src['extra'],
               'local_complete_pct': round(100.0 * n_local_complete / n_local_all, 1),
               'biome_labelled': n_biome_labelled, 'n_biomes': len(biome_order),
               'years_scored': [RY[0], RY[-1]], 'n_years': len(RY)},
}

# ------------------------------------------------------------------ csv appendix
CSVS = [
    ('biodiversity_map.csv', 'The published map: one row per confidently-scored 0.1° cell.'),
    ('grid_cells.csv', 'Core aggregation: every cell, with the 11 index means build_cells computes '
                       '(from this pipeline’s recordings only).'),
    ('grid_cells_yearly.csv', 'The same cells split by year, 1886–2025; effort and richness only, no indices.'),
    ('cell_biome.csv', 'Each cell’s modal biome, and how many biomes it spans.'),
    ('biome_yearly.csv', 'Median richness per biome per year: the source of the heatmap.'),
    ('biome_change.csv', 'First-to-last scored year delta per biome, with the tracked-cell split. Complete file.'),
    ('urban_change.csv', 'City / town / rural change, over the labelled 2015–2022 window and in full. Complete file.'),
    ('cell_change.csv', 'Per-cell first-to-last richness change, cells tracked 2+ years.'),
]
csvdata = []
for fn, desc in CSVS:
    with open(fn, newline='') as f:
        r = csv.reader(f)
        cols = next(r)
        rows, n = [], 0
        for row in r:
            n += 1
            if len(rows) < 100:
                rows.append(row)
    csvdata.append({'file': fn, 'desc': desc, 'cols': cols, 'rows': rows, 'n_rows': n,
                    'n_cols': len(cols), 'bytes': os.path.getsize(fn)})

os.makedirs(OUT, exist_ok=True)
for name, obj in (('payload', payload), ('trend', trend), ('play', play), ('csv', csvdata)):
    with open(os.path.join(OUT, name + '.json'), 'w') as f:
        json.dump(obj, f, separators=(',', ':'), ensure_ascii=False)
    print(f'{name}.json  {os.path.getsize(os.path.join(OUT, name + ".json")):,} bytes')

print('\ncounts', payload['counts'])
print('corrs', [(c['name'], c['rho'], c['n']) for c in corrs])
print('years', RY[0], RY[-1], 'map cells', len(bmap), 'top biome rec share', payload['map_rec_share_top'])
print('cell_biome', payload['cell_biome'], 'map_counts top/bottom3',
      payload['map_counts'].get(biome_order[0]), sum(payload['map_counts'].get(b, 0) for b in biome_order[-3:]))
print('completeness', [(c['name'][:18], c['pct']) for c in completeness])
print('realms', [(r['name'], r['pct']) for r in realms])
print('urban_ix', urban_ix)
print('trend raw_dec', trend['raw_dec'], trend['raw_lo'], trend['raw_hi'], 'tau', trend['raw_tau'], 'p', trend['raw_p'])
print('panel8', trend['panel8']['n'], trend['panel8']['dec'], trend['panel8']['lo'], trend['panel8']['hi'])
print('paired', trend['paired'])
print('pct_new_range', trend['pct_new_range'], 'n_lower', n_lower, '/', n_cmp, 'mw', trend['mw'])
print('forest', [(f_['lab'], f_['n'], f_['est'], f_['lo'], f_['hi']) for f_ in forest])
print('offsets', [(o['name'], o['local_med'], o['extra_med'], o['rank_rho'], o['rho_mixed'], o['rho_local']) for o in offsets])
print('biome deltas', min(b['delta'] for b in biomes), max(b['delta'] for b in biomes),
      'tracked up/down', sum(b['up'] for b in biomes), sum(b['down'] for b in biomes))
print('play scored', len(play_scored), 'effort', len(play_eff), 'n_cap', play['n_cap'])
print('biome tests', payload['biome_tests'])
print('provenance', payload['provenance'])
print('coverage', payload['coverage'])
print('sign_p', sign_p, 'forest paired', forest[-1])
