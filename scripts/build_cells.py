"""Aggregate all continents into 0.1-degree grid cells and test which acoustic
indices track (effort-controlled) species richness.

Per cell we compute:
  - n_rec          : number of recordings
  - S_obs          : total distinct species observed (gen+sp plus the 'also' list)
  - S_rare10       : sample-based (Hurlbert) rarefied richness at 10 recordings
                     -> richness controlled for sampling effort
  - mean of each diversity-relevant acoustic index across the cell's recordings

Then prints Spearman correlation of each acoustic index vs S_rare10 on cells
with n_rec >= 10 (where rarefaction is defined).
"""
import csv
import math
import collections

csv.field_size_limit(10 ** 7)

FILES = {
    'africa': 'score_africa_meta.csv', 'america': 'score_america_meta.csv',
    'asia': 'score_asia_meta.csv', 'australia': 'score_australia_meta.csv',
    'europe': 'score_europe_meta.csv',
}

# Acoustic indices most plausibly tied to biodiversity (cell-mean aggregated).
INDICES = [
    'Acoustic_Complexity_Index__mean', 'Acoustic_Diversity_Index__main_value',
    'Acoustic_Evenness_Index__main_value', 'Bio_acoustic_Index__main_value',
    'Normalized_Difference_Sound_Index__main_value', 'Spectral_Entropy__main_value',
    'Temporal_Entropy__main_value', 'NB_peaks__main_value',
    'Acoustic_Diversity_Index_NR__main_value', 'Bio_acoustic_Index_NR__main_value',
    'Spectral_Entropy_NR__main_value',
]


def fnum(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def cell_key(lat, lon):
    la, lo = fnum(lat), fnum(lon)
    if la is None or lo is None:
        return None
    return (math.floor(la * 10) / 10, math.floor(lo * 10) / 10)


# cell -> aggregation state
cells = collections.defaultdict(lambda: {
    'continent': None, 'n': 0,
    'sums': collections.defaultdict(float), 'cnts': collections.defaultdict(int),
    'sp_incidence': collections.Counter(),  # species -> # recordings containing it
})

total = skipped_coord = 0
for cont, path in FILES.items():
    with open(path, newline='') as f:
        r = csv.reader(f)
        h = next(r)
        ci = {c: i for i, c in enumerate(h)}
        for row in r:
            total += 1
            key = cell_key(row[ci['lat']], row[ci['lon']])
            if key is None:
                skipped_coord += 1
                continue
            c = cells[key]
            c['continent'] = cont
            c['n'] += 1
            # species in this recording: main + co-occurring
            sp = set()
            main = (row[ci['gen']] + ' ' + row[ci['sp']]).strip()
            if main:
                sp.add(main)
            for a in row[ci['also']].split(';'):
                a = a.strip()
                if a:
                    sp.add(a)
            for s in sp:
                c['sp_incidence'][s] += 1
            for idx in INDICES:
                v = fnum(row[ci[idx]])
                if v is not None and not math.isnan(v):
                    c['sums'][idx] += v
                    c['cnts'][idx] += 1

print(f'recordings read: {total}   skipped (no coord): {skipped_coord}')
print(f'grid cells (0.1 deg): {len(cells)}')


def rarefy(incidence, N, m):
    """Hurlbert sample-based rarefaction: expected #species in m of N recordings."""
    if N < m:
        return None
    log_choose_N = math.lgamma(N + 1) - math.lgamma(m + 1) - math.lgamma(N - m + 1)
    est = 0.0
    for ni in incidence.values():
        if N - ni < m:
            est += 1.0
        else:
            log_c = math.lgamma(N - ni + 1) - math.lgamma(m + 1) - math.lgamma(N - ni - m + 1)
            est += 1.0 - math.exp(log_c - log_choose_N)
    return est


# Build per-cell records and write table
M = 10
out_rows = []
with open('grid_cells.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['lat_cell', 'lon_cell', 'continent', 'n_rec', 'S_obs', f'S_rare{M}'] + INDICES)
    for (la, lo), c in cells.items():
        s_obs = len(c['sp_incidence'])
        s_rare = rarefy(c['sp_incidence'], c['n'], M)
        means = [(c['sums'][i] / c['cnts'][i]) if c['cnts'][i] else '' for i in INDICES]
        w.writerow([la, lo, c['continent'], c['n'], s_obs,
                    '' if s_rare is None else round(s_rare, 3)] + means)
        if s_rare is not None:
            out_rows.append((s_rare, means, c['n'], s_obs))

print(f"cells with n_rec >= {M} (usable for validation): {len(out_rows)}")


def spearman(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if x != '' and y != '']
    if len(pairs) < 10:
        return None, 0
    n = len(pairs)

    def ranks(vals):
        order = sorted(range(n), key=lambda i: vals[i])
        rk = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                rk[order[k]] = avg
            i = j + 1
        return rk

    xr = ranks([p[0] for p in pairs])
    yr = ranks([p[1] for p in pairs])
    mx = sum(xr) / n
    my = sum(yr) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xr, yr))
    den = math.sqrt(sum((a - mx) ** 2 for a in xr) * sum((b - my) ** 2 for b in yr))
    return (num / den if den else None), n


print(f"\nSpearman correlation: cell-mean acoustic index vs effort-controlled richness (S_rare{M})")
print(f"{'index':<48} {'rho':>7} {'n_cells':>8}")
ys = [r[0] for r in out_rows]
results = []
for j, idx in enumerate(INDICES):
    xs = [r[1][j] for r in out_rows]
    rho, n = spearman(xs, ys)
    results.append((idx, rho))
    print(f"{idx:<48} {('%+.3f' % rho) if rho is not None else '  n/a':>7} {n:>8}")
