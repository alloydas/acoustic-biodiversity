"""Per-(cell, year) effort-controlled richness, plus year-over-year change.

Same metric as build_cells.py but split by calendar year (from metadata 'date').
Outputs:
  grid_cells_yearly.csv   - one row per (cell, year) with n_rec, S_obs, S_rare10
  cell_change.csv         - cells scored in >=2 years, with first/last richness + trend
Prints global and per-continent richness trend across 2023-2025.
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
M = 10  # rarefaction sample size


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


def rarefy(incidence, N, m):
    if N < m:
        return None
    logN = math.lgamma(N + 1) - math.lgamma(m + 1) - math.lgamma(N - m + 1)
    est = 0.0
    for ni in incidence.values():
        if N - ni < m:
            est += 1.0
        else:
            lc = math.lgamma(N - ni + 1) - math.lgamma(m + 1) - math.lgamma(N - ni - m + 1)
            est += 1.0 - math.exp(lc - logN)
    return est


# (cell, year) -> {continent, n, sp_incidence}
cy = collections.defaultdict(lambda: {'continent': None, 'n': 0, 'sp': collections.Counter()})

for cont, path in FILES.items():
    with open(path, newline='') as f:
        r = csv.reader(f)
        h = next(r); ci = {c: i for i, c in enumerate(h)}
        for row in r:
            key = cell_key(row[ci['lat']], row[ci['lon']])
            if key is None:
                continue
            d = row[ci['date']]
            if len(d) < 4 or not d[:4].isdigit():
                continue
            year = d[:4]
            st = cy[(key, year)]
            st['continent'] = cont
            st['n'] += 1
            sp = set()
            main = (row[ci['gen']] + ' ' + row[ci['sp']]).strip()
            if main:
                sp.add(main)
            for a in row[ci['also']].split(';'):
                a = a.strip()
                if a:
                    sp.add(a)
            for s in sp:
                st['sp'][s] += 1

# write per-(cell,year) table + collect scored
scored = collections.defaultdict(dict)  # cell -> {year: S_rare10}
trend = collections.defaultdict(lambda: collections.defaultdict(list))  # year -> continent -> [richness]
with open('grid_cells_yearly.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['lat_cell', 'lon_cell', 'continent', 'year', 'n_rec', 'S_obs', f'S_rare{M}'])
    for (key, year), st in cy.items():
        s_obs = len(st['sp'])
        s_rare = rarefy(st['sp'], st['n'], M)
        w.writerow([key[0], key[1], st['continent'], year, st['n'], s_obs,
                    '' if s_rare is None else round(s_rare, 3)])
        if s_rare is not None:
            scored[key][year] = s_rare
            trend[year][st['continent']].append(s_rare)

# coverage
cellyears = sum(1 for v in cy.values())
scored_cy = sum(len(v) for v in scored.values())
multi = {k: v for k, v in scored.items() if len(v) >= 2}
all3 = {k: v for k, v in scored.items() if len(v) == 3}
print(f'cell-years total: {cellyears}   scored (n>={M}): {scored_cy}')
print(f'cells scored in >=2 years: {len(multi)}   in all 3 years: {len(all3)}')

# global + continent trend (median richness of scored cells per year)
print('\nMedian effort-controlled richness by year (scored cells):')
years = ['2023', '2024', '2025']
allconts = sorted(FILES)
header = f"{'region':<11}" + ''.join(f'{y:>10}' for y in years)
print(header)
glob = {y: [] for y in years}
for cont in allconts:
    line = f'{cont:<11}'
    for y in years:
        vals = sorted(trend[y][cont]); glob[y].extend(vals)
        med = vals[len(vals) // 2] if vals else float('nan')
        line += f'{med:>10.1f}'
    print(line)
line = f'{"ALL":<11}'
for y in years:
    v = sorted(glob[y]); line += f'{v[len(v)//2]:>10.1f}'
print(line)

# per-cell change (first vs last scored year)
with open('cell_change.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['lat', 'lon', 'first_year', 'last_year', 'first_richness', 'last_richness', 'delta', 'direction'])
    up = down = flat = 0
    for key, yd in multi.items():
        ys = sorted(yd)
        fy, ly = ys[0], ys[-1]
        a, b = yd[fy], yd[ly]
        delta = b - a
        d = 'up' if delta > 0.5 else ('down' if delta < -0.5 else 'flat')
        if d == 'up':
            up += 1
        elif d == 'down':
            down += 1
        else:
            flat += 1
        w.writerow([round(key[0] + 0.05, 4), round(key[1] + 0.05, 4), fy, ly,
                    round(a, 2), round(b, 2), round(delta, 2), d])
print(f'\nCells tracked over time ({len(multi)}): up={up}  down={down}  flat={flat}')
print('wrote grid_cells_yearly.csv and cell_change.csv')
