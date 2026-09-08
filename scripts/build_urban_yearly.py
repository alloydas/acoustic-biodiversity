#!/usr/bin/env python
"""Year-by-year effort-controlled richness, split by URBAN CLASS (city/town/rural).

Mirrors build_biome_yearly.py. A grid cell is assigned the MODAL urban class of
the recordings inside it.

IMPORTANT COVERAGE NOTE. The GCTB built-up polygons stop at 2022, so
classify_urban.py only labels 2015-2022 recordings. Urban class is a property of
the PLACE, not of the recording year, so the cell-level class derived from
2015-2022 is carried forward to 2023-2025 cell-years. Those three years are
flagged `class_carried_forward=1` in the output and drawn differently in the
report, because no GCTB polygon actually covers them.

Inputs
  recordings_urban_class.csv   per-recording lat/lon/year/urban_class (2015-2022)
  grid_cells_yearly.csv        per (cell, year): n_rec, S_obs, S_rare10
  cell_change.csv              per tracked cell: first/last richness, direction

Outputs
  urban_yearly.csv   per (urban_class, year): n_cells, n_rec, median/mean S_rare10
  urban_change.csv   per class: first vs last year median, delta, tracked-cell split
"""
import csv
import math
import collections

csv.field_size_limit(10 ** 9)

MIN_CELLS = 5
LABEL_YEARS = ('2015', '2022')   # window in which GCTB actually labels recordings
CLASSES = ['city', 'town', 'rural']


def ikey(lat, lon):
    return (math.floor(lat * 10), math.floor(lon * 10))


def ikey_from_cell(lat_cell, lon_cell):
    return (round(float(lat_cell) * 10), round(float(lon_cell) * 10))


# ---- cell -> modal urban class ----
votes = collections.defaultdict(collections.Counter)
n_lab = 0
with open('recordings_urban_class.csv', newline='') as f:
    r = csv.reader(f)
    h = next(r); ci = {c: i for i, c in enumerate(h)}
    for row in r:
        u = row[ci['urban_class']]
        if not u:
            continue
        try:
            k = ikey(float(row[ci['lat']]), float(row[ci['lon']]))
        except ValueError:
            continue
        votes[k][u] += 1
        n_lab += 1
cell_class = {k: c.most_common(1)[0][0] for k, c in votes.items()}
mixed = sum(1 for c in votes.values() if len(c) > 1)
pure = sum(c.most_common(1)[0][1] for c in votes.values()) / max(n_lab, 1)
print(f'labelled recordings ({LABEL_YEARS[0]}-{LABEL_YEARS[1]}): {n_lab:,}')
print(f'cells with a class: {len(cell_class):,}   mixed-class cells: {mixed:,} '
      f'({100*mixed/max(len(cell_class),1):.1f}%)')
print(f"recordings in their cell's modal class: {100*pure:.1f}%")
print('class distribution over cells: ' +
      '  '.join(f'{c}={sum(1 for v in cell_class.values() if v == c):,}' for c in CLASSES))

# ---- class x year aggregation ----
by = collections.defaultdict(lambda: {'rich': [], 'n_rec': 0, 'cells': 0})
years, unmatched = set(), 0
with open('grid_cells_yearly.csv', newline='') as f:
    r = csv.reader(f)
    h = next(r); ci = {c: i for i, c in enumerate(h)}
    for row in r:
        k = ikey_from_cell(row[ci['lat_cell']], row[ci['lon_cell']])
        u = cell_class.get(k)
        if u is None:
            unmatched += 1
            continue
        y = row[ci['year']]
        years.add(y)
        st = by[(u, y)]
        st['n_rec'] += int(row[ci['n_rec']])
        st['cells'] += 1
        s = row[ci['S_rare10']]
        if s != '':
            st['rich'].append(float(s))
print(f'cell-years with no class match: {unmatched:,}')

YEARS = sorted(years)


def med(v):
    v = sorted(v)
    if not v:
        return None
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


with open('urban_yearly.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['urban_class', 'year', 'n_cell_years', 'n_recordings',
                'n_scored_cells', 'median_S_rare10', 'mean_S_rare10',
                'class_carried_forward'])
    for u in CLASSES:
        for y in YEARS:
            st = by.get((u, y))
            if not st:
                continue
            m = med(st['rich'])
            w.writerow([u, y, st['cells'], st['n_rec'], len(st['rich']),
                        '' if m is None else round(m, 3),
                        '' if not st['rich'] else round(sum(st['rich']) / len(st['rich']), 3),
                        1 if y > LABEL_YEARS[1] else 0])
print('wrote urban_yearly.csv')

# ---- per-class change + tracked-cell direction split ----
dir_by_class = collections.defaultdict(collections.Counter)
with open('cell_change.csv', newline='') as f:
    r = csv.reader(f)
    h = next(r); ci = {c: i for i, c in enumerate(h)}
    for row in r:
        # cell_change stores CELL CENTRES (lat_cell + 0.05)
        k = ikey(float(row[ci['lat']]) - 0.05, float(row[ci['lon']]) - 0.05)
        u = cell_class.get(k)
        if u:
            dir_by_class[u][row[ci['direction']]] += 1

rows = []
for u in CLASSES:
    ser = [(y, med(by[(u, y)]['rich'])) for y in YEARS
           if (u, y) in by and len(by[(u, y)]['rich']) >= MIN_CELLS]
    if len(ser) < 2:
        continue
    # headline delta over the LABELLED window only (2015-2022), where the class is real
    lab = [(y, v) for y, v in ser if y <= LABEL_YEARS[1]]
    d = dir_by_class[u]
    rows.append({
        'urban_class': u,
        'first_year': ser[0][0], 'last_year': ser[-1][0],
        'first_median': round(ser[0][1], 3), 'last_median': round(ser[-1][1], 3),
        'delta_full': round(ser[-1][1] - ser[0][1], 3),
        'labelled_first_year': lab[0][0], 'labelled_last_year': lab[-1][0],
        'delta_labelled_window': round(lab[-1][1] - lab[0][1], 3),
        'total_recordings': sum(by[(u, y)]['n_rec'] for y in YEARS if (u, y) in by),
        'cells_up': d['up'], 'cells_down': d['down'], 'cells_flat': d['flat'],
    })
with open('urban_change.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print('wrote urban_change.csv')

print(f'\nMedian effort-controlled richness by urban class x year '
      f'(| marks the end of real GCTB labelling)')
print(f"{'class':<8}" + ''.join(f"{y[2:]:>7}" + ('  |' if y == LABEL_YEARS[1] else '')
                                for y in YEARS))
for u in CLASSES:
    line = f'{u:<8}'
    for y in YEARS:
        st = by.get((u, y))
        line += (f"{med(st['rich']):>7.1f}" if st and len(st['rich']) >= MIN_CELLS else f"{'-':>7}")
        if y == LABEL_YEARS[1]:
            line += '  |'
    print(line)
print(f"\n{'class':<8}{'up':>7}{'down':>7}{'flat':>7}   (individually tracked cells)")
for r_ in rows:
    print(f"{r_['urban_class']:<8}{r_['cells_up']:>7}{r_['cells_down']:>7}{r_['cells_flat']:>7}")
