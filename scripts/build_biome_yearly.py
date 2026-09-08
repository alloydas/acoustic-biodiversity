#!/usr/bin/env python
"""Year-by-year effort-controlled richness, split by BIOME.

Joins the ecoregion labels onto the per-(cell, year) richness table. A grid cell
is assigned the MODAL biome of the recordings inside it (cells straddling a biome
boundary go to whichever biome contributed most recordings).

Inputs
  recordings_ecoregion.csv   per-recording lat/lon/year/biome_name
  grid_cells_yearly.csv      per (cell, year): n_rec, S_obs, S_rare10
  cell_change.csv            per tracked cell: first/last richness, direction

Outputs
  biome_yearly.csv   per (biome, year): n_cells, n_rec, median/mean S_rare10
  biome_change.csv   per biome: first vs last year median, delta, and the
                     up/down/flat split of its individually tracked cells
"""
import csv
import math
import collections

csv.field_size_limit(10 ** 9)

MIN_CELLS = 5   # don't report a biome-year median off fewer scored cells than this


def ikey(lat, lon):
    """Integer 0.1-degree cell key, matching build_cells.py's floor() binning."""
    return (math.floor(lat * 10), math.floor(lon * 10))


def ikey_from_cell(lat_cell, lon_cell):
    """Same key from an already-binned lat_cell/lon_cell value."""
    return (round(float(lat_cell) * 10), round(float(lon_cell) * 10))


# ---- cell -> modal biome, from the per-recording ecoregion labels ----
votes = collections.defaultdict(collections.Counter)
with open('recordings_ecoregion.csv', newline='') as f:
    r = csv.reader(f)
    h = next(r); ci = {c: i for i, c in enumerate(h)}
    n_lab = 0
    for row in r:
        b = row[ci['biome_name']]
        if not b or b == 'N/A':
            continue
        try:
            k = ikey(float(row[ci['lat']]), float(row[ci['lon']]))
        except ValueError:
            continue
        votes[k][b] += 1
        n_lab += 1
cell_biome = {k: c.most_common(1)[0][0] for k, c in votes.items()}
# how often is a cell mixed? (honesty check on the modal assignment)
mixed = sum(1 for c in votes.values() if len(c) > 1)
pure_frac = sum(c.most_common(1)[0][1] for c in votes.values()) / max(n_lab, 1)
print(f'labelled recordings: {n_lab:,}')
print(f'cells with a biome: {len(cell_biome):,}   spanning >1 biome: {mixed:,} '
      f'({100*mixed/max(len(cell_biome),1):.1f}%)')
print(f'recordings sitting in their cell\'s modal biome: {100*pure_frac:.1f}%')

# ---- biome x year aggregation over the cell-year richness table ----
by = collections.defaultdict(lambda: {'rich': [], 'n_rec': 0, 'cells': 0})
years = set()
unmatched = 0
with open('grid_cells_yearly.csv', newline='') as f:
    r = csv.reader(f)
    h = next(r); ci = {c: i for i, c in enumerate(h)}
    for row in r:
        k = ikey_from_cell(row[ci['lat_cell']], row[ci['lon_cell']])
        b = cell_biome.get(k)
        if b is None:
            unmatched += 1
            continue
        y = row[ci['year']]
        years.add(y)
        st = by[(b, y)]
        st['n_rec'] += int(row[ci['n_rec']])
        st['cells'] += 1
        s = row[ci['S_rare10']]
        if s != '':
            st['rich'].append(float(s))
print(f'cell-years with no biome match: {unmatched:,}')

YEARS = sorted(years)


def med(v):
    v = sorted(v)
    if not v:
        return None
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


BIOMES = sorted({b for (b, _) in by})

with open('biome_yearly.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['biome_name', 'year', 'n_cell_years', 'n_recordings',
                'n_scored_cells', 'median_S_rare10', 'mean_S_rare10'])
    for b in BIOMES:
        for y in YEARS:
            st = by.get((b, y))
            if not st:
                continue
            m = med(st['rich'])
            w.writerow([b, y, st['cells'], st['n_rec'], len(st['rich']),
                        '' if m is None else round(m, 3),
                        '' if not st['rich'] else round(sum(st['rich']) / len(st['rich']), 3)])
print('wrote biome_yearly.csv')

# ---- per-biome first-vs-last change, plus the tracked-cell up/down split ----
dir_by_biome = collections.defaultdict(collections.Counter)
with open('cell_change.csv', newline='') as f:
    r = csv.reader(f)
    h = next(r); ci = {c: i for i, c in enumerate(h)}
    for row in r:
        # cell_change stores CELL CENTRES (lat_cell + 0.05); undo that before keying
        k = ikey(float(row[ci['lat']]) - 0.05, float(row[ci['lon']]) - 0.05)
        b = cell_biome.get(k)
        if b:
            dir_by_biome[b][row[ci['direction']]] += 1

rows = []
for b in BIOMES:
    ser = [(y, med(by[(b, y)]['rich'])) for y in YEARS
           if (b, y) in by and len(by[(b, y)]['rich']) >= MIN_CELLS]
    if len(ser) < 2:
        continue
    (fy, fv), (ly, lv) = ser[0], ser[-1]
    d = dir_by_biome[b]
    tot_rec = sum(by[(b, y)]['n_rec'] for y in YEARS if (b, y) in by)
    rows.append({
        'biome_name': b, 'first_year': fy, 'last_year': ly,
        'first_median': round(fv, 3), 'last_median': round(lv, 3),
        'delta': round(lv - fv, 3), 'n_years_reported': len(ser),
        'total_recordings': tot_rec,
        'cells_up': d['up'], 'cells_down': d['down'], 'cells_flat': d['flat'],
    })
rows.sort(key=lambda x: -x['total_recordings'])
with open('biome_change.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print('wrote biome_change.csv')

# ---- console summary ----
print(f'\nMedian effort-controlled richness by biome x year '
      f'(blank = fewer than {MIN_CELLS} scored cells)')
print(f"{'biome':<52}" + ''.join(f'{y[2:]:>7}' for y in YEARS) + f"{'delta':>8}")
for r_ in rows:
    b = r_['biome_name']
    line = f"{b[:50]:<52}"
    for y in YEARS:
        st = by.get((b, y))
        if st and len(st['rich']) >= MIN_CELLS:
            line += f"{med(st['rich']):>7.1f}"
        else:
            line += f"{'-':>7}"
    line += f"{r_['delta']:>+8.1f}"
    print(line)

print(f"\nTracked-cell direction split within each biome:")
print(f"{'biome':<52}{'up':>7}{'down':>7}{'flat':>7}")
for r_ in rows:
    print(f"{r_['biome_name'][:50]:<52}{r_['cells_up']:>7}{r_['cells_down']:>7}{r_['cells_flat']:>7}")
