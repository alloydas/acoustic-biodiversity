"""Reconstruct the (previously uncommitted) labeling step of the pipeline.

Consumes grid_cells.csv (from build_cells.py) and derives the good/moderate/bad
biodiversity labels that make_report.py reads. Well-sampled cells (S_rare10
defined, i.e. n_rec >= 10) are split into equal-count tertiles of effort-
controlled richness (good = top third, bad = bottom third); everything else is
'insufficient_data'.

Outputs (schemas reverse-engineered from the Jun-25 artifacts):
  grid_cells_labeled.csv           - grid_cells.csv + biodiversity_label
                                     (GLOBAL tertiles over all usable cells)
  grid_cells_labeled_regional.csv  - grid_cells.csv + lat_band + biodiversity_label_regional
                                     (tertiles WITHIN each 10-deg latitude band; the report's
                                      PRIMARY metric)
  biodiversity_map.csv             - one row per usable cell: lat,lon (cell centre),
                                     continent, n_rec, richness(=S_rare10),
                                     label(=regional), label_code (bad=0/moderate=1/good=2)
"""
import csv
import math
import collections

csv.field_size_limit(10 ** 7)

SRC = 'grid_cells.csv'
LABELS = ['bad', 'moderate', 'good']        # tier 0,1,2 (ascending richness)
CODE = {'bad': 0, 'moderate': 1, 'good': 2}
INSUFF = 'insufficient_data'


def tertiles(usable_idx, richness):
    """Assign equal-count tertiles by rank of richness (ascending -> bad..good).
    Returns {row_index: label}. Ties broken by original order (stable)."""
    order = sorted(usable_idx, key=lambda i: richness[i])
    n = len(order)
    out = {}
    for rank, i in enumerate(order):
        out[i] = LABELS[min(2, rank * 3 // n)] if n else INSUFF
    return out


def lat_band(lat_cell):
    return int(math.floor(lat_cell / 10.0) * 10)


# ---- load grid cells ----
with open(SRC, newline='') as f:
    r = csv.reader(f)
    hdr = next(r)
    rows = [row for row in r]
ci = {c: i for i, c in enumerate(hdr)}
si = ci['S_rare10']

# usable = cells where rarefied richness is defined (n_rec >= 10)
richness = {}
usable = []
for i, row in enumerate(rows):
    v = row[si]
    if v != '':
        richness[i] = float(v)
        usable.append(i)

# ---- GLOBAL tertiles ----
gl = tertiles(usable, richness)          # {i: label} for usable only
global_label = [gl.get(i, INSUFF) for i in range(len(rows))]

# ---- REGIONAL tertiles (within 10-deg latitude band) ----
bands = collections.defaultdict(list)
for i in usable:
    bands[lat_band(float(rows[i][ci['lat_cell']]))].append(i)
reg = {}
for band, idxs in bands.items():
    reg.update(tertiles(idxs, richness))
regional_label = [reg.get(i, INSUFF) for i in range(len(rows))]

# ---- write grid_cells_labeled.csv (global) ----
with open('grid_cells_labeled.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(hdr + ['biodiversity_label'])
    for i, row in enumerate(rows):
        w.writerow(row + [global_label[i]])

# ---- write grid_cells_labeled_regional.csv ----
with open('grid_cells_labeled_regional.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(hdr + ['lat_band', 'biodiversity_label_regional'])
    for i, row in enumerate(rows):
        lb = lat_band(float(row[ci['lat_cell']])) if row[ci['lat_cell']] != '' else ''
        w.writerow(row + [lb, regional_label[i]])

# ---- write biodiversity_map.csv (usable cells; centre coords; regional label) ----
with open('biodiversity_map.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['lat', 'lon', 'continent', 'n_rec', 'richness', 'label', 'label_code'])
    for i in usable:
        row = rows[i]
        lat = round(float(row[ci['lat_cell']]) + 0.05, 4)
        lon = round(float(row[ci['lon_cell']]) + 0.05, 4)
        lab = regional_label[i]
        w.writerow([lat, lon, row[ci['continent']], row[ci['n_rec']],
                    row[si], lab, CODE[lab]])

# ---- report ----
gc = collections.Counter(global_label)
rc = collections.Counter(regional_label)
print(f'cells: {len(rows)}   usable (n_rec>=10): {len(usable)}')
print(f'global  labels: good={gc["good"]} moderate={gc["moderate"]} bad={gc["bad"]} insufficient={gc[INSUFF]}')
print(f'regional labels: good={rc["good"]} moderate={rc["moderate"]} bad={rc["bad"]} insufficient={rc[INSUFF]}')
print(f'latitude bands: {sorted(bands)}')
print('wrote grid_cells_labeled.csv, grid_cells_labeled_regional.csv, biodiversity_map.csv')
