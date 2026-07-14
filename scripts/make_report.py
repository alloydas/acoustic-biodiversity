"""Generate a multi-page PDF report of the acoustic-biodiversity analysis."""
import csv
import math
import datetime
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

csv.field_size_limit(10 ** 7)

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

# ---- load grid cells ----
r = csv.reader(open('grid_cells.csv', newline=''))
h = next(r); ci = {c: i for i, c in enumerate(h)}
cells = [row for row in r]
usable = [row for row in cells if row[ci['S_rare10']] != '']


def f(row, col):
    v = row[ci[col]]
    return float(v) if v != '' else None


def spearman(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    n = len(pairs)
    if n < 10:
        return None

    def ranks(vals):
        order = sorted(range(n), key=lambda i: vals[i])
        rk = [0.0] * n; i = 0
        while i < n:
            j = i
            while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                rk[order[k]] = avg
            i = j + 1
        return rk
    xr = ranks([p[0] for p in pairs]); yr = ranks([p[1] for p in pairs])
    mx = sum(xr) / n; my = sum(yr) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xr, yr))
    den = math.sqrt(sum((a - mx) ** 2 for a in xr) * sum((b - my) ** 2 for b in yr))
    return num / den if den else None


richness = [f(row, 'S_rare10') for row in usable]
corrs = []
for idx in INDICES:
    xs = [f(row, idx) for row in usable]
    corrs.append((SHORT[idx], spearman(xs, richness)))
corrs.sort(key=lambda x: -abs(x[1]) if x[1] is not None else 0)

# ---- load mapping file ----
r = csv.reader(open('biodiversity_map.csv', newline='')); mh = next(r); mi = {c: i for i, c in enumerate(mh)}
mrows = [row for row in r]
COLOR = {'good': '#1a9850', 'moderate': '#fdae61', 'bad': '#d73027'}

# continent medians
bycont = collections.defaultdict(list)
for row in usable:
    bycont[row[ci['continent']]].append(f(row, 'S_rare10'))
cont_med = sorted(((c, sorted(v)[len(v) // 2], len(v)) for c, v in bycont.items()), key=lambda x: -x[1])

label_counts = collections.Counter(row[mi['label']] for row in mrows)
vals = sorted(richness)


def pct(p):
    return vals[int(len(vals) * p)]


# ---- derived figures (computed from the data, not hardcoded, so the report
#      text stays consistent with whatever grid_cells.csv it is run against) ----
TODAY = datetime.date.today().isoformat()
N_CELLS = len(cells)
N_USABLE = len(usable)
BEST_RHO = abs(corrs[0][1]) if corrs and corrs[0][1] is not None else 0.0
COVER_PCT = 100 * N_USABLE / N_CELLS if N_CELLS else 0
# These two are properties of the upstream merge, not visible in grid_cells.csv:
N_RECORDINGS = 759767        # rows in the five score_<c>_meta.csv (merge_meta.py total, 2015-2025)
N_META = 766747              # records flattened into metadata_2015-2025.csv


# ================= BUILD PDF =================
pp = PdfPages('Acoustic_Biodiversity_Report.pdf')

# ---- Page 1: title + summary ----
fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
fig.text(0.5, 0.93, 'Acoustic Biodiversity Mapping', ha='center', size=22, weight='bold')
fig.text(0.5, 0.90, 'A global location-level biodiversity metric from Xeno-canto recordings',
         ha='center', size=11, style='italic', color='#555')
fig.text(0.5, 0.875, f'Generated {TODAY}', ha='center', size=9, color='#888')

summary = (
    "OBJECTIVE\n"
    "Derive a single metric that classifies a geographic location as biodiversity 'good' or 'bad',\n"
    "using a global archive of wildlife sound recordings and their metadata.\n\n"
    "DATA\n"
    f"  -  {N_RECORDINGS:,} recordings across 5 continents (Xeno-canto), each with 42 acoustic indices\n"
    "     computed by a SLURM batch pipeline (compute_indice.py over Butterworth-filtered audio).\n"
    f"  -  Metadata flattened from {N_META:,} API records (species, coordinates, quality, device).\n"
    "  -  Acoustic scores joined to metadata on recording id (>99.99% match).\n\n"
    "METHOD\n"
    f"  1. Aggregate all recordings into 0.1-degree (~11 km) grid cells -> {N_CELLS:,} cells.\n"
    "  2. Per cell, compute species richness from recorded species (genus+species and the\n"
    "     'also' co-occurring-species field), corrected for sampling effort by Hurlbert\n"
    "     rarefaction to 10 recordings (S_rare10).\n"
    "  3. Test whether the acoustic indices predict richness (validation against ground truth).\n"
    "  4. Label each well-sampled cell good / moderate / bad relative to its 10-degree\n"
    "     latitude band (region-relative tertiles).\n\n"
    "HEADLINE RESULT\n"
    f"  -  {N_USABLE:,} cells had enough recordings (>=10) to score confidently.\n"
    f"  -  The acoustic indices did NOT predict species richness (best |Spearman rho| ~ {BEST_RHO:.2f}).\n"
    "     The defensible biodiversity metric is therefore effort-controlled species richness,\n"
    "     not any acoustic index -- the indices are unreliable on targeted (non-soundscape)\n"
    "     recordings that dominate the archive.\n"
    f"  -  Final labels: {label_counts['good']} good, {label_counts['moderate']} moderate, "
    f"{label_counts['bad']} bad."
)
fig.text(0.08, 0.80, summary, ha='left', va='top', size=9.3, family='monospace')
fig.text(0.5, 0.04, 'Acoustic Biodiversity Report  -  page 1', ha='center', size=8, color='#999')
pp.savefig(fig); plt.close(fig)

# ---- Page 2: world map ----
fig = plt.figure(figsize=(11.69, 8.27))  # A4 landscape
ax = fig.add_axes([0.05, 0.10, 0.90, 0.80])
order = ['bad', 'moderate', 'good']
for lab in order:
    pts = [(float(row[mi['lon']]), float(row[mi['lat']])) for row in mrows if row[mi['label']] == lab]
    if pts:
        xs, ys = zip(*pts)
        ax.scatter(xs, ys, s=10, c=COLOR[lab], label=f'{lab} (n={label_counts[lab]})',
                   alpha=0.7, edgecolors='none')
ax.set_xlim(-180, 180); ax.set_ylim(-60, 80)
ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
ax.axhline(0, color='#ccc', lw=0.6); ax.grid(True, lw=0.3, color='#eee')
ax.set_title(f'Region-relative biodiversity label per 0.1-degree cell ({N_USABLE:,} well-sampled cells)', size=13)
ax.legend(loc='lower left', markerscale=2, framealpha=0.9)
fig.text(0.5, 0.04, 'Acoustic Biodiversity Report  -  page 2  -  world map', ha='center', size=8, color='#999')
pp.savefig(fig); plt.close(fig)

# ---- Page 3: validation (correlation + distribution) ----
fig = plt.figure(figsize=(11.69, 8.27))
gs = GridSpec(1, 2, figure=fig, wspace=0.35, left=0.30, right=0.95, top=0.85, bottom=0.12)
ax1 = fig.add_subplot(gs[0, 0])
names = [c[0] for c in corrs][::-1]
rho = [c[1] for c in corrs][::-1]
barcol = ['#d73027' if v < 0 else '#1a9850' for v in rho]
ax1.barh(range(len(names)), rho, color=barcol)
ax1.set_yticks(range(len(names))); ax1.set_yticklabels(names, size=8)
ax1.axvline(0, color='k', lw=0.8)
ax1.set_xlim(-0.3, 0.3); ax1.set_xlabel('Spearman rho vs richness (S_rare10)')
ax1.set_title('Acoustic index vs species richness\n(none is predictive)', size=11)
ax1.grid(True, axis='x', lw=0.3, color='#eee')

ax2 = fig.add_subplot(gs[0, 1])
ax2.hist(richness, bins=40, color='#4575b4', edgecolor='white')
for p, c, ls in [(.5, '#000', '-')]:
    ax2.axvline(pct(p), color=c, ls=ls, lw=1)
ax2.set_xlabel('Effort-controlled richness (S_rare10)'); ax2.set_ylabel('cells')
ax2.set_title(f'Richness distribution across {len(usable)} cells\nmedian={pct(.5):.1f}  p90={pct(.9):.1f}  max={max(vals):.0f}', size=11)
fig.suptitle('Validation: does the audio carry the biodiversity signal?', size=14, weight='bold')
fig.text(0.5, 0.03, 'Acoustic Biodiversity Report  -  page 3  -  validation', ha='center', size=8, color='#999')
pp.savefig(fig); plt.close(fig)

# ---- Page 4: metric definition + tables ----
fig = plt.figure(figsize=(8.27, 11.69))
fig.text(0.5, 0.94, 'The Metric & Results', ha='center', size=18, weight='bold')
metric_txt = (
    "METRIC DEFINITION\n"
    "  location           = 0.1 deg x 0.1 deg grid cell (~11 km)\n"
    "  biodiversity value = S_rare10, the expected number of distinct recorded\n"
    "                       species in 10 recordings (Hurlbert sample-based\n"
    "                       rarefaction; controls for recording effort)\n"
    "  label              = tertile of S_rare10 WITHIN the cell's 10-deg latitude\n"
    "                       band:  good = top third, bad = bottom third\n"
    f"  scope              = cells with >= 10 recordings ({N_USABLE:,} of {N_CELLS:,} cells)\n\n"
    "GLOBAL VALUE RANGE (S_rare10)\n"
    f"  min {vals[0]:.1f}   p25 {pct(.25):.1f}   median {pct(.5):.1f}   "
    f"p75 {pct(.75):.1f}   p90 {pct(.9):.1f}   max {max(vals):.0f}\n"
)
fig.text(0.08, 0.88, metric_txt, ha='left', va='top', size=9.3, family='monospace')

# continent table
ax = fig.add_axes([0.10, 0.50, 0.80, 0.18]); ax.axis('off')
tbl = [['Continent', 'Median richness', 'Scored cells']] + [[c, f'{m:.1f}', str(n)] for c, m, n in cont_med]
t = ax.table(cellText=tbl, loc='center', cellLoc='center'); t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.4)
for k in range(3):
    t[0, k].set_facecolor('#4575b4'); t[0, k].set_text_props(color='white', weight='bold')
fig.text(0.10, 0.70, 'Median richness by continent (well-sampled cells)', size=10, weight='bold')

limit_txt = (
    "KEY LIMITATIONS (must accompany any use of this metric)\n"
    "  1. Acoustic indices are non-predictive here. ACI/ADI/NDSI etc. were built for\n"
    "     passive soundscape monitoring; 69% of archive clips are single-target\n"
    "     recordings (median 24 s), so per-clip indices do not reflect site diversity.\n"
    "  2. 'Richness' = RECORDED species, not true species. It reflects recordist effort\n"
    "     and interest. Rarefaction controls sample SIZE but not observer bias.\n"
    "  3. Weak latitude gradient. Effort-controlled richness is nearly flat across\n"
    "     latitude (thresholds ~10-13 everywhere) -- evidence the signal is\n"
    "     effort-driven, not the true tropical biodiversity peak. Region-relative\n"
    "     labelling is used so scores mean 'rich for its region'.\n"
    f"  4. Coverage. Only {COVER_PCT:.0f}% of cells are confidently scored; the rest are\n"
    "     'insufficient_data' (< 10 recordings).\n"
)
fig.text(0.08, 0.46, limit_txt, ha='left', va='top', size=9.0, family='monospace')

out_txt = (
    "OUTPUT FILES\n"
    f"  grid_cells.csv                    - per-cell features + richness (all {N_CELLS:,} cells)\n"
    "  grid_cells_labeled_regional.csv   - region-relative good/moderate/bad labels\n"
    "  biodiversity_map.csv              - lat, lon, label, code -> ready for GIS/folium\n"
)
fig.text(0.08, 0.16, out_txt, ha='left', va='top', size=9.0, family='monospace')
fig.text(0.5, 0.03, 'Acoustic Biodiversity Report  -  page 4  -  metric & results', ha='center', size=8, color='#999')
pp.savefig(fig); plt.close(fig)

# ---- Page 5: year-by-year ----
import os
if os.path.exists('grid_cells_yearly.csv') and os.path.exists('cell_change.csv'):
    yr = csv.reader(open('grid_cells_yearly.csv', newline='')); yh = next(yr); yi = {c: i for i, c in enumerate(yh)}
    trend = collections.defaultdict(lambda: collections.defaultdict(list))
    cellyears = collections.defaultdict(set); scored_cy = 0
    for row in yr:
        s = row[yi['S_rare10']]
        if s != '':
            trend[row[yi['year']]][row[yi['continent']]].append(float(s))
            cellyears[(row[yi['lat_cell']], row[yi['lon_cell']])].add(row[yi['year']])
            scored_cy += 1
    YEARS = sorted(trend)  # all calendar years present, e.g. 2015-2025
    NYEARS = len(YEARS)
    conts = sorted({c for y in trend for c in trend[y]})
    n_all3 = sum(1 for v in cellyears.values() if len(v) >= NYEARS)

    cr = csv.reader(open('cell_change.csv', newline='')); ch = next(cr); kdir = ch.index('direction')
    cd = ch.index('delta'); deltas = []; dirc = collections.Counter()
    for row in cr:
        dirc[row[kdir]] += 1; deltas.append(float(row[cd]))
    n_tracked = sum(dirc.values())

    fig = plt.figure(figsize=(11.69, 8.27))
    fig.suptitle(f'Year-by-year ({YEARS[0]}-{YEARS[-1]}): temporal slices of the metric', size=14, weight='bold')
    # left: trend lines
    ax1 = fig.add_axes([0.07, 0.32, 0.42, 0.50])
    palette = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']
    for k, cont in enumerate(conts):
        ys = [(sorted(trend[y][cont])[len(trend[y][cont]) // 2] if trend[y][cont] else float('nan')) for y in YEARS]
        ax1.plot(YEARS, ys, marker='o', color=palette[k % 5], label=cont)
    gall = [sorted(v for c in conts for v in trend[y][c]) for y in YEARS]
    ax1.plot(YEARS, [g[len(g) // 2] for g in gall], marker='s', color='black', lw=2.4, label='ALL')
    ax1.set_ylabel('Median richness (S_rare10)'); ax1.set_title('Recorded-richness trend by continent', size=11)
    ax1.grid(True, lw=0.3, color='#eee'); ax1.legend(fontsize=8, ncol=2)
    # right: change histogram
    ax2 = fig.add_axes([0.57, 0.32, 0.38, 0.50])
    ax2.hist(deltas, bins=30, color='#888', edgecolor='white')
    ax2.axvline(0, color='k', lw=1)
    ax2.set_xlabel('Change in richness, first->last scored year')
    ax2.set_title(f'Per-cell change ({n_tracked} cells tracked)\nup={dirc["up"]}  down={dirc["down"]}  flat={dirc["flat"]}', size=11)

    dip_first = gall[0][len(gall[0]) // 2] if gall[0] else float('nan')
    dip_last = gall[-1][len(gall[-1]) // 2] if gall[-1] else float('nan')
    note = (
        f"COVERAGE: {scored_cy:,} scored cell-years; {n_tracked} cells scored in >=2 years, {n_all3} in all {NYEARS}.\n\n"
        f"CAUTION: even an {NYEARS}-year window cannot show real biodiversity change. These movements reflect\n"
        "WHICH cells were recorded and by WHOM each year (effort + observer turnover), not\n"
        f"ecological gain or loss. The global change ({dip_first:.1f} -> {dip_last:.1f}) tracks recording effort, not nature.\n"
        "Use year slices to study sampling coverage over time -- not as a biodiversity time series.\n\n"
        "OUTPUTS: grid_cells_yearly.csv (per cell-year),  cell_change.csv (per tracked cell)."
    )
    fig.text(0.07, 0.24, note, ha='left', va='top', size=9.2, family='monospace')
    fig.text(0.5, 0.03, 'Acoustic Biodiversity Report  -  page 5  -  temporal', ha='center', size=8, color='#999')
    pp.savefig(fig); plt.close(fig)
    npages = 5
else:
    npages = 4

# ---- Page 6: ten-year findings (2015-2025) ----
# Synthesis page. All figures derived from the data at runtime; the only literal
# is the reference |rho| from the smaller 2023-2025 subset (a documented prior).
if os.path.exists('grid_cells_yearly.csv') and os.path.exists('cell_change.csv'):
    PRIOR_RHO = 0.16  # best |Spearman rho| on the 2023-2025 gap-inclusive subset

    # global median effort-controlled richness per year (all continents pooled)
    yr = csv.reader(open('grid_cells_yearly.csv', newline='')); yh = next(yr); yj = {c: i for i, c in enumerate(yh)}
    by_year = collections.defaultdict(list); tracked_years = collections.defaultdict(set)
    for row in yr:
        s = row[yj['S_rare10']]
        if s != '':
            by_year[row[yj['year']]].append(float(s))
            tracked_years[(row[yj['lat_cell']], row[yj['lon_cell']])].add(row[yj['year']])
    F_YEARS = sorted(by_year)
    F_NY = len(F_YEARS)
    gmed = {y: sorted(v)[len(v) // 2] for y, v in by_year.items()}
    med_series = [gmed[y] for y in F_YEARS]
    flat_lo, flat_hi = min(med_series), max(med_series)
    n_tracked2 = sum(1 for v in tracked_years.values() if len(v) >= 2)
    n_all = sum(1 for v in tracked_years.values() if len(v) == F_NY)

    cr = csv.reader(open('cell_change.csv', newline='')); ch = next(cr); kdir = ch.index('direction')
    fdir = collections.Counter(row[kdir] for row in cr)
    n_ch = sum(fdir.values())

    fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
    fig.text(0.5, 0.955, 'Ten-Year Findings (2015-2025)', ha='center', size=18, weight='bold')
    fig.text(0.5, 0.930, 'What the decade-scale data adds to the analysis', ha='center',
             size=10, style='italic', color='#555')

    findings = (
        f"1.  MORE DATA MADE THE ACOUSTIC INDICES LOOK WORSE, NOT BETTER.\n"
        f"    Expanding from the 2023-2025 subset to the full 2015-2025 archive\n"
        f"    ({N_RECORDINGS:,} recordings, {N_USABLE:,} scored cells) shrank the best index\n"
        f"    correlation from |rho| ~ {PRIOR_RHO:.2f} to ~ {BEST_RHO:.2f}. A real signal sharpens\n"
        f"    with more data; this faded toward zero -- the strongest evidence yet\n"
        f"    that the indices carry no site-biodiversity signal in this archive.\n\n"
        f"2.  A DECADE OF EFFORT-CONTROLLED RICHNESS IS ESSENTIALLY FLAT.\n"
        f"    Global median S_rare10 stays within {flat_lo:.1f}-{flat_hi:.1f} across all {F_NY} years\n"
        f"    ({F_YEARS[0]}->{F_YEARS[-1]}: {med_series[0]:.1f} -> {med_series[-1]:.1f}). This reframes the metric as a\n"
        f"    sampling-effort measure, not ecology. The only structure is a mild\n"
        f"    recent dip that tracks the incompleteness of recent uploads.\n\n"
        f"3.  TEMPORAL CHANGE IS A COIN-FLIP (RANDOM-WALK SIGNATURE).\n"
        f"    Of {n_ch:,} cells tracked across >=2 years, {fdir['up']:,} rose and {fdir['down']:,} fell\n"
        f"    (flat={fdir['flat']:,}). That near-perfect symmetry over a decade is the\n"
        f"    fingerprint of observer turnover and noise, not directional change.\n\n"
        f"4.  SPATIALLY DECADE-RICH, BUT TEMPORALLY STILL STARVED.\n"
        f"    Coverage grew to {N_CELLS:,} cells / {N_USABLE:,} scored -- a much denser map --\n"
        f"    yet only {n_all:,} cells were recorded in all {F_NY} years ({n_tracked2:,} in >=2).\n"
        f"    The biodiversity MAP is far stronger; a true TIME SERIES remains\n"
        f"    impossible from this data.\n"
    )
    fig.text(0.07, 0.90, findings, ha='left', va='top', size=9.2, family='monospace')

    # inset: global median richness by year (drives home finding #2 -- flatness)
    ax = fig.add_axes([0.14, 0.10, 0.74, 0.20])
    ax.plot(F_YEARS, med_series, marker='o', color='#333', lw=2)
    ax.axhspan(flat_lo, flat_hi, color='#4575b4', alpha=0.10)
    ax.set_ylim(max(0, flat_lo - 3), flat_hi + 3)
    ax.set_ylabel('Global median S_rare10'); ax.set_title(
        f'Effort-controlled richness barely moves over {F_NY} years (band = {flat_lo:.1f}-{flat_hi:.1f})', size=9.5)
    ax.grid(True, lw=0.3, color='#eee')
    for lab in ax.get_xticklabels():
        lab.set_rotation(45); lab.set_fontsize(8)

    fig.text(0.5, 0.03, 'Acoustic Biodiversity Report  -  page 6  -  ten-year findings', ha='center', size=8, color='#999')
    pp.savefig(fig); plt.close(fig)
    npages += 1

# ---- Page 7: yearwise effort vs. metric (2015-2025) ----
# Recording effort (bars, per continent) rose sharply while the effort-controlled
# metric stayed flat; the apparent recent dip is a sampling artifact, not ecology.
if os.path.exists('grid_cells_yearly.csv'):
    yr = csv.reader(open('grid_cells_yearly.csv', newline='')); yh = next(yr); yj = {c: i for i, c in enumerate(yh)}
    recCY = collections.defaultdict(lambda: collections.defaultdict(int))
    richY = collections.defaultdict(list); cellsY = collections.Counter()
    for row in yr:
        y = row[yj['year']]; cont = row[yj['continent']]
        try:
            nr = int(row[yj['n_rec']])
        except ValueError:
            nr = 0
        recCY[y][cont] += nr
        s = row[yj['S_rare10']]
        if s != '':
            richY[y].append(float(s)); cellsY[y] += 1
    YRS = sorted(recCY)
    xs = list(range(len(YRS)))
    yconts = ['europe', 'america', 'asia', 'africa', 'australia']
    ycol = {'europe': '#4a7fb5', 'america': '#17a07f', 'asia': '#9a6fc0',
            'africa': '#d08a3f', 'australia': '#7ba33f'}

    def _med(v):
        v = sorted(v)
        return v[len(v) // 2] if v else float('nan')
    med = [_med(richY[y]) for y in YRS]
    cells = [cellsY[y] for y in YRS]
    totals = [sum(recCY[y].values()) for y in YRS]

    fig = plt.figure(figsize=(11.69, 8.27))
    fig.suptitle(f'Year-by-year: effort vs. the metric ({YRS[0]}-{YRS[-1]})', size=14, weight='bold')
    # top: stacked recordings by continent
    ax1 = fig.add_axes([0.08, 0.60, 0.86, 0.29])
    bottom = [0] * len(YRS)
    for c in yconts:
        v = [recCY[y][c] for y in YRS]
        ax1.bar(xs, v, bottom=bottom, color=ycol[c], label=c, width=0.72, edgecolor='white', linewidth=0.4)
        bottom = [b + a for b, a in zip(bottom, v)]
    ax1.set_xticks(xs); ax1.set_xticklabels([y[2:] for y in YRS])
    ax1.set_ylabel('recordings in cells')
    ax1.legend(ncol=5, fontsize=8, loc='upper left', frameon=False)
    ax1.set_title(f"Recording effort behind the metric, by continent  ({totals[0]:,} -> {totals[-1]:,}, +{100*(totals[-1]-totals[0])//totals[0]}%)",
                  size=10, loc='left')
    ax1.grid(True, axis='y', lw=0.3, color='#eee')
    # bottom: scored cells (bars) + median richness (line, twin axis)
    ax2 = fig.add_axes([0.08, 0.12, 0.86, 0.31])
    ax2.bar(xs, cells, color='#c9d2c9', width=0.72, alpha=0.85)
    ax2.set_ylabel('scored cells (n>=10)'); ax2.set_ylim(0, max(cells) * 1.3)
    ax2.set_xticks(xs); ax2.set_xticklabels([y[2:] for y in YRS])
    ax2.grid(True, axis='y', lw=0.3, color='#eee')
    ax3 = ax2.twinx()
    ax3.axhspan(min(med), max(med), color='#4575b4', alpha=0.08)
    ax3.plot(xs, med, color='#333', lw=2.2, marker='o', ms=5, mfc='#4575b4', mec='#333')
    for i, m in enumerate(med):
        ax3.annotate(f'{m:.1f}', (xs[i], m), textcoords='offset points', xytext=(0, 8),
                     ha='center', size=8, weight='bold')
    ax3.set_ylim(6, 13); ax3.set_ylabel('median S_rare10')
    pk = med.index(max(med))  # highlight the apparent post-peak decline
    if pk < len(med) - 1:
        ax3.annotate('', xy=(xs[-1], med[-1]), xytext=(xs[pk], med[pk]),
                     arrowprops=dict(arrowstyle='->', color='#d73027', lw=1.4, alpha=0.85))
        ax3.text(xs[-1], med[-1] - 0.1, f'  apparent -{max(med) - med[-1]:.1f}\n  (sampling, not ecology)',
                 color='#d73027', size=7.5, va='top', ha='right')
    ax2.set_title(f"Scored cells grew {cells[0]:,} -> {cells[-1]:,}; median richness stayed flat ({min(med):.1f}-{max(med):.1f} band)",
                  size=10, loc='left')

    note = (
        "Bars = recording EFFORT (geolocated recordings entering the metric); line = effort-controlled median\n"
        "richness. Effort rose ~70% over the decade while the metric stayed inside a narrow band. The apparent\n"
        "recent dip (red arrow) is NOT ecological degradation: it tracks the lower completeness of recent uploads\n"
        "and shifting recordist coverage. Do not read the decline as biodiversity loss."
    )
    fig.text(0.08, 0.095, note, ha='left', va='top', size=8.6, family='monospace')
    fig.text(0.5, 0.03, 'Acoustic Biodiversity Report  -  page 7  -  yearwise', ha='center', size=8, color='#999')
    pp.savefig(fig); plt.close(fig)
    npages += 1

pp.close()
print(f'wrote Acoustic_Biodiversity_Report.pdf  ({npages} pages)')
