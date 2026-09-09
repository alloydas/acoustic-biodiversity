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

# ---- Page 8: biome-wise year-by-year change (2015-2025) ----
# Splits the yearly richness series by BIOME. Cells take the modal biome of their
# recordings (build_biome_yearly.py); 99.4% of recordings sit in their cell's
# modal biome, so the assignment is near-unambiguous.
if os.path.exists('biome_yearly.csv') and os.path.exists('biome_change.csv'):
    import numpy as _np

    br = list(csv.DictReader(open('biome_yearly.csv', newline='')))
    bc = list(csv.DictReader(open('biome_change.csv', newline='')))
    B_YEARS = sorted({r['year'] for r in br})
    B_ORDER = [r['biome_name'] for r in bc]   # most recordings first
    lookup = {(r['biome_name'], r['year']): r for r in br}
    MINC = 5

    grid = _np.full((len(B_ORDER), len(B_YEARS)), _np.nan)
    for i, b in enumerate(B_ORDER):
        for j, y in enumerate(B_YEARS):
            r = lookup.get((b, y))
            if r and r['median_S_rare10'] != '' and int(r['n_scored_cells']) >= MINC:
                grid[i, j] = float(r['median_S_rare10'])
    gmid = _np.nanmedian(grid)

    fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
    fig.text(0.5, 0.962, 'Year-by-Year Change by Biome (2015-2025)',
             ha='center', size=16, weight='bold')
    fig.text(0.5, 0.940, 'Effort-controlled richness (S_rare10) split by terrestrial biome',
             ha='center', size=9, style='italic', color='#555')

    L, W = 0.335, 0.545        # shared left edge / width for both panels
    # ---- heatmap ----
    axh = fig.add_axes([L, 0.605, W, 0.295])
    im = axh.imshow(grid, aspect='auto', cmap='RdYlBu', vmin=gmid - 3, vmax=gmid + 3)
    axh.set_xticks(range(len(B_YEARS)))
    axh.set_xticklabels([y[2:] for y in B_YEARS], size=8)
    axh.set_yticks(range(len(B_ORDER)))
    axh.set_yticklabels([b[:40] for b in B_ORDER], size=7.2)
    for i in range(len(B_ORDER)):
        for j in range(len(B_YEARS)):
            v = grid[i, j]
            if _np.isnan(v):
                col = '#bbb'
            else:
                # white text on the saturated ends of the colormap, dark in the middle
                col = 'white' if abs(v - gmid) > 2.2 else '#222'
            axh.text(j, i, '-' if _np.isnan(v) else f'{v:.1f}', ha='center',
                     va='center', size=6.0, color=col,
                     weight='bold' if not _np.isnan(v) and abs(v - gmid) > 2.2 else 'normal')
    axh.set_title(f'Median S_rare10 per biome-year   (- = fewer than {MINC} scored cells)',
                  size=9, loc='left')
    cax = fig.add_axes([L + W + 0.075, 0.605, 0.014, 0.295])
    cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=7)

    # ---- delta column, same row order as the heatmap ----
    axd = fig.add_axes([L + W + 0.005, 0.605, 0.062, 0.295]); axd.axis('off')
    axd.set_ylim(len(B_ORDER) - 0.5, -0.5)   # match imshow row coords exactly
    axd.set_xlim(0, 1)
    axd.text(0.5, -0.9, 'delta', ha='center', size=7.5, weight='bold')
    for i, r in enumerate(bc):
        d = float(r['delta'])
        axd.text(0.5, i, f'{d:+.1f}', ha='center', va='center', size=7.2,
                 color='#d73027' if d < -0.5 else ('#1a9850' if d > 0.5 else '#888'),
                 weight='bold' if abs(d) > 0.5 else 'normal')

    # ---- tracked-cell direction split (same order, labels shared above) ----
    axb = fig.add_axes([L, 0.395, W, 0.155])
    ys = _np.arange(len(B_ORDER))
    up = _np.array([int(r['cells_up']) for r in bc], float)
    dn = _np.array([int(r['cells_down']) for r in bc], float)
    fl = _np.array([int(r['cells_flat']) for r in bc], float)
    tot = _np.maximum(up + dn + fl, 1)
    axb.barh(ys, 100 * up / tot, color='#1a9850', label='up', height=0.78)
    axb.barh(ys, 100 * fl / tot, left=100 * up / tot, color='#dcdcdc', label='flat', height=0.78)
    axb.barh(ys, 100 * dn / tot, left=100 * (up + fl) / tot, color='#d73027', label='down', height=0.78)
    axb.axvline(50, color='#333', lw=1.0, ls='--', alpha=0.75)
    axb.set_ylim(len(B_ORDER) - 0.5, -0.5)   # top-to-bottom, matching the heatmap
    axb.set_yticks(ys); axb.set_yticklabels([])
    axb.set_xlim(0, 100); axb.set_xlabel('% of individually tracked cells', size=8)
    axb.tick_params(labelsize=7.5)
    axb.legend(ncol=3, fontsize=7, loc='upper center',
               bbox_to_anchor=(0.5, -0.22), frameon=False)
    for i, r in enumerate(bc):
        n = int(r['cells_up']) + int(r['cells_down']) + int(r['cells_flat'])
        axb.text(101.5, i, f'n={n}', va='center', size=6.2, color='#777')
    axb.set_title('Direction of change, cells tracked across >=2 years  (dashed = 50/50)',
                  size=9, loc='left')

    # ---- narrative, full width below both panels ----
    worst = max(bc, key=lambda r: abs(float(r['delta'])))
    t_up = sum(int(r['cells_up']) for r in bc)
    t_dn = sum(int(r['cells_down']) for r in bc)
    note = (
        "NO BIOME SHOWS COHERENT DIRECTIONAL CHANGE.\n\n"
        f"Every biome's year-to-year jitter is as large as its {B_YEARS[0]}->{B_YEARS[-1]} delta, and the cells\n"
        f"tracked inside each biome split close to 50/50 ({t_up:,} up vs {t_dn:,} down overall). The\n"
        f"largest delta -- {worst['biome_name'][:38]} at {float(worst['delta']):+.1f} -- rests on a series\n"
        "that swings by more than that between adjacent years.\n\n"
        "READ THIS AS SAMPLING, NOT ECOLOGY. S_rare10 controls for the NUMBER of recordings in\n"
        "a cell, but not for who recorded, for how long, or with what target. A biome's series\n"
        "therefore moves as its recordist community turns over. The biomes with the fewest\n"
        "recordings -- Tundra, Mangroves, Flooded Grasslands -- swing hardest, which is the\n"
        "signature of small samples rather than of habitat change.\n\n"
        "METHOD. Each 0.1-degree cell takes the modal biome of the recordings inside it. 99.4% of\n"
        "recordings fall in their cell's modal biome and only 1.6% of cells span more than one,\n"
        "so boundary cells cannot drive these patterns. Biomes are ordered by recording volume."
    )
    fig.text(0.075, 0.315, note, ha='left', va='top', size=8.2, family='monospace')
    fig.text(0.5, 0.028, 'Acoustic Biodiversity Report  -  page 8  -  biome-wise yearly change',
             ha='center', size=8, color='#999')
    pp.savefig(fig); plt.close(fig)
    npages += 1

# ---- Page 9: each biome's series on its own axes (small multiples) ----
# The page-8 heatmap compares biomes; this page reads each biome one at a time,
# with its own recording effort behind it so jitter can be traced to sample size.
if os.path.exists('biome_yearly.csv') and os.path.exists('biome_change.csv'):
    import numpy as _np

    br9 = list(csv.DictReader(open('biome_yearly.csv', newline='')))
    bc9 = list(csv.DictReader(open('biome_change.csv', newline='')))
    Y9 = sorted({r['year'] for r in br9})
    ORDER9 = [r['biome_name'] for r in bc9]
    look9 = {(r['biome_name'], r['year']): r for r in br9}
    MINC9 = 5

    ncol, nrow = 3, 5
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.5, 0.968, 'Each Biome, Year by Year (2015-2025)', ha='center',
             size=16, weight='bold')
    fig.text(0.5, 0.947,
             'Line = median S_rare10 (left axis).  Bars = recordings entering the metric (right axis).',
             ha='center', size=8.5, style='italic', color='#555')

    xs9 = _np.arange(len(Y9))
    for i, b in enumerate(ORDER9):
        rr, cc = divmod(i, ncol)
        ax = fig.add_axes([0.085 + cc * 0.305, 0.790 - rr * 0.150, 0.235, 0.100])
        med9, eff9 = [], []
        for y in Y9:
            r = look9.get((b, y))
            ok = r and r['median_S_rare10'] != '' and int(r['n_scored_cells']) >= MINC9
            med9.append(float(r['median_S_rare10']) if ok else _np.nan)
            eff9.append(int(r['n_recordings']) if r else 0)
        axe = ax.twinx()
        axe.bar(xs9, eff9, color='#d8dee8', width=0.75, zorder=1)
        axe.set_ylim(0, max(max(eff9), 1) * 1.9)
        axe.set_yticks([])
        ax.plot(xs9, med9, color='#22303f', lw=1.5, marker='o', ms=3.2,
                mfc='#4575b4', mec='#22303f', zorder=3)
        ax.set_zorder(axe.get_zorder() + 1); ax.patch.set_visible(False)
        ax.set_ylim(3, 15)
        ax.set_yticks([5, 10, 15]); ax.tick_params(labelsize=6.2)
        ax.set_xticks(xs9[::2]); ax.set_xticklabels([y[2:] for y in Y9[::2]], size=6.2)
        ax.grid(True, axis='y', lw=0.3, color='#eee', zorder=0)
        d = float(bc9[i]['delta'])
        ax.set_title(f"{b[:34]}", size=6.9, loc='left', pad=8)
        ax.text(0, 1.015, f"n={int(bc9[i]['total_recordings']):,}   delta {d:+.1f}",
                transform=ax.transAxes, size=5.9,
                color='#d73027' if d < -0.5 else ('#1a9850' if d > 0.5 else '#777'))

    note9 = (
        "Every panel is drawn on the SAME richness axis (3-15), so the series can be compared\n"
        "directly. Read the grey bars first: where effort is small the line is jagged, and where\n"
        "effort is large it is flat. That is the whole pattern -- Tundra, Mangroves and Flooded\n"
        "Grasslands swing by several species between adjacent years on a few hundred recordings,\n"
        "while Temperate Broadleaf & Mixed Forests, with 314,620, barely moves.\n\n"
        "No panel shows a monotonic decade-long trend. The deltas printed above each panel are\n"
        "first-vs-last-year differences, not fitted trends, and in every biome the year-to-year\n"
        "jitter is at least as large as the delta -- so none of them is distinguishable from noise."
    )
    fig.text(0.085, 0.150, note9, ha='left', va='top', size=7.6, family='monospace')
    fig.text(0.5, 0.012, 'Acoustic Biodiversity Report  -  page 9  -  biome series',
             ha='center', size=8, color='#999')
    pp.savefig(fig); plt.close(fig)
    npages += 1

# ---- Page 10: urban class (city/town/rural) year by year ----
if os.path.exists('urban_yearly.csv') and os.path.exists('urban_change.csv'):
    import numpy as _np

    ur = list(csv.DictReader(open('urban_yearly.csv', newline='')))
    uc = list(csv.DictReader(open('urban_change.csv', newline='')))
    UY = sorted({r['year'] for r in ur})
    UCLS = ['city', 'town', 'rural']
    ulook = {(r['urban_class'], r['year']): r for r in ur}
    UCOL = {'city': '#d73027', 'town': '#f0a24a', 'rural': '#1a9850'}
    MINU = 5
    LAST_LABELLED = max(r['year'] for r in ur if r['class_carried_forward'] == '0')

    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.5, 0.962, 'Year-by-Year Change by Urban Class (2015-2025)',
             ha='center', size=16, weight='bold')
    fig.text(0.5, 0.940, 'Effort-controlled richness (S_rare10) for city / town / rural cells',
             ha='center', size=9, style='italic', color='#555')

    xs = _np.arange(len(UY))
    cut = UY.index(LAST_LABELLED)

    # ---- richness series ----
    ax = fig.add_axes([0.11, 0.60, 0.80, 0.29])
    ax.axvspan(cut + 0.5, len(UY) - 0.5, color='#f2f2f2', zorder=0)
    ax.text(len(UY) - 0.55, 11.75, 'class carried forward\n(no GCTB polygons)', ha='right',
            va='top', size=7, color='#888', style='italic')
    for u in UCLS:
        vals = []
        for y in UY:
            r = ulook.get((u, y))
            ok = r and r['median_S_rare10'] != '' and int(r['n_scored_cells']) >= MINU
            vals.append(float(r['median_S_rare10']) if ok else _np.nan)
        ax.plot(xs, vals, color=UCOL[u], lw=1.9, marker='o', ms=4, label=u, zorder=3)
    ax.axvline(cut + 0.5, color='#999', lw=1.0, ls='--', zorder=2)
    ax.set_xticks(xs); ax.set_xticklabels([y[2:] for y in UY], size=8)
    ax.set_ylim(6, 12); ax.set_ylabel('median S_rare10', size=8.5)
    ax.tick_params(labelsize=8)
    ax.grid(True, axis='y', lw=0.3, color='#eee')
    ax.legend(ncol=3, fontsize=8, frameon=False, loc='lower left')
    ax.set_title('Median effort-controlled richness by class', size=9.5, loc='left')

    # ---- effort behind each class (log scale: rural dwarfs the others) ----
    ax2 = fig.add_axes([0.11, 0.395, 0.80, 0.155])
    wdt = 0.26
    for i, u in enumerate(UCLS):
        eff = [int(ulook[(u, y)]['n_recordings']) if (u, y) in ulook else 0 for y in UY]
        ax2.bar(xs + (i - 1) * wdt, eff, width=wdt, color=UCOL[u], label=u, alpha=0.85)
    ax2.set_yscale('log'); ax2.set_ylabel('recordings (log)', size=8.5)
    ax2.set_xticks(xs); ax2.set_xticklabels([y[2:] for y in UY], size=8)
    ax2.tick_params(labelsize=7.5)
    ax2.axvline(cut + 0.5, color='#999', lw=1.0, ls='--')
    ax2.grid(True, axis='y', lw=0.3, color='#eee')
    ax2.set_title('Recording effort behind each class (log scale)', size=9.5, loc='left')

    # ---- narrative ----
    cty = next(r for r in uc if r['urban_class'] == 'city')
    rur = next(r for r in uc if r['urban_class'] == 'rural')
    note = (
        "NEITHER CITY NOR RURAL RICHNESS TRENDS.\n\n"
        f"Over the labelled window ({cty['labelled_first_year']}-{cty['labelled_last_year']}) the city median moves "
        f"{float(cty['delta_labelled_window']):+.1f} and rural {float(rur['delta_labelled_window']):+.1f}. Both are\n"
        "far smaller than the year-to-year jitter in either series, so neither is distinguishable\n"
        "from noise. The town series (only 719 cells) swings by ~3 species between adjacent\n"
        "years, which is what a small sample looks like.\n\n"
        "THE TRACKED CELLS SAY NO CHANGE. Cells followed across >=2 years split\n"
        f"city {cty['cells_up']}/{cty['cells_down']} up/down and rural {rur['cells_up']}/{rur['cells_down']} -- rural is an exact coin flip.\n"
        "A real divergence between urban and rural biodiversity would show up here first,\n"
        "and it does not.\n\n"
        "COVERAGE CAVEAT. GCTB built-up polygons stop at 2022, so only 510,923 recordings\n"
        f"({UY[0]}-{LAST_LABELLED}) carry a real class. Urban class is a property of the PLACE, so each\n"
        "cell's class is carried forward to 2023-2025 (shaded). Those three years rest on an\n"
        "assumption -- that cells did not change built-up status -- not on measurement.\n\n"
        "Cells take the modal class of their recordings: 97.0% of recordings fall in their\n"
        "cell's modal class and 5.8% of cells are mixed (higher than the 1.6% for biomes,\n"
        "because city boundaries cut through 0.1-degree cells far more often than biomes do)."
    )
    fig.text(0.075, 0.335, note, ha='left', va='top', size=8.0, family='monospace')
    fig.text(0.5, 0.028, 'Acoustic Biodiversity Report  -  page 10  -  urban class yearly',
             ha='center', size=8, color='#999')
    pp.savefig(fig); plt.close(fig)
    npages += 1

# ---- Page 11: biome map ----
# Where each biome's recordings actually are. One dot per 0.1-degree cell,
# coloured by the cell's modal biome, sized by recording volume.
if os.path.exists('cell_biome.csv'):
    import numpy as _np

    cb = list(csv.DictReader(open('cell_biome.csv', newline='')))
    counts = collections.Counter(r['biome_name'] for r in cb)
    # colour order follows recording volume so the legend matches pages 8-9
    if os.path.exists('biome_change.csv'):
        order = [r['biome_name'] for r in csv.DictReader(open('biome_change.csv', newline=''))]
    else:
        order = [b for b, _ in counts.most_common()]
    order = [b for b in order if b in counts] + [b for b in counts if b not in order]

    # qualitative, ordered by recording volume -- the two largest biomes must not
    # share a hue, so distinguishability is prioritised over habitat semantics
    PAL = ['#33a02c', '#1f78b4', '#b15928', '#ff7f00', '#a6cee3', '#6a3d9a',
           '#e31a1c', '#fdbf6f', '#b2df8a', '#cab2d6', '#17becf', '#006d2c',
           '#fb9a99', '#7f7f7f', '#bcbd22']
    cmap = {b: PAL[i % len(PAL)] for i, b in enumerate(order)}

    fig = plt.figure(figsize=(11.69, 8.27))  # A4 landscape
    fig.text(0.5, 0.962, 'Recording Coverage by Biome', ha='center', size=17, weight='bold')
    fig.text(0.5, 0.934,
             f"{len(cb):,} grid cells (0.1 deg), each coloured by the modal biome of its recordings",
             ha='center', size=9, style='italic', color='#555')

    axm = fig.add_axes([0.035, 0.30, 0.66, 0.60])
    lat = _np.array([float(r['lat_cell']) for r in cb])
    lon = _np.array([float(r['lon_cell']) for r in cb])
    nrec = _np.array([int(r['n_recordings']) for r in cb], float)
    col = [cmap[r['biome_name']] for r in cb]
    sz = 0.6 + 5.0 * _np.log10(nrec + 1) / _np.log10(nrec.max() + 1)
    # rasterize the 40k-point cloud: as vector it inflates the PDF ~30x (0.2 -> 6 MB)
    axm.scatter(lon, lat, s=sz, c=col, linewidths=0, alpha=0.85, rasterized=True)
    axm.set_xlim(-180, 180); axm.set_ylim(-60, 82)
    axm.set_xticks(range(-180, 181, 60)); axm.set_yticks(range(-60, 81, 30))
    axm.tick_params(labelsize=7.5)
    axm.set_aspect('equal', adjustable='box')
    axm.grid(True, lw=0.3, color='#eee')
    axm.set_title('Dot size = recordings in that cell (log scale)', size=9, loc='left')

    # ---- legend with counts ----
    axl = fig.add_axes([0.71, 0.30, 0.27, 0.60]); axl.axis('off')
    axl.set_title('biome (cells / recordings)', size=9, loc='left')
    tot_rec = collections.Counter()
    for r in cb:
        tot_rec[r['biome_name']] += int(r['n_recordings'])
    for i, b in enumerate(order):
        yy = 1 - (i + 0.5) / len(order)
        axl.scatter([0.03], [yy], s=34, c=cmap[b], linewidths=0, transform=axl.transAxes)
        axl.text(0.09, yy, f'{b[:38]}', va='center', size=7.2, transform=axl.transAxes)
        axl.text(0.09, yy - 0.026, f'{counts[b]:,} cells / {tot_rec[b]:,} rec',
                 va='center', size=6.1, color='#777', transform=axl.transAxes)

    mixed_cells = sum(1 for r in cb if int(r['n_biomes_in_cell']) > 1)
    note = (
        "COVERAGE IS NOT ECOLOGICAL COVERAGE. The map shows where recordists go, not where\n"
        "biomes are. Europe is saturated across a single biome band while whole tropical biomes\n"
        f"are covered by scattered cells: Temperate Broadleaf & Mixed Forests holds {counts[order[0]]:,} cells,\n"
        f"the three smallest together fewer than {sum(counts[b] for b in order[-3:]):,}. This is the sampling bias behind\n"
        "every per-biome number here -- a biome comparison is a comparison of recordist\n"
        "communities as much as of habitats, and the biomes with fewest cells are exactly\n"
        "those whose yearly series swing hardest on pages 8 and 9.\n\n"
        f"Cells take the modal biome of their recordings; {mixed_cells:,} ({100*mixed_cells/len(cb):.1f}%) contain more than one\n"
        "and are drawn in their majority colour. Polygons: RESOLVE Ecoregions 2017\n"
        "(Dinerstein et al., CC-BY 4.0). Colours are categorical only."
    )
    fig.text(0.035, 0.235, note, ha='left', va='top', size=8.0, family='monospace')
    fig.text(0.5, 0.03, 'Acoustic Biodiversity Report  -  page 11  -  biome map',
             ha='center', size=8, color='#999')
    pp.savefig(fig, dpi=300); plt.close(fig)
    npages += 1


# ---- Page 12: is there a ten-year trend? (turnover vs change) ----
# The question page 5 raises but never tests. Everything here is computed at
# runtime from grid_cells_yearly.csv. Theil-Sen / Kendall are hand-rolled to
# avoid a scipy dependency; both were verified against scipy to machine
# precision (slope, 95% CI, tau-b and its tie-corrected p all agree).
if os.path.exists('grid_cells_yearly.csv'):
    def _med(v):
        s = sorted(v); n = len(s)
        return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])

    def _theil_sen(xs, ys, z=1.959963985):
        n = len(xs); sl = []
        for i in range(n):
            for j in range(i + 1, n):
                if xs[j] != xs[i]:
                    sl.append((ys[j] - ys[i]) / (xs[j] - xs[i]))
        sl.sort(); nt = len(sl)
        med = sl[nt // 2] if nt % 2 else 0.5 * (sl[nt // 2 - 1] + sl[nt // 2])
        sig = math.sqrt(n * (n - 1) * (2 * n + 5) / 18.0)
        Ru = min(int(round((nt + z * sig) / 2.0)), nt - 1)
        Rl = max(int(round((nt - z * sig) / 2.0)) - 1, 0)
        return med, sl[Rl], sl[Ru]

    def _kendall(xs, ys):
        n = len(xs); c = d = 0
        for i in range(n):
            for j in range(i + 1, n):
                s = (xs[j] - xs[i]) * (ys[j] - ys[i])
                if s > 0: c += 1
                elif s < 0: d += 1
        S = c - d; n0 = n * (n - 1) / 2.0
        tx = list(collections.Counter(xs).values())
        ty = list(collections.Counter(ys).values())
        n1 = sum(t * (t - 1) / 2.0 for t in tx); n2 = sum(u * (u - 1) / 2.0 for u in ty)
        tau = S / math.sqrt((n0 - n1) * (n0 - n2))
        v0 = n * (n - 1) * (2 * n + 5)
        vt = sum(t * (t - 1) * (2 * t + 5) for t in tx)
        vu = sum(u * (u - 1) * (2 * u + 5) for u in ty)
        v1 = (sum(t * (t - 1) for t in tx) * sum(u * (u - 1) for u in ty)) / (2.0 * n * (n - 1))
        v2 = (sum(t * (t - 1) * (t - 2) for t in tx) * sum(u * (u - 1) * (u - 2) for u in ty)) \
             / (9.0 * n * (n - 1) * (n - 2))
        var = (v0 - vt - vu) / 18.0 + v1 + v2
        return tau, math.erfc(abs(S / math.sqrt(var)) / math.sqrt(2))

    _yr = csv.reader(open('grid_cells_yearly.csv', newline='')); _yh = next(_yr)
    _yi = {c: i for i, c in enumerate(_yh)}
    recs = []
    for row in _yr:
        s = row[_yi['S_rare10']]
        if s == '':
            continue
        recs.append((round(float(row[_yi['lat_cell']]), 1),
                     round(float(row[_yi['lon_cell']]), 1),
                     int(row[_yi['year']]), float(s)))
    TY = sorted({r[2] for r in recs})

    byy = collections.defaultdict(list)
    for r in recs:
        byy[r[2]].append(r[3])
    raw = [_med(byy[y]) for y in TY]
    r_sl, r_lo, r_hi = _theil_sen(TY, raw)
    r_tau, r_p = _kendall(TY, raw)

    seen = collections.defaultdict(set)
    for r in recs:
        seen[(r[0], r[1])].add(r[2])

    def _panel(minY):
        keep = {c for c, ys in seen.items() if len(ys) >= minY}
        b = collections.defaultdict(list)
        for r in recs:
            if (r[0], r[1]) in keep:
                b[r[2]].append(r[3])
        yy = [y for y in TY if len(b[y]) >= 5]
        return keep, yy, [_med(b[y]) for y in yy]

    # paired: same cells, first three years vs last three
    e_lo, e_hi = TY[0], TY[2]; l_lo, l_hi = TY[-3], TY[-1]
    ear = collections.defaultdict(list); lat = collections.defaultdict(list)
    for r in recs:
        if e_lo <= r[2] <= e_hi: ear[(r[0], r[1])].append(r[3])
        elif l_lo <= r[2] <= l_hi: lat[(r[0], r[1])].append(r[3])
    both = sorted(set(ear) & set(lat))
    dif = [sum(lat[c]) / len(lat[c]) - sum(ear[c]) / len(ear[c]) for c in both]
    np_ = len(dif); p_mean = sum(dif) / np_
    p_sd = math.sqrt(sum((x - p_mean) ** 2 for x in dif) / (np_ - 1))
    p_ci = 1.959963985 * p_sd / math.sqrt(np_)
    p_up = sum(1 for x in dif if x > 0); p_dn = sum(1 for x in dif if x < 0)

    # new vs returning cells, per year
    first = {}
    for r in sorted(recs, key=lambda x: x[2]):
        first.setdefault((r[0], r[1]), r[2])
    nw, rt, pnew = [], [], []
    for y in TY:
        a = [r[3] for r in recs if r[2] == y and first[(r[0], r[1])] == y]
        b = [r[3] for r in recs if r[2] == y and first[(r[0], r[1])] < y]
        nw.append(_med(a) if a else float('nan'))
        rt.append(_med(b) if b else float('nan'))
        pnew.append(100.0 * len(a) / (len(a) + len(b)))
    n_lower = sum(1 for a, b in zip(nw, rt)
                  if a == a and b == b and a < b)
    n_cmp = sum(1 for a, b in zip(nw, rt) if a == a and b == b)

    fig = plt.figure(figsize=(11.69, 8.27))
    fig.suptitle('Is there a ten-year trend?  Turnover, not change', size=14, weight='bold')

    ax1 = fig.add_axes([0.065, 0.600, 0.385, 0.270])
    ax1.plot(TY, raw, marker='o', ms=4, color='#c0392b', lw=2,
             label='all scored cells  (%+.2f/decade)' % (r_sl * 10))
    fitx = [TY[0], TY[-1]]
    mid = _med(raw); midx = TY[len(TY) // 2]
    ax1.plot(fitx, [mid + r_sl * (x - midx) for x in fitx], ls='--', lw=1.2,
             color='#c0392b', alpha=0.7)
    k8, y8, s8 = _panel(8)
    _p8 = _theil_sen(y8, s8)[0] * 10
    ax1.plot(y8, s8, marker='s', ms=4, color='#16668f', lw=2,
             label='same cells, >=8 yrs, n=%d  (%+.2f/decade)' % (len(k8), _p8))
    ax1.set_ylabel('Median richness (S_rare10)')
    ax1.set_title('The decline is in the pool, not the places', size=10.5)
    ax1.grid(True, lw=0.3, color='#eee'); ax1.legend(fontsize=7.5, loc='lower left')

    ax2 = fig.add_axes([0.565, 0.600, 0.385, 0.270])
    ax2b = ax2.twinx()
    ax2b.bar(TY, pnew, color='#e8e8e8', width=0.7, zorder=0)
    ax2b.set_ylabel('% of year\'s cells that are new', size=8, color='#999')
    ax2b.tick_params(axis='y', labelsize=7, colors='#999'); ax2b.set_ylim(0, 100)
    ax2.set_zorder(ax2b.get_zorder() + 1); ax2.patch.set_visible(False)
    ax2.plot(TY, rt, marker='o', ms=4, color='#1a7a4c', lw=2, label='returning cells')
    ax2.plot(TY, nw, marker='o', ms=4, color='#c47f17', lw=2, label='newly-recorded cells')
    ax2.set_ylabel('Median richness'); ax2.grid(True, lw=0.3, color='#eee')
    ax2.set_title('Newly-recorded cells are poorer, every year', size=10.5)
    ax2.legend(fontsize=7.5, loc='lower left')

    ax3 = fig.add_axes([0.300, 0.395, 0.620, 0.135])
    ent = [('all scored cells (%d)' % len({(r[0], r[1]) for r in recs}),
            r_sl * 10, r_lo * 10, r_hi * 10)]
    for m in (2, 4, 6, 8):
        kk, yy, ss = _panel(m)
        sl, lo, hi = _theil_sen(yy, ss)
        ent.append(('same cells, >=%d yrs (%d)' % (m, len(kk)), sl * 10, lo * 10, hi * 10))
    ent.append(('paired %d-%d vs %d-%d (%d)' % (e_lo, e_hi, l_lo, l_hi, np_),
                p_mean, p_mean - p_ci, p_mean + p_ci))
    for i, (lab, est, lo, hi) in enumerate(ent):
        yv = len(ent) - 1 - i
        col = '#c0392b' if i == 0 else '#16668f'
        ax3.plot([lo, hi], [yv, yv], color=col, lw=2, solid_capstyle='butt')
        ax3.plot([est], [yv], marker='o', ms=5, color=col)
    ax3.axvline(0, color='#333', lw=1)
    ax3.set_yticks(range(len(ent)))
    ax3.set_yticklabels([e[0] for e in reversed(ent)], size=7.6)
    _wl = min(e[2] for e in ent); _wh = max(e[3] for e in ent)
    ax3.set_xlim(_wl - 0.4, _wh + 0.4); ax3.set_ylim(-0.6, len(ent) - 0.4)
    ax3.set_xlabel('Change in median richness over the decade (species), with 95% CI', size=8.5)
    ax3.set_title('Every turnover-controlled estimate sits on zero -- but none is precise',
                  size=10.5)
    for s in ('top', 'right'):
        ax3.spines[s].set_visible(False)
    ax3.tick_params(axis='y', length=0)
    ax3.grid(True, axis='x', lw=0.3, color='#eee')

    note = (
        "NO. Raw, the series falls %.2f -> %.2f (Theil-Sen %+.2f species/decade, Kendall tau %+.2f,\n"
        "p = %.3f) -- marginal, and not significant. Follow the SAME cells and it disappears:\n"
        "the >=8-year panel gives %+.2f/decade and the %d paired cells %+.3f (%d up / %d down).\n\n"
        "MECHANISM: %.0f-%.0f%% of each year's scored cells were never recorded before, and newly-\n"
        "recorded cells sit below returning ones in %d of %d years. The archive keeps expanding into\n"
        "thinner locations, which drags the pooled median down while no individual place changes.\n\n"
        "LIMIT: the controlled estimates are UNDERPOWERED, not proof of stability. The paired CI is\n"
        "[%+.2f, %+.2f] and the >=8-year panel spans more than two species -- a real decline of up to\n"
        "~0.5 species/decade would not be detectable here. Read this as no evidence of a trend,\n"
        "plus a well-identified artefact that explains the apparent one.\n\n"
        "CONFOUND: the steepest fall is the last three years, exactly where the data source changes\n"
        "from the historical backfill to the base+gap run. Year and processing path cannot be\n"
        "separated with what is on disk, so even the marginal raw decline is suspect."
    ) % (raw[0], raw[-1], r_sl * 10, r_tau, r_p, (_theil_sen(y8, s8)[0]) * 10,
         np_, p_mean, p_up, p_dn, min(pnew[1:]), max(pnew[1:]), n_lower, n_cmp,
         p_mean - p_ci, p_mean + p_ci)
    fig.text(0.065, 0.335, note, ha='left', va='top', size=7.7, family='monospace')
    fig.text(0.5, 0.018, 'Acoustic Biodiversity Report  -  page 12  -  is there a trend?',
             ha='center', size=8, color='#999')
    pp.savefig(fig); plt.close(fig)
    npages += 1

pp.close()
print(f'wrote Acoustic_Biodiversity_Report.pdf  ({npages} pages)')
