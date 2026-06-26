"""Generate a multi-page PDF report of the acoustic-biodiversity analysis."""
import csv
import math
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


# ================= BUILD PDF =================
pp = PdfPages('Acoustic_Biodiversity_Report.pdf')

# ---- Page 1: title + summary ----
fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
fig.text(0.5, 0.93, 'Acoustic Biodiversity Mapping', ha='center', size=22, weight='bold')
fig.text(0.5, 0.90, 'A global location-level biodiversity metric from Xeno-canto recordings',
         ha='center', size=11, style='italic', color='#555')
fig.text(0.5, 0.875, 'Generated 2026-06-26', ha='center', size=9, color='#888')

summary = (
    "OBJECTIVE\n"
    "Derive a single metric that classifies a geographic location as biodiversity 'good' or 'bad',\n"
    "using a global archive of wildlife sound recordings and their metadata.\n\n"
    "DATA\n"
    "  -  192,196 recordings across 5 continents (Xeno-canto), each with 42 acoustic indices\n"
    "     computed by a SLURM batch pipeline (compute_indice.py over Butterworth-filtered audio).\n"
    "  -  Metadata flattened from 274,301 API records (species, coordinates, quality, device).\n"
    "  -  Acoustic scores joined to metadata on recording id (>99.99% match).\n\n"
    "METHOD\n"
    "  1. Aggregate all recordings into 0.1-degree (~11 km) grid cells -> 15,743 cells.\n"
    "  2. Per cell, compute species richness from recorded species (genus+species and the\n"
    "     'also' co-occurring-species field), corrected for sampling effort by Hurlbert\n"
    "     rarefaction to 10 recordings (S_rare10).\n"
    "  3. Test whether the acoustic indices predict richness (validation against ground truth).\n"
    "  4. Label each well-sampled cell good / moderate / bad relative to its 10-degree\n"
    "     latitude band (region-relative tertiles).\n\n"
    "HEADLINE RESULT\n"
    f"  -  3,277 cells had enough recordings (>=10) to score confidently.\n"
    f"  -  The acoustic indices did NOT predict species richness (best |Spearman rho| ~ 0.11).\n"
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
ax.set_title('Region-relative biodiversity label per 0.1-degree cell (3,277 well-sampled cells)', size=13)
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
    "  scope              = cells with >= 10 recordings (3,277 of 15,743 cells)\n\n"
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
    "  4. Coverage. Only 21% of cells are confidently scored; the rest are\n"
    "     'insufficient_data' (< 10 recordings).\n"
)
fig.text(0.08, 0.46, limit_txt, ha='left', va='top', size=9.0, family='monospace')

out_txt = (
    "OUTPUT FILES\n"
    "  grid_cells.csv                    - per-cell features + richness (all 15,743 cells)\n"
    "  grid_cells_labeled_regional.csv   - region-relative good/moderate/bad labels\n"
    "  biodiversity_map.csv              - lat, lon, label, code -> ready for GIS/folium\n"
)
fig.text(0.08, 0.16, out_txt, ha='left', va='top', size=9.0, family='monospace')
fig.text(0.5, 0.03, 'Acoustic Biodiversity Report  -  page 4  -  metric & results', ha='center', size=8, color='#999')
pp.savefig(fig); plt.close(fig)

pp.close()
print('wrote Acoustic_Biodiversity_Report.pdf  (4 pages)')
