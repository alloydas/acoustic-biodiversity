# Acoustic Biodiversity Mapping

A global, location-level biodiversity metric derived from ~760k Xeno-canto wildlife
recordings (2015–2025) and their metadata. See **[report/Acoustic_Biodiversity_Report.pdf](report/Acoustic_Biodiversity_Report.pdf)**
for the full write-up, world map, and results.

## The metric

- **Location** = a 0.1° (~11 km) grid cell.
- **Biodiversity value** = `S_rare10` — the expected number of distinct recorded species
  in 10 recordings (Hurlbert sample-based rarefaction; controls for recording effort).
- **Label** = good / moderate / bad, by tertile of `S_rare10` *within the cell's 10° latitude
  band* (region-relative, so a score means "rich for its region").
- **Scope** = cells with ≥ 10 recordings (11,088 of 40,914 cells scored confidently).

## Headline finding

The 42 acoustic indices (ACI, ADI, NDSI, Bioacoustic Index, …) **do not predict species
richness** on this archive (best Spearman ρ ≈ 0.13). They were designed for passive
soundscape monitoring, but 69% of the archive is single-target recordings (median 24 s).
The defensible biodiversity metric is therefore **effort-controlled species richness from
the metadata**, not any acoustic index. Full limitations are in the report.

## Pipeline (scripts/)

| Script | Purpose |
|--------|---------|
| `download.py` | Download Xeno-canto recordings for one continent (`--continent`, optional `--quality`/`--year`/`--months`) |
| `convert_metadata.py` | Flatten Xeno-canto metadata JSON pages → one CSV (766,747 records) |
| `merge_meta.py` | Left-join metadata onto each `score_<continent>.csv` (base **+ gap + 10-year hist**) by recording `id` |
| `build_cells.py` | Aggregate to 0.1° cells, compute rarefied richness, correlate indices vs richness |
| `build_cells_yearly.py` | Same metric split by calendar year (2015–2025) + per-cell change |
| `build_labels.py` | Assign good/moderate/bad tertile labels (global + region-relative) → `grid_cells_labeled*.csv`, `biodiversity_map.csv` |
| `make_report.py` | Render the multi-page PDF report |

The upstream acoustic-index computation (`compute_indice.py`, `acoustic_index.py`, the SLURM
jobs, and the raw audio/score CSVs) lives in the parent processing project and is **not**
included here — only the analysis layer and its outputs.

## Outputs (outputs/)

| File | Contents |
|------|----------|
| `grid_cells.csv` | Per-cell acoustic features + richness (all 40,914 cells) |
| `grid_cells_labeled.csv` | Global-tertile good/moderate/bad labels |
| `grid_cells_labeled_regional.csv` | Region-relative (latitude-band) labels — **primary** |
| `biodiversity_map.csv` | `lat, lon, continent, n_rec, richness, label, label_code` — ready for GIS/folium |
| `grid_cells_yearly.csv` | Per-(cell, year) richness for 2015–2025 |
| `cell_change.csv` | Cells scored in ≥2 years with first/last richness and trend |

### Temporal note (2015–2025)

The metric is also computed per year (report page 5). **Even an 11-year window cannot reveal
real biodiversity change** — year-to-year movement reflects *which* cells were recorded and by
*whom* each year (effort + observer turnover), not ecological gain/loss. Treat the yearly
slices as sampling-coverage diagnostics, not a biodiversity time series. 2,697 cells are
trackable across ≥2 years (only 10 across all eleven), and per-cell change is a near coin-flip
(1,176 up / 1,188 down / 333 flat) — the signature of noise, not a trend.

## Reproduce

```bash
# download recordings (one continent at a time; matches the existing dataset)
python3 scripts/download.py --continent africa     # ... america asia australia europe
python3 scripts/convert_metadata.py <metadata_dir> metadata.csv
python3 scripts/merge_meta.py            # produces score_<continent>_meta.csv (base + gap + hist)
python3 scripts/build_cells.py           # produces grid_cells.csv + correlations
python3 scripts/build_cells_yearly.py    # produces grid_cells_yearly.csv + cell_change.csv
python3 scripts/build_labels.py          # produces grid_cells_labeled*.csv + biodiversity_map.csv
python3 scripts/make_report.py           # produces the PDF
```
