# Acoustic Biodiversity Mapping

A global, location-level biodiversity metric derived from ~192k Xeno-canto wildlife
recordings and their metadata. See **[report/Acoustic_Biodiversity_Report.pdf](report/Acoustic_Biodiversity_Report.pdf)**
for the full write-up, world map, and results.

## The metric

- **Location** = a 0.1° (~11 km) grid cell.
- **Biodiversity value** = `S_rare10` — the expected number of distinct recorded species
  in 10 recordings (Hurlbert sample-based rarefaction; controls for recording effort).
- **Label** = good / moderate / bad, by tertile of `S_rare10` *within the cell's 10° latitude
  band* (region-relative, so a score means "rich for its region").
- **Scope** = cells with ≥ 10 recordings (3,277 of 15,743 cells scored confidently).

## Headline finding

The 42 acoustic indices (ACI, ADI, NDSI, Bioacoustic Index, …) **do not predict species
richness** on this archive (best Spearman ρ ≈ 0.11). They were designed for passive
soundscape monitoring, but 69% of the archive is single-target recordings (median 24 s).
The defensible biodiversity metric is therefore **effort-controlled species richness from
the metadata**, not any acoustic index. Full limitations are in the report.

## Pipeline (scripts/)

| Script | Purpose |
|--------|---------|
| `convert_metadata.py` | Flatten Xeno-canto metadata JSON pages → one CSV (274,301 records) |
| `merge_meta.py` | Left-join metadata onto each `score_<continent>.csv` by recording `id` |
| `build_cells.py` | Aggregate to 0.1° cells, compute rarefied richness, correlate indices vs richness |
| `make_report.py` | Render the multi-page PDF report |

The upstream acoustic-index computation (`compute_indice.py`, `acoustic_index.py`, the SLURM
jobs, and the raw audio/score CSVs) lives in the parent processing project and is **not**
included here — only the analysis layer and its outputs.

## Outputs (outputs/)

| File | Contents |
|------|----------|
| `grid_cells.csv` | Per-cell acoustic features + richness (all 15,743 cells) |
| `grid_cells_labeled.csv` | Global-tertile good/moderate/bad labels |
| `grid_cells_labeled_regional.csv` | Region-relative (latitude-band) labels — **primary** |
| `biodiversity_map.csv` | `lat, lon, continent, n_rec, richness, label, label_code` — ready for GIS/folium |

## Reproduce

```bash
python3 scripts/convert_metadata.py <metadata_dir> metadata.csv
python3 scripts/merge_meta.py            # produces score_<continent>_meta.csv
python3 scripts/build_cells.py           # produces grid_cells.csv + correlations
python3 scripts/make_report.py           # produces the PDF
```
