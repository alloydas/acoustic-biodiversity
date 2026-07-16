#!/usr/bin/env python
"""Point-in-polygon classify each Xeno-canto recording as city / town / rural
against the GCTB (Global City & Town Boundaries) 30 m polygons.

Restricted to recordings dated 2015-2022, each matched to its EXACT-year GCTB
Cities/Towns shapefiles (GCTB coverage stops at 2022). City takes priority over
town where polygons overlap.

Inputs
  /work/mech-ai-scratch/alloy/Acoustic_Indices/score_<continent>_meta.csv
      (id, lat, lon, date, + 40 acoustic-index columns)
  /work/mech-ai-scratch/alloy/GCTB/Cities_*/Cities_<year>.shp
  /work/mech-ai-scratch/alloy/GCTB/Towns_*/Towns_<year>.shp   (WGS-84)

Outputs (written next to this script)
  recordings_urban_class.csv   one row per classifiable recording:
      id, continent, lat, lon, year, urban_class
  urban_class_summary.csv      per-class count + mean/median of key indices
"""
import os, csv, glob, sys
os.environ.setdefault('PROJ_DATA',
    '/work/mech-ai-scratch/alloy/.conda/envs/geo/share/proj')
import pandas as pd
import geopandas as gpd
from shapely import STRtree

csv.field_size_limit(10**9)

ROOT = '/work/mech-ai-scratch/alloy/Acoustic_Indices'
GCTB = '/work/mech-ai-scratch/alloy/GCTB'
CONTINENTS = ['africa', 'america', 'asia', 'australia', 'europe']
YEARS = list(range(2015, 2023))  # 2015..2022 inclusive

# key indices to summarise by class (subset of the 40 columns)
SUMMARY_COLS = [
    'Acoustic_Complexity_Index__main_value',
    'Bio_acoustic_Index__main_value',
    'Acoustic_Diversity_Index__main_value',
    'Acoustic_Evenness_Index__main_value',
    'Normalized_Difference_Sound_Index__main_value',
    'Spectral_Entropy__main_value',
    'Temporal_Entropy__main_value',
    'NB_peaks__main_value',
    'Wave_SNR__SNR',
    'Spectral_centroid__mean',
    'ZCR__mean',
    'RMS_energy__mean',
]


def find_shp(kind, year):
    """kind in {'Cities','Towns'}; return path to that year's shapefile."""
    hits = glob.glob(f'{GCTB}/{kind}_*/{kind}_{year}.shp')
    if not hits:
        raise FileNotFoundError(f'{kind}_{year}.shp not found under {GCTB}')
    return hits[0]


def load_recordings():
    """Read all continents' meta CSVs -> DataFrame of classifiable recordings."""
    keep = ['id', 'lat', 'lon', 'date'] + SUMMARY_COLS
    frames = []
    for c in CONTINENTS:
        path = f'{ROOT}/score_{c}_meta.csv'
        df = pd.read_csv(path, usecols=lambda x: x in keep, engine='python',
                         on_bad_lines='skip', dtype=str)
        df['continent'] = c
        frames.append(df)
        print(f'  {c}: {len(df)} rows', flush=True)
    df = pd.concat(frames, ignore_index=True)
    df['lat'] = pd.to_numeric(df['lat'], errors='coerce')
    df['lon'] = pd.to_numeric(df['lon'], errors='coerce')
    df['year'] = pd.to_numeric(df['date'].str[:4], errors='coerce')
    for col in SUMMARY_COLS:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    before = len(df)
    df = df.dropna(subset=['lat', 'lon', 'year'])
    df = df[df['year'].between(2015, 2022)].copy()
    df['year'] = df['year'].astype(int)
    print(f'  classifiable (coord + 2015-2022): {len(df)} of {before}', flush=True)
    return df


def in_polygons(pts_gdf, shp_path):
    """Return a boolean Series (index-aligned to pts_gdf) True where the point
    falls inside any polygon of shp_path, via an STRtree + 'intersects'
    (point-in-polygon, boundary inclusive)."""
    polys = gpd.read_file(shp_path, engine='pyogrio', columns=[])  # geometry only
    geoms = polys.geometry.values
    tree = STRtree(geoms)
    pt_geoms = pts_gdf.geometry.values
    # candidate (point, polygon) pairs whose bboxes intersect, exact-confirmed
    # by 'intersects'. shapely applies it as point.intersects(polygon).
    pi, gi = tree.query(pt_geoms, predicate='intersects')
    inside = pd.Series(False, index=pts_gdf.index)
    if len(pi):
        inside.iloc[pd.unique(pi)] = True
    del polys, geoms, tree
    return inside


def main():
    print('Loading recordings...', flush=True)
    df = load_recordings()

    df['urban_class'] = 'rural'
    for year in YEARS:
        sub = df[df['year'] == year]
        if sub.empty:
            continue
        gdf = gpd.GeoDataFrame(
            sub, geometry=gpd.points_from_xy(sub['lon'], sub['lat']),
            crs='EPSG:4326')
        city_shp = find_shp('Cities', year)
        town_shp = find_shp('Towns', year)
        print(f'[{year}] {len(gdf)} pts  cities={os.path.basename(city_shp)} '
              f'towns={os.path.basename(town_shp)}', flush=True)
        in_city = in_polygons(gdf, city_shp)
        in_town = in_polygons(gdf, town_shp)
        cls = pd.Series('rural', index=gdf.index)
        cls[in_town] = 'town'
        cls[in_city] = 'city'   # city wins over town on overlap
        df.loc[gdf.index, 'urban_class'] = cls
        vc = cls.value_counts()
        print(f'[{year}] city={vc.get("city",0)} town={vc.get("town",0)} '
              f'rural={vc.get("rural",0)}', flush=True)

    # ---- per-recording labeled CSV ----
    out = df[['id', 'continent', 'lat', 'lon', 'year', 'urban_class']]
    out_path = f'{ROOT}/recordings_urban_class.csv'
    out.to_csv(out_path, index=False)
    print(f'\nWrote {out_path}  ({len(out)} rows)', flush=True)
    print(df['urban_class'].value_counts().to_string(), flush=True)

    # ---- index-by-class summary ----
    import numpy as np
    # a few recordings carry +/-inf index values (pipeline silent-failure /
    # log artifacts); treat them as missing so they don't poison the means.
    df[SUMMARY_COLS] = df[SUMMARY_COLS].replace([np.inf, -np.inf], np.nan)
    rows = []
    for cls, g in df.groupby('urban_class'):
        rec = {'urban_class': cls, 'n_recordings': len(g)}
        for col in SUMMARY_COLS:
            rec[f'{col}__mean'] = g[col].mean()
            rec[f'{col}__median'] = g[col].median()
        rows.append(rec)
    summ = pd.DataFrame(rows).set_index('urban_class')
    summ_path = f'{ROOT}/urban_class_summary.csv'
    summ.to_csv(summ_path)
    print(f'\nWrote {summ_path}', flush=True)
    with pd.option_context('display.width', 200, 'display.max_columns', 50):
        print(summ[['n_recordings'] + [f'{c}__mean' for c in SUMMARY_COLS]].T.to_string())

    # ---- per-YEAR x class summary ----
    rows = []
    for (yr, cls), g in df.groupby(['year', 'urban_class']):
        rec = {'year': yr, 'urban_class': cls, 'n_recordings': len(g)}
        for col in SUMMARY_COLS:
            rec[f'{col}__mean'] = g[col].mean()
            rec[f'{col}__median'] = g[col].median()
        rows.append(rec)
    yr_summ = pd.DataFrame(rows)
    yr_path = f'{ROOT}/urban_class_by_year.csv'
    yr_summ.to_csv(yr_path, index=False)
    print(f'\nWrote {yr_path}  ({len(yr_summ)} year x class rows)', flush=True)

    # readable counts pivot: rows=year, cols=class
    counts = (df.groupby(['year', 'urban_class']).size()
                .unstack('urban_class').fillna(0).astype(int))
    counts = counts.reindex(columns=['city', 'town', 'rural'])
    counts['total'] = counts.sum(axis=1)
    counts['city_%'] = (100 * counts['city'] / counts['total']).round(1)
    counts['town_%'] = (100 * counts['town'] / counts['total']).round(1)
    print('\nRecordings per year by class:')
    print(counts.to_string())


if __name__ == '__main__':
    main()
