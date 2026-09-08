#!/usr/bin/env python
"""Point-in-polygon label every Xeno-canto recording with its terrestrial
ecoregion / biome / biogeographic realm.

Two reference layers are applied to the same points, so results can be reported
against either:

  RESOLVE Ecoregions 2017  (Dinerstein et al. 2017, CC-BY 4.0)   -- PRIMARY
      847 polygons, one row per ecoregion, human-readable BIOME_NAME/REALM.
  WWF TEOW 2001            (Olson et al. 2001, CC BY-NC 3.0)     -- COMPARISON
      14,458 polygon parts over 827 ecoregions; BIOME/REALM are numeric/2-letter
      codes, mapped to names here. This is the layer published on Data Basin
      (dataset 68635d7c77f1475f9b6c1d1dbe0a4c4c).

Unlike the GCTB urban layer, ecoregions carry NO year dimension, so every
recording with coordinates is classifiable -- no 2015-2022 restriction.

Points that fall in no polygon (offshore GPS, small islands, coastal offsets)
get a nearest-polygon fallback capped at NEAREST_MAX_DEG; beyond that they are
left unassigned.

Inputs
  ROOT/score_<continent>_meta.csv        (id, lat, lon, date, cnt + 40 index cols)
  ECO/resolve_2017/Ecoregions2017.shp    (WGS-84)
  ECO/wwf_teow_2001/wwf_terr_ecos.shp    (WGS-84)

Outputs (written to ROOT)
  recordings_ecoregion.csv       one row per recording with coordinates
  ecoregion_biome_summary.csv    per RESOLVE biome: n + mean/median of key indices
  ecoregion_summary.csv          per RESOLVE ecoregion (n >= MIN_N)
  ecoregion_realm_summary.csv    per realm
"""
import os, csv, sys
os.environ.setdefault('PROJ_DATA',
    '/work/mech-ai-scratch/alloy/.conda/envs/geo/share/proj')
import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from shapely import STRtree

csv.field_size_limit(10**9)

ROOT = '/work/mech-ai-scratch/alloy/Acoustic_Indices'
ECO = '/work/mech-ai-scratch/alloy/ECOREGIONS'
RESOLVE_SHP = f'{ECO}/resolve_2017/Ecoregions2017.shp'
TEOW_SHP = f'{ECO}/wwf_teow_2001/wwf_terr_ecos.shp'
CONTINENTS = ['africa', 'america', 'asia', 'australia', 'europe']

CHUNK = 100_000          # points per STRtree query batch
NEAREST_MAX_DEG = 0.1    # ~11 km at the equator; fallback cap for unmatched pts
MIN_N = 30               # min recordings for a per-ecoregion summary row

# same key indices the urban classifier summarises, for direct comparability
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

# Olson et al. 2001 biome codes (TEOW 'BIOME' field)
TEOW_BIOMES = {
    1: 'Tropical & Subtropical Moist Broadleaf Forests',
    2: 'Tropical & Subtropical Dry Broadleaf Forests',
    3: 'Tropical & Subtropical Coniferous Forests',
    4: 'Temperate Broadleaf & Mixed Forests',
    5: 'Temperate Conifer Forests',
    6: 'Boreal Forests/Taiga',
    7: 'Tropical & Subtropical Grasslands, Savannas & Shrublands',
    8: 'Temperate Grasslands, Savannas & Shrublands',
    9: 'Flooded Grasslands & Savannas',
    10: 'Montane Grasslands & Shrublands',
    11: 'Tundra',
    12: 'Mediterranean Forests, Woodlands & Scrub',
    13: 'Deserts & Xeric Shrublands',
    14: 'Mangroves',
    98: 'Lake',
    99: 'Rock and Ice',
}
TEOW_REALMS = {
    'AT': 'Afrotropic', 'AA': 'Australasia', 'IM': 'Indo-Malay',
    'NA': 'Nearctic', 'NT': 'Neotropic', 'OC': 'Oceania',
    'PA': 'Palearctic', 'AN': 'Antarctic',
}


def load_recordings():
    """Read all continents' meta CSVs -> DataFrame of recordings with coords."""
    keep = ['id', 'lat', 'lon', 'date', 'cnt'] + SUMMARY_COLS
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
    df = df.dropna(subset=['lat', 'lon']).copy()
    # guard against out-of-range coords
    df = df[df['lat'].between(-90, 90) & df['lon'].between(-180, 180)].copy()
    df = df.reset_index(drop=True)
    print(f'  with usable coordinates: {len(df)} of {before} '
          f'({before - len(df)} dropped)', flush=True)
    return df


def assign(points, shp, attr_cols, label):
    """Point-in-polygon assign attr_cols from shp to each point.

    Returns (DataFrame of attributes aligned to points, match-type Series).
    Uses an STRtree with the 'intersects' predicate -- shapely applies the
    predicate as point.<pred>(polygon), so 'intersects' is the correct
    point-in-polygon test ('covers' would be point.covers(polygon) == always
    False for a point vs a polygon).
    """
    print(f'[{label}] reading {os.path.basename(shp)} ...', flush=True)
    polys = gpd.read_file(shp, engine='pyogrio', columns=attr_cols)
    geoms = polys.geometry.values
    print(f'[{label}] {len(polys)} polygons; preparing + indexing ...', flush=True)
    shapely.prepare(geoms)
    tree = STRtree(geoms)

    n = len(points)
    poly_idx = np.full(n, -1, dtype=np.int64)
    for start in range(0, n, CHUNK):
        stop = min(start + CHUNK, n)
        pi, gi = tree.query(points[start:stop], predicate='intersects')
        if len(pi):
            # first polygon wins where a point sits on a shared boundary
            first = np.full(stop - start, -1, dtype=np.int64)
            # iterate in reverse so the lowest gi lands last (== wins)
            order = np.argsort(-gi, kind='stable')
            first[pi[order]] = gi[order]
            poly_idx[start:stop] = first
        print(f'[{label}]   {stop}/{n} points  '
              f'({(poly_idx[:stop] >= 0).sum()} matched)', flush=True)

    match = pd.Series(np.where(poly_idx >= 0, 'exact', 'none'), name='match')

    # ---- nearest-polygon fallback for unmatched points ----
    miss = np.flatnonzero(poly_idx < 0)
    if len(miss):
        print(f'[{label}] {len(miss)} unmatched -> nearest-polygon fallback '
              f'(cap {NEAREST_MAX_DEG} deg)', flush=True)
        near = tree.nearest(points[miss])
        dist = shapely.distance(points[miss], geoms[near])
        ok = dist <= NEAREST_MAX_DEG
        poly_idx[miss[ok]] = near[ok]
        match.iloc[miss[ok]] = 'nearest'
        print(f'[{label}]   {int(ok.sum())} recovered within cap, '
              f'{int((~ok).sum())} left unassigned', flush=True)

    attrs = polys[attr_cols].reset_index(drop=True)
    out = attrs.reindex(poly_idx).reset_index(drop=True)
    out[poly_idx < 0] = np.nan
    del polys, geoms, tree
    return out, match


def summarise(df, by, path, min_n=1):
    """Per-group n + mean/median of the key indices -> CSV."""
    rows = []
    for key, g in df.groupby(by, dropna=True):
        if len(g) < min_n:
            continue
        rec = {}
        if isinstance(by, list):
            for k, v in zip(by, key if isinstance(key, tuple) else (key,)):
                rec[k] = v
        else:
            rec[by] = key
        rec['n_recordings'] = len(g)
        for col in SUMMARY_COLS:
            rec[f'{col}__mean'] = g[col].mean()
            rec[f'{col}__median'] = g[col].median()
        rows.append(rec)
    out = pd.DataFrame(rows).sort_values('n_recordings', ascending=False)
    out.to_csv(path, index=False)
    print(f'\nWrote {path}  ({len(out)} rows)', flush=True)
    return out


def main():
    print('Loading recordings...', flush=True)
    df = load_recordings()

    points = shapely.points(df['lon'].values, df['lat'].values)

    # ---- RESOLVE Ecoregions 2017 (primary) ----
    res, res_match = assign(
        points, RESOLVE_SHP,
        ['ECO_NAME', 'BIOME_NAME', 'BIOME_NUM', 'REALM', 'ECO_ID', 'NNH_NAME'],
        'RESOLVE2017')
    df['eco_name'] = res['ECO_NAME'].values
    df['biome_name'] = res['BIOME_NAME'].values
    df['biome_num'] = res['BIOME_NUM'].values
    df['realm'] = res['REALM'].values
    df['eco_id'] = res['ECO_ID'].values
    df['nnh_name'] = res['NNH_NAME'].values
    df['eco_match'] = res_match.values

    # ---- WWF TEOW 2001 (the Data Basin layer, for comparison) ----
    teow, teow_match = assign(
        points, TEOW_SHP, ['ECO_NAME', 'BIOME', 'REALM', 'eco_code'], 'TEOW2001')
    df['teow_eco_name'] = teow['ECO_NAME'].values
    df['teow_biome_name'] = (pd.to_numeric(teow['BIOME'], errors='coerce')
                             .map(TEOW_BIOMES).values)
    df['teow_realm'] = teow['REALM'].map(TEOW_REALMS).values
    df['teow_eco_code'] = teow['eco_code'].values
    df['teow_match'] = teow_match.values

    # ---- per-recording labeled CSV ----
    out_cols = ['id', 'continent', 'cnt', 'lat', 'lon', 'year',
                'eco_id', 'eco_name', 'biome_name', 'realm', 'nnh_name',
                'eco_match',
                'teow_eco_code', 'teow_eco_name', 'teow_biome_name',
                'teow_realm', 'teow_match']
    out_path = f'{ROOT}/recordings_ecoregion.csv'
    df[out_cols].to_csv(out_path, index=False)
    print(f'\nWrote {out_path}  ({len(df)} rows)', flush=True)

    print('\nRESOLVE 2017 match types:')
    print(df['eco_match'].value_counts().to_string(), flush=True)
    print('\nTEOW 2001 match types:')
    print(df['teow_match'].value_counts().to_string(), flush=True)

    # +/-inf index values are pipeline silent-failure artifacts; drop them so
    # they don't poison the means (same guard as classify_urban.py).
    df[SUMMARY_COLS] = df[SUMMARY_COLS].replace([np.inf, -np.inf], np.nan)

    summarise(df, 'biome_name', f'{ROOT}/ecoregion_biome_summary.csv')
    summarise(df, 'realm', f'{ROOT}/ecoregion_realm_summary.csv')
    summarise(df, 'eco_name', f'{ROOT}/ecoregion_summary.csv', min_n=MIN_N)

    # ---- readable console tables ----
    with pd.option_context('display.width', 220, 'display.max_columns', 40,
                           'display.max_rows', 60):
        counts = (df.groupby('biome_name').size().sort_values(ascending=False)
                    .rename('n_recordings').to_frame())
        counts['pct'] = (100 * counts['n_recordings'] / len(df)).round(2)
        print('\nRecordings per RESOLVE biome:')
        print(counts.to_string())

        print('\nRecordings per realm:')
        rc = (df.groupby('realm').size().sort_values(ascending=False)
                .rename('n_recordings').to_frame())
        rc['pct'] = (100 * rc['n_recordings'] / len(df)).round(2)
        print(rc.to_string())

        print('\nBiome x continent counts:')
        print(df.pivot_table(index='biome_name', columns='continent',
                             values='id', aggfunc='count',
                             fill_value=0).to_string())

    # agreement between the two layers at biome level
    both = df.dropna(subset=['biome_name', 'teow_biome_name'])
    agree = (both['biome_name'] == both['teow_biome_name']).mean()
    print(f'\nRESOLVE-vs-TEOW biome agreement: {100*agree:.1f}% '
          f'of {len(both)} doubly-labeled recordings', flush=True)


if __name__ == '__main__':
    main()
