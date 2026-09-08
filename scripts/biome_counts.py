#!/usr/bin/env python
"""Per-biome / per-realm / per-ecoregion DATA POINT COUNTS.

Not means -- counts. For each group: how many recordings exist, and how many of
them actually carry usable values for the acoustic metrics.

A metric value counts as MISSING if it is blank, NaN, or +/-inf. Blanks are the
dominant failure mode in the gap-fill CSVs (rows written with trailing commas),
and +/-inf are pipeline silent-failure artifacts -- both must be excluded or the
counts overstate how much real data exists.

Inputs   ROOT/recordings_ecoregion.csv     (id, continent, biome_name, realm, eco_name)
         ROOT/score_<continent>_meta.csv   (id + the 40 index columns)
Outputs  ROOT/biome_data_counts.csv        per biome: n + complete + per-metric valid
         ROOT/realm_data_counts.csv
         ROOT/ecoregion_data_counts.csv    (all ecoregions, no n minimum)
"""
import os, csv
os.environ.setdefault('PROJ_DATA',
    '/work/mech-ai-scratch/alloy/.conda/envs/geo/share/proj')
import numpy as np
import pandas as pd

csv.field_size_limit(10**9)
ROOT = '/work/mech-ai-scratch/alloy/Acoustic_Indices'
CONTINENTS = ['africa', 'america', 'asia', 'australia', 'europe']

# the 11 indices build_cells.py actually aggregates downstream -- the set that
# matters for the biodiversity analysis
CORE = [
    'Acoustic_Complexity_Index__mean',
    'Acoustic_Diversity_Index__main_value',
    'Acoustic_Evenness_Index__main_value',
    'Bio_acoustic_Index__main_value',
    'Normalized_Difference_Sound_Index__main_value',
    'Spectral_Entropy__main_value',
    'Temporal_Entropy__main_value',
    'NB_peaks__main_value',
    'Wave_SNR__SNR',
    'Acoustic_Diversity_Index_NR__main_value',
    'Bio_acoustic_Index_NR__main_value',
]


def index_columns():
    """The 40 acoustic-index columns: everything between 'id' and 'gen'."""
    hdr = pd.read_csv(f'{ROOT}/score_africa_meta.csv', nrows=0).columns.tolist()
    return hdr[hdr.index('id') + 1: hdr.index('gen')]


def load_metrics(idx_cols):
    """id + continent + finite-ness of every index column, for all continents."""
    frames = []
    for c in CONTINENTS:
        keep = set(['id'] + idx_cols)
        df = pd.read_csv(f'{ROOT}/score_{c}_meta.csv',
                         usecols=lambda x: x in keep, engine='python',
                         on_bad_lines='skip', dtype=str)
        df['continent'] = c
        frames.append(df)
        print(f'  {c}: {len(df)} rows', flush=True)
    df = pd.concat(frames, ignore_index=True)
    for col in idx_cols:                       # blank -> NaN, inf -> NaN
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df[idx_cols] = df[idx_cols].replace([np.inf, -np.inf], np.nan)
    return df


def counts_by(df, key, idx_cols, path):
    """n, complete-case counts, and per-metric valid counts for one grouping."""
    ok_all = df[idx_cols].notna().all(axis=1)
    ok_core = df[CORE].notna().all(axis=1)
    n_ok = df[idx_cols].notna().sum(axis=1)

    g = df.groupby(key, dropna=False)
    out = pd.DataFrame({
        'n_recordings':        g.size(),
        'n_complete_all40':    ok_all.groupby(df[key], dropna=False).sum(),
        'n_complete_core11':   ok_core.groupby(df[key], dropna=False).sum(),
        'n_zero_metrics':      (n_ok == 0).groupby(df[key], dropna=False).sum(),
        'median_metrics_present': n_ok.groupby(df[key], dropna=False).median(),
    })
    out['pct_complete_all40'] = (100 * out.n_complete_all40 / out.n_recordings).round(1)
    out['pct_complete_core11'] = (100 * out.n_complete_core11 / out.n_recordings).round(1)
    # per-metric valid counts
    for col in idx_cols:
        out[f'valid__{col}'] = df[col].notna().groupby(df[key], dropna=False).sum()
    out = out.sort_values('n_recordings', ascending=False)
    out.to_csv(path)
    print(f'\nWrote {path}  ({len(out)} rows)', flush=True)
    return out


def main():
    idx_cols = index_columns()
    print(f'{len(idx_cols)} index columns; {len(CORE)} core (build_cells) columns\n')

    print('Loading metrics...', flush=True)
    met = load_metrics(idx_cols)

    print('\nLoading ecoregion labels...', flush=True)
    eco = pd.read_csv(f'{ROOT}/recordings_ecoregion.csv',
                      usecols=['id', 'continent', 'biome_name', 'realm',
                               'eco_name', 'eco_match'], dtype=str)
    print(f'  {len(eco)} labeled recordings', flush=True)

    # join on (id, continent): ids are not globally unique (asia/australia overlap)
    df = eco.merge(met, on=['id', 'continent'], how='left')
    print(f'  joined: {len(df)} rows, '
          f'{df[idx_cols[0]].notna().sum()} with metric data', flush=True)

    for key, path in [('biome_name', f'{ROOT}/biome_data_counts.csv'),
                      ('realm',      f'{ROOT}/realm_data_counts.csv'),
                      ('eco_name',   f'{ROOT}/ecoregion_data_counts.csv')]:
        t = counts_by(df, key, idx_cols, path)
        if key != 'eco_name':
            with pd.option_context('display.width', 200):
                print(t[['n_recordings', 'n_complete_all40', 'pct_complete_all40',
                         'n_complete_core11', 'pct_complete_core11',
                         'n_zero_metrics']].to_string())

    # overall totals
    ok_all = df[idx_cols].notna().all(axis=1)
    ok_core = df[CORE].notna().all(axis=1)
    print(f'\n=== TOTAL over {len(df)} labeled recordings ===')
    print(f'  complete on all 40 metrics : {int(ok_all.sum()):>7,} '
          f'({100*ok_all.mean():.1f}%)')
    print(f'  complete on core 11        : {int(ok_core.sum()):>7,} '
          f'({100*ok_core.mean():.1f}%)')
    print(f'  zero metrics present       : {int((df[idx_cols].notna().sum(axis=1)==0).sum()):>7,}')

    print('\nPer-metric valid counts (all labeled recordings):')
    vc = df[idx_cols].notna().sum().sort_values()
    vc = pd.DataFrame({'valid': vc, 'missing': len(df) - vc})
    vc['pct_missing'] = (100 * vc.missing / len(df)).round(2)
    print(vc.to_string())


if __name__ == '__main__':
    main()
