#!/usr/bin/env python
"""India-focused figures for the urban (city/town/rural) classification.

Produces two PNGs under figures/:
  india_recordings_yearwise.png  -- per-year counts (2015-2025) + city/town/rural
                                    split (2015-2022, two panels)
  india_recordings_map.png       -- India map with each georeferenced recording
                                    coloured by its urban class

Inputs (in the processing project, not committed here):
  score_<continent>_meta.csv     -- id, cnt (country), date, lat, lon, indices
  recordings_urban_class.csv     -- id, urban_class  (from classify_urban.py)
  GCTB/ref/ne_10m_admin_0_countries.zip  -- Natural Earth country outlines

Run with the `geo` conda env (geopandas/shapely/pyogrio); see classify_urban.py.
"""
import os, glob
os.environ.setdefault('PROJ_DATA',
    '/work/mech-ai-scratch/alloy/.conda/envs/geo/share/proj')
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROC = '/work/mech-ai-scratch/alloy/Acoustic_Indices'
NE = '/vsizip//work/mech-ai-scratch/alloy/GCTB/ref/ne_10m_admin_0_countries.zip'
FIG = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(FIG, exist_ok=True)

COUNTRY = 'India'
COLORS = {'rural': '#4C9A6B', 'town': '#E0A458', 'city': '#C1503A'}


def load_country():
    """recordings for COUNTRY, merged with urban_class (2015-2022 where present)."""
    frames = []
    for f in glob.glob(f'{PROC}/score_*_meta.csv'):
        d = pd.read_csv(f, usecols=['id', 'cnt', 'date', 'lat', 'lon'],
                        engine='python', on_bad_lines='skip', dtype=str)
        frames.append(d[d['cnt'] == COUNTRY])
    df = pd.concat(frames, ignore_index=True)
    df['lat'] = pd.to_numeric(df['lat'], errors='coerce')
    df['lon'] = pd.to_numeric(df['lon'], errors='coerce')
    df['year'] = pd.to_numeric(df['date'].str[:4], errors='coerce')
    cls = pd.read_csv(f'{PROC}/recordings_urban_class.csv',
                      usecols=['id', 'urban_class'], dtype={'id': str})
    return df.merge(cls, on='id', how='left')


def fig_yearwise(df):
    df = df.dropna(subset=['year']).copy()
    df['year'] = df['year'].astype(int)
    ally, urby = list(range(2015, 2026)), list(range(2015, 2023))
    per_year = df.groupby('year').size().reindex(ally, fill_value=0)
    pivot = (df[df['urban_class'].notna()].groupby(['year', 'urban_class']).size()
             .unstack('urban_class').reindex(index=urby)
             .reindex(columns=['rural', 'town', 'city']).fillna(0).astype(int))

    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6), dpi=140)
    bars = ax[0].bar(ally, per_year.values, color='#3A6EA5', width=0.7)
    ax[0].bar_label(bars, fontsize=8, padding=2)
    ax[0].set_title(f'{COUNTRY} — Xeno-canto recordings per year (2015–2025)',
                    fontsize=12, fontweight='bold')
    ax[0].set_xlabel('year'); ax[0].set_ylabel('number of recordings')
    ax[0].set_xticks(ally); ax[0].tick_params(axis='x', rotation=45)
    ax[0].spines[['top', 'right']].set_visible(False); ax[0].margins(y=0.12)
    ax[0].text(0.5, -0.28, f'total = {len(df):,} recordings', transform=ax[0].transAxes,
               ha='center', fontsize=9, color='#555')

    bottom = np.zeros(len(urby))
    for k in ['rural', 'town', 'city']:
        ax[1].bar(urby, pivot[k].values, bottom=bottom, color=COLORS[k], label=k, width=0.7)
        bottom += pivot[k].values
    totals = pivot.sum(axis=1)
    for x, t in zip(urby, totals):
        ax[1].text(x, t + max(totals) * 0.01, str(int(t)), ha='center', va='bottom', fontsize=8)
    ax[1].set_title(f'{COUNTRY} — city / town / rural split (2015–2022)',
                    fontsize=12, fontweight='bold')
    ax[1].set_xlabel('year'); ax[1].set_ylabel('number of recordings (with coordinates)')
    ax[1].set_xticks(urby); ax[1].tick_params(axis='x', rotation=45)
    ax[1].spines[['top', 'right']].set_visible(False); ax[1].margins(y=0.12)
    ax[1].legend(title='class', frameon=False)
    urb, tot = pivot[['city', 'town']].sum().sum(), pivot.sum().sum()
    ax[1].text(0.5, -0.28, f'{urb:,} urban of {tot:,} classified ({100*urb/tot:.1f}% city+town)',
               transform=ax[1].transAxes, ha='center', fontsize=9, color='#555')

    fig.suptitle(f'{COUNTRY} acoustic recordings — yearly overview',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    out = f'{FIG}/india_recordings_yearwise.png'
    fig.savefig(out, bbox_inches='tight'); plt.close(fig)
    print('wrote', out)


def fig_map(df):
    d = df.dropna(subset=['lat', 'lon', 'urban_class'])
    world = gpd.read_file(NE)
    namecol = [c for c in world.columns if c.upper() in ('ADMIN', 'NAME', 'SOVEREIGNT')][0]
    country = world[world[namecol].str.fullmatch(COUNTRY, case=False, na=False)]

    z = {'rural': 1, 'town': 2, 'city': 3}
    s = {'rural': 6, 'town': 14, 'city': 16}
    fig, ax = plt.subplots(figsize=(9, 10), dpi=150)
    country.plot(ax=ax, color='#F5F5F2', zorder=0)
    country.boundary.plot(ax=ax, color='#333333', linewidth=0.9)
    for k in ['rural', 'town', 'city']:
        g = d[d['urban_class'] == k]
        ax.scatter(g['lon'], g['lat'], s=s[k], c=COLORS[k],
                   alpha=0.65 if k == 'rural' else 0.85, edgecolors='none',
                   zorder=z[k], label=f'{k} ({len(g):,})')
    minx, miny, maxx, maxy = country.total_bounds
    ax.set_xlim(minx - 2, maxx + 2); ax.set_ylim(miny - 1, maxy + 2)
    ax.set_aspect(1 / np.cos(np.radians((miny + maxy) / 2)))
    ax.set_title(f'{COUNTRY} — Xeno-canto recordings by urban class (2015–2022)',
                 fontsize=13, fontweight='bold')
    ax.set_xlabel('longitude'); ax.set_ylabel('latitude')
    ax.legend(title='class (GCTB city/town polygons)', loc='lower left',
              frameon=True, framealpha=0.9)
    ax.grid(True, ls=':', lw=0.4, color='#cccccc')
    fig.text(0.5, 0.055,
             f'{len(d):,} georeferenced recordings, 2015–2022  •  city/town = inside a '
             f'GCTB built-up polygon, rural = outside', ha='center', fontsize=8, color='#555')
    fig.tight_layout()
    out = f'{FIG}/india_recordings_map.png'
    fig.savefig(out, bbox_inches='tight'); plt.close(fig)
    print('wrote', out)


if __name__ == '__main__':
    df = load_country()
    fig_yearwise(df)
    fig_map(df)
