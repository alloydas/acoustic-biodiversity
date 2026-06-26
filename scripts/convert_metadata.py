"""Flatten Xeno-canto metadata JSON pages into a single CSV.

Reads every pageN.json under a metadata folder (each holds a 'recordings' list
from the Xeno-canto API) and writes one row per recording. The nested image-URL
dicts ('sono', 'osci') are dropped; the 'also' list of co-occurring species is
joined with ';'. Usage: python convert_metadata.py <metadata_dir> <out.csv>
"""
import os
import sys
import csv
import glob
import json

# Scalar fields to keep, in output order (from the Xeno-canto v3 record schema).
FIELDS = [
    'id', 'gen', 'sp', 'ssp', 'grp', 'en', 'rec', 'cnt', 'loc', 'lat', 'lon',
    'alt', 'type', 'sex', 'stage', 'method', 'url', 'file', 'file-name', 'lic',
    'q', 'length', 'time', 'date', 'uploaded', 'also', 'rmk', 'animal-seen',
    'playback-used', 'temp', 'regnr', 'auto', 'dvc', 'mic', 'smp',
]


def main(meta_dir, out_path):
    pages = sorted(glob.glob(os.path.join(meta_dir, 'page*.json')),
                   key=lambda p: int(''.join(filter(str.isdigit, os.path.basename(p)))))
    print(f"Found {len(pages)} page files")

    n = 0
    with open(out_path, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=FIELDS, extrasaction='ignore', restval='')
        writer.writeheader()
        for i, page in enumerate(pages, 1):
            try:
                data = json.load(open(page))
            except (json.JSONDecodeError, IOError) as e:
                print(f"  skip {os.path.basename(page)}: {e}")
                continue
            for rec in data.get('recordings', []):
                row = dict(rec)
                also = row.get('also')
                if isinstance(also, list):
                    row['also'] = ';'.join(also)
                writer.writerow(row)
                n += 1
            if i % 200 == 0:
                print(f"  {i}/{len(pages)} pages, {n} rows")

    print(f"Wrote {n} rows to {out_path}")


if __name__ == '__main__':
    meta_dir = sys.argv[1] if len(sys.argv) > 1 else 'dataset/metadata/since:2023-01-01'
    out_path = sys.argv[2] if len(sys.argv) > 2 else 'metadata_since_2023-01-01.csv'
    main(meta_dir, out_path)
