"""Left-join Xeno-canto metadata onto each score_<continent>.csv by 'id'.

Loads the flattened metadata CSV once, then for each continent writes
score_<continent>_meta.csv = score columns + metadata columns (minus the id key).
"""
import sys
import csv

csv.field_size_limit(10 ** 7)

META = 'metadata_since_2023-01-01.csv'
CONTINENTS = ['africa', 'america', 'asia', 'australia']


def load_meta(path):
    with open(path, newline='') as f:
        r = csv.reader(f)
        hdr = next(r)
        idx = hdr.index('id')
        cols = [c for i, c in enumerate(hdr) if i != idx]
        meta = {row[idx]: [v for i, v in enumerate(row) if i != idx] for row in r}
    return meta, cols


def merge(continent, meta, mcols):
    src = f'score_{continent}.csv'
    out = f'score_{continent}_meta.csv'
    blank = [''] * len(mcols)
    matched = total = 0
    with open(src, newline='') as f, open(out, 'w', newline='') as fo:
        r = csv.reader(f)
        w = csv.writer(fo)
        shdr = next(r)
        sid = shdr.index('id')
        w.writerow(shdr + mcols)
        for row in r:
            total += 1
            m = meta.get(row[sid])
            if m is not None:
                matched += 1
            w.writerow(row + (m if m is not None else blank))
    pct = 100 * matched / total if total else 0
    print(f'{continent:10} {total:>7} rows  {matched:>7} matched ({pct:5.1f}%)  unmatched={total - matched}  -> {out}')


if __name__ == '__main__':
    meta, mcols = load_meta(META)
    print(f'loaded {len(meta)} metadata records, {len(mcols)} metadata columns\n')
    for c in (sys.argv[1:] or CONTINENTS):
        merge(c, meta, mcols)
