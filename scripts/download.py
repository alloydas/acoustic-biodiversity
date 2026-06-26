"""Download Xeno-canto recordings for one continent.

Replaces the five copy-pasted per-continent scripts (africa.py ... europe.py),
three of which had been overwritten with Europe settings (area/outdir="europe").
Continent is now an argument, so that class of copy-paste bug can't recur.

Reproduce the existing dataset (all recordings per continent since 2023-01-01,
all quality grades, all years 2023-2025):

    python download.py --continent africa
    python download.py --continent america
    python download.py --continent asia
    python download.py --continent australia
    python download.py --continent europe

Optional filters (omit to match the existing all-grades / all-years data):
    --quality A          restrict to a single quality grade (A best ... E worst)
    --year 2024          restrict to one year
    --months 4-6         restrict to a month range (inclusive)
"""
import os
import argparse

from xenopy import Query

API_KEY = "174b9af76767d05178f2c9ae7aae7010c8121c94"
# Xeno-canto 'area' values; continent name maps 1:1.
AREAS = {"africa", "america", "asia", "australia", "europe"}


def parse_args():
    p = argparse.ArgumentParser(description="Download Xeno-canto recordings for a continent")
    p.add_argument("--continent", required=True, choices=sorted(AREAS),
                   help="continent / Xeno-canto area to download")
    p.add_argument("--outdir", default=None, help="output dir (default: the continent name)")
    p.add_argument("--since", default="2023-01-01", help="earliest upload/record date (default 2023-01-01)")
    p.add_argument("--quality", default=None, help="single quality grade A-E (default: all grades)")
    p.add_argument("--year", default=None, help="restrict to one year (default: all years)")
    p.add_argument("--months", default=None, help="month range 'M' or 'M1-M2' inclusive (default: all)")
    p.add_argument("--nproc", type=int, default=72)
    p.add_argument("--attempts", type=int, default=100)
    return p.parse_args()


def month_range(spec):
    if spec is None:
        return [None]
    if "-" in spec:
        lo, hi = (int(x) for x in spec.split("-"))
    else:
        lo = hi = int(spec)
    return list(range(lo, hi + 1))


def main():
    args = parse_args()
    outdir = args.outdir or args.continent
    os.makedirs(outdir, exist_ok=True)

    for month in month_range(args.months):
        # Only pass filters that are set, so an omitted filter means "no restriction".
        opts = {"api_key": API_KEY, "since": args.since, "area": args.continent}
        if args.quality:
            opts["q"] = args.quality
        if args.year:
            opts["year"] = args.year
        if month is not None:
            opts["month"] = f"{month:02d}"

        label = f"{args.continent}" + (f" month {month:02d}" if month is not None else "")
        print(f"Downloading recordings: {label}  (since {args.since}) -> {outdir}/")
        try:
            q = Query(**opts)
            q.retrieve_recordings(multiprocess=True, nproc=args.nproc,
                                  attempts=args.attempts, outdir=outdir)
        except Exception as e:
            # Log instead of silently swallowing, so failures are visible in job logs.
            print(f"  WARNING: download failed for {label}: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
