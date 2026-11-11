#!/usr/bin/env python3
"""Summarise validate_lemurs.py reports for Technical Validation.

Usage:  python3 summarize_qc.py [directory with lemurs_qc*report*.json]   (default: .)
Prints every file with a flag other than the expected count deficit, then
per-detector totals: files, stored showers, errors, non-finite or negative
voxels, and showers with zero total energy.
"""
import glob, json, os, re, sys
from collections import defaultdict

DETS = ("Par04SiW", "Par04SciPb", "Par04PbWO4", "ODD", "FCCeeCLD", "FCCeeALLEGRO")


def detector(name):
    flat = name.replace("_", "")
    for d in DETS:
        if d in flat:
            return d
    if "Par04" in flat:          # production files of Par04 SiW: ddsim_mesh_Par04_gamma_...
        return "Par04SiW"
    return "?"


folder = sys.argv[1] if len(sys.argv) > 1 else "."
reports = sorted(glob.glob(os.path.join(folder, "lemurs_qc*report*.json")))
if not reports:
    sys.exit(f"no lemurs_qc*report*.json in {folder}")

# Same file name seen more than once (e.g. float64 originals and the released float32
# copies in float32/): keep the record with the fewest errors, i.e. the released copy.
best = {}
for rep in reports:
    for f in json.load(open(rep))["files"]:
        name = os.path.basename(f["path"])
        if name not in best or len(f.get("errors", [])) < len(best[name].get("errors", [])):
            best[name] = f
seen, tot = set(best), defaultdict(lambda: defaultdict(int))
print("Files with flags other than the expected count deficit:")
for name, f in sorted(best.items()):
    if True:
        det = detector(name)
        s = f.get("showers_summary", {})
        t = tot[det]
        t["files"] += 1
        t["showers"] += f.get("event_count", 0) or 0
        t["errors"] += len(f.get("errors", []))
        t["nonfinite"] += s.get("events_with_nonfinite_voxels", 0) or 0
        t["negative"] += s.get("events_with_finite_negative_voxels", 0) or 0
        t["zero"] += s.get("zero_total_energy_events", 0) or 0
        t["testing"] += "1000events" in name
        extra = [w["code"] for w in f.get("warnings", []) if w.get("code") != "nominal_count_differs"]
        if extra or f.get("errors"):
            print(f"  {name[:72]:72s} {extra} errors={len(f.get('errors', []))} "
                  f"zero_energy={s.get('zero_total_energy_events')}")

print("\nPer detector:")
print(f"  {'detector':14s} {'files':>5s} {'(test)':>6s} {'showers':>9s} {'errors':>6s} {'nonfin':>6s} {'neg':>5s} {'zero-E':>6s}")
for det in sorted(tot):
    t = tot[det]
    print(f"  {det:14s} {t['files']:5d} {t['testing']:6d} {t['showers']:9d} {t['errors']:6d} "
          f"{t['nonfinite']:6d} {t['negative']:5d} {t['zero']:6d}")
print(f"\nReports read: {len(reports)}; distinct files: {len(seen)} (LEMURS has 90: 5 detectors x (10 training + 8 testing))")
