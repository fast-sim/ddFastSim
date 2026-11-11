#!/usr/bin/env python3
"""Check saved Geant4 random-engine states for duplicates across batch jobs.

Expected layout (any depth below the paths you pass):
    <DETECTOR>/seeds/<CLUSTER_ID>.<RUN_ID>/currentRun.rndm
'randomSeeds' is accepted as the directory name too (as written by the
ddfastsim v1.0.0 condor wrappers).

currentRun.rndm holds the engine state at the start of each job's run, so two
jobs with an identical state started from the same seed and would have
produced identical showers (given the same configuration).

Usage:
    python3 check_duplicate_seeds.py /eos/.../ddfastsim            # all detectors below
    python3 check_duplicate_seeds.py Par04/seeds ODD/seeds --csv seeds.csv
Exit code: 0 = no duplicates, 1 = duplicates found, 2 = nothing usable found.
"""
import argparse, csv, hashlib, re, sys
from collections import defaultdict
from pathlib import Path

SEED_DIRS = {"seeds", "randomSeeds"}
JOB_RE = re.compile(r"^(\d+)\.(\d+)$")
MIXMAX_RE = re.compile(r"V\[N\]\s*=\s*\{([^}]*)\}\s*;\s*counter\s*=\s*(\d+)", re.S)


def detector_of(path: Path) -> str:
    """Label = path part just above the seeds directory (e.g. 'ODD')."""
    parts = path.parts
    for i in range(len(parts) - 1, -1, -1):
        if parts[i] in SEED_DIRS:
            # EOS layout DETECTOR/<dataset...>/randomSeeds: keep both levels
            if i >= 2 and parts[i - 1].startswith("dataset"):
                return f"{parts[i - 2]}/{parts[i - 1]}"
            return parts[i - 1] if i >= 1 else "?"
    return path.parent.parent.name or "?"


def state_key(text: str):
    """Return (kind, key) identifying the engine state in the file."""
    m = MIXMAX_RE.search(text)
    if m:
        values = tuple(int(v) for v in re.findall(r"-?\d+", m.group(1)))
        return "mixmax", hashlib.sha256(repr((values, int(m.group(2)))).encode()).hexdigest()
    # other engines: compare whitespace-normalised content
    norm = " ".join(text.split())
    return "raw", hashlib.sha256(norm.encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", type=Path, help="directories to search (recursively)")
    ap.add_argument("--name", default="currentRun.rndm", help="state file name (default: currentRun.rndm)")
    ap.add_argument("--csv", help="optional per-file CSV output")
    args = ap.parse_args()

    files = sorted({f for p in args.paths for f in (p.rglob(args.name) if p.is_dir() else [p])})
    if not files:
        print(f"No {args.name} files found below: {', '.join(map(str, args.paths))}")
        sys.exit(2)

    by_key = defaultdict(list)
    per_det = defaultdict(int)
    problems, rows, kinds = [], [], defaultdict(int)
    for f in files:
        det, job = detector_of(f), f.parent.name
        if not JOB_RE.match(job):
            problems.append(f"job directory not <CLUSTER>.<RUN>: {f}")
        try:
            text = f.read_text(errors="replace")
        except OSError as e:
            problems.append(f"unreadable: {f} ({e})")
            continue
        if not text.strip():
            problems.append(f"empty: {f}")
            continue
        kind, key = state_key(text)
        kinds[kind] += 1
        per_det[det] += 1
        by_key[key].append((det, job, f))
        rows.append({"detector": det, "job": job, "engine_format": kind, "state_sha256": key, "file": str(f)})

    print(f"Checked {sum(per_det.values())} state files")
    for det in sorted(per_det):
        print(f"  {det:30s} {per_det[det]:6d} jobs")
    print(f"Engine formats: {dict(kinds)}  (MixMax expected)")

    dups = {k: v for k, v in by_key.items() if len(v) > 1}
    if dups:
        within = sum(1 for v in dups.values() if len({d for d, _, _ in v}) == 1)
        print(f"\nDUPLICATES: {len(dups)} identical state(s) "
              f"({within} within one detector, {len(dups) - within} across detectors)")
        for k, v in dups.items():
            print(f"  state {k[:12]}… in {len(v)} jobs:")
            for det, job, f in v:
                print(f"     {det:24s} {job:16s} {f}")
    else:
        print("\nNo duplicate states: every job started from a different random-engine state.")

    if problems:
        print(f"\n{len(problems)} file(s) needing attention:")
        for p in problems[:50]:
            print("  " + p)
        if len(problems) > 50:
            print(f"  ... and {len(problems) - 50} more")

    if args.csv:
        for r in rows:
            r["duplicate_group_size"] = len(by_key[r["state_sha256"]])
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\nPer-file table written to {args.csv}")

    sys.exit(1 if dups else 0)


if __name__ == "__main__":
    main()
