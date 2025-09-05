#!/usr/bin/env python3
# check_duplicates.py (verbose, multi-dir, cross-dir)
#
# Usage:
#   python check_duplicates.py /path/to/dir1 /path/to/dir2 ...
#
# Notes:
# - Non-recursive: only the files directly inside each given directory are checked.
# - Reports duplicates within each directory AND across all directories.
# - Exit codes: 0 = no duplicates, 1 = duplicates found, 2 = error/no files.

from pathlib import Path
import hashlib
import sys
import argparse
from collections import defaultdict

BUF_SIZE = 1024 * 1024  # 1 MiB

def sha256sum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(BUF_SIZE), b""):
            h.update(chunk)
    return h.hexdigest()

def list_h5(dirpath: Path):
    return sorted(p for p in dirpath.iterdir() if p.is_file() and p.suffix == ".h5")

def check_within_dir(dirpath: Path, global_seen: dict) -> tuple[int, list[tuple[str, str]], list[tuple[Path, str]]]:
    """
    Returns:
      - status (0/1/2),
      - list of (a,b) intra-dir duplicate filename pairs,
      - list of (file_path, digest) for all files (to contribute to cross-dir map)
    """
    if not dirpath.is_dir():
        print(f"[ERROR] {dirpath} is not a directory")
        return 2, [], []

    files = list_h5(dirpath)
    print(f"\n[CHECK] {dirpath}")
    print(f"  Found {len(files)} .h5 files")

    if not files:
        print("  (no .h5 files here)")
        return 2, [], []

    local_seen = {}
    local_dups: list[tuple[str, str]] = []
    items: list[tuple[Path, str]] = []

    for f in files:
        print(f"  -> hashing {f.name} ...", end="", flush=True)
        digest = sha256sum(f)
        print(f" done (sha256={digest[:12]}...)")

        items.append((f, digest))
        if digest in local_seen:
            print(f"     !! duplicate WITHIN DIR of {local_seen[digest].name}")
            local_dups.append((local_seen[digest].name, f.name))
        else:
            print(f"     ++ unique so far in this dir")
            local_seen[digest] = f

        # Track in global_seen for optional live global visibility (not required but handy)
        global_seen.setdefault(digest, []).append(f)

    if local_dups:
        print(f"[RESULT] {len(local_dups)} intra-dir duplicates in {dirpath}:")
        for a, b in local_dups:
            print(f"  - {a} == {b}")
        return 1, local_dups, items
    else:
        print("[RESULT] no duplicates within this directory")
        return 0, [], items

def main():
    ap = argparse.ArgumentParser(
        description="Check for duplicate .h5 files by content within and across provided directories (non-recursive)."
    )
    ap.add_argument(
        "dirs",
        nargs="+",
        help="One or more directories to check.",
    )
    args = ap.parse_args()

    overall_exit = 0
    any_files = False

    # For cross-dir aggregation
    global_seen: dict[str, list[Path]] = {}

    # Per-dir checks (and fill global map)
    for d in args.dirs:
        dirpath = Path(d).resolve()
        status, _, items = check_within_dir(dirpath, global_seen)
        overall_exit = max(overall_exit, 1 if status == 1 else status)
        if items:
            any_files = True

    if not any_files:
        print("\n[GLOBAL] No .h5 files found in any provided directory.")
        sys.exit(2)

    # Cross-dir check (groups of identical files spanning directories)
    print("\n[GLOBAL] Cross-directory duplicate analysis")
    cross_dups = []
    for digest, paths in global_seen.items():
        if len(paths) > 1:
            # If there is more than one file with the same digest, consider it a global duplicate group
            cross_dups.append((digest, paths))

    if cross_dups:
        print(f"  Found {len(cross_dups)} duplicate content group(s) across the provided directories:")
        for digest, paths in cross_dups:
            print(f"  - sha256={digest} present in {len(paths)} files:")
            for p in paths:
                print(f"      * {p}")
        overall_exit = 1
    else:
        print("  No cross-directory duplicates detected.")

    # Final status
    if overall_exit == 0:
        print("\n[SUMMARY] no duplicates")
    elif overall_exit == 1:
        print("\n[SUMMARY] duplicates found")
    else:
        print("\n[SUMMARY] finished with errors or no files")

    sys.exit(overall_exit)

if __name__ == "__main__":
    main()
