#!/usr/bin/env python3
"""Read-only integrity audit of extracted LEMURS HDF5 files.

Usage:
    python3 validate_lemurs.py /path/to/LEMURS --output lemurs_qc_report.json
    python3 validate_lemurs.py file1.h5 file2.h5 --output report.json

Dependencies: numpy, h5py. Directories are searched recursively for .h5/.hdf5.
Archives are not extracted. Every selected file is scanned, in event batches.
Only the JSON report is written; input HDF5 files are opened in mode 'r'.

Schema checked: showers (S, 9, 16, 45), and incident_energy/incident_theta/
incident_phi each (S,), with float32 values. Energies are interpreted as MeV
and angles as radians, following the intended LEMURS schema. Those units
cannot be established from numerical values alone.

Checks: keys, dataset types, dimensions, matching event counts, float32,
readability, NaN/infinity, finite negative energies, nonpositive incident
energy, theta in [0, pi], phi in [-pi, pi], and zero-total-energy showers.
A small angular tolerance accommodates float32 roundoff. Zero-energy showers
are reported for review, not automatically classified as data errors.

Observed ranges, storage/compression metadata, event counts and total-energy
ratios are reported. Ratios above 1 + tolerance are diagnostics for review,
not a universal physics rejection rule. Counts in filenames are reported as
nominal labels, not proof of how many events were actually simulated.

Not checked: exact duplicates, random-seed independence, production history,
geometry/material correctness, distributional agreement, model performance,
full-release completeness, or exact source-energy/acceptance settings. This
script does not generate showers, re-filter data or change the data licence.

This is a newly recreated revision utility, not the historical ddfastsim
plotting script and not evidence that an audit was performed in production.
Exit status: 0 = completed without errors (review warnings may remain);
1 = one or more errors; 2 = invalid invocation/no files; 130 = interrupted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import sys
import tempfile
from datetime import datetime, timezone
from typing import Any

try:
    import h5py
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "Missing dependency: " + str(exc)
        + "\nInstall in your Python environment: python3 -m pip install numpy h5py"
    )

VERSION = "1.0.0-recreated"
KEYS = ("showers", "incident_energy", "incident_theta", "incident_phi")
SHOWER_SHAPE = (9, 16, 45)
UNITS = {"showers": "MeV", "incident_energy": "MeV",
         "incident_theta": "rad", "incident_phi": "rad"}


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def issue(result: dict[str, Any], severity: str, code: str, message: str) -> None:
    result[severity].append({"code": code, "message": message})


def stats_new() -> dict[str, Any]:
    return {"values_read": 0, "finite_count": 0, "nan_count": 0,
            "positive_infinity_count": 0, "negative_infinity_count": 0,
            "finite_negative_count": 0, "zero_count": 0,
            "finite_min": None, "finite_max": None}


def stats_add(stats: dict[str, Any], values: np.ndarray) -> None:
    a = np.asarray(values)
    finite = np.isfinite(a)
    stats["values_read"] += int(a.size)
    stats["finite_count"] += int(np.count_nonzero(finite))
    stats["nan_count"] += int(np.count_nonzero(np.isnan(a)))
    stats["positive_infinity_count"] += int(np.count_nonzero(np.isposinf(a)))
    stats["negative_infinity_count"] += int(np.count_nonzero(np.isneginf(a)))
    stats["finite_negative_count"] += int(np.count_nonzero(finite & (a < 0)))
    stats["zero_count"] += int(np.count_nonzero(a == 0))
    if finite.any():
        lo, hi = float(np.min(a[finite])), float(np.max(a[finite]))
        stats["finite_min"] = lo if stats["finite_min"] is None else min(lo, stats["finite_min"])
        stats["finite_max"] = hi if stats["finite_max"] is None else max(hi, stats["finite_max"])


def nominal_count(filename: str) -> int | None:
    match = re.search(r"_(\d+)(k?)events(?:_|\.)", filename, re.IGNORECASE)
    return None if match is None else int(match[1]) * (1000 if match[2] else 1)


def audit_file(path: Path, batch: int, angle_tol: float, ratio_tol: float) -> dict[str, Any]:
    result: dict[str, Any] = {"path": str(path), "status": "pending", "errors": [],
                              "warnings": [], "datasets": {}, "event_count": None,
                              "nominal_count_from_filename": nominal_count(path.name)}
    try:
        st = path.stat()
        result["size_bytes"] = st.st_size
        result["mtime_ns"] = st.st_mtime_ns
        with h5py.File(path, "r") as f:
            result["root_keys"] = sorted(f.keys())
            datasets: dict[str, Any] = {}
            for key in KEYS:
                if key not in f:
                    issue(result, "errors", "missing_key", f"Missing dataset: {key}")
                    continue
                ds = f[key]
                if not isinstance(ds, h5py.Dataset):
                    issue(result, "errors", "not_a_dataset", f"{key} is not an HDF5 dataset")
                    continue
                info: dict[str, Any] = {
                    "shape": list(ds.shape), "dtype": str(ds.dtype),
                    "interpreted_unit": UNITS[key], "chunks": list(ds.chunks) if ds.chunks else None,
                    "compression": ds.compression, "compression_options": ds.compression_opts,
                    "scan_complete": False, "statistics": stats_new(),
                }
                result["datasets"][key] = info
                expected_rank = 4 if key == "showers" else 1
                if ds.ndim != expected_rank or (key == "showers" and ds.shape[1:] != SHOWER_SHAPE):
                    expected = "(S, 9, 16, 45)" if key == "showers" else "(S,)"
                    issue(result, "errors", "unexpected_shape", f"{key}: {ds.shape}, expected {expected}")
                if ds.dtype.kind != "f" or ds.dtype.itemsize != 4:
                    issue(result, "errors", "unexpected_dtype", f"{key}: {ds.dtype}, expected float32")
                if ds.dtype.kind not in "fiu":
                    issue(result, "errors", "unsupported_numeric_type", f"Cannot numerically scan {key}: {ds.dtype}")
                    continue
                datasets[key] = ds

            showers = datasets.get("showers")
            n = int(showers.shape[0]) if showers is not None and showers.ndim >= 1 else None
            result["event_count"] = n
            if n is not None:
                if n == 0:
                    issue(result, "warnings", "empty_file", "The shower dataset contains zero events")
                for key, ds in datasets.items():
                    if ds.ndim >= 1 and ds.shape[0] != n:
                        issue(result, "errors", "event_count_mismatch", f"{key}: {ds.shape[0]} records; showers: {n}")
                nominal = result["nominal_count_from_filename"]
                if nominal is not None:
                    result["nominal_minus_stored"] = nominal - n
                    if nominal != n:
                        issue(result, "warnings", "nominal_count_differs",
                              f"Filename label: {nominal}; stored: {n}. Cause not inferred by this audit.")

            # Scan incident arrays separately so even length mismatches are reported fully.
            for key in KEYS[1:]:
                if key not in datasets:
                    continue
                ds = datasets[key]
                info = result["datasets"][key]
                out_of_range = 0
                for start in range(0, ds.shape[0] if ds.ndim else 1, batch):
                    a = ds[start:start + batch] if ds.ndim else ds[()]
                    stats_add(info["statistics"], np.asarray(a))
                    finite = np.isfinite(a)
                    if key == "incident_theta":
                        out_of_range += int(np.count_nonzero(finite & ((a < -angle_tol) | (a > math.pi + angle_tol))))
                    elif key == "incident_phi":
                        out_of_range += int(np.count_nonzero(finite & ((a < -math.pi - angle_tol) | (a > math.pi + angle_tol))))
                info["scan_complete"] = True
                if key != "incident_energy":
                    info["values_outside_physical_angle_range"] = out_of_range
                    if out_of_range:
                        issue(result, "errors", "angle_out_of_range", f"{key}: {out_of_range} values outside the specified angle convention")

            if showers is not None:
                info = result["datasets"]["showers"]
                event_summary: dict[str, Any] = {
                    "events_with_nonfinite_voxels": 0, "events_with_finite_negative_voxels": 0,
                    "finite_nonnegative_events": 0, "zero_total_energy_events": 0,
                    "total_energy_MeV": stats_new(), "nonzero_voxel_counts": stats_new(),
                    "deposited_to_incident_ratio": stats_new(),
                    "ratio_events_excluded": 0, "ratio_above_one_plus_tolerance": 0,
                    "derived_nonfinite_totals": 0, "derived_nonfinite_ratios": 0,
                }
                correct_shape = showers.ndim == 4 and showers.shape[1:] == SHOWER_SHAPE
                energy = datasets.get("incident_energy")
                ratio_possible = (correct_shape and energy is not None and energy.ndim == 1 and energy.shape[0] == n)
                event_summary["ratio_computed"] = bool(ratio_possible)
                for start in range(0, showers.shape[0] if showers.ndim else 1, batch):
                    a = showers[start:start + batch] if showers.ndim else showers[()]
                    stats_add(info["statistics"], np.asarray(a))
                    if not correct_shape:
                        continue
                    flat = a.reshape(a.shape[0], -1)
                    row_finite = np.all(np.isfinite(flat), axis=1)
                    row_negative = np.any(np.isfinite(flat) & (flat < 0), axis=1)
                    valid = row_finite & ~row_negative
                    event_summary["events_with_nonfinite_voxels"] += int(np.count_nonzero(~row_finite))
                    event_summary["events_with_finite_negative_voxels"] += int(np.count_nonzero(row_negative))
                    event_summary["finite_nonnegative_events"] += int(np.count_nonzero(valid))
                    with np.errstate(over="ignore", invalid="ignore"):
                        total = np.sum(flat, axis=1, dtype=np.float64)
                    event_summary["derived_nonfinite_totals"] += int(np.count_nonzero(valid & ~np.isfinite(total)))
                    valid_totals = valid & np.isfinite(total)
                    event_summary["zero_total_energy_events"] += int(np.count_nonzero(valid_totals & (total == 0)))
                    stats_add(event_summary["total_energy_MeV"], total[valid_totals])
                    stats_add(event_summary["nonzero_voxel_counts"], np.count_nonzero(flat[valid], axis=1))
                    if ratio_possible:
                        incident = np.asarray(energy[start:start + a.shape[0]], dtype=np.float64)
                        eligible = valid_totals & np.isfinite(incident) & (incident > 0)
                        event_summary["ratio_events_excluded"] += int(np.count_nonzero(~eligible))
                        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                            ratio = total[eligible] / incident[eligible]
                        event_summary["derived_nonfinite_ratios"] += int(np.count_nonzero(~np.isfinite(ratio)))
                        stats_add(event_summary["deposited_to_incident_ratio"], ratio)
                        event_summary["ratio_above_one_plus_tolerance"] += int(np.count_nonzero(np.isfinite(ratio) & (ratio > 1 + ratio_tol)))
                info["scan_complete"] = True
                if correct_shape:
                    result["showers_summary"] = event_summary
                    if event_summary["zero_total_energy_events"]:
                        issue(result, "warnings", "zero_total_energy_showers",
                              f"{event_summary['zero_total_energy_events']} finite nonnegative showers have zero total energy. Review their treatment; no records were removed.")
                    if event_summary["ratio_above_one_plus_tolerance"]:
                        issue(result, "warnings", "ratio_above_one",
                              f"{event_summary['ratio_above_one_plus_tolerance']} ratios exceed 1 + {ratio_tol}; diagnostic only, not a physics pass/fail rule.")
                    if event_summary["derived_nonfinite_totals"] or event_summary["derived_nonfinite_ratios"]:
                        issue(result, "errors", "derived_nonfinite", "Nonfinite values arose while computing totals or ratios")

            for key, info in result["datasets"].items():
                s = info["statistics"]
                nonfinite = s["nan_count"] + s["positive_infinity_count"] + s["negative_infinity_count"]
                if nonfinite:
                    issue(result, "errors", "nonfinite_values", f"{key}: {nonfinite} NaN/infinite values")
                if key in ("showers", "incident_energy") and s["finite_negative_count"]:
                    issue(result, "errors", "negative_energy", f"{key}: {s['finite_negative_count']} finite negative energies")
                if key == "incident_energy" and s["zero_count"]:
                    issue(result, "errors", "zero_incident_energy", f"incident_energy: {s['zero_count']} zeros")
        after = path.stat()
        if after.st_size != st.st_size or after.st_mtime_ns != st.st_mtime_ns:
            issue(result, "warnings", "input_changed_during_scan", "File size/mtime changed during the audit; rerun on a stable copy")
    except Exception as exc:
        # Keep the report for other files and mark any incomplete scan explicitly.
        issue(result, "errors", "read_error", f"{type(exc).__name__}: {exc}")
    result["complete_scan"] = all(
        key in result["datasets"] and result["datasets"][key]["scan_complete"] for key in KEYS
    )
    result["status"] = "error" if result["errors"] else ("review" if result["warnings"] else "ok")
    return result


def discover(inputs: list[str]) -> list[Path]:
    found: set[Path] = set()
    for value in inputs:
        path = Path(value).expanduser()
        if not path.exists():
            raise ValueError(f"Input does not exist: {path}")
        if path.is_dir():
            found.update(p.resolve() for p in path.rglob("*") if p.is_file() and p.suffix.lower() in (".h5", ".hdf5"))
        elif path.suffix.lower() in (".h5", ".hdf5"):
            found.add(path.resolve())
        else:
            raise ValueError(f"Expected an extracted HDF5 file or directory, not: {path}")
    if not found:
        raise ValueError("No .h5/.hdf5 files found. Extract the testing archive first; ZIP/TAR files are not scanned.")
    return sorted(found)


def write_report(path: Path, report: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp: str | None = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, suffix=".tmp", delete=False) as stream:
            temp = stream.name
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write("\n")
        os.replace(temp, path)
    finally:
        if temp is not None and os.path.exists(temp):
            os.unlink(temp)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", nargs="+", help="Extracted HDF5 file(s) or directory/directories")
    parser.add_argument("--output", type=Path, default=Path("lemurs_qc_report.json"))
    parser.add_argument("--chunk-events", type=int, default=256, help="Events read per batch (default: 256)")
    parser.add_argument("--angle-tolerance", type=float, default=1e-6, help="Absolute tolerance in radians (default: 1e-6)")
    parser.add_argument("--ratio-tolerance", type=float, default=1e-5, help="Flag ratios > 1 + tolerance for review (default: 1e-5)")
    args = parser.parse_args()
    if args.chunk_events < 1:
        parser.error("--chunk-events must be positive")
    if any(not math.isfinite(x) or x < 0 for x in (args.angle_tolerance, args.ratio_tolerance)):
        parser.error("Tolerances must be finite and nonnegative")
    try:
        paths = discover(args.inputs)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    output = args.output.expanduser().resolve()
    if output.suffix.lower() != ".json":
        parser.error("--output must end in .json to avoid overwriting a data file")
    if output in paths:
        parser.error("The report must not overwrite an input file")
    report: dict[str, Any] = {
        "audit_version": VERSION, "started_utc": utcnow(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "h5py": h5py.__version__, "hdf5": h5py.version.hdf5_version},
        "parameters": {"chunk_events": args.chunk_events, "angle_tolerance_rad": args.angle_tolerance,
                       "ratio_review_tolerance": args.ratio_tolerance},
        "scope": "All HDF5 files discovered under the supplied inputs, not an independently verified release inventory",
        "units_assumed_from_schema": UNITS,
        "not_checked": ["exact duplicates", "seed independence or replay", "production history",
                        "physics accuracy", "geometry or materials", "units independently of documentation",
                        "generation-specific ranges or distributions", "full-release completeness"],
        "files": [],
    }
    interrupted = False
    try:
        for index, path in enumerate(paths, 1):
            print(f"[{index}/{len(paths)}] {path.name}", flush=True)
            entry = audit_file(path, args.chunk_events, args.angle_tolerance, args.ratio_tolerance)
            report["files"].append(entry)
            print(f"  {entry['status'].upper()}: {entry['event_count']} events; "
                  f"{len(entry['errors'])} errors; {len(entry['warnings'])} review flags", flush=True)
    except KeyboardInterrupt:
        interrupted = True
    entries = report["files"]
    report["finished_utc"] = utcnow()
    report["summary"] = {
        "files_discovered": len(paths), "files_attempted": len(entries),
        "files_with_complete_scans": sum(r["complete_scan"] for r in entries),
        "files_ok": sum(r["status"] == "ok" for r in entries),
        "files_needing_review": sum(r["status"] == "review" for r in entries),
        "files_with_errors": sum(r["status"] == "error" for r in entries),
        "stored_events_in_attempted_files": sum(r["event_count"] or 0 for r in entries),
        "interrupted": interrupted,
        "discovery_list": [str(p) for p in paths],
    }
    try:
        write_report(output, report)
    except (OSError, ValueError) as exc:
        print(f"Could not write report: {exc}", file=sys.stderr)
        return 1
    print(f"\nReport: {output}")
    print(json.dumps({k: v for k, v in report['summary'].items() if k != 'discovery_list'}, indent=2))
    return 130 if interrupted else (1 if report["summary"]["files_with_errors"] else 0)


if __name__ == "__main__":
    raise SystemExit(main())
