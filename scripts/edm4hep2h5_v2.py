#!/usr/bin/env python3
"""Convert existing LEMURS EDM4hep ROOT showers with an RDataFrame filter."""

import argparse
import math
import os
from pathlib import Path

import h5py
import numpy as np


GEOMETRY_COMMIT = "a25ddc8d9bcfd5610e389cda487d3f00172e9443"
GEOMETRIES = {
    "CLD": dict(name="FCCee_o2_v04_noTracker", sides=12,
                barrel_r=2150., barrel_outer_r=2352., barrel_half_z=2210.,
                endcap_z=2307., endcap_back_z=2509., endcap_inner_r=340.,
                endcap_outer_r=2350.),
    "ALLEGRO": dict(name="ALLEGRO_o1_v02_noTracker", sides=0,
                    barrel_r=2149., barrel_outer_r=2721., barrel_half_z=3100.,
                    endcap_z=3260., endcap_back_z=3710., endcap_inner_r=410.,
                    endcap_outer_r=2750.),
}

_NORMALS = [(math.cos(k * math.pi / 6), math.sin(k * math.pi / 6))
            for k in range(12)]


def support(x, y, geometry):
    if geometry["sides"] == 0:
        return math.hypot(x, y)
    return max(x * nx + y * ny for nx, ny in _NORMALS)


CPP_FILTER = r"""
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>

namespace lemurs_v2_rdf {
constexpr double pi = 3.14159265358979323846;
constexpr double c_mm_per_ns = 299.792458;

inline double support(double x, double y, int sides) {
    if (sides == 0) return std::hypot(x, y);
    double result = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < sides; ++k) {
        const double angle = 2. * pi * k / sides;
        result = std::max(result, x * std::cos(angle) + y * std::sin(angle));
    }
    return result;
}

template <class Vector3>
bool finite3(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
           std::isfinite(value.z);
}

inline bool upstream(double x, double y, double z, int sides,
                     double barrel_r, double endcap_z, double margin) {
    return std::isfinite(x) && std::isfinite(y) && std::isfinite(z) &&
           support(x, y, sides) < barrel_r - margin &&
           std::abs(z) < endcap_z - margin;
}

inline double entry_distance(double px, double py, double pz, int sides,
                             double barrel_r, double barrel_half_z,
                             double endcap_z, double endcap_inner_r,
                             double endcap_outer_r) {
    const double norm = std::sqrt(px*px + py*py + pz*pz);
    if (!(norm > 0.) || !std::isfinite(norm))
        return std::numeric_limits<double>::infinity();
    const double dx = px / norm, dy = py / norm, dz = pz / norm;
    double first = std::numeric_limits<double>::infinity();
    const double transverse = support(dx, dy, sides);
    if (transverse > 0.) {
        const double distance = barrel_r / transverse;
        if (std::abs(distance * dz) <= barrel_half_z) first = distance;
    }
    if (std::abs(dz) > 1.e-12) {
        const double distance = endcap_z / std::abs(dz);
        const double radius = support(distance * dx, distance * dy, sides);
        if (endcap_inner_r <= radius && radius <= endcap_outer_r)
            first = std::min(first, distance);
    }
    return first;
}

template <class Particles>
bool keep_event(const Particles& particles, int sides, double barrel_r,
                double barrel_half_z, double endcap_z,
                double endcap_inner_r, double endcap_outer_r,
                double boundary_margin, double time_margin) {
    if (particles.empty()) return false;
    const auto& primary = particles[0];
    if (primary.PDG != 22 || primary.generatorStatus != 1 ||
        !finite3(primary.vertex) || !finite3(primary.endpoint) ||
        !finite3(primary.momentum) || !std::isfinite(primary.time) ||
        std::abs(primary.vertex.x) > 1. ||
        std::abs(primary.vertex.y) > 1. ||
        std::abs(primary.vertex.z) > 1.) return false;

    for (std::size_t i = 1; i < particles.size(); ++i)
        if (particles[i].generatorStatus != 0) return false;

    const double distance = entry_distance(
        primary.momentum.x, primary.momentum.y, primary.momentum.z, sides,
        barrel_r, barrel_half_z, endcap_z, endcap_inner_r, endcap_outer_r);
    if (!std::isfinite(distance)) return false;
    if (upstream(primary.endpoint.x, primary.endpoint.y, primary.endpoint.z,
                 sides, barrel_r, endcap_z, boundary_margin)) return false;

    const double entry_time = primary.time + distance / c_mm_per_ns;
    for (std::size_t i = 1; i < particles.size(); ++i) {
        const auto& secondary = particles[i];
        if (!finite3(secondary.vertex) || !std::isfinite(secondary.time))
            return false;
        if (secondary.time < entry_time - time_margin &&
            upstream(secondary.vertex.x, secondary.vertex.y,
                     secondary.vertex.z, sides, barrel_r, endcap_z,
                     boundary_margin)) return false;
    }
    return true;
}
}
"""

COLUMNS = ("energyCell", "rhoCell", "phiCell", "zCell",
           "energyMC", "thetaMC", "phiMC")


def selection_expression(args):
    geometry = GEOMETRIES[args.detector]
    values = (geometry["sides"], geometry["barrel_r"],
              geometry["barrel_half_z"], geometry["endcap_z"],
              geometry["endcap_inner_r"], geometry["endcap_outer_r"],
              args.boundary_margin_mm, args.time_margin_ns)
    return "lemurs_v2_rdf::keep_event(MCParticles, {})".format(
        ", ".join(f"{value:.17g}" for value in values))


def build_dataframe(ROOT, source, args):
    ROOT.gInterpreter.Declare(CPP_FILTER)
    hits = args.collection
    return (
        ROOT.RDataFrame("events", str(source))
        .Filter(selection_expression(args), "pre-ECAL interaction veto")
        .Define("energyCell",
                f"std::vector<float> out; for (const auto& hit : {hits}) "
                "out.push_back(hit.energy * 1.e3f); return out;")
        .Define("rhoCell",
                f"std::vector<int> out; for (const auto& hit : {hits}) "
                "out.push_back(static_cast<int>(hit.position.x)); return out;")
        .Define("phiCell",
                f"std::vector<int> out; for (const auto& hit : {hits}) "
                "out.push_back(static_cast<int>(hit.position.y)); return out;")
        .Define("zCell",
                f"std::vector<int> out; for (const auto& hit : {hits}) "
                "out.push_back(static_cast<int>(hit.position.z)); return out;")
        .Define("energyMC",
                "const auto& p=MCParticles[0].momentum; "
                "return std::sqrt(p.x*p.x+p.y*p.y+p.z*p.z)*1.e3;")
        .Define("thetaMC",
                "const auto& p=MCParticles[0].momentum; "
                "return std::atan2(std::hypot(p.x,p.y),p.z);")
        .Define("phiMC",
                "const auto& p=MCParticles[0].momentum; "
                "return std::atan2(p.y,p.x);")
    )


def shower_block(arrays, start, stop, shape, sampling_fraction):
    block = np.zeros((stop - start,) + shape, dtype=np.float32)
    skipped = 0
    for output_row, event in enumerate(range(start, stop)):
        hits = zip(arrays["energyCell"][event], arrays["rhoCell"][event],
                   arrays["phiCell"][event], arrays["zCell"][event])
        for energy, r, phi, z in hits:
            index = (int(r), int(phi), int(z))
            if any(i < 0 or i >= limit for i, limit in zip(index, shape)):
                skipped += 1
                continue
            block[(output_row,) + index] = np.float32(
                energy / sampling_fraction)
    return block, skipped


def convert(args):
    source = Path(args.inputFile).expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty ROOT file: {source}")

    output_dir = Path(args.outputDir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    name = args.output_name or source.stem + "_v2.h5"
    if Path(name).name != name or not name.endswith("_v2.h5"):
        raise ValueError("--output-name must be a basename ending in _v2.h5")
    output = output_dir / name
    partial = output.with_name(output.name + ".partial")
    if output.exists():
        if args.skip_existing:
            print(f"Already exists: {output}")
            return output
        raise FileExistsError(f"Output exists: {output}")
    if partial.exists():
        raise FileExistsError(f"Partial output exists: {partial}")

    import ROOT
    if hasattr(ROOT, "IsImplicitMTEnabled") and ROOT.IsImplicitMTEnabled():
        ROOT.DisableImplicitMT()
    arrays = build_dataframe(ROOT, source, args).AsNumpy(columns=list(COLUMNS))
    count = len(arrays["energyMC"])
    shape = (args.numR, args.numPhi, args.numZ)
    skipped = 0

    with h5py.File(partial, "x") as h5:
        h5.create_dataset("incident_energy", data=arrays["energyMC"],
                          dtype="f4", compression="gzip", compression_opts=9)
        h5.create_dataset("incident_theta", data=arrays["thetaMC"],
                          dtype="f4", compression="gzip", compression_opts=9)
        h5.create_dataset("incident_phi", data=arrays["phiMC"],
                          dtype="f4", compression="gzip", compression_opts=9)
        showers = h5.create_dataset(
            "showers", shape=(count,) + shape, maxshape=(None,) + shape,
            chunks=(args.batch_size,) + shape, dtype="f4",
            compression="gzip", compression_opts=9)
        for start in range(0, count, args.batch_size):
            stop = min(start + args.batch_size, count)
            block, outside = shower_block(
                arrays, start, stop, shape, args.samplingFraction)
            showers[start:stop] = block
            skipped += outside

        geometry = GEOMETRIES[args.detector]
        h5.attrs["detector"] = args.detector
        h5.attrs["geometry"] = geometry["name"]
        h5.attrs["geometry_source_commit"] = GEOMETRY_COMMIT
        h5.attrs["selection"] = "pre-ECAL interaction veto"
        h5.attrs["boundary_margin_mm"] = args.boundary_margin_mm
        h5.attrs["time_margin_ns"] = args.time_margin_ns
        h5.attrs["units"] = "energies: MeV; angles: rad"

    os.replace(partial, output)
    print(f"Wrote {count} showers to {output}")
    if skipped:
        print(f"Warning: skipped {skipped} hits outside {shape}")
    return output


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--inputFile", required=True)
    parser.add_argument("-o", "--outputDir", required=True)
    parser.add_argument("--detector", choices=sorted(GEOMETRIES), required=True)
    parser.add_argument("--output-name")
    parser.add_argument("--collection", default="ECalBarrelCollection")
    parser.add_argument("--numR", type=int, default=9)
    parser.add_argument("--numPhi", type=int, default=16)
    parser.add_argument("--numZ", type=int, default=45)
    parser.add_argument("--samplingFraction", type=float, default=1.)
    parser.add_argument("--boundary-margin-mm", type=float, default=1.)
    parser.add_argument("--time-margin-ns", type=float, default=.1)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args(argv)
    if min(args.numR, args.numPhi, args.numZ, args.batch_size) <= 0:
        parser.error("Dimensions and batch size must be positive")
    if not math.isfinite(args.samplingFraction) or args.samplingFraction <= 0:
        parser.error("Sampling fraction must be finite and positive")
    if not 0 <= args.boundary_margin_mm < 100:
        parser.error("Invalid boundary margin")
    if not 0 <= args.time_margin_ns < 100:
        parser.error("Invalid time margin")
    return args


if __name__ == "__main__":
    convert(parse_args())
