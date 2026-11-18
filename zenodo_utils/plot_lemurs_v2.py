#!/usr/bin/env python3
"""Plot LEMURS v2 selection controls (2 columns) and particle locations (3 columns).

The selection page is made from the existing scalar converter caches.  The
particle-location page reads only MCParticles from the original ROOT files,
sequentially, and joins each event by source file plus original ROOT row.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_selection_v2 import binned_mean, hist_density, make_bins, read_inputs


PORTRAIT_PAGE = (8.27, 11.69)
PANEL_LABELS = "abcdef"
COLOURS = ("#333333", "#d95f02", "#1b9e77")
LABELS = ("Inclusive", "numMC==1", "v2")
ROW_NAMES = (
    "Simulated primary photon endpoints",
    "Secondary production vertices",
    "Secondary endpoints",
)


def expand_inputs(inputs):
    paths = []
    for value in inputs:
        path = Path(value).expanduser().resolve()
        paths.extend(sorted(path.glob("*.selection.npz")) if path.is_dir()
                     else [path])
    if not paths:
        raise ValueError("No selection caches found")
    return paths


def input_parts(sources):
    parts = []
    for source in sources:
        name = Path(source["input"]["path"]).name
        match = re.search(r"(?:^|[_-])part(\d+)(?:\D|$)", name,
                          flags=re.IGNORECASE)
        parts.append(int(match.group(1)) if match else name)
    return parts


def parts_text(parts):
    return ", ".join(str(part) for part in parts)


def style_axis(ax):
    ax.tick_params(labelsize=8.5)
    ax.grid(alpha=0.18, linewidth=0.5)
    ax.xaxis.label.set_size(9.5)
    ax.yaxis.label.set_size(9.5)
    ax.title.set_size(10)


def save_page(fig, output_base):
    pdf = output_base.with_suffix(".pdf")
    png = output_base.with_suffix(".png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    return pdf, png


def plot_selection_page(data, output_dir):
    energy = data["incident_energy_MeV"].astype(float) / 1000.
    total = data["etot_MeV"].astype(float) / 1000.
    response = np.divide(total, energy, out=np.full(len(energy), np.nan),
                         where=energy > 0)
    costheta = np.cos(data["incident_theta"].astype(float))
    phi = data["incident_phi"].astype(float)
    valid = (np.isfinite(energy) & (energy > 0) & np.isfinite(costheta)
             & np.isfinite(phi) & np.isfinite(response))
    if not valid.all():
        raise ValueError("Invalid conditioning/response values in selection cache")

    masks = (np.ones(len(energy), dtype=bool), data["nmc"] == 1,
             data["reason"] == 0)
    energy_hist_bins = make_bins(energy, 50)
    energy_curve_bins = make_bins(energy, 20)
    costheta_bins = make_bins(costheta, 50)
    phi_bins = np.linspace(-np.pi, np.nextafter(np.pi, np.inf), 51)
    response_bins = make_bins(response, 70, lower=0.)

    fig, axes = plt.subplots(3, 2, figsize=PORTRAIT_PAGE)
    axes = axes.ravel()
    spectra = (
        (energy, energy_hist_bins, r"$E_{\rm gen}$ (GeV)", "Spectrum of generated energy"),
        (costheta, costheta_bins, r"$\cos\theta$", r"Spectrum of $\cos\theta$"),
        (phi, phi_bins, r"$\phi$ (rad)", r"Spectrum of $\phi$"),
    )
    for panel, (values, bins, xlabel, title) in enumerate(spectra):
        ax = axes[panel]
        for mask, colour in zip(masks, COLOURS):
            ax.stairs(hist_density(values[mask], bins), bins, color=colour,
                      linewidth=1.35)
        ax.set(xlabel=xlabel, ylabel="Density",
               title=f"({PANEL_LABELS[panel]}) {title}")
        style_axis(ax)

    ax = axes[3]
    denominator, _ = np.histogram(energy, energy_curve_bins)
    for mask, colour in zip(masks, COLOURS):
        retained, _ = np.histogram(energy[mask], energy_curve_bins)
        efficiency = np.divide(
            retained, denominator, out=np.full(len(retained), np.nan),
            where=denominator > 0)
        ax.stairs(efficiency, energy_curve_bins, color=colour, linewidth=1.35,
                  baseline=None)
    ax.set(xlabel=r"Generated energy $E_{\rm gen}$ (GeV)",
           ylabel="Retained / original events", ylim=(0., 1.025),
           title="(d) Selection efficiency")
    style_axis(ax)

    ax = axes[4]
    centres = (energy_curve_bins[:-1] + energy_curve_bins[1:]) / 2.
    for mask, colour in zip(masks, COLOURS):
        _, mean, sem = binned_mean(
            energy, response, mask, energy_curve_bins)
        ax.errorbar(centres, mean, yerr=sem, color=colour, marker="o",
                    markersize=2.4, linewidth=1., capsize=1.5)
    ax.set(xlabel=r"Generated energy $E_{\rm gen}$ (GeV)",
           ylabel=r"Mean $E_{\rm voxels}/E_{\rm gen}$",
           title="(e) Mean response")
    style_axis(ax)

    ax = axes[5]
    for mask, colour in zip(masks, COLOURS):
        ax.stairs(hist_density(response[mask], response_bins), response_bins,
                  color=colour, linewidth=1.35)
    ax.set(xlabel=r"$E_{\rm voxels}/E_{\rm gen}$", ylabel="Density",
           title="(f) Response distribution")
    style_axis(ax)

    fig.legend([series.lines[0] for series in axes[4].containers], LABELS,
               loc="upper center", ncol=3, frameon=False, fontsize=9)
    fig.subplots_adjust(left=.12, right=.975, bottom=.06, top=.94,
                        hspace=.38, wspace=.29)
    paths = save_page(fig, output_dir / "selection_validation_combined")

    scenarios = {}
    for name, label, mask in zip(("inclusive", "numMC1", "v2"), LABELS, masks):
        finite_response = response[mask & np.isfinite(response)]
        scenarios[name] = {
            "label": label,
            "retained_events": int(mask.sum()),
            "fraction": float(mask.mean()),
            "zero_energy_showers": int(np.count_nonzero(mask & (total == 0))),
            "mean_response": (float(finite_response.mean())
                              if len(finite_response) else None),
        }
    return paths, scenarios


def _finite_xyz(value):
    xyz = (float(value.x), float(value.y), float(value.z))
    if not all(math.isfinite(component) for component in xyz):
        raise ValueError("Non-finite saved MC-particle coordinates")
    return xyz


def _coord(value):
    x, y, z = _finite_xyz(value)
    return z, math.hypot(x, y)


def _validate_primary(particles, source, row):
    if len(particles) == 0:
        raise ValueError(f"No MCParticles in {source}, ROOT row {row}")
    primary = particles[0]
    if int(primary.PDG) != 22 or int(primary.generatorStatus) != 1:
        raise ValueError(
            f"MCParticles[0] is not the generated photon in {source}, row {row}")
    origin = _finite_xyz(primary.vertex)
    _finite_xyz(primary.endpoint)
    _finite_xyz(primary.momentum)
    if not math.isfinite(float(primary.time)) or max(map(abs, origin)) > 1.:
        raise ValueError(
            f"Primary photon fails the current IP/finite checks in {source}, row {row}")
    for index in range(1, len(particles)):
        if int(particles[index].generatorStatus) != 0:
            raise ValueError(
                f"Additional generated particle in {source}, ROOT row {row}")
    return primary


def source_records(paths, sources, data):
    records = []
    offset = 0
    for cache_path, metadata in zip(paths, sources):
        event_count = int(metadata["total_events"])
        with np.load(cache_path, allow_pickle=False) as cache:
            cache_meta = json.loads(str(cache["metadata_json"].item()))
            if cache_meta["input"]["path"] != metadata["input"]["path"]:
                raise ValueError(f"Cache/source metadata mismatch: {cache_path}")
            if "root_row" not in cache:
                raise ValueError(f"Missing root_row in {cache_path}")
            root_rows = np.asarray(cache["root_row"], dtype=np.int64)
        if root_rows.shape != (event_count,):
            raise ValueError(f"Malformed root_row in {cache_path}")
        if np.any(root_rows < 0) or np.any(root_rows >= event_count):
            raise ValueError(f"Out-of-range root_row in {cache_path}")
        if len(np.unique(root_rows)) != event_count:
            raise ValueError(f"Duplicate root_row in {cache_path}")
        lookup = np.full(event_count, -1, dtype=np.int64)
        lookup[root_rows] = np.arange(event_count, dtype=np.int64)

        source = Path(metadata["input"]["path"])
        if not source.is_file():
            raise FileNotFoundError(f"Original ROOT input not found: {source}")
        stat = source.stat()
        expected = metadata["input"]
        if (stat.st_size != expected["size_bytes"]
                or stat.st_mtime_ns != expected["mtime_ns"]):
            raise ValueError(f"Original ROOT input changed since conversion: {source}")
        records.append({
            "source": source,
            "cache": cache_path,
            "events": event_count,
            "offset": offset,
            "lookup": lookup,
        })
        offset += event_count
    if offset != len(data["reason"]):
        raise ValueError("Combined cache length does not match source metadata")
    return records


def _open_tree(ROOT, record):
    root_file = ROOT.TFile.Open(str(record["source"]), "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"ROOT could not open {record['source']}")
    tree = root_file.Get("events")
    if not tree or not tree.InheritsFrom("TTree"):
        root_file.Close()
        raise ValueError(f"Expected events TTree in {record['source']}")
    if int(tree.GetEntries()) != record["events"]:
        root_file.Close()
        raise ValueError(f"ROOT/cache event-count mismatch: {record['source']}")
    if not tree.GetBranch("MCParticles"):
        root_file.Close()
        raise ValueError(f"Missing MCParticles branch: {record['source']}")
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("MCParticles", 1)
    tree.SetBranchStatus("MCParticles*", 1)
    tree.SetCacheSize(8 * 1024 * 1024)
    return root_file, tree


def _empty_stat():
    return {"count": 0, "min_z": math.inf, "max_z": -math.inf,
            "max_r": 0.}


def _update_stat(stat, z, radius):
    stat["count"] += 1
    stat["min_z"] = min(stat["min_z"], z)
    stat["max_z"] = max(stat["max_z"], z)
    stat["max_r"] = max(stat["max_r"], radius)


def scan_truth_ranges(records, data):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    if hasattr(ROOT, "IsImplicitMTEnabled") and ROOT.IsImplicitMTEnabled():
        ROOT.DisableImplicitMT()

    stats = [_empty_stat() for _ in ROW_NAMES]
    counts = {
        "original_events": 0,
        "v2_events": 0,
        "primary_all": 0,
        "primary_v2": 0,
        "secondary_all": 0,
        "secondary_v2": 0,
        "verified_root_rows": 0,
    }
    for record in records:
        root_file, tree = _open_tree(ROOT, record)
        try:
            for root_row in range(record["events"]):
                if int(tree.GetEntry(root_row)) <= 0:
                    raise OSError(
                        f"ROOT read failed: {record['source']}, row {root_row}")
                local = int(record["lookup"][root_row])
                global_row = record["offset"] + local
                particles = tree.MCParticles
                primary = _validate_primary(particles, record["source"], root_row)
                if int(data["nmc"][global_row]) != len(particles):
                    raise ValueError(
                        f"MC multiplicity mismatch: {record['source']}, row {root_row}")
                endpoint = _finite_xyz(primary.endpoint)
                cached = (float(data["endpoint_x_mm"][global_row]),
                          float(data["endpoint_y_mm"][global_row]),
                          float(data["endpoint_z_mm"][global_row]))
                if not np.allclose(endpoint, cached, rtol=0., atol=1.e-6):
                    raise ValueError(
                        f"Primary endpoint/root_row mismatch: {record['source']}, "
                        f"row {root_row}")
                selected = int(data["reason"][global_row]) == 0
                z, radius = _coord(primary.endpoint)
                _update_stat(stats[0], z, radius)
                counts["original_events"] += 1
                counts["primary_all"] += 1
                counts["v2_events"] += int(selected)
                counts["primary_v2"] += int(selected)
                for particle_index in range(1, len(particles)):
                    secondary = particles[particle_index]
                    production = _coord(secondary.vertex)
                    endpoint_coord = _coord(secondary.endpoint)
                    _update_stat(stats[1], *production)
                    _update_stat(stats[2], *endpoint_coord)
                    counts["secondary_all"] += 1
                    counts["secondary_v2"] += int(selected)
                counts["verified_root_rows"] += 1
        finally:
            root_file.Close()
    return stats, counts


def nice_ceiling(value):
    if not math.isfinite(value) or value <= 0:
        return 1.
    step = 10. ** (math.floor(math.log10(value)) - 1)
    return math.ceil(value / step) * step


def coordinate_edges(stats, geometry, z_bins, r_bins):
    factor = (1. / math.cos(math.pi / geometry["sides"])
              if geometry["sides"] else 1.)
    geometry_z = max(geometry["barrel_half_z"], geometry["endcap_back_z"])
    geometry_r = max(geometry["barrel_outer_r"],
                     geometry["endcap_outer_r"]) * factor
    result = []
    for row, stat in enumerate(stats):
        if stat["count"] == 0:
            z_extent, r_extent = geometry_z, geometry_r
        else:
            z_extent = max(abs(stat["min_z"]), abs(stat["max_z"]), geometry_z)
            r_extent = max(stat["max_r"], geometry_r)
        z_limit = nice_ceiling(z_extent * 1.02)
        r_limit = nice_ceiling(r_extent * 1.02)
        if row == 2:
            z_limit = r_limit = 10000.  # Common display window for secondary endpoints.
        result.append((np.linspace(-z_limit, z_limit, z_bins + 1),
                       np.linspace(0., r_limit, r_bins + 1)))
    return result


class HistogramBuffers:
    def __init__(self, edges, batch_size):
        self.edges = edges
        self.batch_size = batch_size
        self.histograms = [[np.zeros((len(z) - 1, len(r) - 1), dtype=np.int64)
                            for _ in LABELS] for z, r in edges]
        self.overflow = [[0] * len(LABELS) for _ in edges]
        self.buffers = {(row, column): ([], [])
                        for row in range(3) for column in range(len(LABELS))}
        self.buffered = 0

    def add(self, row, column, z, radius):
        z_values, r_values = self.buffers[(row, column)]
        z_values.append(z)
        r_values.append(radius)
        self.buffered += 1
        if self.buffered >= self.batch_size:
            self.flush()

    def flush(self):
        for (row, column), (z_values, r_values) in self.buffers.items():
            if not z_values:
                continue
            z_edges, r_edges = self.edges[row]
            block, _, _ = np.histogram2d(z_values, r_values,
                                         bins=(z_edges, r_edges))
            used = int(block.sum())
            self.histograms[row][column] += block.astype(np.int64)
            self.overflow[row][column] += len(z_values) - used
            z_values.clear()
            r_values.clear()
        self.buffered = 0


def fill_truth_histograms(records, data, edges, batch_size):
    import ROOT
    accumulator = HistogramBuffers(edges, batch_size)
    for record in records:
        root_file, tree = _open_tree(ROOT, record)
        try:
            for root_row in range(record["events"]):
                if int(tree.GetEntry(root_row)) <= 0:
                    raise OSError(
                        f"ROOT read failed: {record['source']}, row {root_row}")
                global_row = record["offset"] + int(record["lookup"][root_row])
                selected = int(data["reason"][global_row]) == 0
                particles = tree.MCParticles
                primary = particles[0]
                z, radius = _coord(primary.endpoint)
                accumulator.add(0, 0, z, radius)
                if data["nmc"][global_row] == 1:
                    accumulator.add(0, 1, z, radius)
                if selected:
                    accumulator.add(0, 2, z, radius)
                for particle_index in range(1, len(particles)):
                    secondary = particles[particle_index]
                    for row, value in ((1, secondary.vertex),
                                       (2, secondary.endpoint)):
                        z, radius = _coord(value)
                        accumulator.add(row, 0, z, radius)
                        if selected:
                            accumulator.add(row, 2, z, radius)
        finally:
            root_file.Close()
    accumulator.flush()
    return accumulator.histograms, accumulator.overflow


def overlay_geometry(ax, geometry):
    colour = "#555555"
    dash = (0, (4, 2))
    half_z = geometry["barrel_half_z"]
    inner = geometry["barrel_r"]
    outer = geometry["barrel_outer_r"]
    sides = geometry["sides"]
    factor = 1. / math.cos(math.pi / sides) if sides else 1.

    if sides:
        for radius in (inner, outer):
            ax.fill_between((-half_z, half_z), radius, radius * factor,
                            facecolor=colour, alpha=.08, edgecolor="none")
            ax.plot((-half_z, half_z), (radius, radius), color=colour,
                    linestyle=dash, linewidth=.8)
            ax.plot((-half_z, half_z), (radius * factor, radius * factor),
                    color=colour, linestyle=dash, linewidth=.8)
        for z in (-half_z, half_z):
            ax.plot((z, z), (inner, outer * factor), color=colour,
                    linestyle=dash, linewidth=.8)
    else:
        ax.add_patch(Rectangle((-half_z, inner), 2 * half_z, outer - inner,
                               fill=False, edgecolor=colour, linestyle=dash,
                               linewidth=.9))

    endcap_inner = geometry["endcap_inner_r"]
    endcap_outer = geometry["endcap_outer_r"] * factor
    front = geometry["endcap_z"]
    back = geometry["endcap_back_z"]
    for sign in (-1, 1):
        left = front if sign > 0 else -back
        ax.add_patch(Rectangle((left, endcap_inner), back - front,
                               endcap_outer - endcap_inner, fill=False,
                               edgecolor=colour, linestyle=dash, linewidth=.9))


def plot_spatial_page(geometry, counts, histograms, overflow, edges, output_dir):
    fig = plt.figure(figsize=(11.7, 13.0))
    grid = fig.add_gridspec(
        3, 4, left=.075, right=.935, bottom=.06, top=.895,
        hspace=.36, wspace=.13, width_ratios=(1, 1, 1, .048))
    axes = np.asarray([[fig.add_subplot(grid[row, column])
                        for column in range(3)] for row in range(3)])
    colourbar_axes = [fig.add_subplot(grid[row, 3]) for row in range(3)]
    primary = (counts["primary_all"], counts["numMC1_events"], counts["primary_v2"])
    secondary = (counts["secondary_all"], 0, counts["secondary_v2"])
    count_rows = (primary, secondary, secondary)
    subtitles = ("All simulated events", "scripts/edm4hep2h5.py",
                 "scripts/edm4hep2h5_v2.py")
    for ax, label, subtitle in zip(axes[0], LABELS, subtitles):
        x = ax.get_position().x0 + ax.get_position().width / 2
        fig.text(x, .975, label, ha="center", fontsize=15, weight="semibold")
        fig.text(x, .952, subtitle, ha="center", fontsize=12)

    for row in range(3):
        norm = LogNorm(vmin=1, vmax=max(2, int(histograms[row][0].max())))
        z_edges, r_edges = edges[row]
        image = None
        axes[row, 1].set_title(ROW_NAMES[row], fontsize=12, pad=7)
        for column in range(3):
            ax = axes[row, column]
            if histograms[row][column].sum() + overflow[row][column] != count_rows[row][column]:
                raise ValueError(f"Truth count mismatch in {ROW_NAMES[row]}, column {column}")
            values = np.ma.masked_less(histograms[row][column].T, 1)
            image = ax.pcolormesh(z_edges, r_edges, values, norm=norm,
                                  cmap="viridis", shading="auto",
                                  rasterized=True)
            overlay_geometry(ax, geometry)
            ax.set_xlim(z_edges[0], z_edges[-1])
            ax.set_ylim(r_edges[0], r_edges[-1])
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3, symmetric=True))
            ax.set_xlabel(r"$z$ (mm)", fontsize=12)
            if column == 0:
                ax.set_ylabel(r"$r$ (mm)", fontsize=12)
            else:
                ax.tick_params(labelleft=False)
            if column == 1 and row > 0:
                ax.text(.5, .5, "No saved secondaries", transform=ax.transAxes,
                        ha="center", va="center", fontsize=11,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=.9))
            ax.tick_params(labelsize=10.5)
            ax.grid(alpha=.18, linewidth=.5)
        colourbar = fig.colorbar(image, cax=colourbar_axes[row])
        colourbar.set_label("primary photons/bin" if row == 0
                            else "saved particles/bin", fontsize=11)
        colourbar.ax.tick_params(labelsize=10)

    paths = save_page(fig, output_dir / "particle_locations_combined")
    return paths, {
        ROW_NAMES[row]: {
            "z_range_mm": [float(edges[row][0][0]),
                           float(edges[row][0][-1])],
            "r_range_mm": [float(edges[row][1][0]),
                           float(edges[row][1][-1])],
            "all_count": int(count_rows[row][0]),
            "v2_count": int(count_rows[row][-1]),
            "overflow_all": int(overflow[row][0]),
            "overflow_v2": int(overflow[row][-1]),
            "numMC1_count": int(count_rows[row][1]),
            "overflow_numMC1": int(overflow[row][1]),
        }
        for row in range(3)
    }


def write_captions(output_dir, detector, parts, scenarios, counts, overflow):
    sample = f"{detector} " + ("part " if len(parts) == 1 else "parts ") + parts_text(parts)
    total, selected = counts["original_events"], counts["v2_events"]
    historical_context = (
        "CLD v1 used Inclusive for parts 1–4 and numMC==1 for parts 5–10."
        if detector == "CLD" else "ALLEGRO v1 used numMC==1.")
    selection = (
        f"Figure 1. {sample}. Inclusive, numMC==1 and v2 retain {total:,}, "
        f"{scenarios['numMC1']['retained_events']:,} and {selected:,} events "
        f"from the same simulation. The samples overlap. {historical_context} "
        "The v2 rule vetoes upstream primary endpoints and early upstream secondary "
        "production relative to nominal ECAL entry. Spectra and response densities "
        "have unit area. Efficiency uses all simulated events in each energy bin. "
        "Response includes empty showers. Error bars show standard errors of the mean."
    )
    spatial = (
        f"Figure 2. Positions of saved MC particles in {sample}, in global signed z "
        f"and cylindrical radius r. Inclusive, numMC==1 and v2 retain {total:,}, "
        f"{counts['numMC1_events']:,} and {selected:,} primary photons (one per event). "
        f"The same events contain {counts['secondary_all']:,}, zero and "
        f"{counts['secondary_v2']:,} saved secondaries, respectively. Each secondary "
        "contributes one vertex and one endpoint. The third row shows |z| and r up to "
        f"10,000 mm. Outside this window lie {overflow[2][0]:,} inclusive and "
        f"{overflow[2][2]:,} v2 endpoints. Each row shares a logarithmic count scale. "
        "Dashed lines mark projected nominal ECAL envelopes. "
        + ("Bands span the face and corner radii of the CLD barrel. "
           if detector == "CLD" else "") +
        "The maps alone cannot establish trajectories, direct parentage or timing."
    )
    path = output_dir / "captions.txt"
    captions = {"selection_validation_combined": selection,
                "particle_locations_combined": spatial}
    path.write_text("\n\n".join(captions.values()) + "\n")
    return path, captions


def run(args):
    paths = expand_inputs(args.inputs)
    data, detector, settings, sources = read_inputs(paths)
    parts = input_parts(sources)
    output_dir = Path(args.outdir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    expected_outputs = [
        output_dir / "selection_validation_combined.pdf",
        output_dir / "selection_validation_combined.png",
        output_dir / "particle_locations_combined.pdf",
        output_dir / "particle_locations_combined.png",
        output_dir / "combined_summary.json",
        output_dir / "captions.txt",
    ]
    if args.three_rules:
        expected_outputs.extend(output_dir / f"particle_locations_three_rules.{ext}"
                                for ext in ("pdf", "png"))
    if not args.replace and any(path.exists() for path in expected_outputs):
        raise FileExistsError("Combined output exists; use --replace to regenerate plots")

    records = source_records(paths, sources, data)
    stats, counts = scan_truth_ranges(records, data)
    if counts["original_events"] != len(data["reason"]):
        raise ValueError("Truth scan did not cover every original ROOT event")
    if counts["v2_events"] != int(np.count_nonzero(data["reason"] == 0)):
        raise ValueError("Truth scan and existing v2 event mask disagree")
    edges = coordinate_edges(
        stats, settings["geometry"], args.z_bins, args.r_bins)
    histograms, overflow = fill_truth_histograms(
        records, data, edges, args.truth_batch_size)
    selection_paths, scenarios = plot_selection_page(data, output_dir)
    counts["numMC1_events"] = scenarios["numMC1"]["retained_events"]
    spatial_paths, spatial_summary = plot_spatial_page(
        settings["geometry"], counts, histograms, overflow,
        edges, output_dir)
    if args.three_rules:
        # Keep the explicit filename used by callers of the preview option.
        aliases = tuple(output_dir / ("particle_locations_three_rules" + path.suffix)
                        for path in spatial_paths)
        for source, alias in zip(spatial_paths, aliases):
            shutil.copyfile(source, alias)
        spatial_paths += aliases
    caption_path, captions = write_captions(
        output_dir, detector, parts, scenarios, counts, overflow)

    summary = {
        "detector": detector,
        "input_parts": parts,
        "original_event_count": counts["original_events"],
        "v2_event_count": counts["v2_events"],
        "sources": [str(record["source"]) for record in records],
        "alignment": {
            "key": "original ROOT file + ROOT row",
            "verified_root_rows": counts["verified_root_rows"],
            "event_counts_match": True,
            "primary_endpoints_match_cache": True,
        },
        "selection_rules": scenarios,
        "particle_counts": {
            "primary_photons_all": counts["primary_all"],
            "primary_photons_numMC1": counts["numMC1_events"],
            "primary_photons_v2": counts["primary_v2"],
            "saved_secondaries_all": counts["secondary_all"],
            "saved_secondaries_numMC1": 0,
            "saved_secondaries_v2": counts["secondary_v2"],
        },
        "spatial_histograms": spatial_summary,
        "geometry": settings["geometry"],
        "outputs": [str(path) for path in
                    (*selection_paths, *spatial_paths, caption_path)],
        "captions": captions,
    }
    summary_path = output_dir / "combined_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({
        "detector": detector,
        "parts": parts,
        "original_events": counts["original_events"],
        "v2_events": counts["v2_events"],
        "saved_secondaries_all": counts["secondary_all"],
        "saved_secondaries_v2": counts["secondary_v2"],
        "overflow": overflow,
        "outputs": [str(path) for path in expected_outputs],
    }, indent=2))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+",
                        help="*.selection.npz files or detector directory")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--three-rules", action="store_true",
                        help="Also save the default 3×3 map under its earlier particle_locations_three_rules name")
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--truth-batch-size", type=int, default=50000,
                        help="Maximum scalar coordinates buffered before histogramming")
    parser.add_argument("--z-bins", type=int, default=220)
    parser.add_argument("--r-bins", type=int, default=160)
    args = parser.parse_args(argv)
    if min(args.truth_batch_size, args.z_bins, args.r_bins) <= 0:
        parser.error("Batch size and bin counts must be positive")
    return args


if __name__ == "__main__":
    try:
        run(parse_args())
    except (ValueError, OSError, KeyError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
