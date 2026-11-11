#!/usr/bin/env python3
"""Local plotting corrections for the LEMURS manuscript revision.

Derived from ddfastsim v1.0.0, https://doi.org/10.5281/zenodo.17076360.
Reads existing data only. See README.txt for scope and known limitations.
"""
import sys
root_not_found=False
try:
    from ROOT import RDataFrame, gSystem, gInterpreter
    gSystem.Load("libedm4hep")
    gInterpreter.Declare("#include <edm4hep/SimCalorimeterHitData.h>")
    gInterpreter.Declare("#include <edm4hep/MCParticleData.h>")
except ModuleNotFoundError:
    root_not_found=True
    pass
import h5py
import random
import numpy as np
import os
import matplotlib
if not os.environ.get("DRAW_SHOW"):   # write files only; set DRAW_SHOW=1 to open windows
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.axes as _maxes

# Legends: skip when there is nothing to label, use a smaller font, always top right,
# and add headroom above the data so the legend box does not cover the histograms.
_orig_legend = _maxes.Axes.legend
def _legend_with_headroom(self, *args, **kwargs):
    labels = [l for l in self.get_legend_handles_labels()[1] if l and not l.startswith("_")]
    if not labels:
        return None
    lo, hi = self.get_ylim()
    n = len(labels)
    if self.get_yscale() == "log":
        self.set_ylim(lo, hi * 10 ** (0.6 + 0.35 * n))
    else:
        self.set_ylim(lo, lo + (hi - lo) * (1.2 + 0.12 * n))
    kwargs["loc"] = "upper right"
    kwargs.setdefault("fontsize", "small")
    kwargs.setdefault("framealpha", 0.9)
    return _orig_legend(self, *args, **kwargs)
_maxes.Axes.legend = _legend_with_headroom
# Larger labels for figures reduced to manuscript width.
plt.rcParams.update({
    "font.size": 24, "axes.labelsize": 26, "axes.titlesize": 26,
    "xtick.labelsize": 22, "ytick.labelsize": 22, "legend.fontsize": 20,
    "savefig.dpi": 180, "savefig.bbox": "tight",
})
from math import pi,log10
import argparse

# Physical voxel widths for the released 6480-voxel representation,
# from Table 1 of the submitted LEMURS manuscript (dr, dz), in mm.
MESH_WIDTHS_MM = {
    "Par04SiW": (4.65, 3.4),
    "Par04SciPb": (8.0, 5.6),
    "ODD": (4.65, 5.05),
    "FCCeeCLD": (4.65, 5.05),
    "FCCeeALLEGRO": (5.0, 11.0),
}


def mesh_parameters(detector, dataset):
    """Return named mesh dimensions; reject unsupported combinations explicitly."""
    detector = "Par04SiW" if detector == "Par04" else detector
    if detector not in MESH_WIDTHS_MM:
        raise ValueError(f"Unknown detector: {detector}")
    dr, dz = MESH_WIDTHS_MM[detector]
    if dataset == 2:
        return detector, 9, 16, 45, dr, dz
    if dataset == 3 and detector in ("Par04SiW", "ODD"):
        # Explicit legacy option, not a claim of high-granularity LEMURS data.
        return detector, 18, 50, 45, dr / 2, dz
    raise ValueError("--dataset 3 is a legacy option supported only for Par04SiW and ODD; use --dataset 2 for released LEMURS files")


parser = argparse.ArgumentParser(description='Plot existing shower files')
parser.add_argument('--detector', '-d', help='Detector whose physical mesh widths apply (Par04 is an alias for Par04SiW)', choices=['Par04'] + list(MESH_WIDTHS_MM), default="Par04SiW")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--plotEratio', action='store_true')
group.add_argument('--plotEdep', action='store_true')
group.add_argument('--plotEcell', action='store_true')
group.add_argument('--plotLogEcell', action='store_true')
group.add_argument('--plotNumcells', action='store_true')
group.add_argument('--plotLong', action='store_true')
group.add_argument('--plotTrans', action='store_true')
group.add_argument('--plotIncident', action='store_true')
parser.add_argument('--noRatio', action='store_true', help="By default for more than 1 file a ratio is plotted. Use this to just draw plots.")
parser.add_argument('--dataset', type=int, choices=[2, 3],required=True)
parser.add_argument('--infileRoot',nargs='+')
parser.add_argument('--infileEdm4hep',nargs='+')
parser.add_argument('--infileH5',nargs='+')
parser.add_argument('--outfile','-o',default="")
parser.add_argument('--legendRoot',nargs='+',default=[""])
parser.add_argument('--legendH5',nargs='+',default=[""])
parser.add_argument('--samplingFraction', type=float, default=1., help="Sampling fraction to use to rescale energy of ROOT or EDM4hep files. No change is done to H5 files. Defined as f=active/(active+absorber). Used only for eratio and logEcell plots.")
parser.add_argument('--energyFilter', type=float, nargs=2, required=False)
parser.add_argument('--axisRange', type=float, nargs=2, required=False)
parser.add_argument('--ncpus', type=int, default=4, help='Number of CPUs to use in analysis')
parser.add_argument('-n', '--numEvents', type=int, default=100, help='Number of events to plot')
args = parser.parse_args()

# files
filename_root = args.infileRoot
if not args.infileRoot:
    filename_root = []
filename_edm4hep = args.infileEdm4hep
if not args.infileEdm4hep:
    filename_edm4hep = []
if len(filename_root)+len(filename_edm4hep) > 0 and root_not_found:
    print("ROOT is needed to analyse ROOT files (plain root or edm4hep root)")
    sys.exit(1)
filename_h5 = args.infileH5
if not args.infileH5:
    filename_h5 = []

num_events = args.numEvents
# The dataset argument is explicit; all released LEMURS files use dataset 2.
try:
    detector_tag, mesh_r, mesh_phi, mesh_z, mesh_dr, mesh_dz = mesh_parameters(args.detector, args.dataset)
except ValueError as exc:
    parser.error(str(exc))
if args.dataset == 2:
    legend_root = [detector_tag + ', ROOT dataset2']
    ifLowGran = True
    num_r_low, num_phi_low, num_z_low = mesh_r, mesh_phi, mesh_z
    size_r_low, size_z_low = mesh_dr, mesh_dz
    size_phi_low = 2*pi/num_phi_low
else:
    legend_root = [detector_tag + ', ROOT dataset3']
    ifLowGran = False
    num_r, num_phi, num_z = mesh_r, mesh_phi, mesh_z
    size_r, size_z = mesh_dr, mesh_dz
    size_phi = 2*pi/num_phi


if args.legendRoot:
    legend_root = args.legendRoot
if args.legendH5:
    legend_h5 = args.legendH5
if len(legend_root) < len(filename_root) + len(filename_edm4hep):
    if legend_root:
        print("Use the same length for legend and files")
        exit(-1)
    else:
        for i in range(len(filename_root)+len(filename_edm4hep)):
            legend_root.append(legend_root[0])
if len(legend_h5) < len(filename_h5):
    if legend_h5:
        print("Use the same length for legend and files")
        exit(-1)
    else:
        for i in range(len(filename_h5)):
            legend_h5.append(legend_h5[0])

colours=['b','r','orange']
# files
df_list = []
events_list = []
found_angle_theta = False
found_angle_phi = False
for file_root in filename_root:
    if args.energyFilter:
        df = RDataFrame("virtualReadout", file_root)
        filtered = df.Filter(f"MCParticles.momentum > {args.energyFilter[0]} && EnergyMC < {args.energyFilter[1]}").Range(num_events)
        df_list.append(filtered)
        events_list.append(df_list[-1].AsNumpy(columns=["EnergyCell", "rhoCell", "zCell", "phiCell"]))
        print(f"Filtered number of events: {filtered.Count().GetValue()}")
    else:
        df = RDataFrame("virtualReadout", file_root).Range(num_events)
        df2 = RDataFrame("global", file_root).Range(num_events)
        df_list.append(df)
        df_list.append(df2)
        events_list.append(df_list[-2].AsNumpy(columns=["EnergyCell", "rhoCell", "zCell", "phiCell"]))
        events_list[-1]["EnergyMC"] = df_list[-1].AsNumpy(columns=["EnergyMC"])["EnergyMC"]
    print(events_list[-1].keys())
    print(sum(events_list[-1]['EnergyCell'][0]))
    print(len(events_list[-1]['EnergyMC']))
    print((events_list[-1]['EnergyMC'][0]))
    # Plain ROOT files come from Par04 which so far has no implementation of phi/theta
    found_angle_theta = False
    found_angle_phi = False

for file_edm4hep in filename_edm4hep:
    df = RDataFrame("events", file_edm4hep).Range(num_events)
    df_list.append(df.Define("EnergyCell","ROOT::VecOps::RVec<float> result; for (auto&p: ECalBarrelCollection) {result.push_back(p.energy);} return result;")\
    .Define("rhoCell","ROOT::VecOps::RVec<float> result; for (auto&p: ECalBarrelCollection) {result.push_back(p.position.x);} return result;")\
    .Define("zCell","ROOT::VecOps::RVec<float> result; for (auto&p: ECalBarrelCollection) {result.push_back(p.position.z);} return result;")\
    .Define("phiCell","ROOT::VecOps::RVec<float> result; for (auto&p: ECalBarrelCollection) {result.push_back(p.position.y);} return result;")\
    .Define("EnergyMC","ROOT::VecOps::RVec<float> result; for (auto&p: MCParticles) {result.push_back(sqrt(p.momentum.x*p.momentum.x+p.momentum.y*p.momentum.y+p.momentum.z*p.momentum.z));} return std::accumulate(result.begin(),result.end(),0.);")
    .Define("ThetaMC","ROOT::Math::XYZVector mom(MCParticles[0].momentum.x,MCParticles[0].momentum.y,MCParticles[0].momentum.z); return mom.Theta();")\
    .Define("PhiMC","ROOT::Math::XYZVector mom(MCParticles[0].momentum.x,MCParticles[0].momentum.y,MCParticles[0].momentum.z); return mom.Phi();"))
    events_list.append(df_list[-1].AsNumpy(columns=["EnergyCell", "rhoCell", "zCell", "phiCell","EnergyMC","ThetaMC","PhiMC"]))
    events_list[-1]["EnergyCell"] = events_list[-1]["EnergyCell"] * 1000.
    events_list[-1]["EnergyMC"] = events_list[-1]["EnergyMC"] * 1000.
    print(events_list[-1].keys())
    print(sum(events_list[-1]['EnergyCell'][0]))
    print(len(events_list[-1]['EnergyMC']))
    print((events_list[-1]['EnergyMC'][0]))
    found_angle_theta = True
    found_angle_phi = True

#hdf5 files
h5_eMC_list = []
h5_thetaMC_list = []
h5_phiMC_list = []
h5_showers_list = []

for file_h5 in filename_h5:
    with h5py.File(file_h5, 'r') as data:
        # Slice the first num_events entries
        h5_eMC_list.append(data['incident_energy'][:num_events])

        if 'incident_theta' in data:
            h5_thetaMC_list.append(data['incident_theta'][:num_events])
            found_angle_theta = True

        if 'incident_phi' in data:
            h5_phiMC_list.append(data['incident_phi'][:num_events])
            found_angle_phi = True

        expected_shape = (mesh_r, mesh_phi, mesh_z)
        if data['showers'].ndim != 4 or data['showers'].shape[1:] != expected_shape:
            parser.error(f"{file_h5}: expected shower shape (S, {mesh_r}, {mesh_phi}, {mesh_z}) for --dataset {args.dataset}, got {data['showers'].shape}")
        # Source-only plots do not need the potentially multi-gigabyte shower array.
        if not args.plotIncident:
            h5_showers_list.append(data['showers'][:num_events])

def energy_ratio_root(fig, events, colour, label):
    ratios = []
    for energies, eMC in zip(events["EnergyCell"], events["EnergyMC"]):
        edep = np.sum(energies)
        ratios.append(edep/eMC/args.samplingFraction)
    return np.array(ratios).mean(), fig.hist(ratios, bins=np.linspace(0.,1.2,200), facecolor=colour, alpha=0.4, label=label, density = True)

def energy_ratio_hdf5(fig, energies, showers, colour, label):
    # Vectorized calculation
    total_energy = np.sum(showers, axis=(1,2,3))
    ratio = total_energy / energies
    print(ratio.shape, total_energy.shape, energies.shape)
    fig.hist(ratio, facecolor=colour, alpha=0.4, label=label, density=True)
    print("mean = ", ratio.mean())
    return ratio.mean(), ratio

def energy_deposited_root(fig, events, colour, label):
    result = []
    for energies in events["EnergyCell"]:
        edep = np.sum(energies)
        result.append(edep/args.samplingFraction)
    return np.array(result).mean(), fig.hist(result, bins=np.linspace(min(result),max(result),200), facecolor=colour, alpha=0.4, label=label, density = True)

def energy_deposited_hdf5(fig, showers, colour, label):
    energy = np.sum(showers, axis=(1,2,3))
    return np.array(energy).mean(), fig.hist(energy, bins=np.linspace(min(energy),max(energy),200), facecolor=colour, alpha=0.4, label=label, density = True)

def num_cells_root(fig,events, colour, label):
    all_num = []
    for energies in events["EnergyCell"]:
        all_num.append(len(energies))
    if args.axisRange:
        mybins=np.linspace(args.axisRange[0], args.axisRange[1],50)
    else:
        mybins=np.linspace(1e2,2e5,50)
    return np.array(all_num).mean(), fig.hist(all_num, bins=mybins, facecolor=colour, alpha=0.4, label=label, density = True)

def num_cells_hdf5(fig, showers, colour, label):
    # Count non-zero cells
    ax = fig.gca()
    ax.hist(counts, bins=100, color=colour, label=label, histtype='step', density=True)
    return fig

def num_cells_hdf5(fig,showers, colour, label):
    all_num = np.count_nonzero(showers, axis=(1,2,3))
    if args.axisRange:
        mybins=np.linspace(args.axisRange[0], args.axisRange[1],50)
    else:
        mybins=np.linspace(1e2,2e5,50)
    return np.array(all_num).mean(), fig.hist(all_num, bins=mybins, facecolor=colour, alpha=0.4, label=label, density = True)

def cell_energy_root(fig,events, colour,label):
    all_en = []
    for energies in events["EnergyCell"]:
        sumG4 = np.sum(energies)
        for val in [en for en in energies]:
            all_en.append(val/args.samplingFraction)
    return fig.hist(all_en, bins=np.linspace(1e-2,1,1000), facecolor=colour, alpha=0.6, label=label, density = True)

def cell_energy_log_root(fig,events, colour, label):
    all_logEn = []
    for energies in events["EnergyCell"]:
        sumG4 = np.sum(energies)
        for val in [log10(en/args.samplingFraction) for en in energies]:
            all_logEn.append(val)
    return fig.hist(all_logEn, bins=np.linspace(-3,4,100), facecolor=colour, alpha=0.4, label=label, density = True)

def cell_energy_log_hdf5(fig, showers, colour, label):
    energies = showers[showers > 0]  # 1D array of all non-zero energies
    all_logEn = np.log10(energies)
    hist_output = fig.hist(all_logEn, bins=np.linspace(-3, 4, 100),
                          facecolor=colour, alpha=0.4, label=label, density=True)
    return hist_output

def phi_distribution_root(fig,events, colour):
    all_phis = []
    all_weights = []
    for phi_ids, energies in zip(events["phiCell"], events["EnergyCell"]):
        sumG4 = np.sum(energies)
        phiMean_G4 = np.sum([phi_ids[i] * size_phi * energies[i] for i in range(len(phi_ids))])
        phiMean_G4 /= sumG4
        for val in [phi_ids[i] * size_phi - phiMean_G4 for i in range(len(phi_ids))]:
            all_phis.append(val)
        for val in [energies[i] for i in range(len(phi_ids))]:
            all_weights.append(val)
    return fig.hist(all_phis, bins=np.linspace((-num_phi+0.5)*size_phi, (num_phi - 0.5)*size_phi,2*num_phi-1), weights=all_weights , facecolor=colour, alpha=0.5, label=label, density = False)

def long_profile_root(fig,events, colour,label, lowGran = False):
    if lowGran:
        numz = num_z_low
        numr = num_r_low
        numphi = num_phi_low
        sizer = size_r_low
        sizephi = size_phi_low
        sizez = size_z_low
    else:
        numz = num_z
        numr = num_r
        numphi = num_phi
        sizer = size_r
        sizephi = size_phi
        sizez = size_z
    return fig.hist(np.concatenate(events["zCell"]), bins=np.linspace(0,numz-1,numz), weights=np.concatenate(events["EnergyCell"]) , facecolor=colour, alpha=0.4, label=label, density = False)

def long_profile_hdf5(fig,showers, colour,label, lowGran = False):
    if lowGran:
        numz = num_z_low
        numr = num_r_low
        numphi = num_phi_low
        sizer = size_r_low
        sizephi = size_phi_low
        sizez = size_z_low
    else:
        numz = num_z
        numr = num_r
        numphi = num_phi
        sizer = size_r
        sizephi = size_phi
        sizez = size_z
    profile = np.sum(showers, axis=(1,2))
    mean_profile = np.mean(profile, axis=0)
    z_ids = np.linspace(0, numz-1, numz)
    # One interval for each of the 45 discrete longitudinal layer indices.
    z_edges = np.arange(numz + 1) - 0.5
    return fig.hist(z_ids, bins=z_edges, weights=mean_profile , facecolor=colour, alpha=0.4, label=label, density = False)

def trans_profile_root(fig,events, colour,label, lowGran = False):
    if lowGran:
        numz = num_z_low
        numr = num_r_low
        numphi = num_phi_low
        sizer = size_r_low
        sizephi = size_phi_low
        sizez = size_z_low
    else:
        numz = num_z
        numr = num_r
        numphi = num_phi
        sizer = size_r
        sizephi = size_phi
        sizez = size_z
    return fig.hist(np.concatenate(events["rhoCell"]), bins=np.linspace(0,numr-1,numr), weights=np.concatenate(events["EnergyCell"]) , facecolor=colour, alpha=0.4, label=label, density = False)

def trans_profile_hdf5(fig,showers, colour,label,lowGran=False):
    if lowGran:
        numz = num_z_low
        numr = num_r_low
        numphi = num_phi_low
        sizer = size_r_low
        sizephi = size_phi_low
        sizez = size_z_low
    else:
        numz = num_z
        numr = num_r
        numphi = num_phi
        sizer = size_r
        sizephi = size_phi
        sizez = size_z

    # Sum over longitudinal axis (z)
    transverse_2D = np.sum(showers, axis=3)  # shape: (N_showers, num_r, num_phi)

    # Sum over phi to get radial profile per shower
    radial_profile_per_shower = np.sum(transverse_2D, axis=2)  # shape: (N_showers, num_r)

    # Mean over all showers
    mean_radial_profile = np.mean(radial_profile_per_shower, axis=0)  # shape: (num_r,)
    r_axis = np.arange(numr) * sizer + sizer/2
    # N radial bins require N+1 edges; retain the outermost radial bin.
    r_axis_bins = np.arange(numr + 1) * sizer

    return fig.hist(r_axis, bins=r_axis_bins, weights=mean_radial_profile , facecolor=colour, alpha=0.4, label=label, density = False)

def long_moments_root(fig_first, fig_second, events, colour, label,lowGran=False):
    if lowGran:
        numz = num_z_low
        numr = num_r_low
        numphi = num_phi_low
        sizer = size_r_low
        sizephi = size_phi_low
        sizez = size_z_low
    else:
        numz = num_z
        numr = num_r
        numphi = num_phi
        sizer = size_r
        sizephi = size_phi
        sizez = size_z
    all_first_moments = []
    all_second_moments = []
    for z_ids, energies in zip(events["zCell"], events["EnergyCell"]):
        sumG4 = np.sum(energies)
        first = np.sum([z_ids[i] * sizez * energies[i] for i in range(len(z_ids))])
        first /= sumG4
        all_first_moments.append(first)
        second = np.sum([pow(z_ids[i] * sizez - first, 2) * energies[i] for i in range(len(z_ids))])
        second /= sumG4
        all_second_moments.append(second)
    return fig_first.hist(all_first_moments, bins=np.linspace(0, 0.4*numz*sizez, 32) , facecolor=colour, alpha=0.4, label=label, density = False), fig_second.hist(all_second_moments, bins=np.linspace(0, pow(numz*sizez,2) / 35., 32) , facecolor=colour, alpha=0.4, label=label, density = False)

def long_moments_hdf5(fig_first, fig_second, showers, colour, label,lowGran=False):
    if lowGran:
        numz = num_z_low
        numr = num_r_low
        numphi = num_phi_low
        sizer = size_r_low
        sizephi = size_phi_low
        sizez = size_z_low
    else:
        numz = num_z
        numr = num_r
        numphi = num_phi
        sizer = size_r
        sizephi = size_phi
        sizez = size_z
    profile = np.sum(showers, axis=(1,2))  # shape: (N_showers, numz)

    # Preserve the archived index-scaled lower-edge convention.
    z_positions = np.linspace(0, numz - 1, numz) * sizez  # shape (numz,)
    # Total energy per shower
    total_energy = np.sum(profile, axis=1)[:, np.newaxis]  # shape: (N_showers, 1)

    # First moment: sum(z * energy) / total_energy
    first_moments = np.sum(profile * z_positions, axis=1) / total_energy[:,0]  # shape: (N_showers,)

    # Second moment: sum((z - mean_z)^2 * energy) / total_energy
    second_moments = np.sum(profile * (z_positions - first_moments[:, np.newaxis])**2, axis=1) / total_energy[:,0]

    # Histograms (same bins as original)
    fig1_hist = fig_first.hist(first_moments, bins=np.linspace(0, 0.5*numz*sizez, 32),
                               facecolor=colour, alpha=0.4, label=label, density=False)
    fig2_hist = fig_second.hist(second_moments, bins=np.linspace(0, (numz*sizez)**2 / 35., 32),
                                facecolor=colour, alpha=0.4, label=label, density=False)

    return fig1_hist, fig2_hist

def trans_moments_root(fig_first, fig_second, events, colour, label, lowGran = False):
    if lowGran:
        numz = num_z_low
        numr = num_r_low
        numphi = num_phi_low
        sizer = size_r_low
        sizephi = size_phi_low
        sizez = size_z_low
    else:
        numz = num_z
        numr = num_r
        numphi = num_phi
        sizer = size_r
        sizephi = size_phi
        sizez = size_z
    all_first_moments = []
    all_second_moments = []
    for r_ids, energies in zip(events["rhoCell"], events["EnergyCell"]):
        sumG4 = np.sum(energies)
        first = np.sum([r_ids[i] * sizer * energies[i] for i in range(len(r_ids))])
        first /= sumG4
        all_first_moments.append(first)
        second = np.sum([pow(r_ids[i] * sizer - first, 2) * energies[i] for i in range(len(r_ids))])
        second /= sumG4
        all_second_moments.append(second)
    return fig_first.hist(all_first_moments, bins=np.linspace(0, 0.75*numr*sizer,128) , facecolor=colour, alpha=0.4, label=label, density = False), fig_second.hist(all_second_moments, bins=np.linspace(0, pow(numr*sizer,2) / 8., 128) , facecolor=colour, alpha=0.4, label=label, density = False)

def trans_moments_hdf5(fig_first, fig_second, showers, colour, label,lowGran=False):
    if lowGran:
        numz = num_z_low
        numr = num_r_low
        numphi = num_phi_low
        sizer = size_r_low
        sizephi = size_phi_low
        sizez = size_z_low
    else:
        numz = num_z
        numr = num_r
        numphi = num_phi
        sizer = size_r
        sizephi = size_phi
        sizez = size_z

    # Radial positions in mm
    # Preserve the archived index-scaled lower-edge convention.
    r_positions = np.linspace(0, numr - 1, numr) * sizer  # shape (numr,)

    # Sum over longitudinal (z) and angular (phi) directions to get transverse profile
    # showers shape: (N_showers, num_r, num_phi, numz)
    profile = np.sum(showers, axis=(2,3))  # shape: (N_showers, numr)

    # Total energy per shower
    total_energy = np.sum(profile, axis=1)[:, np.newaxis]  # shape: (N_showers, 1)

    # First moment: mean radial position weighted by energy
    first_moments = np.sum(profile * r_positions, axis=1) / total_energy[:,0]  # shape: (N_showers,)

    # Second moment: variance of radial distribution weighted by energy
    second_moments = np.sum(profile * (r_positions - first_moments[:, np.newaxis])**2, axis=1) / total_energy[:,0]

    # Histograms (same bins as original transverse calculation)
    fig1_hist = fig_first.hist(first_moments, bins=np.linspace(0, 1.5*sizer, 61),  # finer range (R1, Fig. 13)
                               facecolor=colour, alpha=0.4, label=label, density=False)
    fig2_hist = fig_second.hist(second_moments, bins=np.linspace(0, (numr*sizer)**2 / 35., 32),
                                facecolor=colour, alpha=0.4, label=label, density=False)

    return fig1_hist, fig2_hist

def plot_cell_energy():
    fig1, (ax1) = plt.subplots(1, 1, figsize=(9,7), gridspec_kw={'height_ratios': [1]})
    for it in range(len(events_list)):
        cell_energy_root(ax1,events_list[it],colours[it],legend_root[it])
    ax1.set_xlabel('E (MeV)')
    ax1.set_ylabel('Probability density')
    ax1.set_title('Voxel energy distribution')
    ax1.grid(True, which='both')
    ax1.legend()
    fig1.tight_layout()
    fig1.savefig(args.outfile+'cell_energy_lin.png')

def plot_cell_log_energy():
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        fig1, (ax1, ax1b) = plt.subplots(2, 1, figsize=(9,12), gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig1, (ax1) = plt.subplots(1, 1, figsize=(9,7), gridspec_kw={'height_ratios': [1]})
    val_of_bins_x1 = []
    edges_of_bins_x1 = []
    for it in range(len(events_list)):
        (vals, edges, _) = cell_energy_log_root(ax1,events_list[it],colours[it],legend_root[it])
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
    for it in range(len(h5_eMC_list)):
        (vals, edges, _) = cell_energy_log_hdf5(ax1,h5_showers_list[it],colours[it+len(events_list)],legend_h5[it])
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
    ax1.set_xlabel('log10(E/MeV)')
 #   ax1.set_ylim(bottom=1e-6)
    ax1.set_yscale('log')
    ax1.set_ylabel('Probability density')
    ax1.set_title('Voxel energy distribution')
    ax1.grid(True)
    ax1.legend()
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        ax1b.set_ylabel('ratio')
        bincenter = 0.5 * (edges_of_bins_x1[0][1:] + edges_of_bins_x1[0][:-1])
        if len(events_list) > 0:
            legend_zero = legend_root[0]
        elif len(h5_eMC_list) > 0:
            legend_zero = legend_h5[0]
        for it in range(1,len(events_list)+len(h5_eMC_list)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins_x1[it]))
            np.divide(val_of_bins_x1[it], val_of_bins_x1[0], out=ratio, where=(val_of_bins_x1[0] != 0))
            error = np.divide(val_of_bins_x1[it] * np.sqrt(val_of_bins_x1[0]) + val_of_bins_x1[0] * np.sqrt(val_of_bins_x1[it]),np.power(val_of_bins_x1[0], 2), where=(val_of_bins_x1[0] != 0))
            ax1b.scatter(edges_of_bins_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        ax1b.legend()
        ax1b.set_ylim(0.25,4)
        ax1b.legend()
    fig1.tight_layout()
    fig1.savefig(args.outfile+'cell_energy_log.png')

def plot_num_cells():
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        fig3, (ax3, ax3b) = plt.subplots(2, 1, figsize=(9,12), gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig3, (ax3) = plt.subplots(1, 1, figsize=(9,7), gridspec_kw={'height_ratios': [1]})
    val_of_bins_x1 = []
    edges_of_bins_x1 = []
    for it in range(len(events_list)):
        mean_numcells_root, (vals, edges, _) = num_cells_root(ax3,events_list[it],colours[it],legend_root[it])
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
    for it in range(len(h5_eMC_list)):
        mean_numcells_h5, (vals, edges, _) = num_cells_hdf5(ax3,h5_showers_list[it],colours[it+len(events_list)],legend_h5[it])
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
    ax3.set_xlabel('Number of nonzero voxels')
    ax3.set_ylabel('Probability density')
    ax3.set_title('Number of nonzero voxels')
    ax3.grid(True)
    ax3.legend()
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        ax3b.set_ylabel('ratio')
        bincenter = 0.5 * (edges_of_bins_x1[0][1:] + edges_of_bins_x1[0][:-1])
        if len(events_list) > 0:
            legend_zero = legend_root[0]
        elif len(h5_eMC_list) > 0:
            legend_zero = legend_h5[0]
        for it in range(1,len(events_list)+len(h5_eMC_list)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins_x1[it]))
            np.divide(val_of_bins_x1[it], val_of_bins_x1[0], out=ratio, where=(val_of_bins_x1[0] != 0))
            error = np.divide(val_of_bins_x1[it] * np.sqrt(val_of_bins_x1[0]) + val_of_bins_x1[0] * np.sqrt(val_of_bins_x1[it]),np.power(val_of_bins_x1[0], 2), where=(val_of_bins_x1[0] != 0))
            ax3b.scatter(edges_of_bins_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        ax3b.legend()
        ax3b.set_ylim(0.25,4)
        ax3b.legend()
    fig3.tight_layout()
    fig3.savefig(args.outfile+'num_cells.png')

def plot_energy_deposited():
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        fig4, (ax4, ax4b) = plt.subplots(2, 1,figsize=(9, 12), gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig4, (ax4) = plt.subplots(1, 1,figsize=(9, 7), gridspec_kw={'height_ratios': [1]})
    val_of_bins_x1 = []
    edges_of_bins_x1 = []
    for it in range(len(events_list)):
        mean_edep_root, (vals, edges, _) = energy_deposited_root(ax4,events_list[it],colours[it], legend_root[it])
        ax4.axvline(mean_edep_root, color=colours[it], linestyle='dashed', linewidth=1, label=fr"<E> = {mean_edep_root:.2f}")
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
    for it in range(len(h5_eMC_list)):
        mean_edep_hdf5, (vals, edges, _) = energy_deposited_hdf5(ax4, h5_showers_list[it], colours[it+len(events_list)], legend_h5[it])
        ax4.axvline(mean_edep_hdf5, color=colours[it+len(events_list)], linestyle='dashed', linewidth=1, label=fr"<E> = {mean_edep_hdf5:.2f}")
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
    ax4.set_xlabel(r'$E_{\mathrm{dep}}$ (MeV)')
    ax4.set_ylabel('Probability density')
    ax4.set_title('Energy deposited')
    ax4.grid(True)
    ax4.legend()
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        ax4b.set_ylabel('ratio')
        bincenter = 0.5 * (edges_of_bins_x1[0][1:] + edges_of_bins_x1[0][:-1])
        if len(events_list) > 0:
            legend_zero = legend_root[0]
        elif len(h5_eMC_list) > 0:
            legend_zero = legend_h5[0]
        for it in range(1,len(val_of_bins_x1)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins_x1[it]))
            np.divide(val_of_bins_x1[it], val_of_bins_x1[0], out=ratio, where=(val_of_bins_x1[0] != 0))
            error = np.divide(val_of_bins_x1[it] * np.sqrt(val_of_bins_x1[0]) + val_of_bins_x1[0] * np.sqrt(val_of_bins_x1[it]),np.power(val_of_bins_x1[0], 2), where=(val_of_bins_x1[0] != 0))
            ax4b.scatter(edges_of_bins_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        ax4b.legend()
        # ax4b.set_ylim(0.25,4)
    fig4.tight_layout()
    fig4.savefig(args.outfile+'energy_deposited_small.png')

def plot_energy_ratio():
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        fig4, (ax4, ax4b) = plt.subplots(2, 1,figsize=(9, 12), gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig4, (ax4) = plt.subplots(1, 1,figsize=(9, 7), gridspec_kw={'height_ratios': [1]})
    val_of_bins_x1 = []
    edges_of_bins_x1 = []
    for it in range(len(events_list)):
        mean_eratio_root, (vals, edges, _) = energy_ratio_root(ax4,events_list[it],colours[it], legend_root[it])
        ax4.axvline(mean_eratio_root, color=colours[it], linestyle='dashed', linewidth=1, label=fr"$\langle f\rangle = {mean_eratio_root:.5f}$")
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
    for it in range(len(h5_eMC_list)):
        mean_eratio_hdf5, ratios = energy_ratio_hdf5(ax4, h5_eMC_list[it], h5_showers_list[it], colours[it+len(events_list)], legend_h5[it])
        ax4.axvline(mean_eratio_hdf5, color=colours[it+len(events_list)], linestyle='dashed', linewidth=1, label=fr"$\langle f\rangle = {mean_eratio_hdf5:.5f}$")
        val_of_bins_x1.append(ratios)
        edges_of_bins_x1.append(np.linspace(0., 1.2, 200))
    ax4.set_xlabel(r'$f = E_{\mathrm{dep}} / E_{\mathrm{MC}}$')
    ax4.set_ylabel('Probability density')
    ax4.set_title('Deposited to incident energy ratio')
    ax4.grid(True)
    ax4.legend()
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        ax4b.set_ylabel('ratio')
        bincenter = 0.5 * (edges_of_bins_x1[0][1:] + edges_of_bins_x1[0][:-1])
        if len(events_list) > 0:
            legend_zero = legend_root[0]
        elif len(h5_eMC_list) > 0:
            legend_zero = legend_h5[0]
        for it in range(1,len(val_of_bins_x1)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins_x1[it]))
            np.divide(val_of_bins_x1[it], val_of_bins_x1[0], out=ratio, where=(val_of_bins_x1[0] != 0))
            error = np.divide(val_of_bins_x1[it] * np.sqrt(val_of_bins_x1[0]) + val_of_bins_x1[0] * np.sqrt(val_of_bins_x1[it]),np.power(val_of_bins_x1[0], 2), where=(val_of_bins_x1[0] != 0))
            ax4b.scatter(edges_of_bins_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        ax4b.legend()
        ax4b.set_ylim(0.25,4)
    fig4.tight_layout()
    fig4.savefig(args.outfile+'energy_ratio_small.png')


def plot_incident():
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        fig, ((axE, axT, axP), (axEr, axTr, axPr)) = plt.subplots(2, 3 ,figsize=(21, 9), gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig, ((axE, axT, axP)) = plt.subplots(1, 3 ,figsize=(21, 7), gridspec_kw={'height_ratios': [1]})
    val_of_bins_E = []
    edges_of_bins_E = []
    val_of_bins_T = []
    edges_of_bins_T = []
    val_of_bins_P = []
    edges_of_bins_P = []
    #energy distribution
    if args.axisRange:
        mybins=np.linspace(args.axisRange[0], args.axisRange[1],100)
    else:
        mybins=np.linspace(1e3,1e6,30)
    for it in range(len(events_list)):
        (vals, edges, _) = axE.hist(events_list[it]["EnergyMC"], bins = mybins, facecolor=colours[it], alpha=0.4,label=legend_root[it])
        val_of_bins_E.append(vals)
        edges_of_bins_E.append(edges)
    for it in range(len(h5_eMC_list)):
        (vals, edges, _) = axE.hist(h5_eMC_list[it], bins = mybins, facecolor=colours[it+len(events_list)], alpha=0.4,label=legend_h5[it])
        val_of_bins_E.append(vals)
        edges_of_bins_E.append(edges)
    axE.set_xlabel('$E_{MC}$(MeV)')
    axE.set_ylim(bottom=1)
    axE.set_yscale('log')
    axE.set_ylabel('Entries')
    axE.set_title('Incident energy distribution')
    axE.grid(True)
    axE.legend(loc="lower right")
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        axEr.set_ylabel('ratio')
        bincenter = 0.5 * (edges_of_bins_E[0][1:] + edges_of_bins_E[0][:-1])
        if len(events_list) > 0:
            legend_zero = legend_root[0]
        elif len(h5_eMC_list) > 0:
            legend_zero = legend_h5[0]
        for it in range(1,len(events_list)+len(h5_eMC_list)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins_E[0]))
            ratio = np.divide(val_of_bins_E[it], val_of_bins_E[0], out=ratio,where=(val_of_bins_E[0] != 0))
            error = np.divide(val_of_bins_E[it] * np.sqrt(val_of_bins_E[0]) + val_of_bins_E[0] * np.sqrt(val_of_bins_E[it]),np.power(val_of_bins_E[0], 2), where=(val_of_bins_E[0] != 0))
            axEr.scatter(edges_of_bins_E[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        axEr.legend(loc="lower right")
        axEr.set_ylim(0.25,4)
    if found_angle_theta:
        #theta distribution
        for it in range(len(events_list)):
            (vals, edges, _) = axT.hist(events_list[it]["ThetaMC"], bins = np.linspace(0,pi,100), facecolor=colours[it], alpha=0.4,label=legend_root[it])
            val_of_bins_T.append(vals)
            edges_of_bins_T.append(edges)
        for it in range(len(h5_thetaMC_list)):
            (vals, edges, _) = axT.hist(h5_thetaMC_list[it], bins = np.linspace(0,pi,30), facecolor=colours[it+len(events_list)], alpha=0.4,label=legend_h5[it])
            val_of_bins_T.append(vals)
            edges_of_bins_T.append(edges)
        axT.set_xlabel(r'$\theta_{MC}$ (rad)')
        axT.set_ylim(bottom=1)
        axT.set_yscale('log')
        axT.set_ylabel('Entries')
        axT.set_title('Incident polar angle distribution')
        axT.grid(True)
        axT.legend(loc='lower left')
        if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
            axTr.set_ylabel('ratio')
            bincenter = 0.5 * (edges_of_bins_T[0][1:] + edges_of_bins_T[0][:-1])
            for it in range(1,len(events_list)+len(h5_eMC_list)):
                if it < len(events_list):
                    current_legend = legend_root[it]
                else:
                    current_legend = legend_h5[it-len(events_list)]
                ratio = np.zeros(len(val_of_bins_T[0]))
                ratio = np.divide(val_of_bins_T[it], val_of_bins_T[0], out=ratio,where=(val_of_bins_T[0] != 0))
                error = np.divide(val_of_bins_T[it] * np.sqrt(val_of_bins_T[0]) + val_of_bins_T[0] * np.sqrt(val_of_bins_T[it]),np.power(val_of_bins_T[0], 2), where=(val_of_bins_T[0] != 0))
                axTr.scatter(edges_of_bins_T[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
            axTr.legend()
            axTr.set_ylim(0.25,4)
    if found_angle_phi:
        #phi distribution
        for it in range(len(events_list)):
            (vals, edges, _) = axP.hist(events_list[it]["PhiMC"], bins = np.linspace(-pi,pi,100), facecolor=colours[it], alpha=0.4,label=legend_root[it])
            val_of_bins_P.append(vals)
            edges_of_bins_P.append(edges)
        for it in range(len(h5_phiMC_list)):
            (vals, edges, _) = axP.hist(h5_phiMC_list[it], bins = np.linspace(-pi,pi,30), facecolor=colours[it+len(events_list)], alpha=0.4,label=legend_h5[it])
            val_of_bins_P.append(vals)
            edges_of_bins_P.append(edges)
        axP.set_xlabel(r'$\phi_{MC}$ (rad)')
        axP.set_ylim(bottom=1)
        axP.set_yscale('log')
        axP.set_ylabel('Entries')
        axP.set_title('Incident azimuthal angle distribution')
        axP.grid(True)
        axP.legend(loc='lower left')
        if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
            axPr.set_ylabel('ratio')
            bincenter = 0.5 * (edges_of_bins_P[0][1:] + edges_of_bins_P[0][:-1])
            for it in range(1,len(events_list)+len(h5_eMC_list)):
                if it < len(events_list):
                    current_legend = legend_root[it]
                else:
                    current_legend = legend_h5[it-len(events_list)]
                ratio = np.zeros(len(val_of_bins_P[0]))
                np.divide(val_of_bins_P[it], val_of_bins_P[0], out=ratio, where=(val_of_bins_P[0] != 0))
                error = np.divide(val_of_bins_P[it] * np.sqrt(val_of_bins_P[0]) + val_of_bins_P[0] * np.sqrt(val_of_bins_P[it]),np.power(val_of_bins_P[0], 2), where=(val_of_bins_P[0] != 0))
                axPr.scatter(edges_of_bins_P[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
            axPr.legend()
            axPr.set_ylim(0.25,4)
    fig.tight_layout()
    fig.savefig(args.outfile+'incident_distr.png')

def plot_longitudinal():
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        fig5, ((ax5a, ax5b, ax5c),(ax5aa, ax5bb, ax5cc)) = plt.subplots(2, 3, figsize=(21,9), gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig5, ((ax5a, ax5b, ax5c)) = plt.subplots(1, 3, figsize=(21,7), gridspec_kw={'height_ratios': [1]})
    val_of_bins_x1 = []
    edges_of_bins_x1 = []
    val_of_bins2_x1 = []
    edges_of_bins2_x1 = []
    val_of_bins3_x1 = []
    edges_of_bins3_x1 = []
    for it in range(len(events_list)):
        (vals, edges, _) = long_profile_root(ax5a,events_list[it],colours[it],legend_root[it], lowGran=ifLowGran)
        (vals2, edges2, _), (vals3, edges3, _) =long_moments_root(ax5b,ax5c, events_list[it],colours[it],legend_root[it], lowGran=ifLowGran)
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
        val_of_bins2_x1.append(vals2)
        edges_of_bins2_x1.append(edges2)
        val_of_bins3_x1.append(vals3)
        edges_of_bins3_x1.append(edges3)
    for it in range(len(h5_eMC_list)):
        (vals, edges, _) = long_profile_hdf5(ax5a,h5_showers_list[it],colours[it+len(events_list)],legend_h5[it], lowGran=ifLowGran)
        (vals2, edges2, _), (vals3, edges3, _) =long_moments_hdf5(ax5b,ax5c, h5_showers_list[it],colours[it+len(events_list)],legend_h5[it], lowGran=ifLowGran)
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
        val_of_bins2_x1.append(vals2)
        edges_of_bins2_x1.append(edges2)
        val_of_bins3_x1.append(vals3)
        edges_of_bins3_x1.append(edges3)
    ax5a.set_xlabel('Longitudinal layer index')
    ax5a.set_ylabel(r'$\langle E \rangle$ (MeV)')
    ax5a.set_title('Longitudinal profile')
    ax5a.grid(True)
    ax5a.legend()
    fig5.tight_layout()
    ax5b.set_xlabel(r'$\overline{z}$ (mm)')
    ax5b.set_ylabel('Entries')
    ax5b.set_title('Longitudinal first moment')
    ax5c.set_xlabel(r'$\sigma_z^2$ (mm$^2$)')
    ax5c.set_ylabel('Entries')
    ax5c.set_title('Longitudinal variance')
    ax5b.grid(True)
    ax5c.grid(True)
    ax5b.legend()
    ax5c.legend()
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        ax5aa.set_ylabel('ratio')
        if len(events_list) > 0:
            legend_zero = legend_root[0]
        elif len(h5_eMC_list) > 0:
            legend_zero = legend_h5[0]
        for it in range(1,len(events_list)+len(h5_eMC_list)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins_x1[0]))
            np.divide(val_of_bins_x1[it], val_of_bins_x1[0], out=ratio, where=(val_of_bins_x1[0] != 0))
            error = np.divide(val_of_bins_x1[it] * np.sqrt(val_of_bins_x1[0]) + val_of_bins_x1[0] * np.sqrt(val_of_bins_x1[it]), np.power(val_of_bins_x1[0], 2), where=(val_of_bins_x1[0] != 0))
            ax5aa.scatter(edges_of_bins_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        ax5aa.legend()
        ax5aa.set_ylim(0.25,4)
        ax5bb.set_ylabel('ratio')
        for it in range(1,len(events_list)+len(h5_eMC_list)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins2_x1[0]))
            np.divide(val_of_bins2_x1[it], val_of_bins2_x1[0], out=ratio, where=(val_of_bins2_x1[0] != 0))
            error = np.divide(val_of_bins2_x1[it] * np.sqrt(val_of_bins2_x1[0]) + val_of_bins2_x1[0] * np.sqrt(val_of_bins2_x1[it]),np.power(val_of_bins2_x1[0], 2), where=(val_of_bins2_x1[0] != 0))
            ax5bb.scatter(edges_of_bins2_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        ax5bb.legend()
        ax5bb.set_ylim(0.25,4)
        ax5cc.set_ylabel('ratio')
        for it in range(1,len(events_list)+len(h5_eMC_list)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins3_x1[0]))
            np.divide(val_of_bins3_x1[it], val_of_bins3_x1[0], out=ratio, where=(val_of_bins3_x1[0] != 0))
            error = np.divide(val_of_bins3_x1[it] * np.sqrt(val_of_bins3_x1[0]) + val_of_bins3_x1[0] * np.sqrt(val_of_bins3_x1[it]),np.power(val_of_bins3_x1[0], 2), where=(val_of_bins3_x1[0] != 0))
            ax5cc.scatter(edges_of_bins3_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        ax5cc.legend()
        ax5cc.set_ylim(0.25,4)
    fig5.tight_layout()
    fig5.savefig(args.outfile+'longitudinal_profile.png')

def plot_transverse():
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        fig6, ((ax6a, ax6b, ax6c),(ax6aa, ax6bb, ax6cc)) = plt.subplots(2, 3, figsize=(21,9), gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig6, ((ax6a, ax6b, ax6c)) = plt.subplots(1, 3, figsize=(21,7), gridspec_kw={'height_ratios': [1]})
    val_of_bins_x1 = []
    edges_of_bins_x1 = []
    val_of_bins2_x1 = []
    edges_of_bins2_x1 = []
    val_of_bins3_x1 = []
    edges_of_bins3_x1 = []
    for it in range(len(events_list)):
        (vals, edges, _) =trans_profile_root(ax6a,events_list[it],colours[it],legend_root[it], lowGran=ifLowGran)
        (vals2, edges2, _), (vals3, edges3, _) =trans_moments_root(ax6b,ax6c, events_list[it],colours[it],legend_root[it], lowGran=ifLowGran)
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
        val_of_bins2_x1.append(vals2)
        edges_of_bins2_x1.append(edges2)
        val_of_bins3_x1.append(vals3)
        edges_of_bins3_x1.append(edges3)
    for it in range(len(h5_eMC_list)):
        (vals, edges, _) =trans_profile_hdf5(ax6a,h5_showers_list[it],colours[it+len(events_list)],legend_h5[it], lowGran=ifLowGran)
        (vals2, edges2, _), (vals3, edges3, _) =trans_moments_hdf5(ax6b,ax6c, h5_showers_list[it],colours[it+len(events_list)],legend_h5[it], lowGran=ifLowGran)
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
        val_of_bins2_x1.append(vals2)
        edges_of_bins2_x1.append(edges2)
        val_of_bins3_x1.append(vals3)
        edges_of_bins3_x1.append(edges3)
    ax6a.set_xlabel('r (mm)')
    ax6a.set_ylabel(r'$\langle E \rangle$ (MeV)')
    ax6a.set_title('Transverse profile')
    ax6a.set_yscale('log')
    ax6a.grid(True)
    ax6b.set_xlabel(r'$\overline{r}$ (mm)')
    ax6b.set_ylabel('Entries')
    ax6b.set_title('Transverse first moment')
    ax6c.set_xlabel(r'$\sigma_r^2$ (mm$^2$)')
    ax6c.set_ylabel('Entries')
    ax6c.set_title('Transverse variance')
    ax6b.grid(True)
    ax6c.grid(True)
    ax6a.legend()
    ax6b.legend()
    ax6c.legend()
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        ax6aa.set_ylabel('ratio')
        if len(events_list) > 0:
            legend_zero = legend_root[0]
        elif len(h5_eMC_list) > 0:
            legend_zero = legend_h5[0]
        for it in range(1,len(events_list)+len(h5_eMC_list)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins_x1[0]))
            np.divide(val_of_bins_x1[it], val_of_bins_x1[0], out=ratio,  where=(val_of_bins_x1[0] != 0))
            error = np.divide(val_of_bins_x1[it] * np.sqrt(val_of_bins_x1[0]) + val_of_bins_x1[0] * np.sqrt(val_of_bins_x1[it]),np.power(val_of_bins_x1[0], 2), where=(val_of_bins_x1[0] != 0))
            ax6aa.scatter(edges_of_bins_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        ax6aa.legend()
        ax6aa.set_ylim(0.25,4)
        ax6bb.set_ylabel('ratio')
        for it in range(1,len(events_list)+len(h5_eMC_list)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins2_x1[0]))
            np.divide(val_of_bins2_x1[it], val_of_bins2_x1[0], out=ratio,  where=(val_of_bins2_x1[0] != 0))
            error = np.divide(val_of_bins2_x1[it] * np.sqrt(val_of_bins2_x1[0]) + val_of_bins2_x1[0] * np.sqrt(val_of_bins2_x1[it]),np.power(val_of_bins2_x1[0], 2), where=(val_of_bins2_x1[0] != 0))
            ax6bb.scatter(edges_of_bins2_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        ax6bb.legend()
        ax6bb.set_ylim(0.25,4)
        ax6cc.set_ylabel('ratio')
        for it in range(1,len(events_list)+len(h5_eMC_list)):
            if it < len(events_list):
                current_legend = legend_root[it]
            else:
                current_legend = legend_h5[it-len(events_list)]
            ratio = np.zeros(len(val_of_bins3_x1[0]))
            np.divide(val_of_bins3_x1[it], val_of_bins3_x1[0], out=ratio, where=(val_of_bins3_x1[0] != 0))
            error = np.divide(val_of_bins3_x1[it] * np.sqrt(val_of_bins3_x1[0]) + val_of_bins3_x1[0] * np.sqrt(val_of_bins3_x1[it]),np.power(val_of_bins3_x1[0], 2), where=(val_of_bins3_x1[0] != 0))
            ax6cc.scatter(edges_of_bins3_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=current_legend.split(",")[0]+"/"+legend_zero.split(",")[0])
        ax6cc.legend()
        ax6cc.set_ylim(0.25,4)
    fig6.tight_layout()
    fig6.savefig(args.outfile+'transverse_profile.png')

def plot_incident_energy():
    if len(events_list)+len(h5_eMC_list) > 1:
        fig8, (ax8, ax8b) = plt.subplots(2, 1, figsize=(9,12), gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig8, (ax8) = plt.subplots(1, 1, figsize=(9,7), gridspec_kw={'height_ratios': [1]})
    val_of_bins_x1 = []
    edges_of_bins_x1 = []
    if args.axisRange:
        mybins=np.linspace(args.axisRange[0], args.axisRange[1],100)
    else:
        mybins=np.linspace(1e3,1e6,100)
    for it in range(len(events_list)):
        (vals, edges, _) = ax8.hist(events_list[it]["EnergyMC"], bins = mybins, facecolor=colours[it], alpha=0.6,label=legend_root[it])
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
    ax8.set_xlabel('$E_{MC}$(MeV)')
    ax8.set_ylim(bottom=1)
    ax8.set_yscale('log')
    ax8.set_ylabel('Entries')
    ax8.set_title('Incident energy distribution')
    ax8.grid(True)
    ax8.legend()
    if len(events_list)+len(h5_eMC_list) > 1 and not args.noRatio:
        ax8b.set_ylabel('ratio')
        for it in range(1,len(events_list)):
            ratio = np.divide(val_of_bins_x1[it], val_of_bins_x1[0], where=(val_of_bins_x1[0] != 0))
            error = np.divide(val_of_bins_x1[it] * np.sqrt(val_of_bins_x1[0]) + val_of_bins_x1[0] * np.sqrt(val_of_bins_x1[it]),np.power(val_of_bins_x1[0], 2), where=(val_of_bins_x1[0] != 0))
            ax8b.plot(edges_of_bins_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=legend_root[it].split(",")[0]+"/"+legend_root[0].split(",")[0])
        ax8b.legend()
        ax8b.set_ylim(0.25,4)
    fig8.tight_layout()
    fig8.savefig(args.outfile+'incident_energy.png')

if args.plotEdep:
    plot_energy_deposited()
if args.plotEratio:
    plot_energy_ratio()
if args.plotEcell:
    plot_cell_energy()
if args.plotLogEcell:
    plot_cell_log_energy()
if args.plotNumcells:
    plot_num_cells()
if args.plotLong:
    plot_longitudinal()
if args.plotTrans:
    plot_transverse()
if args.plotIncident:
    plot_incident()
plt.show()
