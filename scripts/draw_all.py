root_not_found=False
try:
    from ROOT import RDataFrame
except ModuleNotFoundError:
    root_not_found=True
    pass
import h5py
import random
import numpy as np
import matplotlib.pyplot as plt
from math import pi,log10,pow
import argparse

parser = argparse.ArgumentParser(description='Validation plots')
parser.add_argument('--detector', '-d', help='Detector to simulate', choices=['Par04', 'ODD'], default="Par04")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--plotEratio', action='store_true')
group.add_argument('--plotEcell', action='store_true')
group.add_argument('--plotLogEcell', action='store_true')
group.add_argument('--plotNumcells', action='store_true')
group.add_argument('--plotLong', action='store_true')
group.add_argument('--plotTrans', action='store_true')
group.add_argument('--plotIncident', action='store_true')
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
# datasets
if args.dataset == 2:
    legend_root = [args.detector+', ROOT dataset2']
    ifLowGran = True
    num_phi_low = 16
    num_z_low = 45
    num_r_low = 9
    size_phi_low = 2*pi/num_phi_low
    if args.detector == 'Par04':
        # Values true for Par04 SiW
        size_z_low = 3.4
        size_r_low = 4.65
    elif args.detector == 'ODD':
        # Values true for ODD
        size_z_low = 5.05
        size_r_low = 4.65
elif args.dataset == 3:
    legend_root = [args.detector+', ROOT dataset3']
    ifLowGran = False
    num_phi = 50
    num_z = 45
    num_r = 18
    size_phi = 2*pi/num_phi
    if args.detector == 'Par04':
        # Values true for Par04 SiW
        size_z = 3.4
        size_r = 2.325
    elif args.detector == 'ODD':
        # Values true for ODD
        size_z_low = 5.05
        size_r_low = 2.325


if args.legendRoot:
    legend_root = args.legendRoot
if args.legendH5:
    legend_h5 = args.legendH5
if len(legend_root) < len(filename_root) + len(filename_edm4hep):
    if args.legend:
        print("Use the same length for legend and files")
        exit(-1)
    else:
        for i in range(len(filename_root)+len(filename_edm4hep)):
            legend_root.append(legend_root[0])
if len(legend_h5) < len(filename_h5):
    if args.legend:
        print("Use the same length for legend and files")
        exit(-1)
    else:
        for i in range(len(filename_h5)):
            legend_h5.append(legend_h5[0])

colours=['g','b','r','o']
# files
df_list = []
events_list = []
for file_root in filename_root:
    if args.energyFilter:
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
    data = h5py.File(file_h5)
    randomIndices=random.sample(range(0,data[list(data.keys())[0]].shape[0]),num_events)
    h5_eMC_list.append(np.array([ data['incident_energy'][i] for i in randomIndices ]))
    found_angle_theta = False
    found_angle_phi = False
    if "incident_theta" in data.keys():
        h5_thetaMC_list.append(np.array([ data['incident_theta'][i] for i in randomIndices ]))
        found_angle_theta = True
    if "incident_phi" in data.keys():
        h5_phiMC_list.append(np.array([ data['incident_phi'][i] for i in randomIndices ]))
        found_angle_phi = True
    h5_showers_list.append(np.array([ data['showers'][i] for i in randomIndices ]))

def energy_ratio_root(fig, events, colour, label):
    ratios = []
    for energies, eMC in zip(events["EnergyCell"], events["EnergyMC"]):
        edep = np.sum(energies)
        print(edep, eMC, args.samplingFraction)
        ratios.append(edep/eMC/args.samplingFraction)
    return np.array(ratios).mean(), fig.hist(ratios, bins=np.linspace(0.,1.2,200), facecolor=colour, alpha=0.4, label=label, density = True)

def energy_ratio_hdf5(fig, energies, showers, colour, label):
    ratios = []
    print(showers.shape)
    for en, sh in zip(energies,showers):
        edep = np.sum(sh)
        ratios.append(edep/en)
    return np.array(ratios).mean(), fig.hist(ratios, bins=np.linspace(0.,1.2,200), facecolor=colour, alpha=0.4, label=label, density = True)

def num_cells_root(fig,events, colour, label):
    all_num = []
    for energies in events["EnergyCell"]:
        all_num.append(len(energies))
    if args.axisRange:
        mybins=np.linspace(args.axisRange[0], args.axisRange[1],50)
    else:
        mybins=np.linspace(1e2,2e5,50)
    return np.array(all_num).mean(), fig.hist(all_num, bins=mybins, facecolor=colour, alpha=0.4, label=label, density = True)

def num_cells_hdf5(fig,showers, colour, label):
    all_num = []
    for sh in showers:
        sh = sh[sh != 0]
        all_num.append(len(sh))
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

def cell_energy_log_hdf5(fig,showers, colour, label):
    all_logEn = []
    for sh in showers:
        sh = sh[sh != 0]
        sumG4 = np.sum(sh)
        for val in [log10(en) for en in sh]:
            all_logEn.append(val)
    return fig.hist(all_logEn, bins=np.linspace(-3,4,100), facecolor=colour, alpha=0.4, label=label, density = True)

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
    return fig.hist(all_phis, bins=np.linspace((-num_phi+0.5)*size_phi, (num_phi - 0.5)*size_phi,2*num_phi-1), weights=all_weights , facecolor=colour, alpha=0.5, label=label, density = True)

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
    return fig.hist(np.concatenate(events["zCell"]), bins=np.linspace(0,numz-1,numz), weights=np.concatenate(events["EnergyCell"]) , facecolor=colour, alpha=0.4, label=label, density = True)

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
    z_ids = np.linspace(0, numz-1, numz)
    energies = [np.sum(showers[:,:,:,int(i)]) for i in z_ids]
    return fig.hist(z_ids, bins=np.linspace(0,numz-1,numz), weights=energies , facecolor=colour, alpha=0.4, label=label, density = True)

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
    return fig.hist(np.concatenate(events["rhoCell"]), bins=np.linspace(0,numr-1,numr), weights=np.concatenate(events["EnergyCell"]) , facecolor=colour, alpha=0.4, label=label, density = True)

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
    r_ids = np.linspace(0, numr-1, numr)
    energies = [np.sum(showers[:,int(i)]) for i in r_ids]
    return fig.hist(r_ids, bins=np.linspace(0,numr-1,numr), weights=energies , facecolor=colour, alpha=0.4, label=label, density = True)

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
    return fig_first.hist(all_first_moments, bins=np.linspace(0, 0.4*numz*sizez, 32) , facecolor=colour, alpha=0.4, label=label, density = True), fig_second.hist(all_second_moments, bins=np.linspace(0, pow(numz*sizez,2) / 35., 32) , facecolor=colour, alpha=0.4, label=label, density = True)

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
    all_first_moments = []
    all_second_moments = []
    for sh in showers:
        z_ids = np.linspace(0, numz-1, numz)
        energies = [np.sum(sh[:,:,int(i)]) for i in np.linspace(0,numz-1,numz)]
        sumG4 = np.sum(energies)
        first = np.sum([z_ids[i] * sizez * energies[i] for i in range(len(z_ids))])
        first /= sumG4
        all_first_moments.append(first)
        second = np.sum([pow(z_ids[i] * sizez - first, 2) * energies[i] for i in range(len(z_ids))])
        second /= sumG4
        all_second_moments.append(second)
    return fig_first.hist(all_first_moments, bins=np.linspace(0, 0.4*numz*sizez, 32) , facecolor=colour, alpha=0.4, label=label, density = True), fig_second.hist(all_second_moments, bins=np.linspace(0, pow(numz*sizez,2) / 35., 32) , facecolor=colour, alpha=0.4, label=label, density = True)

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
    return fig_first.hist(all_first_moments, bins=np.linspace(0, 0.75*numr*sizer,128) , facecolor=colour, alpha=0.4, label=label, density = True), fig_second.hist(all_second_moments, bins=np.linspace(0, pow(numr*sizer,2) / 8., 128) , facecolor=colour, alpha=0.4, label=label, density = True)

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
    all_first_moments = []
    all_second_moments = []
    for sh in showers:
        r_ids = np.linspace(0, numr-1, numr)
        energies = [np.sum(sh[int(i)]) for i in np.linspace(0,numr-1,numr)]
        sumG4 = np.sum(energies)
        first = np.sum([r_ids[i] * sizer * energies[i] for i in range(len(r_ids))])
        first /= sumG4
        all_first_moments.append(first)
        second = np.sum([pow(r_ids[i] * sizer - first, 2) * energies[i] for i in range(len(r_ids))])
        second /= sumG4
        all_second_moments.append(second)
    return fig_first.hist(all_first_moments, bins=np.linspace(0, 0.75*numr*sizer,128) , facecolor=colour, alpha=0.4, label=label, density = True), fig_second.hist(all_second_moments, bins=np.linspace(0, pow(numr*sizer,2) / 8., 128) , facecolor=colour, alpha=0.4, label=label, density = True)

def plot_cell_energy():
    fig1, (ax1) = plt.subplots(1, 1, figsize=(9,7), gridspec_kw={'height_ratios': [1]})
    for it in range(len(events_list)):
        cell_energy_root(ax1,events_list[it],colours[it],legend_root[it])
    ax1.set_xlabel('E (MeV)')
    ax1.set_ylabel('Entries')
    ax1.set_title('Cell energy distribution')
    ax1.grid(True, which='both')
    ax1.legend()
    fig1.tight_layout()
    fig1.savefig(args.outfile+'cell_energy_lin.png')

def plot_cell_log_energy():
    if len(events_list)+len(h5_eMC_list) > 1:
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
    ax1.set_xlabel('log10(E//MeV)')
 #   ax1.set_ylim(bottom=1e-6)
    ax1.set_yscale('log')
    ax1.set_ylabel('Entries')
    ax1.set_title('Cell energy distribution')
    ax1.grid(True)
    ax1.legend()
    if len(events_list)+len(h5_eMC_list) > 1:
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
    if len(events_list)+len(h5_eMC_list) > 1:
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
    ax3.set_xlabel('number of cells')
    ax3.set_ylabel('Entries')
    ax3.set_title('Number of cells with E>0.5 keV')
    ax3.grid(True)
    ax3.legend()
    if len(events_list)+len(h5_eMC_list) > 1:
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
    fig3.savefig(args.outfile+'num_cells')

def plot_energy_ratio():
    if len(events_list)+len(h5_eMC_list) > 1:
        fig4, (ax4, ax4b) = plt.subplots(2, 1,figsize=(9, 12), gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig4, (ax4) = plt.subplots(1, 1,figsize=(9, 7), gridspec_kw={'height_ratios': [1]})
    val_of_bins_x1 = []
    edges_of_bins_x1 = []
    for it in range(len(events_list)):
        mean_eratio_root, (vals, edges, _) = energy_ratio_root(ax4,events_list[it],colours[it], legend_root[it])
        ax4.axvline(mean_eratio_root, color=colours[it], linestyle='dashed', linewidth=1, label=fr"<E> = {mean_eratio_root:.2f}")
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
    for it in range(len(h5_eMC_list)):
        mean_eratio_hdf5, (vals, edges, _) = energy_ratio_hdf5(ax4, h5_eMC_list[it], h5_showers_list[it], colours[it+len(events_list)], legend_h5[it])
        ax4.axvline(mean_eratio_hdf5, color=colours[it+len(events_list)], linestyle='dashed', linewidth=1, label=fr"<E> = {mean_eratio_hdf5:.2f}")
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
    ax4.set_xlabel('$E_{dep} / E_{MC}$')
    ax4.set_ylabel('Normalized entries')
    ax4.set_title('Ratio of energy deposited to incident')
    ax4.grid(True)
    ax4.legend()
    if len(events_list)+len(h5_eMC_list) > 1:
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
    fig4.savefig(args.outfile+'energy_ratio_small')


def plot_incident():
    if len(events_list)+len(h5_eMC_list) > 1:
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
        mybins=np.linspace(1e3,1e6,100)
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
    axE.legend()
    if len(events_list)+len(h5_eMC_list) > 1:
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
        axEr.legend()
        axEr.set_ylim(0.25,4)
    if found_angle_theta:
        #theta distribution
        for it in range(len(events_list)):
            (vals, edges, _) = axT.hist(events_list[it]["ThetaMC"], bins = np.linspace(0,pi,100), facecolor=colours[it], alpha=0.4,label=legend_root[it])
            val_of_bins_T.append(vals)
            edges_of_bins_T.append(edges)
        for it in range(len(h5_thetaMC_list)):
            (vals, edges, _) = axT.hist(h5_thetaMC_list[it], bins = np.linspace(0,pi,100), facecolor=colours[it+len(events_list)], alpha=0.4,label=legend_h5[it])
            val_of_bins_T.append(vals)
            edges_of_bins_T.append(edges)
        axT.set_xlabel('${\theta}_{MC}$ (rad)')
        axT.set_ylim(bottom=1)
        axT.set_yscale('log')
        axT.set_ylabel('Entries')
        axT.set_title('Incident energy distribution')
        axT.grid(True)
        axT.legend()
        if len(events_list)+len(h5_eMC_list) > 1:
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
            (vals, edges, _) = axP.hist(h5_phiMC_list[it], bins = np.linspace(-pi,pi,100), facecolor=colours[it+len(events_list)], alpha=0.4,label=legend_h5[it])
            val_of_bins_P.append(vals)
            edges_of_bins_P.append(edges)
        axP.set_xlabel('$\phi_{MC}$ (rad)')
        axP.set_ylim(bottom=1)
        axP.set_yscale('log')
        axP.set_ylabel('Entries')
        axP.set_title('Incident energy distribution')
        axP.grid(True)
        axP.legend()
        if len(events_list)+len(h5_eMC_list) > 1:
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
    fig.savefig(args.outfile+'incident_distr')

def plot_longitudinal():
    if len(events_list)+len(h5_eMC_list) > 1:
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
    ax5a.set_xlabel('z (ID)')
    ax5a.set_ylabel('$< E>$ (MeV)')
    ax5a.set_title('Longitudinal profile')
    ax5a.grid(True)
    ax5a.legend()
    fig5.tight_layout()
    ax5b.set_xlabel('$<\lambda>$ (mm)')
    ax5b.set_ylabel('Entries')
    ax5b.set_title('Longitudinal first moment')
    ax5c.set_xlabel('$<\lambda^{2}>$ (mm^{2})')
    ax5c.set_ylabel('Entries')
    ax5c.set_title('Longitudinal second moment')
    ax5b.grid(True)
    ax5c.grid(True)
    ax5b.legend()
    ax5c.legend()
    if len(events_list)+len(h5_eMC_list) > 1:
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
    fig5.savefig(args.outfile+'longitudinal_profile')

def plot_transverse():
    if len(events_list)+len(h5_eMC_list) > 1:
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
        (vals, edges, _) =trans_profile_hdf5(ax6a,h5_showers_list[it],colours[it+len(h5_eMC_list)],legend_h5[it], lowGran=ifLowGran)
        (vals2, edges2, _), (vals3, edges3, _) =trans_moments_hdf5(ax6b,ax6c, h5_showers_list[it],colours[it+len(h5_eMC_list)],legend_h5[it], lowGran=ifLowGran)
        val_of_bins_x1.append(vals)
        edges_of_bins_x1.append(edges)
        val_of_bins2_x1.append(vals2)
        edges_of_bins2_x1.append(edges2)
        val_of_bins3_x1.append(vals3)
        edges_of_bins3_x1.append(edges3)
    ax6a.set_xlabel('r (ID)')
    ax6a.set_ylabel('$< E >$ (MeV')
    ax6a.set_title('Transverse profile')
    ax6a.set_yscale('log')
    ax6a.grid(True)
    ax6b.set_xlabel('$<r>$ (mm)')
    ax6b.set_ylabel('Entries')
    ax6b.set_title('Transverse first moment')
    ax6c.set_xlabel('$< r^{2}>$(mm^{2})')
    ax6c.set_ylabel('Entries')
    ax6c.set_title('Transverse second moment')
    ax6b.grid(True)
    ax6c.grid(True)
    ax6a.legend()
    ax6b.legend()
    ax6c.legend()
    if len(events_list)+len(h5_eMC_list) > 1:
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
    fig6.savefig(args.outfile+'transverse_profile')

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
    if len(events_list)+len(h5_eMC_list) > 1:
        ax8b.set_ylabel('ratio')
        for it in range(1,len(events_list)):
            ratio = np.divide(val_of_bins_x1[it], val_of_bins_x1[0], where=(val_of_bins_x1[0] != 0))
            error = np.divide(val_of_bins_x1[it] * np.sqrt(val_of_bins_x1[0]) + val_of_bins_x1[0] * np.sqrt(val_of_bins_x1[it]),np.power(val_of_bins_x1[0], 2), where=(val_of_bins_x1[0] != 0))
            ax8b.plot(edges_of_bins_x1[it][:-1], ratio, alpha=0.4, color=colours[it],label=legend_root[it].split(",")[0]+"/"+legend_root[0].split(",")[0])
        ax8b.legend()
        ax8b.set_ylim(0.25,4)
    fig8.tight_layout()
    fig8.savefig(args.outfile+'incident_energy')

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
