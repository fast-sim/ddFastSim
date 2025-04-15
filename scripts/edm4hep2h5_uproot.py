'''
New conversion for EDM4HEP ROOT file to HDF5 file, for key4hep stack 2025-01-28 (using RDataFrame is broken :( ))
Converts an EDM4HEP ROOT file to an HDF5 file, saving the shower energy in a 3D array and shower data in a 2D array.
Units: Energy values are stored in MeV and angles are stored in radians.
'''

#!/bin/env/python
import sys

##### For appending to Python Path when awkward is broken in Key4hep stack (2025-01-28)
sys.path.append('/usr/lib/python3/dist-packages')

import argparse
import numpy as np
import os
import h5py
import ROOT
import uproot
import awkward as ak

def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--outputDir", '-o', type=str, default="./", help="Path to the output directory")
    p.add_argument("--inputFile", '-i', type=str, required=True, help="Name of the EDM4hep file to translate")
    p.add_argument("--numR", type=int, default=9, help="Number of cells in R")
    p.add_argument("--numPhi", type=int, default=16, help="Number of cells in phi")
    p.add_argument("--numZ", type=int, default=45, help="Number of cells in z")
    p.add_argument("--samplingFraction", type=float, default=1., help="Sampling fraction to use to rescale energy. Defined as f=active/(active+absorber)")
    args = p.parse_args()
    return args

def main(argv):
    # Parse commandline arguments
    args = parse_args(argv)
    input_file = args.inputFile
    output_dir = args.outputDir
    # Number of cells in the r, phi & z directions
    num_cells_R = args.numR
    num_cells_phi = args.numPhi
    num_cells_z = args.numZ
    # Sampling fraction that rescales energy of each cell
    sampling_fraction = args.samplingFraction
    if os.stat(input_file).st_size > 0:
        h5_file = h5py.File(
            f"{output_dir}/{os.path.splitext(os.path.basename(input_file))[0]}.h5", "w"
        )
        print(f"Creating output file {output_dir}/{os.path.splitext(os.path.basename(input_file))[0]}.h5")
        # Read EDM4HEP root file
        events = uproot.open(input_file)["events"]
        events.show()
        energy_cell = events["ECalBarrelCollection/ECalBarrelCollection.energy"].array(library="ak")*1e3
        rhoCell = events["ECalBarrelCollection/ECalBarrelCollection.position.x"].array(library="ak")
        zCell = events["ECalBarrelCollection/ECalBarrelCollection.position.z"].array(library="ak")
        phiCell = events["ECalBarrelCollection/ECalBarrelCollection.position.y"].array(library="ak")
        
        #### Assuming one MC Particle!!!
        MCParticles_momx = events["MCParticles/MCParticles.momentum.x"].array(library="ak")
        MCParticles_momy = events["MCParticles/MCParticles.momentum.y"].array(library="ak")
        MCParticles_momz = events["MCParticles/MCParticles.momentum.z"].array(library="ak")

        print(len(MCParticles_momx))     

        energy_MC = np.sqrt(MCParticles_momx**2 + MCParticles_momy**2 + MCParticles_momz**2)
        #### NOTE: SKIPPED OTHER CONDITIONTS FOR NOW- will probably want to use local direction anyhow!!!

        ### This is the energy of the first MC Particle (i.e. the one you want)
        print(energy_MC[:, 0])
        print(len(energy_cell[1]))

        print(len(rhoCell[1]))

        #energy_MC = events["MCParticles/MCParticles.energy"].array(library="ak")*1.e3

        num_showers = len(energy_MC)
        print('num_showers: ', num_showers)

        all_showers = []
        # loop over events
        for event in range(num_showers):
            # Initialize a 3D array with shape nb_events, nb_cells in x,y,z
            shower = np.zeros((num_cells_R, num_cells_phi, num_cells_z))
            for cell in range(len(rhoCell[event])):
                if (
                    (int(rhoCell[event][cell]) < num_cells_R)
                    and (int(phiCell[event][cell]) < num_cells_phi)
                    and (int(zCell[event][cell]) < num_cells_z)
                ):
                    shower[int(rhoCell[event][cell])][
                        int(phiCell[event][cell])][
                        int(zCell[event][cell])
                    ] = energy_cell[event][cell] / sampling_fraction
            all_showers.append(shower)
        # Save dataset
        print(f"Creating datasets with shape 3x{len(energy_MC)} and {np.shape(all_showers)} ")
        h5_file.create_dataset("incident_energy", data=energy_MC[:, 0], compression="gzip", compression_opts=9, dtype='f4')
        #h5_file.create_dataset("incident_phi", data=phiMC, compression="gzip", compression_opts=9, dtype='f4')
        #h5_file.create_dataset("incident_theta", data=thetaMC, compression="gzip", compression_opts=9, dtype='f4')
        h5_file.create_dataset("showers", data=all_showers, compression="gzip", compression_opts=9, dtype='f4')
        h5_file.close()


if __name__ == "__main__":
    exit(main(sys.argv[1:]))