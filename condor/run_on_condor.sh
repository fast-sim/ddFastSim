#!/bin/bash

### First some configuration settings that need to changed for each user
export FCCcore=/afs/cern.ch/work/a/azaborow/FCC/
OUTPATH=/eos/geant4/fastSim/ddodd/Par04/test/
###

source $FCCcore/ddfastsim/install/bin/thisDDFastSim.sh
source /cvmfs/sw.hsf.org/key4hep/setup.sh

echo "Steering file ================"
cat $FCCcore/ddfastsim/options/Par04/Par04_ddsim_steer_pgun.py
cat job.$4.$1GeV.N$2.mac

echo "Simulating for energy" $1 " GeV and for " $2 " showers of " $3

ddsim --steeringFile $FCCcore/ddfastsim/options/Par04/Par04_ddsim_steer_pgun.py --gun.energy $1*GeV --gun.part $3 --numberOfEvents $2 --macroFile job.$4.$1GeV.N$2.mac --outputFile ddsim_mesh_Par04_$4_$2events_$1GeV_edm4hep.root

echo "Copying the output file to " $OUTPATH
mkdir -p $OUTPATH/
cp  ddsim_mesh_Par04_$4_$2events_$1GeV_edm4hep.root $OUTPATH/ddsim_mesh_Par04_$4_$2events_$1GeV_edm4hep_$5.$6.root
cp analysis.root $OUTPATH/analysis_ddsim_mesh_Par04_$4_$2events_$1GeV_edm4hep_$5.$6.root
# Copy also the random seeds to be able to reproduce the run
mkdir -p $OUTPATH/randomSeeds/$5.$6/
cp *.rndm $OUTPATH/randomSeeds/$5.$6/.
