#!/bin/bash

### First some configuration settings that need to changed for each user
export FCCcore=/afs/cern.ch/work/a/azaborow/FCC/
OUTPATH=/eos/geant4/fastSim/ddfastsim/ODD/dataset2_1GeV1TeVFlat_theta0p87to2p27_phiFull/
ENERGY=$1
NUMEVENTS=$2
PARTICLE=$4
CLUSTERID=$5
NODEID=$6
###

source $FCCcore/ddfastsim/install/bin/thisDDFastSim.sh
source $FCCcore/OpenDataDetector/install/bin/this_odd.sh
source /cvmfs/sw.hsf.org/key4hep/releases/2023-11-23/x86_64-almalinux9-gcc11.3.1-opt/key4hep-stack/2023-11-30-gyuooo/setup.sh

echo "Steering file ================"
cat $FCCcore/ddfastsim/options/ODD/ODD_ddsim_steer_mesh_gps.py
cat job.${PARTICLE}.N${NUMEVENTS}.mac

echo "Simulating for energy" $ENERGY " GeV and for " ${NUMEVENTS} " showers of " ${PARTICLE}

ddsim --steeringFile $FCCcore/ddfastsim/options/ODD/ODD_ddsim_steer_mesh_gps.py --numberOfEvents ${NUMEVENTS} --macroFile job.${PARTICLE}.${ENERGY}GeV.N${NUMEVENTS}.mac --outputFile ddsim_mesh_ODD_${PARTICLE}_${NUMEVENTS}events_1GeV1TeV_GPSflat_edm4hep.root

echo "Copying the output file to " $OUTPATH
mkdir -p $OUTPATH/
cp  ddsim_mesh_ODD_${PARTICLE}_${NUMEVENTS}events_1GeV1TeV_GPSflat_edm4hep.root $OUTPATH/ddsim_mesh_ODD_${PARTICLE}_${NUMEVENTS}events_1GeV1TeV_GPSflat_edm4hep_${CLUSTERID}.${NODEID}.root
cp analysis.root $OUTPATH/analysis_ddsim_mesh_ODD_${PARTICLE}_${NUMEVENTS}events_1GeV1TeV_edm4hep_${CLUSTERID}.${NODEID}.root
# Copy also the random seeds to be able to reproduce the run
mkdir -p $OUTPATH/randomSeeds/${CLUSTERID}.${NODEID}/
cp *.rndm $OUTPATH/randomSeeds/${CLUSTERID}.${NODEID}/.
