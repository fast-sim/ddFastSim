#!/bin/bash

### First some configuration settings that need to changed for each user
export FCCcore=/afs/cern.ch/work/a/azaborow/FCC/
OUTPATH=/eos/geant4/fastSim/ddfastsim/Par04PbWO4/dataset2_PbWO4_discrete/
ENERGY=$1
NUMEVENTS=$2
NONESCAPEDPARTICLE=$3
PARTICLE=$4
CLUSTERID=$5
NODEID=$6
PHI=$7
THETA=$8
###

source $FCCcore/ddfastsim/install/bin/thisDDFastSim.sh
source /cvmfs/sw.hsf.org/key4hep/releases/2023-11-23/x86_64-almalinux9-gcc11.3.1-opt/key4hep-stack/2023-11-30-gyuooo/setup.sh

echo "Steering file ================"
cat $FCCcore/ddfastsim/options/Par04_PbWO4/Par04_ddsim_steer_dataset2_mesh_pgun.py
cat job.${PARTICLE}.${ENERGY}GeV.N${NUMEVENTS}.mac

echo "Simulating for energy" $ENERGY " GeV and for " ${NUMEVENTS} " showers of " ${PARTICLE} " for phi = " ${PHI} " and theta = " ${THETA}

ddsim --steeringFile $FCCcore/ddfastsim/options/Par04_PbWO4/Par04_ddsim_steer_dataset2_mesh_pgun.py --gun.energy ${ENERGY}*GeV --gun.part ${NONESCAPEDPARTICLE} --gun.phiMin $PHI --gun.phiMax $PHI --gun.thetaMin $THETA --gun.thetaMax $THETA --numberOfEvents ${NUMEVENTS} --macroFile job.${PARTICLE}.${ENERGY}GeV.N${NUMEVENTS}.mac --outputFile ddsim_mesh_Par04_${PARTICLE}_${NUMEVENTS}events_${ENERGY}GeV_phi${PHI}_theta${THETA}_edm4hep.root

echo "Copying the output file to " $OUTPATH
mkdir -p $OUTPATH/
cp  ddsim_mesh_Par04_${PARTICLE}_${NUMEVENTS}events_${ENERGY}GeV_phi${PHI}_theta${THETA}_edm4hep.root $OUTPATH/ddsim_mesh_Par04_${PARTICLE}_${NUMEVENTS}events_${ENERGY}GeV_phi${PHI}_theta${THETA}_edm4hep_${CLUSTERID}.${NODEID}.root
cp analysis.root $OUTPATH/analysis_ddsim_mesh_Par04_${PARTICLE}_${NUMEVENTS}events_${ENERGY}GeV_phi${PHI}_theta${THETA}_edm4hep_${CLUSTERID}.${NODEID}.root
# Copy also the random seeds to be able to reproduce the run
mkdir -p $OUTPATH/randomSeeds/${CLUSTERID}.${NODEID}/
cp *.rndm $OUTPATH/randomSeeds/${CLUSTERID}.${NODEID}/.
