#!/bin/bash

### First some configuration settings that need to changed for each user
export FCCcore=/afs/cern.ch/work/a/azaborow/FCC/
OUTPATH=/eos/geant4/fastSim/ddodd/Par04/test/
ENERGY=$1
NUMEVENTS=$2
PARTICLE=$4
CLUSTERID=$5
NODEID=$6
###

source $FCCcore/ddfastsim/install/bin/thisDDFastSim.sh
source /cvmfs/sw.hsf.org/key4hep/setup.sh

echo "Steering file ================"
cat $FCCcore/ddfastsim/options/Par04/Par04_ddsim_steer_pgun.py
cat job.${PARTICLE}.${ENERGY}GeV.N${NUMEVENTS}.mac

echo "Simulating for energy" ${ENERGY} " GeV and for " ${NUMEVENTS} " showers of " ${PARTICLE}

ddsim --steeringFile $FCCcore/ddfastsim/options/Par04/Par04_ddsim_steer_pgun.py --gun.energy ${ENERGY}*GeV --gun.part $3 --numberOfEvents ${NUMEVENTS} --macroFile job.${PARTICLE}.${ENERGY}GeV.N${NUMEVENTS}.mac --outputFile ddsim_mesh_Par04_${PARTICLE}_${NUMEVENTS}events_${ENERGY}GeV_edm4hep.root

echo "Copying the output file to " $OUTPATH
mkdir -p $OUTPATH/
cp  ddsim_mesh_Par04_${PARTICLE}_${NUMEVENTS}events_${ENERGY}GeV_edm4hep.root $OUTPATH/ddsim_mesh_Par04_${PARTICLE}_${NUMEVENTS}events_${ENERGY}GeV_edm4hep_${CLUSTERID}.${NODEID}.root
cp analysis.root $OUTPATH/analysis_ddsim_mesh_Par04_${PARTICLE}_${NUMEVENTS}events_${ENERGY}GeV_edm4hep_${CLUSTERID}.${NODEID}.root
# Copy also the random seeds to be able to reproduce the run
mkdir -p $OUTPATH/randomSeeds/${CLUSTERID}.${NODEID}/
cp *.rndm $OUTPATH/randomSeeds/${CLUSTERID}.${NODEID}/.
