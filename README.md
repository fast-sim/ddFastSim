# ddFastSim

This repository contains plugins for fast simulation with DDG4.

To start with, a (so far partial) replication of Geant4 extended/parameterisations/Par04 example is presented.
It demonstrates how to define a very simple fast shower model and how to attach it to the region of parameterisation.

## Getting started

Get all needed packages (Geant4, DD4hep), e.g. from:
```
source /cvmfs/sw.hsf.org/key4hep/setup.sh
```

Then compile this library that contains few plugins for DDG4 simulation:
```
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install
source ../install/bin/thisDDFastSim.sh
```