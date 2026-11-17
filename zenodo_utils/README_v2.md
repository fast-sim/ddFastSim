# LEMURS v2

## Why this revision was needed

The CLD v1 training files were not made with one consistent event selection.
Parts 1--4 had no `MCParticles` cut, while parts 5--10 required `numMC == 1`.
That cut can also remove photons which reach the ECAL normally but have extra
stored activity later in the event.

LEMURS v2 replaces this with the same selection for every part. The change is
also applied to ALLEGRO, although there is no claim that ALLEGRO had the same
part-to-part problem as CLD.

## What changed

The v2 converter reads the existing EDM4hep ROOT files; no showers are
resimulated. It requires the first stored particle to be the primary photon
from the interaction point, but it does not require `numMC == 1`.

An event is rejected when either:

1. the saved primary endpoint is still clearly upstream of the nominal ECAL
   entry; or
2. a stored secondary is produced in that upstream region before the nominal
   entry time.

The entry estimate uses the original photon direction and detector-specific
ECAL region envelopes. The geometry is frozen from `k4geo-FastSim.zip` at
commit `a25ddc8d9bcfd5610e389cda487d3f00172e9443`.

The converter in this repository now follows the usual ROOT `RDataFrame`
pattern: `Filter`, `Define`, and `AsNumpy`. It writes only the four HDF5
datasets: `showers`, `incident_energy`, `incident_theta`, and `incident_phi`,
stored as float32 with gzip level 9.

The completed production also contains JSON reports and small `.selection.npz`
files made during the production run. Those caches contain scalar diagnostics
for all original ROOT rows and were used for the plots; they do not contain
copies of the shower tensors. The JSON reports record the exact converter hash
used for production.

The selection is still a geometry/time approximation, not a measurement of the
track at ECAL entry. The converter currently reads `ECalBarrelCollection` only.

## Validation plots

The validation includes generated-energy and angular distributions, deposited
energy and response, primary endpoints, and selection efficiency. The old
`numMC == 1` rule and the v2 rule are applied to the same ROOT events so that
the effect of the selection can be seen directly. There are plots for each part
and for all parts combined.

The plots and converted files are available in the
[LEMURS Zenodo record](https://doi.org/10.5281/zenodo.17045562).

`zenodo_utils/plot_selection_v2.py` reproduces the selection plots from the
`.selection.npz` files. The usual `draw_all.py` and
`zenodo_utils/validate_lemurs.py` checks were also run on every HDF5 part.

## Reproducing the conversion

Run in the key4hep environment used for the ROOT files:

```bash
python3 scripts/edm4hep2h5_v2.py --detector CLD \
  -i "$INPUT" -o "$OUTPUT/CLD" --output-name part1_v2.h5 --skip-existing
```

Use `--detector ALLEGRO` for ALLEGRO. Tests which do not require a ROOT
installation can be run with:

```bash
python3 -m unittest discover -s tests -p 'test_lemurs_v2.py'
```
