
# LEMURS v2

Favaro, Giammanco and Krause reported a non-uniform incident-energy spectrum in
the smaller CLD v1 training files ([Appendix D, Fig. D1](https://arxiv.org/pdf/2601.05289#page=32)).
We traced it to a non uniform translation change in selection of showers:
parts 1–4 included every event, while parts 5–10 required exactly one saved
MC particle (`numMC==1`). That cut can discard a shower even when the photon
reaches the ECAL just because leakage produces a secondary elsewhere afterwards, and in case
of CLD that leakage was rather significant because of the detector placed next to calorimetric endcap.
In case of ALLEGRO the change is minimal. It also missed some early interactions that should flag the event
for removal, but since DD4hep did not store all the secondaries, it passed the filter of just one MC particle.

V2 applies a new rule to every file in order to regain the maximum available statistics and remove rare air interactions
that could have altered the depositions in the shower grid.
It rejects events if the primary photon (originated at IP) ends upstream of the ECAL, or if a saved secondary
is produced upstream before the photon's estimated arrival. The estimate uses the initial photon direction
and the CLD or ALLEGRO ECAL envelope from `k4geo-FastSim`
(`a25ddc8d9bcfd5610e389cda487d3f00172e9443`). The Par04 samples do not need this veto.

## Conversion

From the repository root, in the key4hep environment for the input ROOT file:

```bash
python3 scripts/edm4hep2h5_v2.py --detector CLD \
  -i "$INPUT" -o "$OUTPUT/CLD" --output-name part1_v2.h5
```

Use `--detector ALLEGRO` for ALLEGRO. The original
[edm4hep2h5.py](../scripts/edm4hep2h5.py) implements the multiplicity cut.

## Validation plots

The release note `v2_log.pdf` illustrates the change for CLD and ALLEGRO and
will accompany the v2 files on [Zenodo](https://doi.org/10.5281/zenodo.17045561).
Reproducing its figures requires the production `.selection.npz` cache and
the original ROOT file named in it as it iterates over MC particles. The converter writes only HDF5.

```bash
python3 zenodo_utils/plot_lemurs_v2.py "$DATA/CLD/part1_v2.selection.npz" \
  --outdir "$PLOTS/CLD/combined_part1"
```

Both figures compare Inclusive (filter applied in LEMURS-v1 for CLD part1-4), `numMC==1`
(filter applied in LEMURS-v1 for CLD part5-10) and LEMURS-v2 on the same simulated events.
`particle_locations_combined` maps primary endpoints, secondary vertices and
secondary endpoints (3×3). `selection_validation_combined` shows generated
spectra, selection efficiency and shower response (3×2). Both are saved as PDF and PNG.
