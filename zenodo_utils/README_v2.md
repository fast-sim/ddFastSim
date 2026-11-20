
# LEMURS v2

Favaro, Giammanco and Krause reported a non-uniform incident-energy spectrum in
the smaller CLD v1 training files ([Appendix D, Fig. D1](https://arxiv.org/pdf/2601.05289#page=32)).
The selection changed between training parts:
parts 1–4 included every event, while parts 5–10 required exactly one saved
MC particle (`numMC==1`). In CLD, many high-energy events contained a saved
secondary produced in the forward region after the photon's expected ECAL
arrival. This is consistent with shower leakage but does not identify the
material where the interaction occurred. The multiplicity cut rejected these
events regardless of timing. It removed far fewer events in ALLEGRO.
Conversely, an upstream interaction could pass if no secondary was saved.

V2 uses the same selection for all CLD and ALLEGRO files. It rejects an event
if the primary photon ends clearly upstream of the ECAL, or if a saved secondary
is produced upstream more than 0.1 ns before the photon's estimated arrival.
The estimate uses the initial photon direction
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
For the testing files, use the same v2 converter and keep the energy and angles
in the output name. Each detector has eight points at 5 and 50 GeV,
φ = 0.0 and 0.2 rad, and θ = 1.57 and 2.1 rad.

## Validation plots

The release note `v2_log.pdf` illustrates the change for CLD and ALLEGRO and
will accompany the v2 files on [Zenodo](https://doi.org/10.5281/zenodo.17045561).
The repository plotter makes the two combined figures below. The release
note's forward-activity plot versus energy and angle is a separate diagnostic.
For the combined figures, use the production `.selection.npz` cache and the
original ROOT file at the path recorded in it. The current converter writes
only HDF5; it cannot recover rejected events from that file.

If you do not have the cache, regenerate it from the original ROOT input
with the archived converter. Use a fresh directory: the archived script
writes its own HDF5 and JSON report alongside the cache.

```bash
cache_dir=$(mktemp -d)
git show 798947d:scripts/edm4hep2h5_v2.py > "$cache_dir/converter.py"
python3 "$cache_dir/converter.py" --detector CLD -i "$INPUT" \
  -o "$cache_dir" --output-name part1_v2.h5
python3 zenodo_utils/plot_lemurs_v2.py "$cache_dir/part1_v2.selection.npz" \
  --outdir "$PLOTS/CLD/combined_part1"
```

Both figures compare Inclusive (filter applied in LEMURS-v1 for CLD part1-4), `numMC==1`
(filter applied in LEMURS-v1 for CLD part5-10) and LEMURS-v2 on the same simulated events.
`particle_locations_combined` maps primary endpoints, secondary vertices and
secondary endpoints (3×3). `selection_validation_combined` shows generated
spectra, selection efficiency and shower response (3×2). Grey shading marks the
Inclusive density; the vertical axes follow the data. Both figures are saved as 300-dpi PNGs.
