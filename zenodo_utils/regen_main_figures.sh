#!/usr/bin/env bash
# Regenerate the main-text figures (Figs. 9-14) of the LEMURS paper with the
# patched draw_all.py: larger fonts, fixed profile binning, nonzero-voxel title,
# per-detector voxel widths, finer Fig. 13 moment range.
#
# Usage:  bash regen_main_figures.sh [h5_dir] [output_dir]
#         defaults: ~/LEMURS/data and ~/LEMURS/paper_lemurs/audit/images_new
# Output file names match the \includegraphics names in main.tex, so the new
# PNGs can replace the ones in images/ directly.
set -euo pipefail

DATA=${1:-$HOME/LEMURS/data}
OUT=${2:-$HOME/LEMURS/paper_lemurs/audit/images_new}
DRAW="python3 $(dirname "$0")/draw_all.py"
mkdir -p "$OUT"

# find each file anywhere below DATA (part1 files at the top, testing files in <DETECTOR>/testing/)
locate() {
  local hit; hit=$(find "$DATA" -name "$1" -print -quit)
  [[ -n "$hit" ]] || { echo "missing below $DATA: $1" >&2; exit 1; }
  echo "$hit"
}
SIW_TRAIN=$(locate LEMURS_Par04SiW_gamma_100kEvents_1GeV1TeV_GPSflat_part1.h5)
SIW_500=$(locate LEMURS_Par04SiW_gamma_1000events_500GeV_phi0.0_theta1.57.h5)
SIW_50=$(locate LEMURS_Par04SiW_gamma_1000events_50GeV_phi0.0_theta1.57.h5)
ALG_50=$(locate LEMURS_FCCeeALLEGRO_gamma_1000events_50GeV_phi0.0_theta1.57.h5)
ALG_5=$(locate LEMURS_FCCeeALLEGRO_gamma_1000events_5GeV_phi0.0_theta1.57.h5)
echo "inputs:"; printf "  %s\n" "$SIW_TRAIN" "$SIW_500" "$SIW_50" "$ALG_50" "$ALG_5"

run() { echo "+ $*"; eval "$@"; }

# Figs. 9 and 10: first training file of Par04 SiW, all 100k showers (one sample: no legend,
# the caption names the file); nonzero-voxel axis limited to the 9 x 16 x 45 = 6480 voxels
for flag in --plotIncident --plotLogEcell --plotNumcells; do
  extra=""; [[ $flag == --plotNumcells ]] && extra="--axisRange 0 6480"
  run $DRAW -d Par04SiW $flag --dataset 2 --noRatio -n 100000 $extra \
      --infileH5 "\"$SIW_TRAIN\"" --legendH5 "\"\"" \
      -o "\"$OUT/Par04_SiW_part1\""
done

# Figs. 11-13 and 14a: Par04 SiW testing files, 500 vs 50 GeV, theta = 1.57, phi = 0
for flag in --plotLogEcell --plotNumcells --plotLong --plotTrans --plotEratio; do
  extra=""; [[ $flag == --plotNumcells ]] && extra="--axisRange 0 6480"
  run $DRAW -d Par04SiW $flag --dataset 2 --noRatio -n 1000 $extra \
      --infileH5 "\"$SIW_500\"" "\"$SIW_50\"" \
      --legendH5 "\"500 GeV\"" "\"50 GeV\"" \
      -o "\"$OUT/Par04_SiW_50vs500GeV\""
done

# Fig. 14b: FCC-ee ALLEGRO testing files, 50 vs 5 GeV
run $DRAW -d FCCeeALLEGRO --plotEratio --dataset 2 --noRatio -n 1000 \
    --infileH5 "\"$ALG_50\"" "\"$ALG_5\"" \
    --legendH5 "\"50 GeV\"" "\"5 GeV\"" \
    -o "\"$OUT/FCCee_ALLEGRO_5vs50GeV\""

echo
echo "Done. Main text uses these files from $OUT/:"
for n in Par04_SiW_part1incident_distr Par04_SiW_part1cell_energy_log Par04_SiW_part1num_cells \
         Par04_SiW_50vs500GeVcell_energy_log Par04_SiW_50vs500GeVnum_cells \
         Par04_SiW_50vs500GeVlongitudinal_profile Par04_SiW_50vs500GeVtransverse_profile \
         Par04_SiW_50vs500GeVenergy_ratio_small FCCee_ALLEGRO_5vs50GeVenergy_ratio_small; do
  if [[ -f "$OUT/$n.png" ]]; then echo "  ok      $n.png"; else echo "  MISSING $n.png"; fi
done
echo "Copy them over images/ in the paper source and send them to me for the diff."
