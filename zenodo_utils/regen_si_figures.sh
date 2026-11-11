#!/usr/bin/env bash
# Regenerate the Supplementary Information plots (voxel energy and nonzero-voxel
# distributions for every testing file) with the patched draw_all.py, under the
# file names supplementary.tex expects, e.g. images_si/ODD_50vs500GeV_phi0.0_theta1.57.cell_energy_log.png
#
# Usage:  bash regen_si_figures.sh [h5_dir] [output_dir]
#         defaults: ~/LEMURS/data and ~/LEMURS/paper_lemurs/audit/images_si
set -euo pipefail
DATA=${1:-$HOME/LEMURS/data}
OUT=${2:-$HOME/LEMURS/paper_lemurs/audit/images_si}
DRAW="python3 $(dirname "$0")/draw_all.py"
mkdir -p "$OUT"

locate() {
  local hit; hit=$(find "$DATA" -name "$1" -print -quit)
  [[ -n "$hit" ]] || { echo "missing below $DATA: $1" >&2; exit 1; }
  echo "$hit"
}

made=0
for det in Par04SiW Par04SciPb ODD FCCeeCLD FCCeeALLEGRO; do
  if [[ $det == FCCee* ]]; then hi=50; lo=5; tag="5vs50GeV"; else hi=500; lo=50; tag="50vs500GeV"; fi
  for phi in 0.0 0.2; do
    for th in 1.57 2.1; do
      fhi=$(locate "LEMURS_${det}_gamma_1000events_${hi}GeV_phi${phi}_theta${th}.h5")
      flo=$(locate "LEMURS_${det}_gamma_1000events_${lo}GeV_phi${phi}_theta${th}.h5")
      for flag in --plotLogEcell --plotNumcells; do
        extra=""; [[ $flag == --plotNumcells ]] && extra="--axisRange 0 6480"
        echo "+ $det phi=$phi theta=$th $flag"
        $DRAW -d "$det" $flag --dataset 2 --noRatio -n 1000 $extra \
              --infileH5 "$fhi" "$flo" --legendH5 "${hi} GeV" "${lo} GeV" \
              -o "$OUT/${det}_${tag}_phi${phi}_theta${th}." > /dev/null
        made=$((made + 1))
      done
    done
  done
done
echo "Done: $made plots in $OUT (expected 40). Copy them over images/ and send them to me."
