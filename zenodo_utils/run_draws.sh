#!/usr/bin/env bash
# run_draws.sh
# Hardcoded runner for training & validation/testing plots.
# Supports --dry-run flag to only print commands without executing.

set -euo pipefail

# -----------------------------
# CONFIG (hardcoded)
# -----------------------------
declare -a DETS=("Par04SiW" "Par04SciPb" "ODD" "FCCeeCLD" "FCCeeALLEGRO")

# -d flag for draw_all.py (same for all, per your example)
DFLAG="Par04"

# Legends (pretty names per detector)
declare -A LEGEND_PREFIX=(
  ["Par04SiW"]="Par04 SiW"
  ["Par04SciPb"]="Par04 SciPb"
  ["ODD"]="ODD"
  ["FCCeeCLD"]="FCCee CLD"
  ["FCCeeALLEGRO"]="FCCee ALLEGRO"
)

# Output directory for figures
FIGDIR="fig"

# draw_all.py location (assumed current dir)
DRAW="python draw_all.py"

# Dry run mode? (default: off)
DRYRUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
  DRYRUN=1
  shift
fi

# -----------------------------
# helpers
# -----------------------------
echo_run () {
  local cmd="$*"
  echo ">> $cmd"
  if [[ $DRYRUN -eq 0 ]]; then
    eval "$cmd"
  fi
}

ensure_figdir () {
  mkdir -p "$FIGDIR"
}

legend_for_training () {
  local det="$1"
  local file="$2"
  local part
  part=$(grep -oE 'part[0-9]+' <<<"$(basename "$file")" || true)
  if [[ -z "$part" ]]; then part="part?"; fi
  part="${part/part/part }"
  echo "${LEGEND_PREFIX[$det]}, ${part}"
}

outfile_for_training () {
  local det="$1"
  local file="$2"
  local part
  part=$(grep -oE 'part[0-9]+' <<<"$(basename "$file")" || true)
  if [[ -z "$part" ]]; then part="part"; fi
  echo "${FIGDIR}/${det}_${part}_"
}

phi_theta_token () {
  local file="$1"
  grep -oE 'phi[0-9.]+_theta[0-9.]+' <<<"$(basename "$file")" || true
}

# -----------------------------
# training
# -----------------------------
run_training () {
  local det
  for det in "${DETS[@]}"; do
    local dir="to_upload/${det}/training"
    [[ -d "$dir" ]] || continue

    shopt -s nullglob
    local files=("$dir"/*.h5)
    shopt -u nullglob

    for f in "${files[@]}"; do
      local legend outfile
      legend="$(legend_for_training "$det" "$f")"
      outfile="$(outfile_for_training "$det" "$f")"

      echo_run "$DRAW -d \"$DFLAG\" --plotIncident --dataset 2 \
        --infileH5 \"$f\" \
        --legendH5 \"$legend\" \
        --outfile \"$outfile\" \
        -n 100000"
    done
  done
}

# -----------------------------
# testing
# - FCCee* detectors:      5 vs 50 GeV
# - other detectors:       50 vs 500 GeV (original behavior)
# -----------------------------
run_testing () {
  local det
  for det in "${DETS[@]}"; do
    local dir="to_upload/${det}/testing/"
    [[ -d "$dir" ]] || continue

    # Decide energy pair based on detector name
    local LO_GE V HI_GEV OUTPAIR
    if [[ "$det" == FCCee* ]]; then
      LO_GEV=5
      HI_GEV=50
      OUTPAIR="5vs50GeV"
    else
      LO_GEV=50
      HI_GEV=500
      OUTPAIR="50vs500GeV"
    fi

    shopt -s nullglob
    # Anchor iteration on the HIGH-energy files so we can derive the matching LOW-energy filename
    local fHI
    for fHI in "$dir"/LEMURS_"$det"_gamma_1000events_"${HI_GEV}"GeV_phi*_theta*.h5; do
      local phi_theta fLO outname
      phi_theta="$(phi_theta_token "$fHI")"
      [[ -n "$phi_theta" ]] || continue

      local base_HI base_LO
      base_HI="$(basename "$fHI")"
      base_LO="${base_HI/${HI_GEV}GeV_/${LO_GEV}GeV_}"

      fLO="${dir}/${base_LO}"
      if [[ ! -f "$fLO" ]]; then
        # allow extra IDs/suffixes after the nominal name
        local prefix="${base_LO%.h5}"
        local candidates=("$dir"/"$prefix"*.h5)
        if [[ ${#candidates[@]} -gt 0 && -f "${candidates[0]}" ]]; then
          fLO="${candidates[0]}"
        else
          echo "[WARN] No matching ${LO_GEV} GeV for: $fHI"
          continue
        fi
      fi

      outname="${FIGDIR}/${det}_${OUTPAIR}_${phi_theta}"

      echo_run "$DRAW -d \"$DFLAG\" --plotLogEcell --dataset 2 \
        --infileH5 \"$fHI\" \"$fLO\" \
        --legendH5 \"${LEGEND_PREFIX[$det]}, ${HI_GEV} GeV photons\" \"${LEGEND_PREFIX[$det]}, ${LO_GEV} GeV photons\" \
        --noRatio \
        --outfile \"$outname\" \
        -n 1000 \
        --axisRange 0 6480"
    done
    shopt -u nullglob
  done
}

# -----------------------------
# main
# -----------------------------
ensure_figdir
run_training
run_testing

echo "All done. (dry run = $DRYRUN)"
