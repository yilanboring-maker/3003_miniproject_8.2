#!/usr/bin/env bash
set -Eeuo pipefail

RUN_ID="${RUN_ID:?RUN_ID is required}"
CANDIDATE_ID="${CANDIDATE_ID:-$RUN_ID}"
STAGE_NS="${STAGE_NS:?STAGE_NS is required, e.g. 20}"
BASE_DIR="${BASE_DIR:-/tmp/ica2_runs}"
SCRIPT_DIR="${SCRIPT_DIR:-$BASE_DIR/scripts}"
WORK_DIR="${WORK_DIR:-$BASE_DIR/output/$RUN_ID}"
RESULTS_DIR="${RESULTS_DIR:-$BASE_DIR/results/$RUN_ID}"
LOG_DIR="${LOG_DIR:-$BASE_DIR/logs}"
GMX="${GMX:-gmx}"
LOG_FILE="$LOG_DIR/${RUN_ID}_stage_${STAGE_NS}ns_analysis.log"

TPR="md_${STAGE_NS}ns.tpr"
XTC="md_${STAGE_NS}ns_center.xtc"
PREFIX="stage_${STAGE_NS}ns"

mkdir -p "$RESULTS_DIR/analysis/$PREFIX" "$RESULTS_DIR/figures/$PREFIX" "$LOG_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

cd "$WORK_DIR"
echo "[analysis] started_at=$(date)"
echo "[analysis] run_id=$RUN_ID candidate_id=$CANDIDATE_ID stage=${STAGE_NS}ns"
echo "[analysis] work_dir=$WORK_DIR"
echo "[analysis] results_dir=$RESULTS_DIR"
echo "[analysis] tpr=$TPR xtc=$XTC index=index_ab_ag.ndx"

for f in "$TPR" "$XTC" index_ab_ag.ndx; do
  if [[ ! -s "$f" ]]; then
    echo "[analysis:error] missing $f" >&2
    exit 2
  fi
done

ANALYSIS_WORK="$WORK_DIR/analysis_$PREFIX"
mkdir -p "$ANALYSIS_WORK"

echo "[analysis] RMSD"
printf "5\n5\n" | "$GMX" rms -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$ANALYSIS_WORK/rmsd_complex.xvg" -tu ns
printf "3\n3\n" | "$GMX" rms -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$ANALYSIS_WORK/rmsd_antibody.xvg" -tu ns
printf "4\n4\n" | "$GMX" rms -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$ANALYSIS_WORK/rmsd_antigen.xvg" -tu ns

echo "[analysis] RMSF"
echo "8" | "$GMX" rmsf -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -res -o "$ANALYSIS_WORK/rmsf_complex_ca.xvg"
echo "6" | "$GMX" rmsf -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -res -o "$ANALYSIS_WORK/rmsf_antibody_ca.xvg"
echo "7" | "$GMX" rmsf -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -res -o "$ANALYSIS_WORK/rmsf_antigen_ca.xvg"

echo "[analysis] Rg"
echo "2" | "$GMX" gyrate -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$ANALYSIS_WORK/rg_complex.xvg"
echo "0" | "$GMX" gyrate -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$ANALYSIS_WORK/rg_antibody.xvg"
echo "1" | "$GMX" gyrate -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$ANALYSIS_WORK/rg_antigen.xvg"

echo "[analysis] SASA"
"$GMX" sasa -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$ANALYSIS_WORK/sasa_complex.xvg" -surface 'group "Complex"' -output 'group "Complex"'
"$GMX" sasa -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$ANALYSIS_WORK/sasa_antibody.xvg" -surface 'group "Antibody"' -output 'group "Antibody"'
"$GMX" sasa -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$ANALYSIS_WORK/sasa_antigen.xvg" -surface 'group "Antigen"' -output 'group "Antigen"'

echo "[analysis] H-bonds"
printf "0\n1\n" | "$GMX" hbond -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -num "$ANALYSIS_WORK/hbonds_antibody_antigen.xvg" -dist "$ANALYSIS_WORK/hbonds_distance.xvg" -ang "$ANALYSIS_WORK/hbonds_angle.xvg" || true

echo "[analysis] Interface contacts"
printf "0\n1\n" | "$GMX" mindist -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -od "$ANALYSIS_WORK/interface_mindist.xvg" -on "$ANALYSIS_WORK/interface_contacts.xvg" -d 0.45 -group || true

cp -f "$ANALYSIS_WORK"/*.xvg "$RESULTS_DIR/analysis/$PREFIX/" 2>/dev/null || true
if ! python "$SCRIPT_DIR/72_plot_md_analysis_generic.py" \
  --analysis-dir "$RESULTS_DIR/analysis/$PREFIX" \
  --output-dir "$RESULTS_DIR/figures/$PREFIX" \
  --summary "$RESULTS_DIR/${PREFIX}_teacher_summary.md" \
  --candidate-id "$CANDIDATE_ID" \
  --run-id "$RUN_ID" \
  --production-label "${STAGE_NS} ns total continuation stage"; then
  echo "[analysis:warn] plotting failed on remote; raw XVG outputs are preserved and figures can be regenerated locally"
fi

if [[ "$STAGE_NS" == "50" && "$CANDIDATE_ID" == "antibody_b005_18" ]]; then
  RECHECK_PREFIX="rechecked_40ns_from_50ns"
  RECHECK_WORK="$WORK_DIR/analysis_$RECHECK_PREFIX"
  RECHECK_ANALYSIS="$RESULTS_DIR/analysis/$RECHECK_PREFIX"
  RECHECK_FIGURES="$RESULTS_DIR/figures/$RECHECK_PREFIX"
  RECHECK_START_PS="${RECHECK_START_PS:-38623}"
  RECHECK_END_PS="${RECHECK_END_PS:-40000}"
  RECHECK_START_NS="${RECHECK_START_NS:-38.623}"
  RECHECK_END_NS="${RECHECK_END_NS:-40.000}"
  mkdir -p "$RECHECK_WORK" "$RECHECK_ANALYSIS" "$RECHECK_FIGURES"
  echo "[analysis] Rechecked 40 ns from new 50 ns trajectory"
  echo "[analysis] Recheck window: ${RECHECK_START_PS}-${RECHECK_END_PS} ps"
  echo "[analysis] Note: this is a checkpoint-window recheck, not a full 0-40 ns trajectory average."

  printf "5\n5\n" | "$GMX" rms -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$RECHECK_WORK/rmsd_complex.xvg" -tu ns -b "$RECHECK_START_NS" -e "$RECHECK_END_NS"
  printf "3\n3\n" | "$GMX" rms -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$RECHECK_WORK/rmsd_antibody.xvg" -tu ns -b "$RECHECK_START_NS" -e "$RECHECK_END_NS"
  printf "4\n4\n" | "$GMX" rms -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$RECHECK_WORK/rmsd_antigen.xvg" -tu ns -b "$RECHECK_START_NS" -e "$RECHECK_END_NS"

  echo "8" | "$GMX" rmsf -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -res -o "$RECHECK_WORK/rmsf_complex_ca.xvg" -b "$RECHECK_START_PS" -e "$RECHECK_END_PS"
  echo "6" | "$GMX" rmsf -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -res -o "$RECHECK_WORK/rmsf_antibody_ca.xvg" -b "$RECHECK_START_PS" -e "$RECHECK_END_PS"
  echo "7" | "$GMX" rmsf -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -res -o "$RECHECK_WORK/rmsf_antigen_ca.xvg" -b "$RECHECK_START_PS" -e "$RECHECK_END_PS"

  echo "2" | "$GMX" gyrate -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$RECHECK_WORK/rg_complex.xvg" -b "$RECHECK_START_PS" -e "$RECHECK_END_PS"
  echo "0" | "$GMX" gyrate -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$RECHECK_WORK/rg_antibody.xvg" -b "$RECHECK_START_PS" -e "$RECHECK_END_PS"
  echo "1" | "$GMX" gyrate -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$RECHECK_WORK/rg_antigen.xvg" -b "$RECHECK_START_PS" -e "$RECHECK_END_PS"

  "$GMX" sasa -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$RECHECK_WORK/sasa_complex.xvg" -surface 'group "Complex"' -output 'group "Complex"' -b "$RECHECK_START_PS" -e "$RECHECK_END_PS"
  "$GMX" sasa -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$RECHECK_WORK/sasa_antibody.xvg" -surface 'group "Antibody"' -output 'group "Antibody"' -b "$RECHECK_START_PS" -e "$RECHECK_END_PS"
  "$GMX" sasa -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -o "$RECHECK_WORK/sasa_antigen.xvg" -surface 'group "Antigen"' -output 'group "Antigen"' -b "$RECHECK_START_PS" -e "$RECHECK_END_PS"

  printf "0\n1\n" | "$GMX" hbond -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -num "$RECHECK_WORK/hbonds_antibody_antigen.xvg" -dist "$RECHECK_WORK/hbonds_distance.xvg" -ang "$RECHECK_WORK/hbonds_angle.xvg" -b "$RECHECK_START_PS" -e "$RECHECK_END_PS" || true
  printf "0\n1\n" | "$GMX" mindist -s "$TPR" -f "$XTC" -n index_ab_ag.ndx -od "$RECHECK_WORK/interface_mindist.xvg" -on "$RECHECK_WORK/interface_contacts.xvg" -d 0.45 -group -b "$RECHECK_START_PS" -e "$RECHECK_END_PS" || true

  cp -f "$RECHECK_WORK"/*.xvg "$RECHECK_ANALYSIS/" 2>/dev/null || true
  if ! python "$SCRIPT_DIR/72_plot_md_analysis_generic.py" \
    --analysis-dir "$RECHECK_ANALYSIS" \
    --output-dir "$RECHECK_FIGURES" \
    --summary "$RESULTS_DIR/${RECHECK_PREFIX}_teacher_summary.md" \
    --candidate-id "$CANDIDATE_ID" \
    --run-id "$RUN_ID" \
    --production-label "40 ns checkpoint-window recheck from 50 ns continuation (${RECHECK_START_PS}-${RECHECK_END_PS} ps)"; then
    echo "[analysis:warn] rechecked 40 ns plotting failed; raw XVG outputs are preserved"
  fi
fi

echo "[analysis] done_at=$(date)"
