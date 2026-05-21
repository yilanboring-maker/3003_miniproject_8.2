#!/usr/bin/env bash
set -Eeuo pipefail

RUN_ID="${RUN_ID:?RUN_ID is required}"
CANDIDATE_ID="${CANDIDATE_ID:-$RUN_ID}"
BASE_DIR="${BASE_DIR:-/mnt/PPIFlow/ica2_runs}"
SCRIPT_DIR="${SCRIPT_DIR:-$BASE_DIR/scripts}"
WORK_DIR="${WORK_DIR:-$BASE_DIR/output/$RUN_ID}"
RESULTS_DIR="${RESULTS_DIR:-$BASE_DIR/results/$RUN_ID}"
LOG_DIR="${LOG_DIR:-$BASE_DIR/logs}"
GMX="${GMX:-gmx}"
LOG_FILE="$LOG_DIR/${RUN_ID}_analysis.log"

mkdir -p "$RESULTS_DIR/analysis" "$RESULTS_DIR/figures" "$LOG_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

cd "$WORK_DIR"
echo "[analysis] started_at=$(date)"
echo "[analysis] run_id=$RUN_ID candidate_id=$CANDIDATE_ID"
echo "[analysis] work_dir=$WORK_DIR"
echo "[analysis] results_dir=$RESULTS_DIR"

for f in md_10ns.tpr md_center.xtc index_ab_ag.ndx; do
  if [[ ! -s "$f" ]]; then
    echo "[analysis:error] missing $f" >&2
    exit 2
  fi
done

mkdir -p analysis

echo "[analysis] RMSD"
printf "5\n5\n" | "$GMX" rms -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -o analysis/rmsd_complex.xvg -tu ns
printf "3\n3\n" | "$GMX" rms -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -o analysis/rmsd_antibody.xvg -tu ns
printf "4\n4\n" | "$GMX" rms -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -o analysis/rmsd_antigen.xvg -tu ns

echo "[analysis] RMSF"
echo "8" | "$GMX" rmsf -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -res -o analysis/rmsf_complex_ca.xvg
echo "6" | "$GMX" rmsf -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -res -o analysis/rmsf_antibody_ca.xvg
echo "7" | "$GMX" rmsf -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -res -o analysis/rmsf_antigen_ca.xvg

echo "[analysis] Rg"
echo "2" | "$GMX" gyrate -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -o analysis/rg_complex.xvg
echo "0" | "$GMX" gyrate -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -o analysis/rg_antibody.xvg
echo "1" | "$GMX" gyrate -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -o analysis/rg_antigen.xvg

echo "[analysis] SASA"
"$GMX" sasa -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -o analysis/sasa_complex.xvg -surface 'group "Complex"' -output 'group "Complex"'
"$GMX" sasa -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -o analysis/sasa_antibody.xvg -surface 'group "Antibody"' -output 'group "Antibody"'
"$GMX" sasa -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -o analysis/sasa_antigen.xvg -surface 'group "Antigen"' -output 'group "Antigen"'

echo "[analysis] H-bonds"
printf "0\n1\n" | "$GMX" hbond -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -num analysis/hbonds_antibody_antigen.xvg -dist analysis/hbonds_distance.xvg -ang analysis/hbonds_angle.xvg || true

echo "[analysis] Interface contacts"
printf "0\n1\n" | "$GMX" mindist -s md_10ns.tpr -f md_center.xtc -n index_ab_ag.ndx -od analysis/interface_mindist.xvg -on analysis/interface_contacts.xvg -d 0.45 -group || true

cp -f analysis/*.xvg "$RESULTS_DIR/analysis/" 2>/dev/null || true
if ! python "$SCRIPT_DIR/72_plot_md_analysis_generic.py" \
  --analysis-dir "$RESULTS_DIR/analysis" \
  --output-dir "$RESULTS_DIR/figures" \
  --summary "$RESULTS_DIR/md_teacher_summary.md" \
  --candidate-id "$CANDIDATE_ID" \
  --run-id "$RUN_ID"; then
  echo "[analysis:warn] plotting failed on remote; raw XVG outputs are preserved and figures can be regenerated locally"
fi

echo "[analysis] done_at=$(date)"
