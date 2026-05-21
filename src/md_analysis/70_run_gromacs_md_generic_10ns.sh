#!/usr/bin/env bash
set -Eeuo pipefail

RUN_ID="${RUN_ID:?RUN_ID is required, e.g. md_antibody_b005_18_af3_10ns}"
INPUT_PDB="${INPUT_PDB:?INPUT_PDB is required}"
CANDIDATE_ID="${CANDIDATE_ID:-$RUN_ID}"
BASE_DIR="${BASE_DIR:-/mnt/PPIFlow/ica2_runs}"
SCRIPT_DIR="${SCRIPT_DIR:-$BASE_DIR/scripts}"
WORK_DIR="${WORK_DIR:-$BASE_DIR/output/$RUN_ID}"
RESULTS_DIR="${RESULTS_DIR:-$BASE_DIR/results/$RUN_ID}"
LOG_DIR="${LOG_DIR:-$BASE_DIR/logs}"
GMX="${GMX:-gmx}"
THREADS="${THREADS:-16}"
NTOMP="${NTOMP:-$THREADS}"
NTMPI="${NTMPI:-1}"
PROD_NS="${PROD_NS:-10}"
PROD_STEPS="${PROD_STEPS:-5000000}"
REQUIRE_GPU="${REQUIRE_GPU:-1}"
NPT_REFCORD_SCALING_COM="${NPT_REFCORD_SCALING_COM:-0}"
GMX_GPU_ARGS="${GMX_GPU_ARGS:-auto}"
CHECKPOINT="$RESULTS_DIR/checkpoint.json"
LOG_FILE="$LOG_DIR/${RUN_ID}.log"

mkdir -p "$WORK_DIR" "$RESULTS_DIR" "$LOG_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

now_epoch() { date +%s; }

write_checkpoint() {
  local status="$1"
  local stage="$2"
  local detail="${3:-}"
  python - "$CHECKPOINT" "$status" "$stage" "$detail" "$RUN_ID" "$CANDIDATE_ID" "$WORK_DIR" "$RESULTS_DIR" <<'PY'
import json, sys, time
path, status, stage, detail, run_id, candidate_id, work_dir, results_dir = sys.argv[1:9]
payload = {
    "script": "70_run_gromacs_md_generic_10ns.sh",
    "status": status,
    "stage": stage,
    "detail": detail,
    "run_id": run_id,
    "candidate_id": candidate_id,
    "work_dir": work_dir,
    "results_dir": results_dir,
    "updated_at_epoch": int(time.time()),
}
with open(path, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2)
PY
}

run_stage() {
  local stage="$1"
  shift
  echo "[stage:$stage] started_at=$(date)"
  write_checkpoint "running" "$stage" "$*"
  "$@"
  echo "[stage:$stage] done_at=$(date)"
  write_checkpoint "completed" "$stage" "$*"
}

gmx_supports_cuda() {
  "$GMX" --version 2>/dev/null | grep -Eiq "GPU support[[:space:]]*:[[:space:]]*(CUDA|enabled|yes)"
}

resolve_gpu_args() {
  if [[ "$GMX_GPU_ARGS" != "auto" ]]; then
    # shellcheck disable=SC2206
    GPU_ARGS_ARRAY=($GMX_GPU_ARGS)
    return
  fi
  GPU_ARGS_ARRAY=()
  if gmx_supports_cuda; then
    GPU_ARGS_ARRAY=(-nb gpu -pme gpu -bonded gpu)
    if "$GMX" mdrun -h 2>&1 | grep -q -- "-update"; then
      GPU_ARGS_ARRAY+=(-update gpu)
    fi
  fi
}

make_mdp_files() {
  cat > "$WORK_DIR/minim.mdp" <<'EOF'
integrator              = steep
emtol                   = 1000.0
emstep                  = 0.01
nsteps                  = 50000
cutoff-scheme           = Verlet
nstlist                 = 10
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
pbc                     = xyz
EOF

  cat > "$WORK_DIR/nvt.mdp" <<EOF
define                  = -DPOSRES
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstenergy               = 500
nstlog                  = 500
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1 0.1
ref_t                   = 300 300
pcoupl                  = no
pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
EOF

  cat > "$WORK_DIR/npt.mdp" <<EOF
define                  = -DPOSRES
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstenergy               = 500
nstlog                  = 500
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1 0.1
ref_t                   = 300 300
pcoupl                  = Berendsen
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = no
EOF
  if [[ "$NPT_REFCORD_SCALING_COM" == "1" ]]; then
    echo "refcoord_scaling       = com" >> "$WORK_DIR/npt.mdp"
  fi

  cat > "$WORK_DIR/md.mdp" <<EOF
integrator              = md
nsteps                  = $PROD_STEPS
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 5000
compressed-x-grps       = System
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1 0.1
ref_t                   = 300 300
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = no
EOF
}

mdrun_cmd() {
  local deffnm="$1"
  shift || true
  local resume_args=()
  local parallel_args=(-ntomp "$NTOMP")
  if [[ -n "$NTMPI" ]]; then
    parallel_args=(-ntmpi "$NTMPI" -ntomp "$NTOMP")
  fi
  if [[ -s "$deffnm.cpt" && ! -s "$deffnm.gro" && " $* " != *" -cpi "* ]]; then
    resume_args=(-cpi "$deffnm.cpt" -append)
    echo "[mdrun] resuming $deffnm from $deffnm.cpt"
  fi
  "$GMX" mdrun -deffnm "$deffnm" -v "${parallel_args[@]}" "${GPU_ARGS_ARRAY[@]}" "${resume_args[@]}" "$@"
}

echo "[md] run_id=$RUN_ID candidate_id=$CANDIDATE_ID"
echo "[md] started_at=$(date)"
echo "[md] work_dir=$WORK_DIR"
echo "[md] results_dir=$RESULTS_DIR"
echo "[md] input_pdb=$INPUT_PDB"
echo "[md] protocol=amber99sb-ildn/tip3p/dodecahedron_1.0nm/0.15M_NaCl/NVT_100ps/NPT_100ps/production_${PROD_NS}ns_dt2fs"
write_checkpoint "running" "start" "initializing"

if [[ ! -s "$INPUT_PDB" ]]; then
  echo "[error] input PDB not found: $INPUT_PDB" >&2
  write_checkpoint "failed" "input" "missing input PDB"
  exit 2
fi

command -v "$GMX"
"$GMX" --version
resolve_gpu_args
echo "[md] gpu_args=${GPU_ARGS_ARRAY[*]:-none}"
if [[ "$REQUIRE_GPU" == "1" ]] && ! gmx_supports_cuda; then
  echo "[error] REQUIRE_GPU=1 but $GMX does not report CUDA/enabled GPU support; not starting long CPU MD." >&2
  write_checkpoint "failed" "gpu_check" "CUDA-enabled GROMACS not available"
  exit 10
fi

make_mdp_files
cp "$INPUT_PDB" "$WORK_DIR/input_raw.pdb"
cd "$WORK_DIR"

python - <<'PY'
import json
from pathlib import Path

raw = Path("input_raw.pdb")
clean = Path("input_clean.pdb")
removed = []
kept = []
chains = {}
for line in raw.read_text(errors="ignore").splitlines():
    if line.startswith(("ATOM", "HETATM")):
        atom = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip()
        chains.setdefault(chain, set()).add((line[22:26].strip(), line[26].strip()))
        if resname == "GLY" and atom == "CB":
            removed.append({
                "atom_serial": line[6:11].strip(),
                "chain": chain,
                "residue": line[22:26].strip(),
                "icode": line[26].strip(),
            })
            continue
    kept.append(line)
clean.write_text("\n".join(kept) + "\n", encoding="utf-8")
summary = {
    "raw_pdb": str(raw),
    "cleaned_pdb": str(clean),
    "operation": "removed GLY CB atoms only; no atoms were invented",
    "removed_atom_count": len(removed),
    "removed_atoms": removed,
    "chain_residue_counts_raw": {k: len(v) for k, v in sorted(chains.items())},
}
Path("pdb_cleaning_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print("[md] chain_residue_counts_raw", summary["chain_residue_counts_raw"], flush=True)
print(f"[md] pdb cleaning removed_gly_cb={len(removed)}", flush=True)
missing = {"A", "B", "C"} - set(chains)
if missing:
    raise SystemExit(f"Missing expected chains A/B/C in input: {sorted(missing)}")
PY
cp -f pdb_cleaning_summary.json input_raw.pdb input_clean.pdb "$RESULTS_DIR/" 2>/dev/null || true

PDB2GMX_INPUT="input_clean.pdb"
if python - <<'PY'
import openmm  # noqa: F401
import pdbfixer  # noqa: F401
PY
then
  echo "[md] PDBFixer available; rebuilding missing heavy atoms"
  python "$SCRIPT_DIR/23_rebuild_missing_atoms_pdbfixer.py" \
    --input input_clean.pdb \
    --output input_rebuilt_heavy.pdb \
    --summary "$RESULTS_DIR/pdbfixer_rebuild_summary.json"
  cp -f input_rebuilt_heavy.pdb "$RESULTS_DIR/" 2>/dev/null || true
  PDB2GMX_INPUT="input_rebuilt_heavy.pdb"
else
  echo "[md] PDBFixer not available; pdb2gmx will use GLY-CB-cleaned input only"
fi
echo "[md] pdb2gmx_input=$PDB2GMX_INPUT"

if [[ ! -s processed.pdb || ! -s topol.top ]]; then
  run_stage "pdb2gmx" "$GMX" pdb2gmx -f "$PDB2GMX_INPUT" -o processed.pdb -p topol.top -ff amber99sb-ildn -water tip3p -ignh
fi

python "$SCRIPT_DIR/20_make_md_index.py" --pdb processed.pdb --output index_ab_ag.ndx --summary "$RESULTS_DIR/index_summary.json"

if [[ ! -s boxed.gro ]]; then
  run_stage "editconf" "$GMX" editconf -f processed.pdb -o boxed.gro -bt dodecahedron -d 1.0
fi
if [[ ! -s solv.gro ]]; then
  run_stage "solvate" "$GMX" solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top
fi
if [[ ! -s ions.tpr ]]; then
  run_stage "grompp_ions" "$GMX" grompp -f minim.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
fi
if [[ ! -s solv_ions.gro ]]; then
  echo "[stage:genion] started_at=$(date)"
  write_checkpoint "running" "genion" "neutralize and add 0.15 M NaCl"
  echo "SOL" | "$GMX" genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
  echo "[stage:genion] done_at=$(date)"
  write_checkpoint "completed" "genion" "solv_ions.gro"
fi
if [[ ! -s em.tpr ]]; then
  run_stage "grompp_em" "$GMX" grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
fi
if [[ ! -s em.gro ]]; then
  run_stage "mdrun_em" mdrun_cmd em
fi
if [[ ! -s nvt.tpr ]]; then
  run_stage "grompp_nvt" "$GMX" grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
fi
if [[ ! -s nvt.gro ]]; then
  run_stage "mdrun_nvt" mdrun_cmd nvt
fi
if [[ ! -s npt.tpr ]]; then
  run_stage "grompp_npt" "$GMX" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
fi
if [[ ! -s npt.gro ]]; then
  run_stage "mdrun_npt" mdrun_cmd npt
fi
if [[ ! -s md_10ns.tpr ]]; then
  run_stage "grompp_md" "$GMX" grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_10ns.tpr -maxwarn 1
fi

echo "[stage:mdrun_10ns] started_at=$(date)"
write_checkpoint "running" "mdrun_10ns" "10 ns production MD"
if [[ -s md_10ns.cpt && ! -s md_10ns.gro ]]; then
  mdrun_cmd md_10ns -cpi md_10ns.cpt -append
elif [[ ! -s md_10ns.gro ]]; then
  mdrun_cmd md_10ns
else
  echo "[stage:mdrun_10ns] md_10ns.gro exists; skipping production rerun"
fi
echo "[stage:mdrun_10ns] done_at=$(date)"
write_checkpoint "completed" "mdrun_10ns" "production complete or already present"

echo "[stage:trjconv] started_at=$(date)"
write_checkpoint "running" "trjconv" "remove PBC and center"
if [[ ! -s md_nojump.xtc ]]; then
  echo "0" | "$GMX" trjconv -s md_10ns.tpr -f md_10ns.xtc -o md_nojump.xtc -pbc nojump
fi
if [[ ! -s md_center.xtc ]]; then
  printf "1\n0\n" | "$GMX" trjconv -s md_10ns.tpr -f md_nojump.xtc -o md_center.xtc -center -pbc mol -ur compact
fi
if [[ ! -s md_final.pdb ]]; then
  echo "2" | "$GMX" trjconv -s md_10ns.tpr -f md_10ns.gro -o md_final.pdb -n index_ab_ag.ndx
fi
echo "[stage:trjconv] done_at=$(date)"
write_checkpoint "completed" "trjconv" "centered trajectory and final PDB"

cp -f md_10ns.tpr md_10ns.xtc md_10ns.gro md_10ns.edr md_10ns.log md_center.xtc md_final.pdb index_ab_ag.ndx topol.top "$RESULTS_DIR/" 2>/dev/null || true

echo "[stage:analysis] started_at=$(date)"
write_checkpoint "running" "analysis" "GROMACS analysis and plotting"
bash "$SCRIPT_DIR/71_run_gromacs_md_analysis_generic.sh"
echo "[stage:analysis] done_at=$(date)"
write_checkpoint "completed" "analysis" "analysis complete"

echo "[md] completed_at=$(date)"
write_checkpoint "completed" "all_md" "MD run and analysis complete"
