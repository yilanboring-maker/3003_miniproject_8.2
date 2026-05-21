#!/usr/bin/env python3
"""Run corrected AbMPNN + FlowPacker optimization for every IgG candidate.

This script is intended to run on the Matpool/JupyterLab machine under
/mnt/PPIFlow. It keeps checkpoints and partial CSVs so the run can be audited
or resumed.

Corrections vs the first all-IgG run:
- constrain AbMPNN output back to the original IgG framework, mutating only
  PPIFlow variable Ala-placeholder runs in antibody chains A/B;
- keep conserved framework residues, including Cys, from the original IgG;
- copy the actual FlowPacker model output from run_1/best_run, not the
  intermediate residue-renamed input PDB;
- validate antibody side-chain completeness before accepting a packed PDB.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path


PPI_ROOT = Path("/mnt/PPIFlow/PPIFlow-main")
RUN_ROOT = Path("/mnt/PPIFlow/ica2_runs")
REMOTE_PRODIGY = RUN_ROOT / "results/prodigy_all_results.csv"
WORK_DIR = Path(os.environ.get("ABMPNN_FLOWPACKER_WORK_DIR", "/tmp/ica2_abmpnn_flowpacker/output/abmpnn_flowpacker_all_igg_corrected"))
RESULTS_DIR = Path(os.environ.get("ABMPNN_FLOWPACKER_RESULTS_DIR", "/tmp/ica2_abmpnn_flowpacker/results/abmpnn_flowpacker_all_igg_corrected"))
LOG_DIR = RUN_ROOT / "logs"
SCRIPTS_DIR = RUN_ROOT / "scripts"
FLOWPACKER_ZIP = PPI_ROOT / "flowpacker-main.zip"
FLOWPACKER_TOOLS_PARENT = Path(os.environ.get("FLOWPACKER_TOOLS_PARENT", "/tmp/ica2_abmpnn_flowpacker/tools"))
FLOWPACKER_DIR = FLOWPACKER_TOOLS_PARENT / "flowpacker-main"
FLOWPACKER_CHECKPOINT_DIR = Path(os.environ.get("FLOWPACKER_CHECKPOINT_DIR", "/tmp/flowpacker_checkpoints"))
FLOWPACKER_CLUSTER_CKPT = FLOWPACKER_CHECKPOINT_DIR / "cluster.pth"
FLOWPACKER_CONF_CKPT = FLOWPACKER_CHECKPOINT_DIR / "confidence.pth"

HOTSPOTS = "C87,C88,C89,C90,C91,C114,C115,C116,C117"
DESIGNED_CHAINS = ("A", "B")
FIXED_ANTIGEN_CHAIN = "C"
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
MUTABLE_ALA_RUN_MIN = int(os.environ.get("ABMPNN_MUTABLE_ALA_RUN_MIN", "3"))

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AA1_TO_3 = {v: k for k, v in AA3_TO_1.items()}


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[{now()}] {msg}", flush=True)


def ensure_dirs() -> None:
    for path in [
        WORK_DIR,
        RESULTS_DIR,
        LOG_DIR,
        SCRIPTS_DIR,
        FLOWPACKER_TOOLS_PARENT,
        WORK_DIR / "input_original_pdbs",
        WORK_DIR / "abmpnn",
        WORK_DIR / "mutated_pdbs",
        WORK_DIR / "flowpacker_batches",
        WORK_DIR / "flowpacker_raw",
        RESULTS_DIR / "packed_pdbs",
    ]:
        path.mkdir(parents=True, exist_ok=True)


def run_cmd(cmd: list[str], cwd: Path | None, log_path: Path, env: dict[str, str] | None = None) -> int:
    log(f"RUN cwd={cwd or Path.cwd()} cmd={' '.join(cmd)}")
    with log_path.open("a", encoding="utf-8") as lf:
        lf.write(f"\n[{now()}] RUN cwd={cwd or Path.cwd()} cmd={' '.join(cmd)}\n")
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd else None,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="", flush=True)
            lf.write(line)
        rc = proc.wait()
        lf.write(f"[{now()}] EXIT {rc}\n")
    log(f"EXIT {rc}")
    return rc


def read_manifest_rows() -> list[dict[str, str]]:
    if not REMOTE_PRODIGY.exists():
        raise FileNotFoundError(f"Missing remote PRODIGY table: {REMOTE_PRODIGY}")
    rows = list(csv.DictReader(REMOTE_PRODIGY.open(newline="", encoding="utf-8")))
    antibodies = [r for r in rows if r.get("modality") == "antibody"]

    def rank_key(row: dict[str, str]) -> int:
        try:
            return int(float(row.get("prodigy_rank") or row.get("rank") or "999999"))
        except ValueError:
            return 999999

    antibodies.sort(key=rank_key)
    return antibodies


def residue_records(pdb_path: Path) -> dict[str, list[tuple[tuple[str, str], str]]]:
    chains: dict[str, list[tuple[tuple[str, str], str]]] = {}
    seen: set[tuple[str, str, str]] = set()
    with pdb_path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            chain = line[21]
            resseq = line[22:26].strip()
            icode = line[26].strip()
            key = (chain, resseq, icode)
            if key in seen:
                continue
            seen.add(key)
            chains.setdefault(chain, []).append(((resseq, icode), line[17:20].strip()))
    return chains


def chain_sequence(pdb_path: Path, chain: str) -> str:
    chains = residue_records(pdb_path)
    seq = []
    for _, resname in chains.get(chain, []):
        seq.append(AA3_TO_1.get(resname.upper(), "X"))
    return "".join(seq)


def create_input_manifest(limit: int | None, smoke: bool) -> list[dict[str, str]]:
    antibodies = read_manifest_rows()
    if smoke:
        picks = []
        if antibodies:
            picks.append(antibodies[0])
            picks.append(antibodies[len(antibodies) // 2])
            picks.append(antibodies[-1])
        antibodies = picks
    elif limit:
        antibodies = antibodies[:limit]

    manifest = []
    missing = []
    input_dir = WORK_DIR / ("smoke_input_original_pdbs" if smoke else "input_original_pdbs")
    input_dir.mkdir(parents=True, exist_ok=True)
    for row in antibodies:
        design_id = row["design_id"]
        src = Path(row["pdb_path"])
        if not src.exists():
            missing.append({"design_id": design_id, "pdb_path": str(src)})
            continue
        dst = input_dir / f"{design_id}.pdb"
        shutil.copy2(src, dst)
        chain_a = chain_sequence(dst, "A")
        chain_b = chain_sequence(dst, "B")
        chain_c = chain_sequence(dst, "C")
        manifest.append({
            "design_id": design_id,
            "prodigy_rank": row.get("prodigy_rank", ""),
            "original_prodigy_delta_g": row.get("prodigy_delta_g", ""),
            "original_prodigy_kd": row.get("prodigy_kd", row.get("prodigy_kd_m", "")),
            "source_pdb_path": str(src),
            "input_pdb_path": str(dst),
            "original_seq_A": chain_a,
            "original_seq_B": chain_b,
            "original_seq_C": chain_c,
            "len_A": str(len(chain_a)),
            "len_B": str(len(chain_b)),
            "len_C": str(len(chain_c)),
            "hotspots": HOTSPOTS,
        })

    manifest_path = RESULTS_DIR / ("smoke_igg_input_manifest.csv" if smoke else "igg_input_manifest.csv")
    write_csv(manifest_path, manifest)
    if missing:
        write_csv(RESULTS_DIR / ("smoke_missing_inputs.csv" if smoke else "missing_inputs.csv"), missing)
    log(f"Prepared manifest rows={len(manifest)} missing={len(missing)} smoke={smoke}")
    return manifest


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        keys: list[str] = []
        for row in rows:
            for key in row:
                if key not in keys:
                    keys.append(key)
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=keys)
            writer.writeheader()
            writer.writerows(rows)
    else:
        path.write_text("", encoding="utf-8")


def mutable_placeholder_mask(seq: str, min_run: int = MUTABLE_ALA_RUN_MIN) -> list[bool]:
    """Treat only long Ala runs as mutable PPIFlow design placeholders."""
    mask = [False] * len(seq)
    i = 0
    while i < len(seq):
        if seq[i] != "A":
            i += 1
            continue
        j = i
        while j < len(seq) and seq[j] == "A":
            j += 1
        if j - i >= min_run:
            for k in range(i, j):
                mask[k] = True
        i = j
    return mask


def constrained_sequence(original: str, sampled: str, mask: list[bool]) -> str:
    out: list[str] = []
    for idx, aa in enumerate(original):
        sampled_aa = sampled[idx] if idx < len(sampled) else ""
        if idx < len(mask) and mask[idx] and sampled_aa in AA1_TO_3:
            out.append(sampled_aa)
        else:
            out.append(aa)
    return "".join(out)


def mask_to_ranges(mask: list[bool]) -> str:
    ranges: list[str] = []
    start: int | None = None
    for idx, value in enumerate(mask + [False], 1):
        if value and start is None:
            start = idx
        elif not value and start is not None:
            end = idx - 1
            ranges.append(str(start) if start == end else f"{start}-{end}")
            start = None
    return ",".join(ranges)


def update_checkpoint(stage: str, rows_done: int, rows_total: int, extra: dict[str, object] | None = None) -> None:
    payload = {
        "timestamp": now(),
        "stage": stage,
        "done": rows_done,
        "total": rows_total,
        "percent": round(100.0 * rows_done / rows_total, 2) if rows_total else 0,
    }
    if extra:
        payload.update(extra)
    (RESULTS_DIR / "checkpoint.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def run_abmpnn(manifest: list[dict[str, str]], smoke: bool, force: bool) -> Path:
    input_dir = Path(manifest[0]["input_pdb_path"]).parent
    out_dir = WORK_DIR / ("abmpnn_smoke" if smoke else "abmpnn")
    seq_dir = out_dir / "seqs"
    expected = [seq_dir / f"{row['design_id']}.fa" for row in manifest]
    if not force and expected and all(p.exists() for p in expected):
        log(f"AbMPNN outputs already exist for {len(expected)} targets; skipping generation.")
        return out_dir

    if force and out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"45_abmpnn_all_igg_corrected_{'smoke' if smoke else 'full'}_{time.strftime('%Y%m%d_%H%M%S')}.log"
    # The workshop copy of protein_mpnn_run.py references fixed_positions_dict
    # even when no --position_list is provided. Use an auditable patched copy
    # that only initializes the missing local variable; design settings stay the
    # same and the original PPIFlow file is not modified.
    safe_runner = SCRIPTS_DIR / "protein_mpnn_run_abmpnn_safe.py"
    original_runner = PPI_ROOT / "ProteinMPNN/protein_mpnn_run.py"
    runner_text = original_runner.read_text(encoding="utf-8")
    if "fixed_positions_dict = None" not in runner_text[:2000]:
        runner_text = runner_text.replace(
            "def main(args):\n",
            "def main(args):\n    fixed_positions_dict = None\n",
            1,
        )
    safe_runner.write_text(runner_text, encoding="utf-8")
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{PPI_ROOT / 'ProteinMPNN'}:{env.get('PYTHONPATH', '')}"
    cmd = [
        sys.executable,
        str(safe_runner),
        "--path_to_model_weights", "ProteinMPNN/model_weights/",
        "--model_name", "abmpnn",
        "--folder_with_pdbs_path", str(input_dir),
        "--out_folder", str(out_dir),
        "--chain_list", "A B",
        "--num_seq_per_target", "1",
        "--sampling_temp", "0.5",
        "--seed", "37",
        "--batch_size", "1",
    ]
    rc = run_cmd(cmd, PPI_ROOT, log_path, env=env)
    if rc != 0:
        raise RuntimeError(f"AbMPNN failed with exit code {rc}; see {log_path}")
    return out_dir


def parse_fasta_design(fasta_path: Path) -> str:
    lines = []
    for line in fasta_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if line and not line.startswith(">"):
            lines.append(line)
    if not lines:
        raise ValueError(f"No sequence lines in {fasta_path}")
    return lines[-1].replace(" ", "")


def split_design_sequence(raw_seq: str, len_a: int, len_b: int, len_c: int) -> tuple[str, str]:
    seq = raw_seq.strip()
    parts = [p for p in seq.split("/") if p]
    if len(parts) >= 3:
        for i in range(len(parts) - 2):
            if len(parts[i]) == len_a and len(parts[i + 1]) == len_b and len(parts[i + 2]) == len_c:
                return parts[i], parts[i + 1]
    if len(parts) >= 2:
        for i in range(len(parts) - 1):
            if len(parts[i]) == len_a and len(parts[i + 1]) == len_b:
                return parts[i], parts[i + 1]
    compact = "".join(parts) if parts else seq
    if len(compact) >= len_a + len_b:
        return compact[:len_a], compact[len_a:len_a + len_b]
    raise ValueError(f"Could not split AbMPNN sequence lengths A={len_a}, B={len_b}, C={len_c}: {raw_seq[:120]}")


def mutate_pdb_constrained_input(
    src_pdb: Path,
    dst_pdb: Path,
    seq_a: str,
    seq_b: str,
    mutable_a: list[bool],
    mutable_b: list[bool],
) -> None:
    chains = residue_records(src_pdb)
    index_maps: dict[str, dict[tuple[str, str], str]] = {}
    mutable_maps: dict[str, set[tuple[str, str]]] = {}
    for chain, seq in [("A", seq_a), ("B", seq_b)]:
        records = chains.get(chain, [])
        if len(records) != len(seq):
            raise ValueError(f"{src_pdb.name} chain {chain} residue/sequence length mismatch {len(records)} vs {len(seq)}")
        index_maps[chain] = {records[i][0]: AA1_TO_3.get(seq[i], "ALA") for i in range(len(records))}
        mask = mutable_a if chain == "A" else mutable_b
        mutable_maps[chain] = {records[i][0] for i, is_mutable in enumerate(mask) if is_mutable}

    with src_pdb.open(encoding="utf-8", errors="ignore") as inp, dst_pdb.open("w", encoding="utf-8") as out:
        for line in inp:
            if line.startswith("ATOM") and line[21] in DESIGNED_CHAINS:
                atom = line[12:16].strip()
                chain = line[21]
                key = (line[22:26].strip(), line[26].strip())
                if key in mutable_maps[chain] and atom not in BACKBONE_ATOMS:
                    continue
                new_resname = index_maps[chain].get(key)
                if new_resname:
                    line = line[:17] + new_resname.rjust(3) + line[20:]
            out.write(line)


def parse_abmpnn_outputs(manifest: list[dict[str, str]], abmpnn_dir: Path, smoke: bool) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    fasta_lines: list[str] = []
    mutated_dir = WORK_DIR / ("mutated_pdbs_smoke" if smoke else "mutated_pdbs")
    mutated_dir.mkdir(parents=True, exist_ok=True)

    for idx, row in enumerate(manifest, 1):
        design_id = row["design_id"]
        fasta_path = abmpnn_dir / "seqs" / f"{design_id}.fa"
        status = "ok"
        error = ""
        opt_a = opt_b = ""
        raw_opt_a = raw_opt_b = ""
        mutable_a = mutable_placeholder_mask(row["original_seq_A"])
        mutable_b = mutable_placeholder_mask(row["original_seq_B"])
        mutation_count = 0
        mutated_path = mutated_dir / f"{design_id}.pdb"
        try:
            raw = parse_fasta_design(fasta_path)
            raw_opt_a, raw_opt_b = split_design_sequence(raw, int(row["len_A"]), int(row["len_B"]), int(row["len_C"]))
            opt_a = constrained_sequence(row["original_seq_A"], raw_opt_a, mutable_a)
            opt_b = constrained_sequence(row["original_seq_B"], raw_opt_b, mutable_b)
            mutation_count = sum(a != b for a, b in zip(row["original_seq_A"], opt_a)) + sum(a != b for a, b in zip(row["original_seq_B"], opt_b))
            mutate_pdb_constrained_input(Path(row["input_pdb_path"]), mutated_path, opt_a, opt_b, mutable_a, mutable_b)
            fasta_lines.append(f">{design_id}|chains=A+B|source={row['source_pdb_path']}")
            fasta_lines.append(f"{opt_a}/{opt_b}")
        except Exception as exc:  # keep partial results
            status = "failed_abmpnn_parse_or_mutation"
            error = str(exc)

        out = dict(row)
        out.update({
            "raw_abmpnn_seq_A": raw_opt_a,
            "raw_abmpnn_seq_B": raw_opt_b,
            "optimized_seq_A": opt_a,
            "optimized_seq_B": opt_b,
            "mutable_rule": f"original Ala runs length >= {MUTABLE_ALA_RUN_MIN}",
            "mutable_positions_A": mask_to_ranges(mutable_a),
            "mutable_positions_B": mask_to_ranges(mutable_b),
            "mutable_count_A": sum(mutable_a),
            "mutable_count_B": sum(mutable_b),
            "mutation_count_A": sum(a != b for a, b in zip(row["original_seq_A"], opt_a)) if opt_a else "",
            "mutation_count_B": sum(a != b for a, b in zip(row["original_seq_B"], opt_b)) if opt_b else "",
            "mutation_count_total": mutation_count if opt_a and opt_b else "",
            "abmpnn_fasta": str(fasta_path),
            "mutated_backbone_pdb": str(mutated_path) if mutated_path.exists() else "",
            "status": status,
            "error": error,
        })
        rows.append(out)
        if idx % 20 == 0 or idx == len(manifest):
            write_csv(RESULTS_DIR / "partial_results.csv", rows)
            update_checkpoint("abmpnn_parsed", idx, len(manifest), {"ok": sum(r["status"] == "ok" for r in rows)})
            log(f"AbMPNN parse progress {idx}/{len(manifest)} ok={sum(r['status'] == 'ok' for r in rows)}")

    fasta_path = RESULTS_DIR / "abmpnn_sequences.fasta"
    fasta_path.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    write_csv(RESULTS_DIR / "abmpnn_variants.csv", rows)
    return rows


def ensure_flowpacker() -> None:
    if not FLOWPACKER_DIR.exists():
        if not FLOWPACKER_ZIP.exists():
            raise FileNotFoundError(f"Missing FlowPacker zip: {FLOWPACKER_ZIP}")
        log(f"Extracting FlowPacker to {FLOWPACKER_TOOLS_PARENT}")
        FLOWPACKER_TOOLS_PARENT.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(FLOWPACKER_ZIP) as zf:
            zf.extractall(FLOWPACKER_TOOLS_PARENT)
    base_yaml = FLOWPACKER_DIR / "config/inference/base.yaml"
    if not base_yaml.exists():
        raise FileNotFoundError(f"Missing FlowPacker base config: {base_yaml}")
    if not FLOWPACKER_CLUSTER_CKPT.exists():
        raise FileNotFoundError(f"Missing FlowPacker cluster checkpoint: {FLOWPACKER_CLUSTER_CKPT}")
    if not FLOWPACKER_CONF_CKPT.exists():
        raise FileNotFoundError(f"Missing FlowPacker confidence checkpoint: {FLOWPACKER_CONF_CKPT}")


def write_flowpacker_config(base_yaml: Path, dest_yaml: Path, test_path: Path) -> None:
    lines = base_yaml.read_text(encoding="utf-8").splitlines()
    out = []
    replaced = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("test_path:"):
            indent = line[: len(line) - len(line.lstrip())]
            out.append(f"{indent}test_path: '{test_path}'")
            replaced = True
        elif stripped.startswith("ckpt:"):
            indent = line[: len(line) - len(line.lstrip())]
            out.append(f"{indent}ckpt: '{FLOWPACKER_CLUSTER_CKPT}'")
        elif stripped.startswith("conf_ckpt:"):
            indent = line[: len(line) - len(line.lstrip())]
            out.append(f"{indent}conf_ckpt: '{FLOWPACKER_CONF_CKPT}'")
        else:
            out.append(line)
    if not replaced:
        raise ValueError(f"Could not patch test_path in {base_yaml}")
    dest_yaml.write_text("\n".join(out) + "\n", encoding="utf-8")


def run_flowpacker(rows: list[dict[str, object]], smoke: bool, batch_size: int, force: bool) -> list[dict[str, object]]:
    ensure_flowpacker()
    valid = [r for r in rows if r.get("status") == "ok" and r.get("mutated_backbone_pdb")]
    base_yaml = FLOWPACKER_DIR / "config/inference/base.yaml"
    log_path = LOG_DIR / f"45_flowpacker_all_igg_corrected_{'smoke' if smoke else 'full'}_{time.strftime('%Y%m%d_%H%M%S')}.log"
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{FLOWPACKER_DIR}:{PPI_ROOT}:{env.get('PYTHONPATH', '')}"
    env["CUDA_VISIBLE_DEVICES"] = env.get("CUDA_VISIBLE_DEVICES", "0")

    done = 0
    for start in range(0, len(valid), batch_size):
        batch = valid[start:start + batch_size]
        batch_name = f"{'smoke_' if smoke else ''}batch_{start // batch_size + 1:03d}"
        batch_input = WORK_DIR / "flowpacker_batches" / batch_name
        batch_input.mkdir(parents=True, exist_ok=True)
        batch_csv = WORK_DIR / "flowpacker_batches" / f"{batch_name}.csv"
        batch_cfg = WORK_DIR / "flowpacker_batches" / f"{batch_name}.yaml"
        batch_save = WORK_DIR / "flowpacker_raw" / batch_name / "run"
        if force and batch_save.parent.exists():
            shutil.rmtree(batch_save.parent)
        batch_save.mkdir(parents=True, exist_ok=True)

        csv_rows = []
        for row in batch:
            design_id = str(row["design_id"])
            dst = batch_input / f"{design_id}.pdb"
            shutil.copy2(str(row["mutated_backbone_pdb"]), dst)
            csv_rows.append({"link_name": f"{design_id}.pdb", "seq": row["optimized_seq_A"], "seq_idx": "opt_001"})
        write_csv(batch_csv, csv_rows)
        write_flowpacker_config(base_yaml, batch_cfg, batch_input)

        expected_names = [f"{row['design_id']}_opt_001.pdb" for row in batch]
        final_dir = RESULTS_DIR / "packed_pdbs"
        if not force and all((final_dir / name).exists() for name in expected_names):
            log(f"FlowPacker batch {batch_name} already complete; skipping.")
        else:
            cmd = [
                sys.executable,
                str(PPI_ROOT / "demo_scripts/flowpacker_af3score/sampler_pdb_pipe.py"),
                str(batch_cfg),
                "--save_dir", str(batch_save),
                "--use_gt_masks", "True",
                "--csv_file", str(batch_csv),
            ]
            rc = run_cmd(cmd, FLOWPACKER_TOOLS_PARENT, log_path, env=env)
            if rc != 0:
                for row in batch:
                    row["status"] = "failed_flowpacker"
                    row["error"] = f"FlowPacker batch {batch_name} exited {rc}; see {log_path}"
                write_csv(RESULTS_DIR / "partial_results.csv", rows)
                continue

            for row in batch:
                name = f"{row['design_id']}_opt_001.pdb"
                candidates = [
                    batch_save / "best_run" / name,
                    batch_save / "run_1" / name,
                ]
                src = next((p for p in candidates if p.exists()), candidates[0])
                dst = final_dir / name
                if src.exists():
                    shutil.copy2(src, dst)
                    row["packed_pdb"] = str(dst)
                    row["packed_source"] = str(src)
                    row["status"] = "ok"
                    row["error"] = ""
                else:
                    row["status"] = "failed_flowpacker_missing_output"
                    row["error"] = "Missing packed output; tried " + "; ".join(str(p) for p in candidates)

        done += len(batch)
        ok_count = sum(1 for r in rows if r.get("packed_pdb"))
        write_csv(RESULTS_DIR / "partial_results.csv", rows)
        update_checkpoint("flowpacker", done, len(valid), {"packed_ok": ok_count})
        elapsed = ""
        log(f"FlowPacker progress {done}/{len(valid)} packed_ok={ok_count} batch={batch_name} {elapsed}")

    return rows


def sidechain_stats(pdb: Path, chain: str) -> dict[str, int]:
    residues: dict[tuple[str, str], dict[str, object]] = {}
    with pdb.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM") or line[21] != chain:
                continue
            key = (line[22:26].strip(), line[26].strip())
            atom = line[12:16].strip()
            resname = line[17:20].strip()
            rec = residues.setdefault(key, {"resname": resname, "sidechain": False})
            if atom not in BACKBONE_ATOMS:
                rec["sidechain"] = True
    total = len(residues)
    expected = sum(1 for rec in residues.values() if rec["resname"] != "GLY")
    present = sum(1 for rec in residues.values() if rec["resname"] == "GLY" or rec["sidechain"])
    missing = expected - sum(1 for rec in residues.values() if rec["resname"] != "GLY" and rec["sidechain"])
    return {"total": total, "expected_non_gly": expected, "present_or_gly": present, "missing_non_gly_sidechain": missing}


def validate_outputs(rows: list[dict[str, object]]) -> None:
    for row in rows:
        packed = row.get("packed_pdb")
        if not packed:
            continue
        pdb = Path(str(packed))
        chains = residue_records(pdb)
        row["packed_has_chain_A"] = "A" in chains
        row["packed_has_chain_B"] = "B" in chains
        row["packed_has_chain_C"] = "C" in chains
        row["packed_chain_A_residues"] = len(chains.get("A", []))
        row["packed_chain_B_residues"] = len(chains.get("B", []))
        row["packed_chain_C_residues"] = len(chains.get("C", []))
        row["antigen_C_unchanged"] = chain_sequence(pdb, "C") == row.get("original_seq_C")
        for chain in ("A", "B", "C"):
            stats = sidechain_stats(pdb, chain)
            row[f"chain_{chain}_expected_non_gly_sidechains"] = stats["expected_non_gly"]
            row[f"chain_{chain}_missing_non_gly_sidechains"] = stats["missing_non_gly_sidechain"]
        def missing_count(chain_id: str) -> int:
            value = row.get(f"chain_{chain_id}_missing_non_gly_sidechains")
            try:
                return int(value)
            except (TypeError, ValueError):
                return 999999

        sidechains_ok = missing_count("A") == 0 and missing_count("B") == 0 and missing_count("C") == 0
        if not (row["packed_has_chain_A"] and row["packed_has_chain_B"] and row["packed_has_chain_C"] and row["antigen_C_unchanged"] and sidechains_ok):
            row["status"] = "failed_validation"
            row["error"] = "Packed PDB failed chain/antigen/sidechain validation"


def write_summary(rows: list[dict[str, object]], smoke: bool) -> None:
    validate_outputs(rows)
    all_results = RESULTS_DIR / ("smoke_corrected_igg_abmpnn_flowpacker_all_results.csv" if smoke else "corrected_igg_abmpnn_flowpacker_all_results.csv")
    write_csv(all_results, rows)

    total = len(rows)
    packed_ok = sum(1 for r in rows if r.get("packed_pdb") and r.get("status") == "ok")
    abmpnn_ok = sum(1 for r in rows if r.get("optimized_seq_A") and r.get("optimized_seq_B"))
    failed = total - packed_ok
    summary = {
        "timestamp": now(),
        "mode": "smoke" if smoke else "full",
        "input_table": str(REMOTE_PRODIGY),
        "target": "all IgG antibody candidates",
        "total_candidates": total,
        "abmpnn_success": abmpnn_ok,
        "flowpacker_packed_success": packed_ok,
        "failed": failed,
        "designed_chains": "A+B",
        "fixed_antigen_chain": "C",
        "hotspots": HOTSPOTS,
        "mutable_rule": f"Only original antibody Ala-runs with length >= {MUTABLE_ALA_RUN_MIN} are allowed to change; framework and conserved residues are restored.",
        "note": "Corrected AbMPNN + FlowPacker optimization outputs only; not AF3/MD validated.",
        "all_results_csv": str(all_results),
        "packed_pdb_dir": str(RESULTS_DIR / "packed_pdbs"),
    }
    (RESULTS_DIR / "corrected_igg_abmpnn_flowpacker_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    teacher = RESULTS_DIR / "corrected_igg_abmpnn_flowpacker_teacher_summary.md"
    teacher.write_text(
        "\n".join([
            "# Corrected AbMPNN + FlowPacker IgG Optimization Summary",
            "",
            f"- Scope: all existing IgG antibody candidates (`{total}` attempted).",
            "- Method: one constrained AbMPNN optimized sequence per IgG, followed by one FlowPacker model-output side-chain packed complex PDB.",
            "- Mutable rule: only original antibody Ala-placeholder runs are allowed to change; framework/conserved residues are restored.",
            "- Designed chains: antibody chains `A+B`, constrained to mutable placeholder positions.",
            "- Fixed chain: antigen chain `C`.",
            f"- Hotspot annotation: `{HOTSPOTS}`.",
            f"- AbMPNN parsed successfully: `{abmpnn_ok}/{total}`.",
            f"- FlowPacker packed PDB success: `{packed_ok}/{total}`.",
            f"- Failed or incomplete: `{failed}/{total}`.",
            "",
            "Important: these are corrected optimization outputs only. They have not yet been verified by AF3, PRODIGY, or MD.",
            "",
            f"- Full table: `{all_results.name}`",
            "- Packed PDB folder: `packed_pdbs/`",
        ]) + "\n",
        encoding="utf-8",
    )
    update_checkpoint("done_smoke" if smoke else "done", total, total, summary)


def make_package() -> Path:
    zip_path = RESULTS_DIR / "corrected_igg_abmpnn_flowpacker_result_package.zip"
    if zip_path.exists():
        zip_path.unlink()
    include_names = {
        "corrected_igg_abmpnn_flowpacker_all_results.csv",
        "corrected_igg_abmpnn_flowpacker_summary.json",
        "corrected_igg_abmpnn_flowpacker_teacher_summary.md",
        "abmpnn_sequences.fasta",
        "abmpnn_variants.csv",
        "checkpoint.json",
        "partial_results.csv",
        "igg_input_manifest.csv",
    }
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name in include_names:
            p = RESULTS_DIR / name
            if p.exists():
                zf.write(p, arcname=name)
        packed_dir = RESULTS_DIR / "packed_pdbs"
        for pdb in sorted(packed_dir.glob("*.pdb")):
            zf.write(pdb, arcname=f"packed_pdbs/{pdb.name}")
        script_path = SCRIPTS_DIR / Path(__file__).name
        if script_path.exists():
            zf.write(script_path, arcname=f"scripts/{script_path.name}")
        for log_file in sorted(LOG_DIR.glob("45_*")) + sorted(LOG_DIR.glob("46_*")):
            zf.write(log_file, arcname=f"logs/{log_file.name}")
    return zip_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full"], required=True)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-flowpacker", action="store_true")
    args = parser.parse_args()

    smoke = args.mode == "smoke"
    ensure_dirs()
    log(f"START mode={args.mode} limit={args.limit} batch_size={args.batch_size} force={args.force}")
    log(f"Expected remote work dir: {WORK_DIR}")
    log(f"Expected remote results dir: {RESULTS_DIR}")
    manifest = create_input_manifest(args.limit, smoke=smoke)
    if not manifest:
        raise RuntimeError("No valid IgG input PDBs found.")
    update_checkpoint("manifest", 0, len(manifest), {"mode": args.mode})

    abmpnn_dir = run_abmpnn(manifest, smoke=smoke, force=args.force)
    rows = parse_abmpnn_outputs(manifest, abmpnn_dir, smoke=smoke)
    if not args.skip_flowpacker:
        rows = run_flowpacker(rows, smoke=smoke, batch_size=args.batch_size, force=args.force)
    write_summary(rows, smoke=smoke)
    if not smoke:
        zip_path = make_package()
        log(f"Package written: {zip_path}")
    log("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
