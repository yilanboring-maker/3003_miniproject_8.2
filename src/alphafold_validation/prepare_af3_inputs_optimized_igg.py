#!/usr/bin/env python3
"""Prepare AlphaFold3 inputs for optimized IgG candidates.

This script creates AF3 JSON/FASTA files from optimized packed PDBs without
adding templates, coordinates, MSAs, ligands, or constraints to the AF3 inputs.
The design PDBs are copied only for later post-prediction pose validation.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import time
import zipfile
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CANDIDATES = ROOT / "results" / "af3_followup_candidates" / "af3_candidates_delta_g_le_minus9.csv"
DEFAULT_OUTPUT_PARENT = ROOT / "results"
HOTSPOTS = [87, 88, 89, 90, 91, 114, 115, 116, 117]
SMOKE_DESIGNS = {"antibody_b005_18", "antibody_b006_3", "antibody_b009_19"}

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-csv", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument("--run-id", default=datetime.now().strftime("%Y%m%d_%H%M%S"))
    parser.add_argument("--limit", type=int, default=None)
    return parser.parse_args()


def read_candidates(path: Path, limit: int | None) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"candidate CSV not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows = [row for row in rows if row.get("status", "ok") == "ok"]
    if limit is not None:
        rows = rows[:limit]
    return rows


def parse_pdb_sequences(path: Path) -> tuple[dict[str, str], dict[str, int], list[str]]:
    residues_by_chain: dict[str, dict[tuple[int, str], str]] = {}
    with path.open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            chain = line[21].strip()
            if not chain:
                continue
            resn = line[17:20].strip()
            try:
                resi = int(line[22:26])
            except ValueError:
                continue
            icode = line[26].strip()
            residues_by_chain.setdefault(chain, {}).setdefault((resi, icode), resn)

    sequences: dict[str, str] = {}
    unknowns: dict[str, int] = {}
    warnings: list[str] = []
    for chain, residues in residues_by_chain.items():
        seq_chars = []
        unknown_count = 0
        for key in sorted(residues):
            aa = AA3_TO_1.get(residues[key], "X")
            if aa == "X":
                unknown_count += 1
            seq_chars.append(aa)
        sequences[chain] = "".join(seq_chars)
        unknowns[chain] = unknown_count
        if unknown_count:
            warnings.append(f"chain {chain} has {unknown_count} unknown residues")
    return sequences, unknowns, warnings


def write_fasta(path: Path, design_id: str, sequences: dict[str, str]) -> None:
    lines: list[str] = []
    for chain, role in [("A", "antibody_heavy"), ("B", "antibody_light"), ("C", "antigen")]:
        lines.append(f">{design_id}_{role}_chain_{chain}")
        lines.append(sequences[chain])
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_af3_json(path: Path, design_id: str, sequences: dict[str, str]) -> None:
    payload = {
        "name": design_id,
        "modelSeeds": [1],
        "sequences": [
            {"protein": {"id": "A", "sequence": sequences["A"]}},
            {"protein": {"id": "B", "sequence": sequences["B"]}},
            {"protein": {"id": "C", "sequence": sequences["C"]}},
        ],
        "dialect": "alphafold3",
        "version": 1,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="ascii")


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save_checkpoint(path: Path, status: str, done: int, total: int, output_root: Path, message: str = "") -> None:
    payload = {
        "status": status,
        "done": done,
        "total": total,
        "updated_at": datetime.now().isoformat(timespec="seconds"),
        "output_root": str(output_root),
        "message": message,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def make_zip(zip_path: Path, output_root: Path) -> None:
    include_names = {
        "inputs",
        "design_pdbs",
        "af3_input_manifest.csv",
        "af3_smoke_top3_manifest.csv",
        "checkpoint.json",
        "partial_results.csv",
        "af3_input_summary.json",
        "af3_teacher_summary.md",
    }
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for item in output_root.rglob("*"):
            if item == zip_path or item.is_dir():
                continue
            rel = item.relative_to(output_root)
            if rel.parts[0] in include_names or str(rel) in include_names:
                zf.write(item, rel.as_posix())


def main() -> int:
    args = parse_args()
    output_root = args.output_root
    if output_root is None:
        output_root = DEFAULT_OUTPUT_PARENT / f"af3_optimized_igg_delta_g_le_minus9_{args.run_id}"
    if output_root.exists():
        raise FileExistsError(f"output root already exists, refusing to overwrite: {output_root}")

    json_dir = output_root / "inputs" / "json"
    fasta_dir = output_root / "inputs" / "fasta"
    design_dir = output_root / "design_pdbs"
    logs_dir = output_root / "logs"
    for directory in [json_dir, fasta_dir, design_dir, logs_dir]:
        directory.mkdir(parents=True, exist_ok=False)

    candidates = read_candidates(args.candidate_csv, args.limit)
    total = len(candidates)
    start = time.time()
    manifest_rows: list[dict[str, object]] = []
    fieldnames = [
        "af3_priority_rank",
        "design_id",
        "optimized_prodigy_delta_g",
        "optimized_prodigy_kd_m",
        "original_prodigy_delta_g",
        "delta_g_change",
        "source_pdb",
        "design_pdb_copy",
        "af3_json",
        "af3_fasta",
        "chain_A_length",
        "chain_B_length",
        "chain_C_length",
        "chain_A_unknowns",
        "chain_B_unknowns",
        "chain_C_unknowns",
        "hotspots",
        "is_smoke",
        "status",
        "message",
    ]

    partial_csv = output_root / "partial_results.csv"
    checkpoint = output_root / "checkpoint.json"
    save_checkpoint(checkpoint, "running", 0, total, output_root, "starting AF3 input preparation")

    print(f"[START] preparing AF3 input package", flush=True)
    print(f"[INFO] candidates={total}", flush=True)
    print(f"[INFO] output_root={output_root}", flush=True)
    print(f"[INFO] checkpoint={checkpoint}", flush=True)
    print(f"[INFO] partial_results={partial_csv}", flush=True)

    for index, row in enumerate(candidates, start=1):
        design_id = row["design_id"]
        source_pdb = Path(row["local_packed_pdb_path"])
        message = ""
        status = "ok"
        json_path = json_dir / f"{design_id}_af3_input.json"
        fasta_path = fasta_dir / f"{design_id}_af3_sequences.fasta"
        design_copy = design_dir / f"{design_id}_design.pdb"

        try:
            if not source_pdb.exists():
                raise FileNotFoundError(f"source PDB missing: {source_pdb}")
            sequences, unknowns, warnings = parse_pdb_sequences(source_pdb)
            missing_chains = [chain for chain in ["A", "B", "C"] if not sequences.get(chain)]
            if missing_chains:
                raise ValueError(f"missing or empty chains: {','.join(missing_chains)}")
            if warnings:
                message = "; ".join(warnings)
            write_af3_json(json_path, design_id, sequences)
            write_fasta(fasta_path, design_id, sequences)
            shutil.copy2(source_pdb, design_copy)
            lengths = {chain: len(sequences[chain]) for chain in ["A", "B", "C"]}
        except Exception as exc:  # noqa: BLE001
            status = "failed"
            message = str(exc)
            lengths = {"A": 0, "B": 0, "C": 0}
            unknowns = {"A": 0, "B": 0, "C": 0}

        manifest_rows.append(
            {
                "af3_priority_rank": row.get("af3_priority_rank", index),
                "design_id": design_id,
                "optimized_prodigy_delta_g": row.get("optimized_prodigy_delta_g", ""),
                "optimized_prodigy_kd_m": row.get("optimized_prodigy_kd_m", ""),
                "original_prodigy_delta_g": row.get("original_prodigy_delta_g", ""),
                "delta_g_change": row.get("delta_g_change", ""),
                "source_pdb": str(source_pdb),
                "design_pdb_copy": str(design_copy) if design_copy.exists() else "",
                "af3_json": str(json_path) if json_path.exists() else "",
                "af3_fasta": str(fasta_path) if fasta_path.exists() else "",
                "chain_A_length": lengths["A"],
                "chain_B_length": lengths["B"],
                "chain_C_length": lengths["C"],
                "chain_A_unknowns": unknowns.get("A", 0),
                "chain_B_unknowns": unknowns.get("B", 0),
                "chain_C_unknowns": unknowns.get("C", 0),
                "hotspots": "C" + ",C".join(str(item) for item in HOTSPOTS),
                "is_smoke": design_id in SMOKE_DESIGNS,
                "status": status,
                "message": message,
            }
        )

        write_csv(partial_csv, manifest_rows, fieldnames)
        elapsed = time.time() - start
        rate = index / elapsed if elapsed > 0 else 0
        eta = (total - index) / rate if rate > 0 else 0
        print(
            f"[PROGRESS] done={index}/{total} percent={index / total * 100:.1f} "
            f"elapsed={elapsed:.1f}s ETA={eta:.1f}s current={design_id} status={status}",
            flush=True,
        )
        save_checkpoint(checkpoint, "running", index, total, output_root, f"last={design_id} status={status}")

    manifest_csv = output_root / "af3_input_manifest.csv"
    write_csv(manifest_csv, manifest_rows, fieldnames)
    smoke_rows = [row for row in manifest_rows if row["is_smoke"]]
    write_csv(output_root / "af3_smoke_top3_manifest.csv", smoke_rows, fieldnames)

    ok_count = sum(1 for row in manifest_rows if row["status"] == "ok")
    failed_count = total - ok_count
    summary = {
        "status": "done" if failed_count == 0 else "done_with_failures",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "candidate_csv": str(args.candidate_csv),
        "output_root": str(output_root),
        "total_candidates": total,
        "ok_count": ok_count,
        "failed_count": failed_count,
        "json_count": len(list(json_dir.glob("*.json"))),
        "fasta_count": len(list(fasta_dir.glob("*.fasta"))),
        "design_pdb_count": len(list(design_dir.glob("*.pdb"))),
        "smoke_designs": sorted(SMOKE_DESIGNS),
        "hotspots": [f"C{item}" for item in HOTSPOTS],
        "af3_input_policy": "sequence-only AF3 JSON; no templates/MSA/ligands/coordinates/constraints",
    }
    (output_root / "af3_input_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    teacher_summary = [
        "# AF3 optimized IgG input package",
        "",
        f"- Candidates prepared: {ok_count}/{total}",
        "- Selection rule: corrected optimized IgG with PRODIGY delta G <= -9 kcal/mol",
        "- Chains: antibody heavy A, antibody light B, antigen C",
        "- AF3 JSON policy: sequences only; no design complex coordinates/templates were provided to AF3",
        f"- Smoke candidates: {', '.join(sorted(SMOKE_DESIGNS))}",
        f"- Input manifest: {manifest_csv.name}",
        "",
        "This package prepares the AlphaFold3 complex prediction step. It is not an AF3 result yet.",
    ]
    (output_root / "af3_teacher_summary.md").write_text("\n".join(teacher_summary) + "\n", encoding="utf-8")
    zip_path = output_root / "af3_input_package.zip"
    make_zip(zip_path, output_root)
    save_checkpoint(checkpoint, summary["status"], total, total, output_root, "AF3 inputs prepared")

    print(f"[DONE] status={summary['status']} ok={ok_count}/{total} failed={failed_count}", flush=True)
    print(f"[DONE] manifest={manifest_csv}", flush=True)
    print(f"[DONE] zip={zip_path}", flush=True)
    return 0 if failed_count == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
