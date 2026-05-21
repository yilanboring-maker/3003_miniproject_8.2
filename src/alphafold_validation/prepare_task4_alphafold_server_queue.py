#!/usr/bin/env python3
"""Prepare a Task 4 AlphaFold Server submission queue for optimized IgG."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CANDIDATES = ROOT / "results" / "af3_followup_candidates" / "af3_candidates_delta_g_le_minus9.csv"
DEFAULT_AF3_INPUT_ROOT = ROOT / "results" / "af3_optimized_igg_delta_g_le_minus9_20260517_105733"
DEFAULT_OUTPUT_PARENT = ROOT / "results"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-csv", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--af3-input-root", type=Path, default=DEFAULT_AF3_INPUT_ROOT)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--run-id", default=datetime.now().strftime("%Y%m%d_%H%M%S"))
    return parser.parse_args()


def read_candidates(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows = [row for row in rows if row.get("status", "ok") == "ok"]
    rows.sort(key=lambda row: float(row["optimized_prodigy_delta_g"]))
    return rows


def priority_batch(delta_g: float, rank: int) -> str:
    if rank <= 3:
        return "01_smoke_top3"
    if delta_g <= -12.0:
        return "02_remaining_le_minus12"
    if delta_g <= -9.55:
        return "03_minus12_to_minus9p55"
    return "04_minus9p55_to_minus9"


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    output_root = args.output_root or (
        DEFAULT_OUTPUT_PARENT / f"task4_alphafold_server_validation_{args.run_id}"
    )
    if output_root.exists():
        raise FileExistsError(f"refusing to overwrite existing output root: {output_root}")

    dirs = [
        output_root / "server_inputs" / "fasta",
        output_root / "server_inputs" / "json",
        output_root / "design_pdbs",
        output_root / "af3_raw_downloads",
        output_root / "pymol_sessions",
        output_root / "pymol_screenshots",
        output_root / "logs",
        output_root / "review_tables",
    ]
    for directory in dirs:
        directory.mkdir(parents=True, exist_ok=False)

    candidates = read_candidates(args.candidate_csv)
    manifest_rows: list[dict[str, object]] = []
    for rank, row in enumerate(candidates, start=1):
        design_id = row["design_id"]
        delta_g = float(row["optimized_prodigy_delta_g"])
        fasta_src = args.af3_input_root / "inputs" / "fasta" / f"{design_id}_af3_sequences.fasta"
        json_src = args.af3_input_root / "inputs" / "json" / f"{design_id}_af3_input.json"
        pdb_src = args.af3_input_root / "design_pdbs" / f"{design_id}_design.pdb"
        if not fasta_src.exists():
            raise FileNotFoundError(f"missing FASTA for {design_id}: {fasta_src}")
        if not json_src.exists():
            raise FileNotFoundError(f"missing JSON for {design_id}: {json_src}")
        if not pdb_src.exists():
            raise FileNotFoundError(f"missing design PDB for {design_id}: {pdb_src}")
        fasta_dst = output_root / "server_inputs" / "fasta" / fasta_src.name
        json_dst = output_root / "server_inputs" / "json" / json_src.name
        pdb_dst = output_root / "design_pdbs" / pdb_src.name
        shutil.copy2(fasta_src, fasta_dst)
        shutil.copy2(json_src, json_dst)
        shutil.copy2(pdb_src, pdb_dst)
        manifest_rows.append(
            {
                "submission_rank": rank,
                "priority_batch": priority_batch(delta_g, rank),
                "design_id": design_id,
                "optimized_prodigy_delta_g": f"{delta_g:.3f}",
                "optimized_prodigy_kd_m": row.get("optimized_prodigy_kd_m", ""),
                "original_prodigy_delta_g": row.get("original_prodigy_delta_g", ""),
                "delta_g_change": row.get("delta_g_change", ""),
                "heavy_chain_id": "A",
                "light_chain_id": "B",
                "antigen_chain_id": "C",
                "hotspots": "C87,C88,C89,C90,C91,C114,C115,C116,C117",
                "fasta_path": str(fasta_dst),
                "af3_json_path": str(json_dst),
                "design_pdb_path": str(pdb_dst),
                "alphafold_server_job_name": "",
                "alphafold_server_status": "not_submitted",
                "submitted_at": "",
                "completed_at": "",
                "download_path": "",
                "pose_review_status": "not_reviewed",
                "notes": "",
            }
        )

    fields = [
        "submission_rank",
        "priority_batch",
        "design_id",
        "optimized_prodigy_delta_g",
        "optimized_prodigy_kd_m",
        "original_prodigy_delta_g",
        "delta_g_change",
        "heavy_chain_id",
        "light_chain_id",
        "antigen_chain_id",
        "hotspots",
        "fasta_path",
        "af3_json_path",
        "design_pdb_path",
        "alphafold_server_job_name",
        "alphafold_server_status",
        "submitted_at",
        "completed_at",
        "download_path",
        "pose_review_status",
        "notes",
    ]
    manifest = output_root / "af3_server_job_manifest.csv"
    write_csv(manifest, manifest_rows, fields)
    write_csv(output_root / "partial_results.csv", manifest_rows, fields)

    review_fields = [
        "design_id",
        "optimized_prodigy_delta_g",
        "alphafold_server_job_name",
        "pose_review_status",
        "pymol_session",
        "screenshot_design",
        "screenshot_af3",
        "screenshot_overlay",
        "confidence_summary",
        "manual_notes",
    ]
    for filename in [
        "af3_pose_review_all_results.csv",
        "af3_pose_review_acceptable.csv",
        "af3_pose_review_needs_review.csv",
        "af3_pose_review_reject.csv",
    ]:
        write_csv(output_root / filename, [], review_fields)

    batch_counts: dict[str, int] = {}
    for row in manifest_rows:
        batch_counts[str(row["priority_batch"])] = batch_counts.get(str(row["priority_batch"]), 0) + 1

    checkpoint = {
        "status": "queue_prepared",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "output_root": str(output_root),
        "total_candidates": len(manifest_rows),
        "submitted": 0,
        "completed": 0,
        "reviewed": 0,
        "batch_counts": batch_counts,
        "next_action": "Open AlphaFold Server and submit Top 3 smoke jobs.",
    }
    (output_root / "checkpoint.json").write_text(json.dumps(checkpoint, indent=2), encoding="utf-8")

    summary = [
        "# Task 4 AlphaFold Server validation",
        "",
        f"- Queue prepared: {len(manifest_rows)} optimized IgG candidates with delta G <= -9 kcal/mol.",
        "- Submission order: Top 3 smoke, remaining <= -12, then <= -9.55, then the rest down to -9.",
        "- Input policy: sequence-only heavy A, light B, antigen C; no design coordinates/templates.",
        "- Review policy: PyMOL/manual pose review is primary; contact/distance metrics are auxiliary.",
        f"- Manifest: {manifest.name}",
    ]
    (output_root / "task4_teacher_summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    print(f"[DONE] queue_prepared={len(manifest_rows)}", flush=True)
    print(f"[DONE] output_root={output_root}", flush=True)
    print(f"[DONE] manifest={manifest}", flush=True)
    print(f"[DONE] batch_counts={json.dumps(batch_counts, sort_keys=True)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
