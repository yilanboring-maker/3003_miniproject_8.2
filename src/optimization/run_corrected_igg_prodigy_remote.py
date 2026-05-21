#!/usr/bin/env python3
"""Run PRODIGY on corrected optimized IgG AbMPNN+FlowPacker complexes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path


RT_KCAL_25C = 0.0019872041 * 298.15
DEFAULT_INPUT_DIR = Path("/tmp/ica2_abmpnn_flowpacker/results/abmpnn_flowpacker_all_igg_corrected")
DEFAULT_RESULTS_DIR = Path("/tmp/ica2_abmpnn_flowpacker/results/abmpnn_flowpacker_all_igg_corrected_prodigy")
DEFAULT_LOG_DIR = Path("/mnt/PPIFlow/ica2_runs/logs")
SCRIPT_NAME = "46_run_corrected_igg_prodigy_remote.py"
AA1_TO_3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(message: str) -> None:
    print(f"[{now()}] {message}", flush=True)


def parse_float(value: object, default: float = math.nan) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return default


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    tmp.replace(path)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)


def output_fields(rows: list[dict[str, object]]) -> list[str]:
    if not rows:
        return []
    return [field for field in rows[0].keys() if not field.startswith("_")]


def parse_delta_g(stdout: str) -> float:
    for line in stdout.splitlines():
        parts = line.strip().split()
        if not parts:
            continue
        for token in reversed(parts):
            try:
                return float(token)
            except ValueError:
                continue
    raise ValueError(f"No numeric PRODIGY value found in output: {stdout!r}")


def kd_from_delta_g(delta_g: float) -> float:
    return math.exp(delta_g / RT_KCAL_25C)


def kd_label(kd_m: float) -> str:
    if not math.isfinite(kd_m):
        return ""
    if kd_m < 1e-9:
        return f"{kd_m * 1e12:.2f} pM"
    if kd_m < 1e-6:
        return f"{kd_m * 1e9:.2f} nM"
    if kd_m < 1e-3:
        return f"{kd_m * 1e6:.2f} uM"
    return f"{kd_m:.2e} M"


def unique_by_design(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    seen: set[str] = set()
    out: list[dict[str, str]] = []
    for row in rows:
        design_id = row.get("design_id", "")
        if not design_id or design_id in seen:
            continue
        seen.add(design_id)
        out.append(row)
    return out


def chains_in_pdb(path: Path) -> set[str]:
    chains: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                chains.add(line[21])
    return chains


def build_rows(input_dir: Path, mode: str) -> list[dict[str, object]]:
    source_csv = input_dir / "corrected_igg_abmpnn_flowpacker_all_results.csv"
    if mode == "smoke" and not source_csv.exists():
        source_csv = input_dir / "smoke_corrected_igg_abmpnn_flowpacker_all_results.csv"
    packed_dir = input_dir / "packed_pdbs"
    if not source_csv.exists():
        raise FileNotFoundError(f"Missing optimized IgG result CSV: {source_csv}")
    if not packed_dir.exists():
        raise FileNotFoundError(f"Missing packed PDB directory: {packed_dir}")

    source_rows = unique_by_design(load_csv(source_csv))
    if len(source_rows) != 200 and mode == "full":
        raise RuntimeError(f"Expected 200 optimized IgG rows, found {len(source_rows)}")
    if mode == "smoke":
        source_rows = [source_rows[0], source_rows[len(source_rows) // 2], source_rows[-1]]

    rows: list[dict[str, object]] = []
    for row in source_rows:
        design_id = row["design_id"]
        pdb_path = packed_dir / f"{design_id}_opt_001.pdb"
        original = parse_float(row.get("original_prodigy_delta_g", ""))
        original_rank = row.get("original_prodigy_rank") or row.get("prodigy_rank") or ""
        if not original_rank:
            original_rank = row.get("rank", "")
        rows.append(
            {
                "design_id": design_id,
                "original_prodigy_delta_g": f"{original:.3f}" if math.isfinite(original) else "",
                "original_prodigy_rank": original_rank,
                "mutation_count_total": row.get("mutation_count_total", ""),
                "optimized_prodigy_delta_g": "",
                "optimized_prodigy_kd_m": "",
                "optimized_prodigy_kd_label": "",
                "optimized_prodigy_rank": "",
                "delta_g_change": "",
                "status": "pending",
                "pdb_path": str(pdb_path),
                "prodigy_input_pdb_path": "",
                "prodigy_command": "",
                "pdb_normalization_status": "not_checked",
                "pdb_normalization_notes": "",
                "prodigy_stdout": "",
                "prodigy_stderr": "",
                "_expected_seq_A": row.get("optimized_seq_A", ""),
                "_expected_seq_B": row.get("optimized_seq_B", ""),
                "_expected_seq_C": row.get("original_seq_C", ""),
                "_fallback_seq_A": row.get("original_seq_A", ""),
                "_fallback_seq_B": row.get("original_seq_B", ""),
                "_fallback_seq_C": row.get("original_seq_C", ""),
            }
        )
    return rows


def add_ranks(rows: list[dict[str, object]]) -> None:
    scored = [r for r in rows if r.get("status") == "ok"]
    scored.sort(key=lambda r: parse_float(r.get("optimized_prodigy_delta_g", ""), 9999.0))
    for rank, row in enumerate(scored, start=1):
        row["optimized_prodigy_rank"] = rank


def residue_from_sequence(seq: str, residue_number: str) -> str:
    try:
        index = int(residue_number.strip()) - 1
    except ValueError:
        return ""
    if index < 0 or index >= len(seq):
        return ""
    return AA1_TO_3.get(seq[index].upper(), "")


def expected_residue_name(row: dict[str, object], chain: str, residue_number: str) -> tuple[str, str]:
    expected = residue_from_sequence(str(row.get(f"_expected_seq_{chain}", "")), residue_number)
    if expected:
        return expected, "expected_sequence"
    fallback = residue_from_sequence(str(row.get(f"_fallback_seq_{chain}", "")), residue_number)
    if fallback:
        return fallback, "original_sequence_fallback"
    return "", ""


def normalize_unk_residues_for_prodigy(row: dict[str, object], results_dir: Path) -> Path:
    source = Path(str(row["pdb_path"]))
    row["prodigy_input_pdb_path"] = str(source)
    if not source.exists():
        row["pdb_normalization_status"] = "missing_pdb"
        row["pdb_normalization_notes"] = f"Missing PDB: {source}"
        return source

    has_unk = False
    with source.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 20 and line[17:20] == "UNK":
                has_unk = True
                break
    if not has_unk:
        row["pdb_normalization_status"] = "raw_pdb"
        row["pdb_normalization_notes"] = "No UNK residues detected."
        return source

    out_dir = results_dir / "prodigy_input_pdbs"
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / f"{row.get('design_id')}_prodigy_input.pdb"
    changes: dict[tuple[str, str], str] = {}
    used_fallback = False
    unresolved: set[tuple[str, str]] = set()
    out_lines: list[str] = []
    with source.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 26 and line[17:20] == "UNK":
                chain = line[21]
                residue_number = line[22:26]
                mapped, source_name = expected_residue_name(row, chain, residue_number)
                key = (chain, residue_number.strip())
                if mapped:
                    line = line[:17] + f"{mapped:>3}" + line[20:]
                    if source_name == "original_sequence_fallback":
                        used_fallback = True
                        changes[key] = f"{mapped}(original_sequence_fallback)"
                    else:
                        changes[key] = mapped
                else:
                    unresolved.add(key)
            out_lines.append(line)
    target.write_text("".join(out_lines), encoding="utf-8")
    row["prodigy_input_pdb_path"] = str(target)
    if unresolved:
        row["pdb_normalization_status"] = "unresolved_unk"
    elif used_fallback:
        row["pdb_normalization_status"] = "renamed_unk_from_original_sequence_fallback"
    else:
        row["pdb_normalization_status"] = "renamed_unk_from_expected_sequence"
    changed_notes = [f"{chain}{resnum}->{resname}" for (chain, resnum), resname in sorted(changes.items())]
    unresolved_notes = [f"{chain}{resnum}" for chain, resnum in sorted(unresolved)]
    notes = []
    if changed_notes:
        notes.append("renamed " + ";".join(changed_notes))
    if unresolved_notes:
        notes.append("unresolved " + ";".join(unresolved_notes))
    row["pdb_normalization_notes"] = " | ".join(notes)
    return target


def score_one(row: dict[str, object], timeout: int, results_dir: Path, normalize_unk: bool) -> None:
    pdb_path = Path(str(row["pdb_path"]))
    if not pdb_path.exists():
        row["status"] = "missing_pdb"
        row["prodigy_stderr"] = f"Missing PDB: {pdb_path}"
        return
    if normalize_unk:
        pdb_path = normalize_unk_residues_for_prodigy(row, results_dir)
    else:
        row["prodigy_input_pdb_path"] = str(pdb_path)
        row["pdb_normalization_status"] = "disabled"
        row["pdb_normalization_notes"] = ""
    chains = chains_in_pdb(pdb_path)
    if not {"A", "B", "C"}.issubset(chains):
        row["status"] = "missing_required_chain"
        row["prodigy_stderr"] = f"Required chains A/B/C not all present; found {sorted(chains)}"
        return

    cmd = ["prodigy", "-q", str(pdb_path), "--selection", "A,B", "C"]
    row["prodigy_command"] = " ".join(cmd)
    try:
        completed = subprocess.run(
            cmd,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
        row["prodigy_stdout"] = completed.stdout.strip()
        row["prodigy_stderr"] = completed.stderr.strip()
        if completed.returncode != 0:
            row["status"] = f"failed_exit_{completed.returncode}"
            return
        delta_g = parse_delta_g(completed.stdout)
        kd_m = kd_from_delta_g(delta_g)
        original = parse_float(row.get("original_prodigy_delta_g", ""))
        row["optimized_prodigy_delta_g"] = f"{delta_g:.3f}"
        row["optimized_prodigy_kd_m"] = f"{kd_m:.6e}"
        row["optimized_prodigy_kd_label"] = kd_label(kd_m)
        row["delta_g_change"] = f"{(delta_g - original):.3f}" if math.isfinite(original) else ""
        row["status"] = "ok"
    except Exception as exc:
        row["status"] = "exception"
        row["prodigy_stderr"] = repr(exc)


def write_summary(results_dir: Path, rows: list[dict[str, object]], mode: str, input_dir: Path) -> None:
    scored = [r for r in rows if r.get("status") == "ok"]
    best = min(scored, key=lambda r: parse_float(r.get("optimized_prodigy_delta_g", ""), 9999.0)) if scored else None
    improved = [
        r for r in scored
        if parse_float(r.get("delta_g_change", ""), 0.0) < 0
    ]
    payload = {
        "timestamp": now(),
        "mode": mode,
        "input_dir": str(input_dir),
        "total_candidates": len(rows),
        "prodigy_success": len(scored),
        "failed": len(rows) - len(scored),
        "selection": "A,B C",
        "pdb_normalization": "UNK residue names are restored from optimized AbMPNN sequences when possible; if the optimized sequence contains X, the original residue name at the same chain/residue index is used as a transparent fallback. No atoms or coordinates are added or changed.",
        "best_design_id": best.get("design_id") if best else "",
        "best_optimized_delta_g": best.get("optimized_prodigy_delta_g") if best else "",
        "best_optimized_kd_m": best.get("optimized_prodigy_kd_m") if best else "",
        "improved_vs_original_count": len(improved),
        "note": "PRODIGY rescoring only; not AF3 or MD validated.",
    }
    write_json(results_dir / "corrected_igg_prodigy_summary.json", payload)

    top20_path = results_dir / "corrected_igg_prodigy_top20.csv"
    scored_sorted = sorted(scored, key=lambda r: parse_float(r.get("optimized_prodigy_delta_g", ""), 9999.0))
    fields = output_fields(rows)
    write_csv(top20_path, scored_sorted[:20], fields)

    lines = [
        "# Corrected Optimized IgG PRODIGY Rescoring Summary",
        "",
        f"- Candidates scored: {len(scored)}/{len(rows)}",
        "- Selection: antibody chains `A+B` vs antigen chain `C`",
        "- Command: `prodigy -q <packed_pdb> --selection A,B C`",
        "- PDB input note: `UNK` residue names, when present, were restored from the AbMPNN optimized sequence when possible; if the optimized sequence contained `X`, the original residue name at the same chain/residue index was used as a transparent fallback. No atoms or coordinates were added or moved.",
        "- Scope: corrected AbMPNN + FlowPacker optimized IgG structures only; no AF3/MD validation in this step.",
    ]
    if best:
        lines.extend(
            [
                f"- Top optimized IgG: `{best.get('design_id')}`",
                f"- Top optimized PRODIGY delta G: `{best.get('optimized_prodigy_delta_g')} kcal/mol`",
                f"- Top optimized predicted Kd: `{best.get('optimized_prodigy_kd_label')}`",
            ]
        )
    lines.append(f"- Improved vs original PRODIGY count: `{len(improved)}`")
    (results_dir / "corrected_igg_prodigy_teacher_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def package_results(results_dir: Path) -> Path:
    package = results_dir / "corrected_igg_prodigy_result_package.zip"
    members = [
        "corrected_igg_prodigy_all_results.csv",
        "corrected_igg_prodigy_top20.csv",
        "corrected_igg_prodigy_summary.json",
        "corrected_igg_prodigy_teacher_summary.md",
        "checkpoint.json",
        "partial_results.csv",
    ]
    with zipfile.ZipFile(package, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name in members:
            path = results_dir / name
            if path.exists():
                zf.write(path, arcname=name)
        logs_dir = results_dir / "logs"
        if logs_dir.exists():
            for path in logs_dir.rglob("*"):
                if path.is_file():
                    zf.write(path, arcname=str(path.relative_to(results_dir)))
        scripts_dir = results_dir / "scripts"
        if scripts_dir.exists():
            for path in scripts_dir.rglob("*"):
                if path.is_file():
                    zf.write(path, arcname=str(path.relative_to(results_dir)))
        pdb_input_dir = results_dir / "prodigy_input_pdbs"
        if pdb_input_dir.exists():
            for path in pdb_input_dir.rglob("*"):
                if path.is_file():
                    zf.write(path, arcname=str(path.relative_to(results_dir)))
    return package


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full"], default="full")
    parser.add_argument("--input-dir", default=str(DEFAULT_INPUT_DIR))
    parser.add_argument("--results-dir", default=str(DEFAULT_RESULTS_DIR))
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--normalize-unk-from-sequence", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    results_dir = Path(args.results_dir)
    if args.force and results_dir.exists():
        shutil.rmtree(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "logs").mkdir(exist_ok=True)
    (results_dir / "scripts").mkdir(exist_ok=True)
    DEFAULT_LOG_DIR.mkdir(parents=True, exist_ok=True)

    log(f"START mode={args.mode}")
    log(f"Input dir: {input_dir}")
    log(f"Results dir: {results_dir}")
    if shutil.which("prodigy") is None:
        raise RuntimeError("PRODIGY command not found in PATH")
    log(f"PRODIGY: {shutil.which('prodigy')}")

    script_src = Path(__file__)
    if script_src.exists():
        shutil.copy2(script_src, results_dir / "scripts" / script_src.name)

    rows = build_rows(input_dir, args.mode)
    fields = output_fields(rows)
    all_csv = results_dir / ("smoke_corrected_igg_prodigy_all_results.csv" if args.mode == "smoke" else "corrected_igg_prodigy_all_results.csv")
    partial_csv = results_dir / "partial_results.csv"
    checkpoint = results_dir / "checkpoint.json"

    total = len(rows)
    start = time.time()
    best_delta = math.inf
    best_id = ""
    write_json(checkpoint, {"script": SCRIPT_NAME, "status": "running", "mode": args.mode, "done": 0, "total": total})

    for index, row in enumerate(rows, start=1):
        elapsed = time.time() - start
        rate = (index - 1) / elapsed if index > 1 and elapsed > 0 else 0.0
        eta = (total - index + 1) / rate if rate > 0 else -1.0
        log(
            f"PRODIGY {index}/{total} {100.0 * index / max(total, 1):.1f}% "
            f"elapsed={elapsed:.1f}s eta={eta:.1f}s best={best_id}:{best_delta if math.isfinite(best_delta) else 'NA'} "
            f"design={row.get('design_id')}"
        )

        score_one(row, args.timeout, results_dir, args.normalize_unk_from_sequence)
        delta = parse_float(row.get("optimized_prodigy_delta_g", ""), math.inf)
        if row.get("status") == "ok" and delta < best_delta:
            best_delta = delta
            best_id = str(row.get("design_id", ""))
        add_ranks(rows[:index])
        write_csv(partial_csv, rows[:index], fields)
        write_json(
            checkpoint,
            {
                "script": SCRIPT_NAME,
                "status": "running",
                "mode": args.mode,
                "done": index,
                "total": total,
                "ok": sum(1 for r in rows[:index] if r.get("status") == "ok"),
                "best_design_id": best_id,
                "best_optimized_delta_g": best_delta if math.isfinite(best_delta) else "",
                "updated_at": now(),
            },
        )

    add_ranks(rows)
    scored_sorted = sorted(
        rows,
        key=lambda r: parse_float(r.get("optimized_prodigy_delta_g", ""), 9999.0),
    )
    write_csv(all_csv, scored_sorted, fields)
    if args.mode == "smoke":
        write_csv(results_dir / "corrected_igg_prodigy_all_results.csv", scored_sorted, fields)
    write_summary(results_dir, scored_sorted, args.mode, input_dir)
    write_json(
        checkpoint,
        {
            "script": SCRIPT_NAME,
            "status": "done",
            "mode": args.mode,
            "done": total,
            "total": total,
            "ok": sum(1 for r in rows if r.get("status") == "ok"),
            "best_design_id": best_id,
            "best_optimized_delta_g": best_delta if math.isfinite(best_delta) else "",
            "updated_at": now(),
        },
    )
    package = package_results(results_dir)
    log(f"Package written: {package}")
    log("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
