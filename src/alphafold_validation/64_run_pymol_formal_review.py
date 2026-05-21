#!/usr/bin/env python
"""Run real PyMOL rendering for Task 4 acceptable AF3 candidates.

This script intentionally uses real PyMOL through the conda environment and
stores per-candidate PNG/PSE/log outputs plus checkpointed review tables.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
import zipfile
from datetime import datetime
from pathlib import Path

import pandas as pd


HOTSPOTS = "87+88+89+90+91+114+115+116+117"


def now_text() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def pymol_path(path: Path) -> str:
    return path.resolve().as_posix()


def status_call(row: pd.Series, pymol_ok: bool) -> str:
    if not pymol_ok:
        return "pymol_failed"
    auto = str(row.get("auto_pose_call", "")).strip()
    hs6 = int(row.get("hotspot_contact_6a_count", 0))
    centroid = float(row.get("interface_centroid_distance_to_design_a", 999.0))
    clash = float(row.get("has_clash", 0.0))
    if auto == "acceptable" and hs6 >= 7 and centroid <= 8.5 and clash == 0.0:
        return "confirmed_accept"
    if auto == "acceptable":
        return "uncertain_needs_teacher_review"
    return "reject_after_pymol"


def write_pml(row: pd.Series, pml_path: Path, png_path: Path, pse_path: Path) -> None:
    design_pdb = Path(str(row["design_pdb_path"]))
    af3_pdb = Path(str(row["aligned_pdb_path"]))
    design_id = str(row["design_id"])
    rank = int(row["submission_rank"])
    dg = float(row["optimized_prodigy_delta_g"])
    hs6 = int(row.get("hotspot_contact_6a_count", 0))
    centroid = float(row.get("interface_centroid_distance_to_design_a", 999.0))

    pml = f"""# PyMOL formal visual check for {design_id}
# Design antibody is cyan. AF3 antibody is magenta. Antigen is grey. Hotspots are orange/red.
reinitialize
bg_color white
set ray_opaque_background, off
set cartoon_fancy_helices, 1
set antialias, 2
set transparency_mode, 1
set depth_cue, off

load {pymol_path(design_pdb)}, design
load {pymol_path(af3_pdb)}, af3_aligned

hide everything
show cartoon, design
show cartoon, af3_aligned

color grey80, design and chain C
color cyan, design and chain A+B
color grey60, af3_aligned and chain C
color magenta, af3_aligned and chain A+B
set cartoon_transparency, 0.55, design and chain A+B
set cartoon_transparency, 0.15, af3_aligned and chain A+B
set cartoon_transparency, 0.65, af3_aligned and chain C

select design_hotspot, design and chain C and resi {HOTSPOTS}
select af3_hotspot, af3_aligned and chain C and resi {HOTSPOTS}
show spheres, design_hotspot
show spheres, af3_hotspot
color orange, design_hotspot
color red, af3_hotspot
set sphere_scale, 0.7, design_hotspot or af3_hotspot

select design_interface_ab, design and chain A+B within 6 of design_hotspot
select af3_interface_ab, af3_aligned and chain A+B within 6 of af3_hotspot
show sticks, design_interface_ab or af3_interface_ab
color blue, design_interface_ab
color magenta, af3_interface_ab

distance af3_hotspot_contacts, af3_interface_ab, af3_hotspot, 6
hide labels, af3_hotspot_contacts
color yellow, af3_hotspot_contacts

zoom design_hotspot or af3_hotspot or design_interface_ab or af3_interface_ab, 12
orient design and chain C
png {pymol_path(png_path)}, width=1800, height=1200, dpi=180, ray=1
save {pymol_path(pse_path)}
quit
"""
    pml_path.write_text(pml, encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def package_outputs(root: Path, out_dir: Path, csv_paths: list[Path], summary_path: Path) -> Path:
    zip_path = out_dir / f"pymol_formal_review_result_package_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"
    items: list[Path] = []
    for sub in ["pml", "png", "pse", "logs"]:
        items.extend(sorted((out_dir / sub).glob("*")))
    items.extend([p for p in csv_paths + [summary_path, out_dir / "checkpoint.json"] if p.exists()])
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for item in items:
            zf.write(item, item.relative_to(out_dir.parent).as_posix())
    return zip_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--conda", default=str(Path.home() / "miniforge3" / "Scripts" / "conda.exe"))
    parser.add_argument("--env", default="pymol_review")
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    root = Path(args.root).resolve()
    conda = Path(args.conda).resolve()
    input_csv = root / "af3_pose_review_acceptable.csv"
    out_dir = root / "pymol_manual_review"
    pml_dir = out_dir / "pml"
    png_dir = out_dir / "png"
    pse_dir = out_dir / "pse"
    log_dir = out_dir / "logs"
    for d in [out_dir, pml_dir, png_dir, pse_dir, log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    checkpoint_path = out_dir / "checkpoint.json"
    all_csv = out_dir / "pymol_manual_review_acceptable49.csv"
    confirmed_csv = out_dir / "pymol_confirmed_accept.csv"
    uncertain_csv = out_dir / "pymol_uncertain_needs_teacher_review.csv"
    reject_csv = out_dir / "pymol_reject_after_review.csv"
    summary_path = out_dir / "pymol_teacher_summary.md"

    df = pd.read_csv(input_csv).sort_values("optimized_prodigy_delta_g", ascending=True).reset_index(drop=True)
    if args.limit:
        df = df.head(args.limit).copy()
    total = len(df)
    print(f"[{now_text()}] PyMOL formal review start")
    print(f"root={root}")
    print(f"input_csv={input_csv}")
    print(f"out_dir={out_dir}")
    print(f"conda={conda}")
    print(f"env={args.env}")
    print(f"total={total}", flush=True)

    results: list[dict] = []
    start = time.time()
    for idx, row in df.iterrows():
        done = idx + 1
        rank = int(row["submission_rank"])
        design_id = str(row["design_id"])
        stem = f"{rank:03d}_{design_id}_pymol_formal"
        pml_path = pml_dir / f"{stem}.pml"
        png_path = png_dir / f"{stem}.png"
        pse_path = pse_dir / f"{stem}.pse"
        log_path = log_dir / f"{stem}.log"
        write_pml(row, pml_path, png_path, pse_path)

        cmd = [str(conda), "run", "-n", args.env, "pymol", "-cq", str(pml_path)]
        elapsed = time.time() - start
        pct = 100.0 * (done - 1) / total if total else 100.0
        print(f"[{now_text()}] running {done}/{total} ({pct:.1f}%) rank={rank} design={design_id}", flush=True)
        proc = subprocess.run(cmd, text=True, capture_output=True)
        log_path.write_text(
            "COMMAND: " + " ".join(cmd) + "\n\nSTDOUT:\n" + proc.stdout + "\n\nSTDERR:\n" + proc.stderr,
            encoding="utf-8",
        )
        pymol_ok = proc.returncode == 0 and png_path.exists() and png_path.stat().st_size > 1000 and pse_path.exists()
        call = status_call(row, pymol_ok)
        result = {
            "submission_rank": rank,
            "design_id": design_id,
            "optimized_prodigy_delta_g": row.get("optimized_prodigy_delta_g", ""),
            "optimized_prodigy_kd_m": row.get("optimized_prodigy_kd_m", ""),
            "auto_pose_call": row.get("auto_pose_call", ""),
            "hotspot_contact_6a_count": row.get("hotspot_contact_6a_count", ""),
            "hotspot_contact_5a_count": row.get("hotspot_contact_5a_count", ""),
            "min_hotspot_distance_a": row.get("min_hotspot_distance_a", ""),
            "interface_centroid_distance_to_design_a": row.get("interface_centroid_distance_to_design_a", ""),
            "has_clash": row.get("has_clash", ""),
            "pymol_run_status": "success" if pymol_ok else "failed",
            "pymol_formal_call": call,
            "pml_path": str(pml_path),
            "pymol_png_path": str(png_path),
            "pymol_pse_path": str(pse_path),
            "pymol_log_path": str(log_path),
            "review_note": "Real PyMOL rendering completed; final call combines PyMOL-rendered view availability with hotspot/contact triage metrics.",
        }
        results.append(result)
        write_csv(all_csv, results)
        checkpoint = {
            "updated_at": now_text(),
            "stage": "pymol_formal_review_running",
            "done": done,
            "total": total,
            "last_design_id": design_id,
            "last_call": call,
            "all_results": str(all_csv),
        }
        checkpoint_path.write_text(json.dumps(checkpoint, indent=2), encoding="utf-8")
        elapsed = time.time() - start
        rate = done / elapsed if elapsed > 0 else 0.0
        eta = (total - done) / rate if rate > 0 else 0.0
        print(
            f"[{now_text()}] done={done}/{total} elapsed={elapsed:.1f}s ETA={eta:.1f}s "
            f"status={'success' if pymol_ok else 'failed'} call={call}",
            flush=True,
        )

    confirmed = [r for r in results if r["pymol_formal_call"] == "confirmed_accept"]
    uncertain = [r for r in results if r["pymol_formal_call"] == "uncertain_needs_teacher_review"]
    reject = [r for r in results if r["pymol_formal_call"] == "reject_after_pymol"]
    failed = [r for r in results if r["pymol_formal_call"] == "pymol_failed"]
    for path, rows in [(confirmed_csv, confirmed), (uncertain_csv, uncertain), (reject_csv, reject)]:
        write_csv(path, rows if rows else results[:1])
        if not rows:
            path.write_text(",".join(results[0].keys()) + "\n", encoding="utf-8")

    top_lines = []
    for r in confirmed[:12]:
        top_lines.append(
            f"- rank {r['submission_rank']}: `{r['design_id']}`, "
            f"dG={float(r['optimized_prodigy_delta_g']):.3f} kcal/mol, "
            f"HS6A={int(r['hotspot_contact_6a_count'])}/9, "
            f"centroid={float(r['interface_centroid_distance_to_design_a']):.2f} A"
        )
    summary = "\n".join(
        [
            "# Task 4 PyMOL Formal Review Summary",
            "",
            f"Updated: {now_text()}",
            "",
            "Scope: 49 AF3 candidates previously classified as automated `acceptable`.",
            "",
            "PyMOL environment:",
            f"- conda env: `{args.env}`",
            "- PyMOL: conda-forge `pymol-open-source`",
            "",
            "Results:",
            f"- PyMOL rendered successfully: {sum(1 for r in results if r['pymol_run_status'] == 'success')}/{len(results)}",
            f"- confirmed_accept: {len(confirmed)}",
            f"- uncertain_needs_teacher_review: {len(uncertain)}",
            f"- reject_after_pymol: {len(reject)}",
            f"- pymol_failed: {len(failed)}",
            "",
            "Important caveat:",
            "- The screenshots and `.pse` sessions are generated by real PyMOL. The conclusion table still uses objective hotspot/contact metrics to standardize calls; borderline visual cases should be shown to the teacher.",
            "",
            "Top confirmed candidates:",
            "",
            *top_lines,
            "",
            "Core files:",
            f"- all results: `{all_csv}`",
            f"- confirmed: `{confirmed_csv}`",
            f"- screenshots: `{png_dir}`",
            f"- PyMOL sessions: `{pse_dir}`",
        ]
    )
    summary_path.write_text(summary, encoding="utf-8")

    zip_path = package_outputs(root, out_dir, [all_csv, confirmed_csv, uncertain_csv, reject_csv], summary_path)
    final_checkpoint = {
        "updated_at": now_text(),
        "stage": "pymol_formal_review_complete",
        "done": len(results),
        "total": total,
        "confirmed_accept": len(confirmed),
        "uncertain_needs_teacher_review": len(uncertain),
        "reject_after_pymol": len(reject),
        "pymol_failed": len(failed),
        "package": str(zip_path),
    }
    checkpoint_path.write_text(json.dumps(final_checkpoint, indent=2), encoding="utf-8")
    print(f"[{now_text()}] complete package={zip_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
