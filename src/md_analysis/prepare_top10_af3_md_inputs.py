#!/usr/bin/env python
"""Prepare AF3/PyMOL-confirmed strict candidates for Top10 MD.

The script does not run PRODIGY or MD. It creates a durable input package with
the strict 9/9 AF3 structures, a preliminary Top10 list, and a manifest that a
remote Matpool/GROMACS workflow can consume.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import time
import zipfile
from datetime import datetime
from pathlib import Path

import pandas as pd


def now_text() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def safe_name(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=r"C:\Users\lyl\Desktop\CMML2\ica2_antibody_design")
    parser.add_argument("--timestamp", default=datetime.now().strftime("%Y%m%d_%H%M%S"))
    args = parser.parse_args()

    project = Path(args.project_root).resolve()
    task4 = project / "results" / "task4_alphafold_server_validation_20260517_112024"
    manual_csv = task4 / "pymol_manual_review" / "pymol_manual_confirmation_notes.csv"
    af3_all_csv = task4 / "af3_pose_review_all_results.csv"
    result_root = project / "results" / f"md_top10_af3_confirmed_10ns_{args.timestamp}"
    input_dir = result_root / "inputs"
    af3_pdb_dir = input_dir / "af3_confirmed_pdbs"
    refs_dir = input_dir / "reference_design_pdbs"
    scripts_dir = result_root / "remote_scripts"
    for directory in [result_root, input_dir, af3_pdb_dir, refs_dir, scripts_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    checkpoint_path = result_root / "checkpoint.json"
    print(f"[{now_text()}] prepare Top10 AF3/PyMOL confirmed MD inputs")
    print(f"project={project}")
    print(f"task4={task4}")
    print(f"result_root={result_root}")
    print(f"checkpoint={checkpoint_path}")

    manual = pd.read_csv(manual_csv)
    af3 = pd.read_csv(af3_all_csv)
    strict = manual[manual["manual_confidence_group"].eq("high_confidence_9of9")].copy()
    strict["optimized_prodigy_delta_g"] = pd.to_numeric(strict["optimized_prodigy_delta_g"], errors="coerce")
    strict = strict.sort_values("optimized_prodigy_delta_g", ascending=True).reset_index(drop=True)
    merged = strict.merge(
        af3[
            [
                "design_id",
                "aligned_pdb_path",
                "af3_model_cif_path",
                "design_pdb_path",
                "ranking_score",
                "iptm",
                "ptm",
                "fraction_disordered",
                "has_clash",
                "raw_download_path",
            ]
        ],
        on="design_id",
        how="left",
        suffixes=("", "_af3"),
    )

    rows: list[dict] = []
    start = time.time()
    total = len(merged)
    for idx, row in merged.iterrows():
        done = idx + 1
        design_id = str(row["design_id"])
        submission_rank = int(row["submission_rank"])
        local_name = f"{submission_rank:03d}_{safe_name(design_id)}_af3_confirmed.pdb"
        ref_name = f"{submission_rank:03d}_{safe_name(design_id)}_design_reference.pdb"
        src = Path(str(row["aligned_pdb_path"]))
        ref_src = Path(str(row["design_pdb_path"]))
        dst = af3_pdb_dir / local_name
        ref_dst = refs_dir / ref_name
        if not src.exists():
            raise FileNotFoundError(f"Missing aligned AF3 PDB for {design_id}: {src}")
        shutil.copy2(src, dst)
        if ref_src.exists():
            shutil.copy2(ref_src, ref_dst)
        remote_input_pdb = f"/mnt/PPIFlow/ica2_runs/input/md_top10_af3_confirmed_10ns/{local_name}"
        rows.append(
            {
                "strict_rank_by_optimized_prodigy": idx + 1,
                "submission_rank": submission_rank,
                "design_id": design_id,
                "optimized_prodigy_delta_g": f"{float(row['optimized_prodigy_delta_g']):.3f}",
                "optimized_prodigy_kd_m": row.get("optimized_prodigy_kd_m", ""),
                "hotspot_contact_6a_count": int(row["hotspot_contact_6a_count"]),
                "hotspot_contact_5a_count": int(row["hotspot_contact_5a_count"]),
                "interface_centroid_distance_to_design_a": f"{float(row['interface_centroid_distance_to_design_a']):.6f}",
                "min_hotspot_distance_a": f"{float(row['min_hotspot_distance_a']):.6f}",
                "ranking_score": row.get("ranking_score", ""),
                "iptm": row.get("iptm", ""),
                "ptm": row.get("ptm", ""),
                "fraction_disordered": row.get("fraction_disordered", ""),
                "has_clash": row.get("has_clash", ""),
                "local_af3_pdb": str(dst),
                "local_design_reference_pdb": str(ref_dst) if ref_dst.exists() else "",
                "remote_input_pdb": remote_input_pdb,
                "remote_run_id": f"md_{safe_name(design_id)}_af3_10ns",
                "af3_model_cif_path": row.get("af3_model_cif_path", ""),
                "raw_download_path": row.get("raw_download_path", ""),
                "md_selected_pre_af3_prodigy": "yes" if idx < 10 else "reserve",
            }
        )
        elapsed = time.time() - start
        rate = done / elapsed if elapsed > 0 else 0.0
        eta = (total - done) / rate if rate > 0 else 0.0
        print(
            f"[{now_text()}] done={done}/{total} ({100*done/total:.1f}%) "
            f"design={design_id} elapsed={elapsed:.1f}s ETA={eta:.1f}s",
            flush=True,
        )
        checkpoint_path.write_text(
            json.dumps(
                {
                    "updated_at": now_text(),
                    "stage": "preparing_inputs",
                    "done": done,
                    "total": total,
                    "current_design_id": design_id,
                    "result_root": str(result_root),
                },
                indent=2,
            ),
            encoding="utf-8",
        )

    all_manifest = result_root / "strict_9of9_af3_confirmed_manifest.csv"
    top10_prelim = result_root / "top10_preliminary_by_optimized_prodigy.csv"
    write_csv(all_manifest, rows)
    write_csv(top10_prelim, rows[:10])

    readme = result_root / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# Top10 AF3/PyMOL Confirmed IgG MD Input Package",
                "",
                f"Created: {now_text()}",
                "",
                "This package prepares the 38 strict high-confidence candidates for AF3-structure PRODIGY and Top10 MD.",
                "",
                "Selection basis:",
                "- Source candidates: PyMOL confirmed IgG complexes.",
                "- Strict subset: 9/9 hotspot contacts within 6 A and interface centroid distance <= 8 A.",
                "- Preliminary Top10: sorted by corrected optimized PRODIGY Delta G; final Top10 should be re-ranked by AF3-structure PRODIGY when available.",
                "",
                "Inputs:",
                f"- AF3 confirmed PDBs: `{af3_pdb_dir}`",
                f"- Design reference PDBs: `{refs_dir}`",
                "",
                "Core tables:",
                f"- strict manifest: `{all_manifest}`",
                f"- preliminary Top10: `{top10_prelim}`",
                "",
                "Remote defaults:",
                "- Remote input directory: `/mnt/PPIFlow/ica2_runs/input/md_top10_af3_confirmed_10ns/`",
                "- Remote output root: `/mnt/PPIFlow/ica2_runs/output/md_top10_af3_confirmed_10ns_<timestamp>/`",
                "- Remote results root: `/mnt/PPIFlow/ica2_runs/results/md_top10_af3_confirmed_10ns_<timestamp>/`",
                "",
            ]
        ),
        encoding="utf-8",
    )

    package_path = result_root / f"md_top10_af3_confirmed_input_package_{args.timestamp}.zip"
    package_items = [all_manifest, top10_prelim, readme]
    package_items.extend(sorted(af3_pdb_dir.glob("*.pdb")))
    package_items.extend(sorted(refs_dir.glob("*.pdb")))
    with zipfile.ZipFile(package_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for item in package_items:
            zf.write(item, item.relative_to(result_root).as_posix())

    summary = {
        "updated_at": now_text(),
        "stage": "input_preparation_complete",
        "strict_candidate_count": len(rows),
        "preliminary_top10_count": 10,
        "result_root": str(result_root),
        "strict_manifest": str(all_manifest),
        "top10_preliminary": str(top10_prelim),
        "input_package": str(package_path),
        "next_step": "Run AF3-structure PRODIGY for the 38 strict candidates, then select Top10 for 10 ns MD.",
    }
    (result_root / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    checkpoint_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[{now_text()}] complete strict_count={len(rows)} package={package_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
