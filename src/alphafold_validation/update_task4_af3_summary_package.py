#!/usr/bin/env python
"""Update Task 4 AlphaFold Server status files and build a current package."""

from __future__ import annotations

import argparse
import json
import time
import zipfile
from datetime import datetime
from pathlib import Path

import pandas as pd


def now_text() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def safe_int(value) -> int:
    try:
        return int(value)
    except Exception:
        return 0


def add_file(zf: zipfile.ZipFile, path: Path, root: Path, written: set[str]) -> None:
    if not path.exists() or not path.is_file():
        return
    arcname = str(path.relative_to(root.parent))
    if arcname in written:
        return
    zf.write(path, arcname)
    written.add(arcname)


def add_tree(
    zf: zipfile.ZipFile,
    tree: Path,
    root: Path,
    written: set[str],
    *,
    include_extracted_af3_dirs: bool,
) -> int:
    if not tree.exists():
        return 0
    count = 0
    for path in sorted(tree.rglob("*")):
        if not path.is_file():
            continue
        if tree.name == "af3_raw_downloads" and not include_extracted_af3_dirs:
            if path.suffix.lower() != ".zip":
                continue
        add_file(zf, path, root, written)
        count += 1
    return count


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Task 4 result root")
    ap.add_argument(
        "--include-extracted-af3-dirs",
        action="store_true",
        help="Also include extracted AlphaFold Server result directories, not only original downloaded zips.",
    )
    args = ap.parse_args()

    root = Path(args.root).resolve()
    manifest_path = root / "af3_server_job_manifest.csv"
    all_results_path = root / "af3_pose_review_all_results.csv"
    summary_path = root / "review_tables" / "task4_status_summary.json"
    checkpoint_path = root / "checkpoint.json"
    teacher_summary_path = root / "task4_teacher_summary.md"
    package_path = root / f"task4_af3_current_result_package_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"

    print(f"[{now_text()}] update start root={root}", flush=True)
    print(
        f"[{now_text()}] expected outputs: summary={summary_path}; checkpoint={checkpoint_path}; "
        f"teacher_summary={teacher_summary_path}; package={package_path}",
        flush=True,
    )

    manifest = pd.read_csv(manifest_path)
    all_results = pd.read_csv(all_results_path)

    reviewed_mask = all_results["review_status"].eq("reviewed")
    completed_designs = set(all_results.loc[reviewed_mask, "design_id"].astype(str))
    downloaded_by_design = dict(
        zip(
            all_results.loc[reviewed_mask, "design_id"].astype(str),
            all_results.loc[reviewed_mask, "raw_download_path"].fillna("").astype(str),
        )
    )

    for df in (manifest, all_results):
        design_str = df["design_id"].astype(str)
        completed_mask = design_str.isin(completed_designs)
        if "alphafold_server_status" in df.columns:
            df.loc[completed_mask, "alphafold_server_status"] = "completed"
        if "pose_review_status" in df.columns:
            df.loc[completed_mask, "pose_review_status"] = "reviewed"
        if "download_path" in df.columns:
            for idx, design_id in df.loc[completed_mask, "design_id"].items():
                raw_path = downloaded_by_design.get(str(design_id), "")
                if raw_path:
                    df.at[idx, "download_path"] = raw_path

    manifest.to_csv(manifest_path, index=False)
    all_results.to_csv(all_results_path, index=False)

    total = len(all_results)
    reviewed_count = int(reviewed_mask.sum())
    submitted_count = int(all_results["alphafold_server_status"].ne("not_submitted").sum())
    in_progress_rows = all_results[
        all_results["alphafold_server_status"].eq("in_progress")
        & all_results["review_status"].ne("reviewed")
    ].copy()
    not_submitted_count = int(all_results["alphafold_server_status"].eq("not_submitted").sum())
    pose_counts = {
        str(k): safe_int(v)
        for k, v in all_results["auto_pose_call"].value_counts(dropna=False).to_dict().items()
    }
    status_counts = {
        str(k): safe_int(v)
        for k, v in all_results["alphafold_server_status"].value_counts(dropna=False).to_dict().items()
    }
    if reviewed_count == total:
        next_step = (
            "All downloaded AlphaFold Server jobs have been parsed. Use the "
            "acceptable/needs_review/reject tables and PyMOL-aligned models for the "
            "teacher-facing Task 4 decision; continue with PRODIGY/MD only for "
            "teacher-approved candidates."
        )
    else:
        next_step = "Download/review remaining in-progress jobs once complete."

    summary = {
        "updated_at": now_text(),
        "task": "Task 4 AlphaFold Server validation for 150 optimized IgG candidates",
        "output_root": str(root),
        "total_candidates": total,
        "submitted_to_server": submitted_count,
        "downloaded_reviewed": reviewed_count,
        "status_counts": status_counts,
        "pose_review_counts": pose_counts,
        "still_in_progress_designs": sorted(in_progress_rows["design_id"].astype(str).tolist()),
        "still_in_progress_jobs": sorted(in_progress_rows["alphafold_server_job_name"].astype(str).tolist()),
        "not_submitted_count": not_submitted_count,
        "all_results": str(all_results_path),
        "manifest": str(manifest_path),
        "current_package": str(package_path),
        "next_step": next_step,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    checkpoint_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    acceptable = pose_counts.get("acceptable", 0)
    needs_review = pose_counts.get("needs_review", 0)
    reject = pose_counts.get("reject", 0)
    not_reviewed = pose_counts.get("not_reviewed", 0)
    pending_names = ", ".join(summary["still_in_progress_jobs"]) or "none"

    teacher_summary = f"""# Task 4 AlphaFold Server validation current summary

Updated: {summary["updated_at"]}

Scope: 150 optimized IgG candidates with corrected PRODIGY Delta G <= -9 kcal/mol.

Current progress:
- Submitted to AlphaFold Server: {submitted_count}/150
- Downloaded and parsed/reviewed locally: {reviewed_count}/150
- Still running on server: {len(summary["still_in_progress_jobs"])} ({pending_names})
- Not submitted yet due to quota / next account batch: {not_submitted_count}

Automated pose triage among parsed results:
- acceptable: {acceptable}
- needs_review: {needs_review}
- reject: {reject}
- not_reviewed: {not_reviewed}

Important note: the labels above are automated antigen-aligned pose metrics for triage. Final Task 4 wording should still say that PyMOL visual checking is required/used for the teacher-facing decision.

Key local files:
- Full result table: `{all_results_path}`
- Accepted triage table: `{root / "af3_pose_review_acceptable.csv"}`
- Needs-review triage table: `{root / "af3_pose_review_needs_review.csv"}`
- Reject triage table: `{root / "af3_pose_review_reject.csv"}`
- Current result package: `{package_path}`
"""
    teacher_summary_path.write_text(teacher_summary, encoding="utf-8")

    package_files = [
        manifest_path,
        all_results_path,
        root / "af3_pose_review_acceptable.csv",
        root / "af3_pose_review_needs_review.csv",
        root / "af3_pose_review_reject.csv",
        root / "partial_results.csv",
        checkpoint_path,
        summary_path,
        root / "review_tables" / "af3_pose_review_summary.json",
        root / "review_tables" / "af3_pose_review_model_metrics.csv",
        teacher_summary_path,
    ]
    package_dirs = [
        root / "review_tables",
        root / "server_inputs",
        root / "design_pdbs",
        root / "af3_raw_downloads",
    ]

    written: set[str] = set()
    start = time.time()
    with zipfile.ZipFile(package_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        total_units = len(package_files) + len(package_dirs)
        done = 0
        for path in package_files:
            add_file(zf, path, root, written)
            done += 1
            elapsed = time.time() - start
            print(
                f"[{now_text()}] package_progress done={done}/{total_units} "
                f"elapsed={elapsed:.1f}s current={path.name}",
                flush=True,
            )
        for tree in package_dirs:
            n = add_tree(
                zf,
                tree,
                root,
                written,
                include_extracted_af3_dirs=args.include_extracted_af3_dirs,
            )
            done += 1
            elapsed = time.time() - start
            print(
                f"[{now_text()}] package_progress done={done}/{total_units} "
                f"elapsed={elapsed:.1f}s tree={tree.name} files_added={n}",
                flush=True,
            )

    summary["current_package_size_bytes"] = package_path.stat().st_size
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    checkpoint_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[{now_text()}] update complete package={package_path} size={package_path.stat().st_size}", flush=True)
    print(f"[{now_text()}] summary={summary}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
