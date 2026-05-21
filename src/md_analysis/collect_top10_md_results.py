#!/usr/bin/env python
"""Collect Top10 AF3-confirmed MD status, analysis statistics, and package results."""

from __future__ import annotations

import argparse
import csv
import json
import time
import zipfile
from datetime import datetime
from pathlib import Path


def now_text() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    if not fields:
        fields = ["status"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def candidate_result_dir(results_root: Path, design_id: str) -> Path:
    return results_root / "candidates" / f"md_{design_id}_af3_10ns"


def add_if_exists(zf: zipfile.ZipFile, path: Path, root: Path) -> int:
    if not path.exists():
        return 0
    count = 0
    if path.is_file():
        zf.write(path, path.relative_to(root).as_posix())
        return 1
    for item in path.rglob("*"):
        if item.is_file():
            zf.write(item, item.relative_to(root).as_posix())
            count += 1
            if count % 200 == 0:
                print(f"[package] added={count} current={item}", flush=True)
    return count


def archive_name(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.name


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", required=True)
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--make-zip", action="store_true")
    parser.add_argument("--no-zip", action="store_true")
    args = parser.parse_args()

    selected = Path(args.selected).resolve()
    results_root = Path(args.results_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    selected_rows = read_csv(selected)
    rows: list[dict[str, str]] = []
    for row in selected_rows:
        design_id = row.get("design_id", "")
        run_id = f"md_{design_id}_af3_10ns"
        cand_dir = candidate_result_dir(results_root, design_id)
        checkpoint = read_json(cand_dir / "checkpoint.json")
        stats = read_json(cand_dir / "analysis" / "summary_stats.json")
        out = dict(row)
        out.update(
            {
                "run_id": run_id,
                "candidate_results_dir": str(cand_dir),
                "md_status": checkpoint.get("status", "not_started" if not cand_dir.exists() else "unknown"),
                "md_stage": checkpoint.get("stage", ""),
                "rmsd_complex_last_nm": stats.get("rmsd_complex", {}).get("last"),
                "rmsd_complex_mean_nm": stats.get("rmsd_complex", {}).get("mean"),
                "hbonds_mean": stats.get("hbonds", {}).get("mean"),
                "interface_contacts_mean": stats.get("interface_contacts", {}).get("mean"),
                "interface_mindist_last_nm": stats.get("interface_mindist", {}).get("last"),
                "rg_complex_mean_nm": stats.get("rg_complex", {}).get("mean"),
                "sasa_complex_mean_nm2": stats.get("sasa_complex", {}).get("mean"),
                "buried_sasa_mean_nm2": stats.get("buried_sasa", {}).get("mean"),
                "nis_proxy_mean_nm2": stats.get("nis_proxy", {}).get("mean"),
            }
        )
        rows.append(out)

    summary_csv = output_dir / "md_top10_summary.csv"
    stability_csv = output_dir / "top10_md_stability_rank.csv"
    write_csv(summary_csv, rows)

    completed = [r for r in rows if r.get("md_status") == "completed"]
    completed.sort(
        key=lambda r: (
            float(r["rmsd_complex_last_nm"]) if r.get("rmsd_complex_last_nm") not in (None, "", "None") else 999.0,
            -(float(r["interface_contacts_mean"]) if r.get("interface_contacts_mean") not in (None, "", "None") else -1.0),
        )
    )
    for idx, row in enumerate(completed, start=1):
        row["md_stability_rank"] = str(idx)
    write_csv(stability_csv, completed)

    teacher_summary = output_dir / "md_top10_teacher_summary.md"
    lines = [
        "# Top10 AF3/PyMOL confirmed IgG 10 ns MD summary",
        "",
        f"- Updated: {now_text()}",
        "- Input structures: AF3/PyMOL-confirmed complexes.",
        "- Selection: strict 9/9 hotspot contact candidates ranked by AF3-structure PRODIGY before MD.",
        "- Protocol: amber99sb-ildn, tip3p, dodecahedron 1.0 nm, 0.15 M NaCl, EM, NVT 100 ps, NPT 100 ps, production 10 ns.",
        "- Analyses: ICs, NIS proxy, RMSD, H-bonds, Rg, SASA, RMSF.",
        "",
        f"- Selected candidates: {len(rows)}",
        f"- Completed MD: {len(completed)}",
        f"- Not completed/failed/running: {len(rows) - len(completed)}",
        "",
        "## Completed candidates",
    ]
    for row in completed:
        lines.append(
            f"- {row.get('design_id')}: RMSD last={row.get('rmsd_complex_last_nm')} nm, "
            f"H-bonds mean={row.get('hbonds_mean')}, contacts mean={row.get('interface_contacts_mean')}"
        )
    teacher_summary.write_text("\n".join(lines) + "\n", encoding="utf-8")

    summary_json = output_dir / "summary.json"
    summary = {
        "updated_at": now_text(),
        "selected_count": len(rows),
        "completed_count": len(completed),
        "summary_csv": str(summary_csv),
        "stability_csv": str(stability_csv),
        "teacher_summary": str(teacher_summary),
    }
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)

    if args.make_zip and not args.no_zip:
        package = output_dir / "md_top10_result_package.zip"
        start = time.time()
        print(f"[package] start {package}", flush=True)
        with zipfile.ZipFile(package, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=3) as zf:
            for file_path in [selected, summary_csv, stability_csv, teacher_summary, summary_json]:
                if file_path.exists():
                    zf.write(file_path, archive_name(file_path, output_dir))
            add_if_exists(zf, output_dir / "af3_structure_prodigy", output_dir)
            add_if_exists(zf, output_dir / "candidates", output_dir)
        elapsed = time.time() - start
        print(f"[package] complete {package} elapsed={elapsed:.1f}s size={package.stat().st_size}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
