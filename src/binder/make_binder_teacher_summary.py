#!/usr/bin/env python3
"""Create teacher-facing binder PRODIGY summaries and copy top PDBs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import time
from pathlib import Path


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def parse_float(value: object, default: float = math.nan) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return default


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def teacher_rows(rows: list[dict[str, str]], limit: int) -> list[dict[str, str]]:
    scored = [row for row in rows if row.get("prodigy_status") == "ok"]
    scored.sort(key=lambda item: parse_float(item.get("prodigy_delta_g", ""), 9999.0))
    out_rows: list[dict[str, str]] = []
    for rank, row in enumerate(scored[:limit], start=1):
        clash = parse_bool(row.get("has_clash", "False"))
        out_rows.append(
            {
                "rank": str(rank),
                "design_id": row.get("design_id", ""),
                "prodigy_delta_g_kcal_mol": row.get("prodigy_delta_g", ""),
                "predicted_kd_m": row.get("prodigy_kd_m", row.get("prodigy_kd", "")),
                "predicted_kd_label": row.get("prodigy_kd_label", ""),
                "hotspot_coverage_6a": row.get("hotspot_coverage_6a", ""),
                "hotspot_contacted_6a": row.get("hotspot_contacted_6a", ""),
                "hotspot_total": row.get("hotspot_total", ""),
                "coverage_hotspot_list_6a": row.get("coverage_hotspot_list_6a", ""),
                "clash_status": "clash" if clash else "no_clash",
                "min_interchain_distance": row.get("min_interchain_distance", ""),
                "binder_length": row.get("binder_length", ""),
                "pdb_path": row.get("pdb_path", ""),
            }
        )
    return out_rows


def copy_top_pdbs(rows: list[dict[str, str]], dest: Path, limit: int) -> list[str]:
    dest.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    for row in rows[:limit]:
        source = Path(row.get("pdb_path", ""))
        if not source.exists():
            continue
        target = dest / f"binder_rank{int(row['rank']):02d}_{row['design_id']}.pdb"
        shutil.copy2(source, target)
        copied.append(str(target))
    return copied


def make_top10_png(rows: list[dict[str, str]], png_path: Path) -> str:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        return f"matplotlib unavailable: {exc!r}"

    top10 = rows[:10]
    labels = [row["design_id"] for row in top10]
    values = [parse_float(row["prodigy_delta_g_kcal_mol"], 0.0) for row in top10]
    colors = ["#235789" if row["clash_status"] == "no_clash" else "#b23a48" for row in top10]

    fig, ax = plt.subplots(figsize=(11, 6.5))
    positions = list(range(len(top10)))
    ax.barh(positions, values, color=colors)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("PRODIGY predicted binding free energy, delta G (kcal/mol)")
    ax.set_title("Top 10 protein binders ranked by PRODIGY")
    ax.axvline(0, color="#444444", linewidth=0.8)
    for y, row in zip(positions, top10):
        dg = row["prodigy_delta_g_kcal_mol"]
        kd = row["predicted_kd_label"] or row["predicted_kd_m"]
        coverage = row["hotspot_coverage_6a"]
        ax.text(
            values[y] - 0.05,
            y,
            f"{dg} kcal/mol | {kd} | hotspot {coverage}",
            va="center",
            ha="right",
            color="white",
            fontsize=8,
        )
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=220)
    plt.close(fig)
    return "ok"


def write_markdown(path: Path, summary: dict[str, object], top_rows: list[dict[str, str]]) -> None:
    top = top_rows[0] if top_rows else {}
    lines = [
        "# Binder PRODIGY summary",
        "",
        f"- Generated binder records: {summary['generated_count']}",
        f"- PRODIGY success: {summary['prodigy_success_count']}/{summary['generated_count']}",
        f"- PRODIGY failures: {summary['prodigy_failure_count']}",
        f"- Top binder: {top.get('design_id', 'NA')}",
        f"- Top delta G: {top.get('prodigy_delta_g_kcal_mol', 'NA')} kcal/mol",
        f"- Top predicted Kd: {top.get('predicted_kd_label', top.get('predicted_kd_m', 'NA'))}",
        f"- Top hotspot coverage: {top.get('hotspot_contacted_6a', 'NA')}/{top.get('hotspot_total', 'NA')} at 6 A",
        f"- Top clash status: {top.get('clash_status', 'NA')}",
        "",
        "Files:",
        "",
        f"- Top 10 image: `{summary['top10_png']}`",
        f"- Top 20 table: `{summary['top20_csv']}`",
        f"- All PRODIGY results: `{summary['all_prodigy_csv']}`",
        f"- Top PDB directory: `{summary['top_pdb_dir']}`",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prodigy-csv", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--top-pdb-dir", required=True)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--copy-pdb-n", type=int, default=10)
    args = parser.parse_args()

    prodigy_csv = Path(args.prodigy_csv)
    results_dir = Path(args.results_dir)
    top_pdb_dir = Path(args.top_pdb_dir)
    rows = load_rows(prodigy_csv)
    top_rows = teacher_rows(rows, args.top_n)

    top20_csv = results_dir / "binder_teacher_summary_top20.csv"
    top10_png = results_dir / "binder_teacher_summary_top10.png"
    final_top = results_dir / "binder_final_top_candidates.csv"
    summary_json = results_dir / "binder_teacher_summary.json"
    summary_md = results_dir / "binder_teacher_summary.md"

    fields = [
        "rank",
        "design_id",
        "prodigy_delta_g_kcal_mol",
        "predicted_kd_m",
        "predicted_kd_label",
        "hotspot_coverage_6a",
        "hotspot_contacted_6a",
        "hotspot_total",
        "coverage_hotspot_list_6a",
        "clash_status",
        "min_interchain_distance",
        "binder_length",
        "pdb_path",
    ]
    write_csv(top20_csv, top_rows, fields)
    write_csv(final_top, top_rows, fields)
    copied = copy_top_pdbs(top_rows, top_pdb_dir, args.copy_pdb_n)
    png_status = make_top10_png(top_rows, top10_png)

    success_count = sum(1 for row in rows if row.get("prodigy_status") == "ok")
    failure_count = len(rows) - success_count
    summary = {
        "script": "17_make_binder_teacher_summary.py",
        "generated_count": len(rows),
        "prodigy_success_count": success_count,
        "prodigy_failure_count": failure_count,
        "top_design_id": top_rows[0]["design_id"] if top_rows else None,
        "top_delta_g_kcal_mol": top_rows[0]["prodigy_delta_g_kcal_mol"] if top_rows else None,
        "top_predicted_kd": top_rows[0]["predicted_kd_label"] if top_rows else None,
        "top_hotspot_coverage_6a": top_rows[0]["hotspot_coverage_6a"] if top_rows else None,
        "top_clash_status": top_rows[0]["clash_status"] if top_rows else None,
        "all_prodigy_csv": str(prodigy_csv),
        "top20_csv": str(top20_csv),
        "top10_png": str(top10_png),
        "final_top_candidates": str(final_top),
        "top_pdb_dir": str(top_pdb_dir),
        "copied_pdbs": copied,
        "png_status": png_status,
        "updated_at_epoch": int(time.time()),
    }
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_markdown(summary_md, summary, top_rows)

    print(f"[done] generated_count={len(rows)} prodigy_success={success_count} failures={failure_count}", flush=True)
    print(f"[done] top20_csv={top20_csv}", flush=True)
    print(f"[done] top10_png={top10_png} status={png_status}", flush=True)
    print(f"[done] top_pdb_dir={top_pdb_dir} copied={len(copied)}", flush=True)
    return 0 if success_count > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
