#!/usr/bin/env python
"""Run PRODIGY on AF3/PyMOL-confirmed strict IgG complexes.

Designed for the remote Matpool environment where PRODIGY is installed. The
script is also safe to run locally: it will checkpoint failures instead of
fabricating scores.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import time
from datetime import datetime
from pathlib import Path


R_KCAL = 0.00198720425864083
TEMP_K = 298.15


def now_text() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def kd_from_delta_g(delta_g: float) -> float:
    return math.exp(delta_g / (R_KCAL * TEMP_K))


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_delta_g(stdout: str) -> float:
    for line in stdout.splitlines():
        if "Predicted binding affinity" in line or "kcal" in line:
            for token in line.replace("=", " ").replace(":", " ").split():
                try:
                    value = float(token)
                except ValueError:
                    continue
                if -100.0 < value < 100.0:
                    return value
    for token in stdout.replace("=", " ").replace(":", " ").replace(",", " ").split():
        try:
            value = float(token)
        except ValueError:
            continue
        if -100.0 < value < 0.0:
            return value
    raise ValueError(f"No numeric PRODIGY delta G found in output: {stdout[:500]!r}")


def portable_basename(path_text: str) -> str:
    """Return basename for either POSIX or Windows-style paths."""
    return re.split(r"[\\/]+", path_text.strip())[-1]


def resolve_pdb(row: dict[str, str], input_base: Path | None) -> Path:
    if input_base is not None:
        candidate = input_base / portable_basename(row["local_af3_pdb"])
        if candidate.exists():
            return candidate
    local = Path(row["local_af3_pdb"])
    if local.exists():
        return local
    remote = Path(row.get("remote_input_pdb", ""))
    return remote


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--input-base", default="")
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--selection-left", default="A,B")
    parser.add_argument("--selection-right", default="C")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    manifest = Path(args.manifest).resolve()
    results_dir = Path(args.results_dir).resolve()
    input_base = Path(args.input_base).resolve() if args.input_base else None
    results_dir.mkdir(parents=True, exist_ok=True)
    all_csv = results_dir / "af3_confirmed_prodigy_all_results.csv"
    top_csv = results_dir / "af3_confirmed_prodigy_top10.csv"
    checkpoint = results_dir / "checkpoint.json"
    partial = results_dir / "partial_results.csv"
    summary_json = results_dir / "summary.json"

    rows = load_rows(manifest)
    total = len(rows)
    prodigy_path = shutil.which("prodigy")
    print(f"[{now_text()}] AF3 confirmed PRODIGY start")
    print(f"manifest={manifest}")
    print(f"results_dir={results_dir}")
    print(f"input_base={input_base or ''}")
    print(f"total={total}")
    print(f"prodigy={prodigy_path or 'NOT_FOUND'}", flush=True)

    completed: list[dict[str, str]] = []
    if all_csv.exists() and not args.force:
        completed = load_rows(all_csv)
    done_ids = {row.get("design_id") for row in completed if row.get("status") == "ok"}

    start = time.time()
    for idx, row in enumerate(rows, start=1):
        design_id = row["design_id"]
        out = dict(row)
        if design_id in done_ids and not args.force:
            print(f"[{now_text()}] skip existing {idx}/{total} design={design_id}", flush=True)
            continue
        pdb = resolve_pdb(row, input_base)
        cmd = ["prodigy", "-q", str(pdb), "--selection", args.selection_left, args.selection_right]
        out["af3_prodigy_command"] = " ".join(cmd)
        out["af3_prodigy_pdb_path"] = str(pdb)
        elapsed = time.time() - start
        rate = (idx - 1) / elapsed if elapsed > 0 and idx > 1 else 0.0
        eta = (total - idx + 1) / rate if rate > 0 else 0.0
        print(
            f"[{now_text()}] running {idx}/{total} ({100*idx/total:.1f}%) "
            f"design={design_id} elapsed={elapsed:.1f}s ETA={eta:.1f}s",
            flush=True,
        )
        if not prodigy_path:
            out["status"] = "failed_no_prodigy"
            out["af3_prodigy_delta_g"] = ""
            out["af3_prodigy_kd_m"] = ""
            out["af3_prodigy_stdout"] = ""
            out["af3_prodigy_stderr"] = "prodigy executable not found on PATH"
        elif not pdb.exists():
            out["status"] = "failed_missing_pdb"
            out["af3_prodigy_delta_g"] = ""
            out["af3_prodigy_kd_m"] = ""
            out["af3_prodigy_stdout"] = ""
            out["af3_prodigy_stderr"] = f"PDB not found: {pdb}"
        else:
            proc = subprocess.run(cmd, text=True, capture_output=True)
            out["af3_prodigy_stdout"] = proc.stdout.strip()
            out["af3_prodigy_stderr"] = proc.stderr.strip()
            if proc.returncode != 0:
                out["status"] = f"failed_exit_{proc.returncode}"
                out["af3_prodigy_delta_g"] = ""
                out["af3_prodigy_kd_m"] = ""
            else:
                try:
                    delta_g = parse_delta_g(proc.stdout)
                    kd = kd_from_delta_g(delta_g)
                    out["status"] = "ok"
                    out["af3_prodigy_delta_g"] = f"{delta_g:.3f}"
                    out["af3_prodigy_kd_m"] = f"{kd:.6e}"
                except Exception as exc:
                    out["status"] = "failed_parse"
                    out["af3_prodigy_delta_g"] = ""
                    out["af3_prodigy_kd_m"] = ""
                    out["af3_prodigy_stderr"] = (out["af3_prodigy_stderr"] + "\n" + repr(exc)).strip()

        completed = [r for r in completed if r.get("design_id") != design_id] + [out]
        scored = [r for r in completed if r.get("status") == "ok" and r.get("af3_prodigy_delta_g")]
        scored.sort(key=lambda r: float(r["af3_prodigy_delta_g"]))
        for rank, scored_row in enumerate(scored, start=1):
            scored_row["af3_prodigy_rank"] = str(rank)
            scored_row["md_selected_by_af3_prodigy"] = "yes" if rank <= args.top_n else "reserve"
        write_rows(partial, completed)
        write_rows(all_csv, completed)
        write_rows(top_csv, scored[: args.top_n])
        current_best = scored[0]["af3_prodigy_delta_g"] if scored else ""
        checkpoint.write_text(
            json.dumps(
                {
                    "updated_at": now_text(),
                    "stage": "running",
                    "done": idx,
                    "total": total,
                    "current_design_id": design_id,
                    "current_status": out["status"],
                    "current_best_delta_g": current_best,
                    "all_results": str(all_csv),
                    "top10": str(top_csv),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        print(
            f"[{now_text()}] done={idx}/{total} status={out['status']} "
            f"delta_g={out.get('af3_prodigy_delta_g','')} best={current_best}",
            flush=True,
        )

    completed = load_rows(all_csv) if all_csv.exists() else completed
    scored = [r for r in completed if r.get("status") == "ok" and r.get("af3_prodigy_delta_g")]
    scored.sort(key=lambda r: float(r["af3_prodigy_delta_g"]))
    for rank, row in enumerate(scored, start=1):
        row["af3_prodigy_rank"] = str(rank)
        row["md_selected_by_af3_prodigy"] = "yes" if rank <= args.top_n else "reserve"
    write_rows(all_csv, completed)
    write_rows(top_csv, scored[: args.top_n])
    summary = {
        "updated_at": now_text(),
        "stage": "complete",
        "total": total,
        "ok": len(scored),
        "failed": total - len(scored),
        "top_n": args.top_n,
        "all_results": str(all_csv),
        "top10": str(top_csv),
    }
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    checkpoint.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[{now_text()}] complete ok={len(scored)}/{total} top10={top_csv}")
    return 0 if len(scored) else 2


if __name__ == "__main__":
    raise SystemExit(main())
