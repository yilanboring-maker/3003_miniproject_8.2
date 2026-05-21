#!/usr/bin/env python3
"""Run PRODIGY on ranked PPIFlow candidates with checkpointed progress."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path


RT_KCAL_25C = 0.0019872041 * 298.15


def parse_float(value: str, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_checkpoint(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp_path.replace(path)


def select_rows(rows: list[dict[str, str]], max_per_modality: int | None) -> list[dict[str, str]]:
    if max_per_modality is None:
        return rows
    selected: list[dict[str, str]] = []
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        modality = row.get("modality", "")
        if counts[modality] < max_per_modality:
            selected.append(row)
            counts[modality] += 1
    return selected


def prodigy_command(row: dict[str, str]) -> list[str]:
    pdb_path = row["pdb_path"]
    modality = row.get("modality", "")
    if modality == "antibody":
        return ["prodigy", "-q", pdb_path, "--selection", "A,B", "C"]
    if modality == "nanobody":
        return ["prodigy", "-q", pdb_path, "--selection", "G", "A"]
    raise ValueError(f"Unknown modality: {modality!r}")


def parse_quiet_output(stdout: str) -> float:
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


def add_prodigy_ranks(rows: list[dict[str, str]]) -> None:
    by_modality: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_modality[row.get("modality", "")].append(row)

    for modality_rows in by_modality.values():
        scored = [row for row in modality_rows if row.get("prodigy_status") == "ok"]
        scored.sort(key=lambda item: parse_float(item.get("prodigy_delta_g", ""), 9999.0))
        for rank, row in enumerate(scored, start=1):
            row["prodigy_rank"] = str(rank)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--max-per-modality", type=int, default=20)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()

    input_csv = Path(args.input_csv)
    output_csv = Path(args.output_csv)
    checkpoint = Path(args.checkpoint)
    rows = select_rows(load_rows(input_csv), args.max_per_modality)

    fields = list(rows[0].keys()) if rows else []
    for field in ["prodigy_delta_g", "prodigy_kd", "prodigy_kd_m", "prodigy_rank", "prodigy_status", "prodigy_stdout", "prodigy_stderr"]:
        if field not in fields:
            fields.append(field)

    scored_rows: list[dict[str, str]] = []
    start = time.time()
    total = len(rows)
    write_checkpoint(
        checkpoint,
        {
            "script": "08_run_prodigy_scores.py",
            "status": "running",
            "input_csv": str(input_csv),
            "output_csv": str(output_csv),
            "done": 0,
            "total": total,
            "updated_at_epoch": int(time.time()),
        },
    )

    for index, row in enumerate(rows, start=1):
        elapsed = time.time() - start
        rate = (index - 1) / elapsed if elapsed > 0 and index > 1 else 0.0
        eta = (total - index + 1) / rate if rate > 0 else -1.0
        print(
            f"[prodigy] {index}/{total} {100.0 * index / max(total, 1):.1f}% "
            f"elapsed={elapsed:.1f}s eta={eta:.1f}s sample={row.get('sample')}",
            flush=True,
        )

        out = dict(row)
        try:
            completed = subprocess.run(
                prodigy_command(row),
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=args.timeout,
            )
            out["prodigy_stdout"] = completed.stdout.strip()
            out["prodigy_stderr"] = completed.stderr.strip()
            if completed.returncode != 0:
                out["prodigy_status"] = f"failed_exit_{completed.returncode}"
                out["prodigy_delta_g"] = ""
                out["prodigy_kd_m"] = ""
            else:
                delta_g = parse_quiet_output(completed.stdout)
                kd_m = kd_from_delta_g(delta_g)
                out["prodigy_delta_g"] = f"{delta_g:.3f}"
                out["prodigy_kd"] = f"{kd_m:.6e}"
                out["prodigy_kd_m"] = f"{kd_m:.6e}"
                out["prodigy_status"] = "ok"
        except Exception as exc:
            out["prodigy_status"] = "exception"
            out["prodigy_delta_g"] = ""
            out["prodigy_kd_m"] = ""
            out["prodigy_stdout"] = ""
            out["prodigy_stderr"] = repr(exc)

        scored_rows.append(out)
        add_prodigy_ranks(scored_rows)
        write_rows(output_csv, scored_rows, fields)
        write_checkpoint(
            checkpoint,
            {
                "script": "08_run_prodigy_scores.py",
                "status": "running",
                "input_csv": str(input_csv),
                "output_csv": str(output_csv),
                "done": index,
                "total": total,
                "updated_at_epoch": int(time.time()),
            },
        )

    add_prodigy_ranks(scored_rows)
    scored_rows.sort(
        key=lambda item: (
            item.get("modality", ""),
            parse_float(item.get("prodigy_rank", ""), 9999.0),
            parse_float(item.get("rank", ""), 9999.0),
        )
    )
    write_rows(output_csv, scored_rows, fields)
    write_checkpoint(
        checkpoint,
        {
            "script": "08_run_prodigy_scores.py",
            "status": "completed",
            "input_csv": str(input_csv),
            "output_csv": str(output_csv),
            "done": total,
            "total": total,
            "elapsed_seconds": round(time.time() - start, 2),
            "updated_at_epoch": int(time.time()),
        },
    )
    print(f"[done] output_csv={output_csv}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
