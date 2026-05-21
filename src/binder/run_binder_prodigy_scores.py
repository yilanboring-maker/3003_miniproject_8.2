#!/usr/bin/env python3
"""Run PRODIGY for binder-antigen complexes with checkpointed progress."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import time
from pathlib import Path


RT_KCAL_25C = 0.0019872041 * 298.15


def parse_float(value: object, default: float = math.nan) -> float:
    try:
        return float(str(value))
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
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)


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


def kd_label(kd_m: float) -> str:
    if kd_m < 1e-9:
        return f"{kd_m * 1e12:.2f} pM"
    if kd_m < 1e-6:
        return f"{kd_m * 1e9:.2f} nM"
    if kd_m < 1e-3:
        return f"{kd_m * 1e6:.2f} uM"
    return f"{kd_m:.2e} M"


def add_ranks(rows: list[dict[str, str]]) -> None:
    scored = [row for row in rows if row.get("prodigy_status") == "ok"]
    scored.sort(key=lambda item: parse_float(item.get("prodigy_delta_g", ""), 9999.0))
    for rank, row in enumerate(scored, start=1):
        row["prodigy_rank"] = str(rank)


def command_for(row: dict[str, str], binder_chain: str, target_chain: str) -> list[str]:
    return ["prodigy", "-q", row["pdb_path"], "--selection", binder_chain, target_chain]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--binder-chain", default="T")
    parser.add_argument("--target-chain", default="A")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    input_csv = Path(args.input_csv)
    output_csv = Path(args.output_csv)
    checkpoint = Path(args.checkpoint)
    input_rows = load_rows(input_csv)

    existing_by_id: dict[str, dict[str, str]] = {}
    if output_csv.exists() and not args.force:
        for row in load_rows(output_csv):
            if row.get("prodigy_status") == "ok":
                existing_by_id[row.get("design_id", "")] = row

    fields = list(input_rows[0].keys()) if input_rows else []
    for field in [
        "prodigy_delta_g",
        "prodigy_kd",
        "prodigy_kd_m",
        "prodigy_kd_label",
        "prodigy_rank",
        "prodigy_status",
        "prodigy_stdout",
        "prodigy_stderr",
    ]:
        if field not in fields:
            fields.append(field)

    total = len(input_rows)
    scored_rows: list[dict[str, str]] = []
    start = time.time()
    best_delta = math.inf

    write_checkpoint(
        checkpoint,
        {
            "script": "16_run_binder_prodigy_scores.py",
            "status": "running",
            "input_csv": str(input_csv),
            "output_csv": str(output_csv),
            "done": 0,
            "total": total,
            "updated_at_epoch": int(time.time()),
        },
    )

    for index, row in enumerate(input_rows, start=1):
        elapsed = time.time() - start
        rate = (index - 1) / elapsed if elapsed > 0 and index > 1 else 0.0
        eta = (total - index + 1) / rate if rate > 0 else -1.0
        design_id = row.get("design_id", "")

        if design_id in existing_by_id:
            out = dict(existing_by_id[design_id])
            print(
                f"[prodigy-binder] {index}/{total} {100.0 * index / max(total, 1):.1f}% "
                f"elapsed={elapsed:.1f}s eta={eta:.1f}s sample={design_id} status=reused",
                flush=True,
            )
        else:
            print(
                f"[prodigy-binder] {index}/{total} {100.0 * index / max(total, 1):.1f}% "
                f"elapsed={elapsed:.1f}s eta={eta:.1f}s best_delta_g={best_delta if best_delta < math.inf else 'NA'} "
                f"sample={design_id}",
                flush=True,
            )
            out = dict(row)
            try:
                completed = subprocess.run(
                    command_for(row, args.binder_chain, args.target_chain),
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
                    out["prodigy_kd"] = ""
                    out["prodigy_kd_m"] = ""
                    out["prodigy_kd_label"] = ""
                else:
                    delta_g = parse_quiet_output(completed.stdout)
                    kd_m = kd_from_delta_g(delta_g)
                    out["prodigy_delta_g"] = f"{delta_g:.3f}"
                    out["prodigy_kd"] = f"{kd_m:.6e}"
                    out["prodigy_kd_m"] = f"{kd_m:.6e}"
                    out["prodigy_kd_label"] = kd_label(kd_m)
                    out["prodigy_status"] = "ok"
                    best_delta = min(best_delta, delta_g)
            except Exception as exc:
                out["prodigy_status"] = "exception"
                out["prodigy_delta_g"] = ""
                out["prodigy_kd"] = ""
                out["prodigy_kd_m"] = ""
                out["prodigy_kd_label"] = ""
                out["prodigy_stdout"] = ""
                out["prodigy_stderr"] = repr(exc)

        scored_rows.append(out)
        add_ranks(scored_rows)
        scored_rows.sort(
            key=lambda item: (
                parse_float(item.get("prodigy_rank", ""), 9999.0),
                parse_float(item.get("rank", ""), 9999.0),
            )
        )
        write_rows(output_csv, scored_rows, fields)
        ok_count = sum(1 for item in scored_rows if item.get("prodigy_status") == "ok")
        fail_count = len(scored_rows) - ok_count
        write_checkpoint(
            checkpoint,
            {
                "script": "16_run_binder_prodigy_scores.py",
                "status": "running",
                "input_csv": str(input_csv),
                "output_csv": str(output_csv),
                "done": index,
                "total": total,
                "ok_count": ok_count,
                "fail_count": fail_count,
                "best_delta_g": None if best_delta == math.inf else best_delta,
                "updated_at_epoch": int(time.time()),
            },
        )

    add_ranks(scored_rows)
    scored_rows.sort(
        key=lambda item: (
            parse_float(item.get("prodigy_rank", ""), 9999.0),
            parse_float(item.get("rank", ""), 9999.0),
        )
    )
    write_rows(output_csv, scored_rows, fields)
    ok_count = sum(1 for item in scored_rows if item.get("prodigy_status") == "ok")
    fail_count = len(scored_rows) - ok_count
    write_checkpoint(
        checkpoint,
        {
            "script": "16_run_binder_prodigy_scores.py",
            "status": "completed" if fail_count == 0 else "completed_with_failures",
            "input_csv": str(input_csv),
            "output_csv": str(output_csv),
            "done": total,
            "total": total,
            "ok_count": ok_count,
            "fail_count": fail_count,
            "elapsed_seconds": round(time.time() - start, 2),
            "updated_at_epoch": int(time.time()),
        },
    )
    print(f"[done] ok={ok_count}/{total} fail={fail_count} output_csv={output_csv}", flush=True)
    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
