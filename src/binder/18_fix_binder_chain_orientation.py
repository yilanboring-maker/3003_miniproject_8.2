#!/usr/bin/env python3
"""Normalize generated binder complexes to antigen chain A and binder chain T."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ATOM_PREFIXES = ("ATOM", "HETATM", "ANISOU", "TER")


def chain_lengths(lines: list[str]) -> dict[str, int]:
    residues: dict[str, set[tuple[str, str]]] = {}
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")) or len(line) <= 26:
            continue
        chain = line[21]
        res_id = (line[22:26].strip(), line[26].strip())
        residues.setdefault(chain, set()).add(res_id)
    return {chain: len(items) for chain, items in residues.items()}


def rewrite_chains(lines: list[str], mapping: dict[str, str]) -> list[str]:
    rewritten: list[str] = []
    for line in lines:
        if line.startswith(ATOM_PREFIXES) and len(line) > 21 and line[21] in mapping:
            line = line[:21] + mapping[line[21]] + line[22:]
        rewritten.append(line)
    return rewritten


def desired_mapping(lengths: dict[str, int], target_chain: str, binder_chain: str) -> dict[str, str] | None:
    if target_chain in lengths and binder_chain in lengths:
        if lengths[target_chain] > lengths[binder_chain]:
            return None
        return {target_chain: binder_chain, binder_chain: target_chain}

    if len(lengths) != 2:
        return None

    ordered = sorted(lengths.items(), key=lambda item: item[1])
    binder_raw, binder_len = ordered[0]
    target_raw, target_len = ordered[1]
    if binder_len >= target_len:
        return None
    return {binder_raw: binder_chain, target_raw: target_chain}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-root", required=True)
    parser.add_argument("--target-chain", default="A")
    parser.add_argument("--binder-chain", default="T")
    parser.add_argument("--summary-json", required=True)
    args = parser.parse_args()

    out_root = Path(args.out_root)
    pdbs = sorted(out_root.glob("smoke_test/*.pdb")) + sorted(out_root.glob("binder/**/*.pdb"))
    fixed = 0
    already_ok = 0
    skipped: list[dict[str, object]] = []

    for index, pdb in enumerate(pdbs, start=1):
        lines = pdb.read_text().splitlines(True)
        lengths = chain_lengths(lines)
        mapping = desired_mapping(lengths, args.target_chain, args.binder_chain)
        if mapping is None:
            if args.target_chain in lengths and args.binder_chain in lengths and lengths[args.target_chain] > lengths[args.binder_chain]:
                already_ok += 1
            else:
                skipped.append({"pdb": str(pdb), "lengths": lengths})
            continue
        pdb.write_text("".join(rewrite_chains(lines, mapping)))
        fixed += 1
        if index % 25 == 0 or index == len(pdbs):
            print(f"[fix-chains] {index}/{len(pdbs)} fixed={fixed} already_ok={already_ok} current={pdb.name}", flush=True)

    payload = {
        "out_root": str(out_root),
        "pdb_count": len(pdbs),
        "fixed": fixed,
        "already_ok": already_ok,
        "skipped_count": len(skipped),
        "skipped": skipped[:20],
        "target_chain": args.target_chain,
        "binder_chain": args.binder_chain,
    }
    Path(args.summary_json).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2), flush=True)
    return 0 if not skipped else 2


if __name__ == "__main__":
    raise SystemExit(main())
