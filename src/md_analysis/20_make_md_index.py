#!/usr/bin/env python3
"""Create antibody/antigen GROMACS index groups from processed PDB chain IDs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


BB_ATOMS = {"N", "CA", "C"}


def wrap_numbers(numbers: list[int]) -> str:
    lines = []
    for start in range(0, len(numbers), 15):
        lines.append(" ".join(str(item) for item in numbers[start : start + 15]))
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    pdb = Path(args.pdb)
    groups = {
        "Antibody": [],
        "Antigen": [],
        "Complex": [],
        "Antibody_BB": [],
        "Antigen_BB": [],
        "Complex_BB": [],
        "Antibody_CA": [],
        "Antigen_CA": [],
        "Complex_CA": [],
    }
    residues: dict[str, set[tuple[str, str]]] = {}

    for line in pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 26:
            continue
        atom_id = int(line[6:11])
        atom_name = line[12:16].strip()
        chain = line[21]
        residues.setdefault(chain, set()).add((line[22:26].strip(), line[26].strip()))
        if chain in {"A", "B", "C"}:
            groups["Complex"].append(atom_id)
            if atom_name in BB_ATOMS:
                groups["Complex_BB"].append(atom_id)
            if atom_name == "CA":
                groups["Complex_CA"].append(atom_id)
        if chain in {"A", "B"}:
            groups["Antibody"].append(atom_id)
            if atom_name in BB_ATOMS:
                groups["Antibody_BB"].append(atom_id)
            if atom_name == "CA":
                groups["Antibody_CA"].append(atom_id)
        elif chain == "C":
            groups["Antigen"].append(atom_id)
            if atom_name in BB_ATOMS:
                groups["Antigen_BB"].append(atom_id)
            if atom_name == "CA":
                groups["Antigen_CA"].append(atom_id)

    missing = [name for name, atoms in groups.items() if not atoms]
    if missing:
        raise SystemExit(f"Empty index groups: {missing}")

    out = Path(args.output)
    with out.open("w", encoding="utf-8") as handle:
        for name, atoms in groups.items():
            handle.write(f"[ {name} ]\n")
            handle.write(wrap_numbers(atoms))

    summary = {
        "pdb": str(pdb),
        "output": str(out),
        "chain_residue_counts": {chain: len(items) for chain, items in sorted(residues.items())},
        "group_atom_counts": {name: len(atoms) for name, atoms in groups.items()},
        "group_order": list(groups.keys()),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
