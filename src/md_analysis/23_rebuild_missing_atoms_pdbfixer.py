#!/usr/bin/env python3
"""Rebuild missing heavy atoms with PDBFixer before all-atom GROMACS setup."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def count_atoms(path: Path) -> dict[str, object]:
    atoms = 0
    residues: set[tuple[str, str, str, str]] = set()
    atom_names: dict[str, int] = {}
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms += 1
        atom = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip()
        resseq = line[22:26].strip()
        icode = line[26].strip()
        residues.add((chain, resseq, icode, resname))
        atom_names[atom] = atom_names.get(atom, 0) + 1
    return {
        "atom_count": atoms,
        "residue_count": len(residues),
        "atom_name_counts": atom_names,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    input_path = Path(args.input)
    output_path = Path(args.output)
    summary_path = Path(args.summary)

    before = count_atoms(input_path)
    fixer = PDBFixer(filename=str(input_path))
    fixer.findMissingResidues()
    detected_missing_residues = {str(key): value for key, value in fixer.missingResidues.items()}
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    missing_atoms = {str(key): [str(atom) for atom in value] for key, value in fixer.missingAtoms.items()}
    missing_terminals = {str(key): [str(atom) for atom in value] for key, value in fixer.missingTerminals.items()}
    fixer.addMissingAtoms()
    with output_path.open("w", encoding="utf-8") as handle:
        try:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
        except TypeError:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle)
    after = count_atoms(output_path)
    summary = {
        "input": str(input_path),
        "output": str(output_path),
        "operation": "PDBFixer addMissingAtoms only; missing residues were detected but not added",
        "before": before,
        "after": after,
        "added_atom_count": int(after["atom_count"]) - int(before["atom_count"]),
        "detected_missing_residues_not_added": detected_missing_residues,
        "missing_atoms_added_by_template": missing_atoms,
        "missing_terminal_atoms_added_by_template": missing_terminals,
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
