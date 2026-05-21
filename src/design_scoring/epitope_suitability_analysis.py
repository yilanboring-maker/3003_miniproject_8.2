#!/usr/bin/env python3
"""Assess whether the selected Antigen1 hotspot patch is structurally suitable."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley


ROOT = Path(__file__).resolve().parents[1]
PDB_PATH = ROOT / "data" / "Antigen1_Tetramer.pdb"
OUT_CSV = ROOT / "results" / "epitope_suitability_summary.csv"
OUT_JSON = ROOT / "results" / "epitope_suitability_summary.json"
OUT_MD = ROOT / "report" / "epitope_suitability_note.md"
HOTSPOTS = [87, 88, 89, 90, 91, 114, 115, 116, 117]


def ca_coord(residue) -> np.ndarray | None:
    if "CA" not in residue:
        return None
    atom = residue["CA"]
    return np.array(atom.coord, dtype=float)


def pairwise_distances(coords: list[np.ndarray]) -> list[float]:
    values = []
    for i, a in enumerate(coords):
        for b in coords[i + 1 :]:
            values.append(float(np.linalg.norm(a - b)))
    return values


def analyse_chain(chain) -> dict[str, object]:
    standard_residues = [res for res in chain if res.id[0] == " " and "CA" in res]
    residue_sasa = []
    hotspot_rows = []
    hotspot_coords = []

    for res in standard_residues:
        resseq = int(res.id[1])
        sasa = float(getattr(res, "sasa", 0.0))
        residue_sasa.append(sasa)
        if resseq in HOTSPOTS:
            coord = ca_coord(res)
            if coord is not None:
                hotspot_coords.append(coord)
            hotspot_rows.append(
                {
                    "resi": resseq,
                    "resname": res.resname,
                    "sasa": round(sasa, 3),
                    "above_chain_median": False,
                }
            )

    median_sasa = float(np.median(residue_sasa))
    mean_sasa = float(np.mean(residue_sasa))
    for row in hotspot_rows:
        row["above_chain_median"] = row["sasa"] >= median_sasa

    dists = pairwise_distances(hotspot_coords)
    centroid = np.mean(np.vstack(hotspot_coords), axis=0)
    centroid_distances = [float(np.linalg.norm(coord - centroid)) for coord in hotspot_coords]

    return {
        "chain": chain.id,
        "chain_residue_count": len(standard_residues),
        "chain_mean_sasa": round(mean_sasa, 3),
        "chain_median_sasa": round(median_sasa, 3),
        "hotspot_count": len(hotspot_rows),
        "hotspot_mean_sasa": round(float(np.mean([row["sasa"] for row in hotspot_rows])), 3),
        "hotspot_median_sasa": round(float(np.median([row["sasa"] for row in hotspot_rows])), 3),
        "hotspots_above_chain_median": sum(1 for row in hotspot_rows if row["above_chain_median"]),
        "hotspot_max_ca_distance": round(max(dists), 3),
        "hotspot_mean_ca_distance": round(float(np.mean(dists)), 3),
        "hotspot_max_centroid_distance": round(max(centroid_distances), 3),
        "hotspots": hotspot_rows,
    }


def main() -> int:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("Antigen1", str(PDB_PATH))
    ShrakeRupley(n_points=100).compute(structure, level="R")

    model = structure[0]
    summaries = [analyse_chain(model[chain_id]) for chain_id in ["A", "C"]]

    OUT_JSON.write_text(json.dumps(summaries, indent=2), encoding="ascii")
    with OUT_CSV.open("w", encoding="ascii", newline="") as handle:
        fields = [
            "chain",
            "chain_residue_count",
            "chain_mean_sasa",
            "chain_median_sasa",
            "hotspot_count",
            "hotspot_mean_sasa",
            "hotspot_median_sasa",
            "hotspots_above_chain_median",
            "hotspot_max_ca_distance",
            "hotspot_mean_ca_distance",
            "hotspot_max_centroid_distance",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in summaries:
            writer.writerow({field: row[field] for field in fields})

    lines = [
        "# Epitope suitability note",
        "",
        "The selected hotspot patch is residues 87-91 and 114-117 on Antigen1.",
        "Biopython Shrake-Rupley SASA was used as a reproducible local surface-accessibility check.",
        "",
    ]
    for row in summaries:
        lines.extend(
            [
                f"- Chain {row['chain']}: hotspot mean SASA {row['hotspot_mean_sasa']} A^2 "
                f"versus chain mean SASA {row['chain_mean_sasa']} A^2; "
                f"{row['hotspots_above_chain_median']}/{row['hotspot_count']} hotspots are above the chain median SASA. "
                f"The hotspot C-alpha cluster has max pairwise distance {row['hotspot_max_ca_distance']} A "
                f"and max centroid distance {row['hotspot_max_centroid_distance']} A.",
            ]
        )
    lines.extend(
        [
            "",
            "Interpretation: this does not prove global optimality, but it supports the patch as an accessible, spatially defined epitope candidate consistent with the course workshop PPIFlow setup.",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="ascii")

    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_JSON}")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
