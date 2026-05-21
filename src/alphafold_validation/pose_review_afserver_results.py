#!/usr/bin/env python
"""Review AlphaFold Server complex poses against the original designed epitope.

This script is deliberately checkpointed and verbose because Task 4 is a
long-running web/server workflow. It does not submit jobs. It only parses
downloaded AlphaFold Server result packages already saved under af3_raw_downloads.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Superimposer


HOTSPOT_RESIDUES = [87, 88, 89, 90, 91, 114, 115, 116, 117]
ANTIBODY_CHAINS = ["A", "B"]
ANTIGEN_CHAIN = "C"


def now_text() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def protein_residues(chain):
    return [res for res in chain if res.id[0] == " "]


def chain_residue_map(structure, chain_id: str):
    chain = structure[0][chain_id]
    return {res.id[1]: res for res in protein_residues(chain)}


def atom_coord_array(residues) -> np.ndarray:
    coords = []
    for res in residues:
        for atom in res:
            if atom.element != "H":
                coords.append(atom.coord)
    if not coords:
        return np.empty((0, 3), dtype=float)
    return np.asarray(coords, dtype=float)


def residue_min_distance(residue, other_coords: np.ndarray) -> float:
    own = atom_coord_array([residue])
    if own.size == 0 or other_coords.size == 0:
        return math.inf
    diff = own[:, None, :] - other_coords[None, :, :]
    return float(np.sqrt(np.sum(diff * diff, axis=2)).min())


def ca_centroid(residue_map: dict[int, object], residue_numbers: Iterable[int]) -> np.ndarray | None:
    coords = []
    for resnum in residue_numbers:
        res = residue_map.get(int(resnum))
        if res is not None and "CA" in res:
            coords.append(res["CA"].coord)
    if not coords:
        return None
    return np.asarray(coords, dtype=float).mean(axis=0)


def interface_antigen_residues(structure, cutoff: float = 6.0) -> list[int]:
    antibody_residues = []
    for chain_id in ANTIBODY_CHAINS:
        antibody_residues.extend(protein_residues(structure[0][chain_id]))
    antibody_coords = atom_coord_array(antibody_residues)
    c_map = chain_residue_map(structure, ANTIGEN_CHAIN)
    interface = []
    for resnum, res in c_map.items():
        if residue_min_distance(res, antibody_coords) <= cutoff:
            interface.append(int(resnum))
    return sorted(interface)


def align_af3_to_design_antigen(design_structure, af3_structure) -> float:
    design_c = chain_residue_map(design_structure, ANTIGEN_CHAIN)
    af3_c = chain_residue_map(af3_structure, ANTIGEN_CHAIN)
    fixed = []
    moving = []
    for resnum in sorted(set(design_c) & set(af3_c)):
        if "CA" in design_c[resnum] and "CA" in af3_c[resnum]:
            fixed.append(design_c[resnum]["CA"])
            moving.append(af3_c[resnum]["CA"])
    if len(fixed) < 3:
        raise ValueError(f"Need at least 3 matched antigen CA atoms for alignment, got {len(fixed)}")
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    sup.apply(list(af3_structure.get_atoms()))
    return float(sup.rms)


def confidence_for_model(job_dir: Path, job_name: str, model_index: int) -> dict:
    path = job_dir / f"fold_{job_name}_summary_confidences_{model_index}.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def pick_best_model(job_dir: Path, job_name: str) -> tuple[int | None, list[dict]]:
    rows = []
    for summary_path in sorted(job_dir.glob(f"fold_{job_name}_summary_confidences_*.json")):
        model_index = int(summary_path.stem.rsplit("_", 1)[-1])
        data = json.loads(summary_path.read_text(encoding="utf-8"))
        rows.append(
            {
                "model_index": model_index,
                "ranking_score": data.get("ranking_score"),
                "iptm": data.get("iptm"),
                "ptm": data.get("ptm"),
                "fraction_disordered": data.get("fraction_disordered"),
                "has_clash": data.get("has_clash"),
            }
        )
    if not rows:
        return None, rows
    best = max(rows, key=lambda r: (float(r.get("ranking_score") or -999), -int(r["model_index"])))
    return int(best["model_index"]), rows


@dataclass
class ReviewResult:
    design_id: str
    job_name: str
    status: str
    auto_pose_call: str
    best_model_index: int | None = None
    antigen_alignment_rmsd: float | None = None
    hotspot_contact_6a_count: int | None = None
    hotspot_contact_5a_count: int | None = None
    min_hotspot_distance_a: float | None = None
    design_interface_antigen_residue_count: int | None = None
    af3_interface_antigen_residue_count: int | None = None
    interface_centroid_distance_to_design_a: float | None = None
    ranking_score: float | None = None
    iptm: float | None = None
    ptm: float | None = None
    fraction_disordered: float | None = None
    has_clash: float | None = None
    aligned_pdb_path: str = ""
    af3_model_cif_path: str = ""
    raw_download_path: str = ""
    notes: str = ""


def review_one(row: pd.Series, root: Path, parser_pdb: PDBParser, parser_cif: MMCIFParser) -> tuple[ReviewResult, list[dict]]:
    design_id = str(row["design_id"])
    job_name = str(row.get("alphafold_server_job_name") or "")
    if not job_name or job_name == "nan":
        return ReviewResult(design_id, "", "not_submitted", "not_reviewed", notes="No AlphaFold Server job name recorded."), []

    raw_root = root / "af3_raw_downloads"
    extract_dir = raw_root / f"{job_name}_afserver_raw" / job_name
    if not extract_dir.exists():
        design_prefixed = sorted(raw_root.glob(f"*{design_id}*_afserver_raw/{job_name}"))
        if design_prefixed:
            extract_dir = design_prefixed[0]
    if not extract_dir.exists():
        any_download_parent = sorted(raw_root.glob(f"*/{job_name}"))
        if any_download_parent:
            extract_dir = any_download_parent[0]
    if not extract_dir.exists():
        return ReviewResult(design_id, job_name, "not_downloaded", "not_reviewed", notes="Downloaded/extracted AF3 result directory not found."), []

    best_model, model_rows = pick_best_model(extract_dir, job_name)
    if best_model is None:
        return ReviewResult(design_id, job_name, "no_confidence_json", "not_reviewed", raw_download_path=str(extract_dir), notes="No summary_confidences JSON found."), model_rows

    cif_path = extract_dir / f"fold_{job_name}_model_{best_model}.cif"
    if not cif_path.exists():
        return ReviewResult(design_id, job_name, "no_model_cif", "not_reviewed", best_model_index=best_model, raw_download_path=str(extract_dir), notes=f"Missing model CIF for best model {best_model}."), model_rows

    design_pdb = Path(str(row["design_pdb_path"]))
    if not design_pdb.exists():
        return ReviewResult(design_id, job_name, "missing_design_pdb", "not_reviewed", best_model_index=best_model, af3_model_cif_path=str(cif_path), raw_download_path=str(extract_dir), notes=f"Missing design PDB: {design_pdb}"), model_rows

    design = parser_pdb.get_structure(f"{design_id}_design", str(design_pdb))
    af3 = parser_cif.get_structure(f"{design_id}_af3", str(cif_path))
    for chain_id in ["A", "B", "C"]:
        if chain_id not in design[0] or chain_id not in af3[0]:
            return ReviewResult(design_id, job_name, "missing_chain", "not_reviewed", best_model_index=best_model, af3_model_cif_path=str(cif_path), raw_download_path=str(extract_dir), notes=f"Missing chain {chain_id} in design or AF3 model."), model_rows

    antigen_rmsd = align_af3_to_design_antigen(design, af3)
    design_c_map = chain_residue_map(design, ANTIGEN_CHAIN)
    af3_c_map = chain_residue_map(af3, ANTIGEN_CHAIN)
    af3_antibody_res = []
    for chain_id in ANTIBODY_CHAINS:
        af3_antibody_res.extend(protein_residues(af3[0][chain_id]))
    af3_antibody_coords = atom_coord_array(af3_antibody_res)

    hotspot_distances = {}
    for resnum in HOTSPOT_RESIDUES:
        res = af3_c_map.get(resnum)
        hotspot_distances[resnum] = residue_min_distance(res, af3_antibody_coords) if res is not None else math.inf
    count6 = sum(1 for d in hotspot_distances.values() if d <= 6.0)
    count5 = sum(1 for d in hotspot_distances.values() if d <= 5.0)
    min_hotspot = min(hotspot_distances.values()) if hotspot_distances else math.inf

    design_interface = interface_antigen_residues(design, cutoff=6.0)
    af3_interface = interface_antigen_residues(af3, cutoff=6.0)
    design_centroid = ca_centroid(design_c_map, design_interface)
    af3_centroid = ca_centroid(af3_c_map, af3_interface)
    centroid_dist = math.inf
    if design_centroid is not None and af3_centroid is not None:
        centroid_dist = float(np.linalg.norm(af3_centroid - design_centroid))

    if count6 >= 7 and centroid_dist <= 8.0:
        auto_call = "acceptable"
    elif (4 <= count6 <= 6) or (8.0 < centroid_dist <= 15.0):
        auto_call = "needs_review"
    else:
        auto_call = "reject"

    aligned_dir = root / "review_tables" / "aligned_af3_models"
    aligned_dir.mkdir(parents=True, exist_ok=True)
    aligned_pdb = aligned_dir / f"{job_name}_model_{best_model}_aligned_to_design_antigen.pdb"
    io = PDBIO()
    io.set_structure(af3)
    io.save(str(aligned_pdb))

    conf = confidence_for_model(extract_dir, job_name, best_model)
    return (
        ReviewResult(
            design_id=design_id,
            job_name=job_name,
            status="reviewed",
            auto_pose_call=auto_call,
            best_model_index=best_model,
            antigen_alignment_rmsd=antigen_rmsd,
            hotspot_contact_6a_count=count6,
            hotspot_contact_5a_count=count5,
            min_hotspot_distance_a=float(min_hotspot),
            design_interface_antigen_residue_count=len(design_interface),
            af3_interface_antigen_residue_count=len(af3_interface),
            interface_centroid_distance_to_design_a=float(centroid_dist),
            ranking_score=conf.get("ranking_score"),
            iptm=conf.get("iptm"),
            ptm=conf.get("ptm"),
            fraction_disordered=conf.get("fraction_disordered"),
            has_clash=conf.get("has_clash"),
            aligned_pdb_path=str(aligned_pdb),
            af3_model_cif_path=str(cif_path),
            raw_download_path=str(extract_dir),
            notes="Automated antigen-aligned pose metrics; PyMOL visual review can use aligned_pdb_path plus design_pdb_path.",
        ),
        model_rows,
    )


def result_to_dict(result: ReviewResult, row: pd.Series) -> dict:
    base = row.to_dict()
    base.update(
        {
            "review_status": result.status,
            "auto_pose_call": result.auto_pose_call,
            "best_model_index": result.best_model_index,
            "antigen_alignment_rmsd": result.antigen_alignment_rmsd,
            "hotspot_contact_6a_count": result.hotspot_contact_6a_count,
            "hotspot_contact_5a_count": result.hotspot_contact_5a_count,
            "min_hotspot_distance_a": result.min_hotspot_distance_a,
            "design_interface_antigen_residue_count": result.design_interface_antigen_residue_count,
            "af3_interface_antigen_residue_count": result.af3_interface_antigen_residue_count,
            "interface_centroid_distance_to_design_a": result.interface_centroid_distance_to_design_a,
            "ranking_score": result.ranking_score,
            "iptm": result.iptm,
            "ptm": result.ptm,
            "fraction_disordered": result.fraction_disordered,
            "has_clash": result.has_clash,
            "aligned_pdb_path": result.aligned_pdb_path,
            "af3_model_cif_path": result.af3_model_cif_path,
            "raw_download_path": result.raw_download_path,
            "review_notes": result.notes,
        }
    )
    return base


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Task 4 output root")
    ap.add_argument("--only-design-id", action="append", default=[], help="Review only selected design_id; can repeat")
    args = ap.parse_args()

    root = Path(args.root)
    manifest_path = root / "af3_server_job_manifest.csv"
    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path)
    manifest = pd.read_csv(manifest_path)
    if args.only_design_id:
        manifest_to_review = manifest[manifest["design_id"].isin(args.only_design_id)].copy()
    else:
        manifest_to_review = manifest.copy()

    out_all = root / "af3_pose_review_all_results.csv"
    out_models = root / "review_tables" / "af3_pose_review_model_metrics.csv"
    out_models.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_path = root / "checkpoint.json"

    parser_pdb = PDBParser(QUIET=True)
    parser_cif = MMCIFParser(QUIET=True)
    results = []
    model_metric_rows = []
    start = time.time()
    total = len(manifest_to_review)
    print(f"[{now_text()}] Task4 pose review start total={total} root={root}", flush=True)
    print(f"[{now_text()}] outputs: {out_all}; checkpoint={checkpoint_path}", flush=True)

    for idx, (_, row) in enumerate(manifest_to_review.iterrows(), start=1):
        design_id = row["design_id"]
        try:
            result, model_rows = review_one(row, root, parser_pdb, parser_cif)
        except Exception as exc:
            result = ReviewResult(str(design_id), str(row.get("alphafold_server_job_name") or ""), "error", "not_reviewed", notes=repr(exc))
            model_rows = []
        results.append(result_to_dict(result, row))
        for mr in model_rows:
            mr.update({"design_id": design_id, "job_name": row.get("alphafold_server_job_name")})
            model_metric_rows.append(mr)

        reviewed = sum(1 for r in results if r["review_status"] == "reviewed")
        acceptable = sum(1 for r in results if r["auto_pose_call"] == "acceptable")
        elapsed = time.time() - start
        rate = idx / elapsed if elapsed > 0 else 0
        eta = (total - idx) / rate if rate > 0 else 0
        print(
            f"[{now_text()}] done={idx}/{total} ({idx/total*100:.1f}%) "
            f"reviewed={reviewed} acceptable={acceptable} "
            f"current={design_id} status={result.status} call={result.auto_pose_call} "
            f"elapsed={elapsed:.1f}s eta={eta:.1f}s",
            flush=True,
        )

        pd.DataFrame(results).to_csv(root / "partial_results.csv", index=False)
        checkpoint = {
            "updated_at": now_text(),
            "stage": "pose_review",
            "root": str(root),
            "done": idx,
            "total": total,
            "reviewed": reviewed,
            "acceptable": acceptable,
            "current_design_id": str(design_id),
            "all_results": str(out_all),
            "partial_results": str(root / "partial_results.csv"),
        }
        checkpoint_path.write_text(json.dumps(checkpoint, indent=2, ensure_ascii=False), encoding="utf-8")

    all_results = pd.DataFrame(results)
    if not args.only_design_id and len(all_results) < len(manifest):
        pass
    out_all.write_text(all_results.to_csv(index=False), encoding="utf-8")
    if model_metric_rows:
        pd.DataFrame(model_metric_rows).to_csv(out_models, index=False)
    else:
        out_models.write_text("", encoding="utf-8")

    for call, name in [
        ("acceptable", "af3_pose_review_acceptable.csv"),
        ("needs_review", "af3_pose_review_needs_review.csv"),
        ("reject", "af3_pose_review_reject.csv"),
    ]:
        subset = all_results[all_results["auto_pose_call"] == call]
        subset.to_csv(root / name, index=False)

    summary = {
        "updated_at": now_text(),
        "reviewed_count": int((all_results["review_status"] == "reviewed").sum()) if len(all_results) else 0,
        "auto_pose_call_counts": all_results["auto_pose_call"].value_counts(dropna=False).to_dict() if len(all_results) else {},
        "all_results": str(out_all),
        "model_metrics": str(out_models),
    }
    (root / "review_tables" / "af3_pose_review_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"[{now_text()}] Task4 pose review complete summary={summary}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
