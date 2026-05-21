#!/usr/bin/env python
"""Prepare PyMOL-compatible visual review artifacts for Task 4 AF3 poses.

The local machine may not have PyMOL installed. This script therefore creates:
- per-candidate PyMOL .pml scripts for formal review in PyMOL;
- static PNG inspection panels from the same aligned design/AF3 structures;
- review CSV/checkpoint files with progress output.
"""

from __future__ import annotations

import argparse
import json
import math
import textwrap
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser


HOTSPOTS = [87, 88, 89, 90, 91, 114, 115, 116, 117]
ANTIBODY_CHAINS = ("A", "B")
ANTIGEN_CHAIN = "C"


def now_text() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def as_path(value: object) -> Path:
    return Path(str(value)).expanduser()


def protein_residues(chain):
    return [res for res in chain if res.id[0] == " "]


def ca_coords(structure, chain_ids: tuple[str, ...] | list[str]) -> np.ndarray:
    coords = []
    for chain_id in chain_ids:
        if chain_id not in structure[0]:
            continue
        for res in protein_residues(structure[0][chain_id]):
            if "CA" in res:
                coords.append(res["CA"].coord)
    if not coords:
        return np.empty((0, 3), dtype=float)
    return np.asarray(coords, dtype=float)


def hotspot_coords(structure) -> np.ndarray:
    coords = []
    chain = structure[0][ANTIGEN_CHAIN]
    for res in protein_residues(chain):
        if res.id[1] in HOTSPOTS and "CA" in res:
            coords.append(res["CA"].coord)
    if not coords:
        return np.empty((0, 3), dtype=float)
    return np.asarray(coords, dtype=float)


def atom_coords(structure, chain_ids: tuple[str, ...] | list[str]) -> np.ndarray:
    coords = []
    for chain_id in chain_ids:
        if chain_id not in structure[0]:
            continue
        for res in protein_residues(structure[0][chain_id]):
            for atom in res:
                if atom.element != "H":
                    coords.append(atom.coord)
    if not coords:
        return np.empty((0, 3), dtype=float)
    return np.asarray(coords, dtype=float)


def centroid(coords: np.ndarray) -> np.ndarray:
    if coords.size == 0:
        return np.zeros(3, dtype=float)
    return coords.mean(axis=0)


def set_equal_axes(ax, coords: np.ndarray, pad: float = 4.0) -> None:
    if coords.size == 0:
        return
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    center = (mins + maxs) / 2.0
    radius = max(float((maxs - mins).max()) / 2.0 + pad, 8.0)
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)


def plot_trace(ax, coords: np.ndarray, color: str, label: str, alpha: float = 1.0, lw: float = 1.5) -> None:
    if coords.size == 0:
        return
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color=color, alpha=alpha, lw=lw, label=label)


def plot_points(ax, coords: np.ndarray, color: str, label: str, size: float = 36, alpha: float = 1.0, marker: str = "o") -> None:
    if coords.size == 0:
        return
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=color, s=size, alpha=alpha, marker=marker, label=label, depthshade=True)


def write_pymol_script(row: pd.Series, out_path: Path) -> None:
    design_pdb = Path(str(row["design_pdb_path"])).as_posix()
    af3_pdb = Path(str(row["aligned_pdb_path"])).as_posix()
    design_id = str(row["design_id"])
    hotspot_expr = "+".join(str(x) for x in HOTSPOTS)
    session_path = out_path.with_suffix(".pse").as_posix()
    png_path = out_path.with_suffix(".png").as_posix()
    text = f"""\
    # PyMOL visual check for {design_id}
    # Design antibody: cyan; AF3 antibody: magenta; antigen: grey; hotspots: orange/red.
    reinitialize
    bg_color white
    set ray_opaque_background, off
    set cartoon_fancy_helices, 1
    set antialias, 2
    set transparency_mode, 1

    load {design_pdb}, design
    load {af3_pdb}, af3_aligned

    hide everything
    show cartoon, design
    show cartoon, af3_aligned

    color grey80, design and chain C
    color cyan, design and chain A+B
    color grey60, af3_aligned and chain C
    color magenta, af3_aligned and chain A+B
    set cartoon_transparency, 0.55, design and chain A+B
    set cartoon_transparency, 0.15, af3_aligned and chain A+B
    set cartoon_transparency, 0.65, af3_aligned and chain C

    select design_hotspot, design and chain C and resi {hotspot_expr}
    select af3_hotspot, af3_aligned and chain C and resi {hotspot_expr}
    show spheres, design_hotspot
    show spheres, af3_hotspot
    color orange, design_hotspot
    color red, af3_hotspot
    set sphere_scale, 0.7, design_hotspot or af3_hotspot

    select design_interface_ab, design and chain A+B within 6 of design_hotspot
    select af3_interface_ab, af3_aligned and chain A+B within 6 of af3_hotspot
    show sticks, design_interface_ab or af3_interface_ab
    color blue, design_interface_ab
    color tv_magenta, af3_interface_ab

    distance af3_hotspot_contacts, af3_interface_ab, af3_hotspot, 6
    hide labels, af3_hotspot_contacts
    color yellow, af3_hotspot_contacts

    zoom design_hotspot or af3_hotspot or design_interface_ab or af3_interface_ab, 12
    orient design and chain C
    png {png_path}, width=1800, height=1200, dpi=180, ray=1
    save {session_path}
    """
    out_path.write_text(textwrap.dedent(text), encoding="utf-8")


@dataclass
class ReviewArtifact:
    submission_rank: int
    design_id: str
    optimized_prodigy_delta_g: float
    hotspot_contact_6a_count: int
    centroid_distance: float
    ranking_score: float | None
    iptm: float | None
    pml_path: str
    png_path: str
    visual_priority: str
    notes: str


def make_png(row: pd.Series, out_path: Path, parser: PDBParser) -> str:
    design = parser.get_structure(f"{row['design_id']}_design", str(row["design_pdb_path"]))
    af3 = parser.get_structure(f"{row['design_id']}_af3", str(row["aligned_pdb_path"]))

    antigen = ca_coords(design, [ANTIGEN_CHAIN])
    design_ab = ca_coords(design, ANTIBODY_CHAINS)
    af3_ab = ca_coords(af3, ANTIBODY_CHAINS)
    hs = hotspot_coords(design)
    all_coords = np.vstack([x for x in [antigen, design_ab, af3_ab, hs] if x.size])

    hs_center = centroid(hs)
    design_ab_center = centroid(design_ab)
    af3_ab_center = centroid(af3_ab)
    title = (
        f"{row['submission_rank']}. {row['design_id']} | "
        f"dG={float(row['optimized_prodigy_delta_g']):.3f} | "
        f"HS6A={int(float(row['hotspot_contact_6a_count']))}/9 | "
        f"centroid={float(row['interface_centroid_distance_to_design_a']):.2f} A"
    )

    fig = plt.figure(figsize=(16, 10))
    views = [
        (25, 35, "overlay view"),
        (10, 115, "side view"),
        (70, -75, "top view"),
    ]
    for i, (elev, azim, label) in enumerate(views, start=1):
        ax = fig.add_subplot(2, 2, i, projection="3d")
        plot_trace(ax, antigen, "#9a9a9a", "antigen C", alpha=0.9, lw=1.6)
        plot_trace(ax, design_ab, "#3f8fd2", "design antibody A+B", alpha=0.35, lw=1.2)
        plot_trace(ax, af3_ab, "#b21fb3", "AF3 antibody A+B", alpha=0.9, lw=1.5)
        plot_points(ax, hs, "#ff9f1c", "hotspots C87-91/C114-117", size=54, alpha=1.0)
        plot_points(ax, np.asarray([design_ab_center]), "#3f8fd2", "design Ab centroid", size=80, marker="^")
        plot_points(ax, np.asarray([af3_ab_center]), "#b21fb3", "AF3 Ab centroid", size=80, marker="^")
        ax.plot(
            [hs_center[0], af3_ab_center[0]],
            [hs_center[1], af3_ab_center[1]],
            [hs_center[2], af3_ab_center[2]],
            color="#f4d35e",
            lw=2.0,
            alpha=0.9,
        )
        set_equal_axes(ax, all_coords, pad=5.0)
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(label, fontsize=11)
        ax.set_axis_off()

    ax = fig.add_subplot(2, 2, 4)
    ax.axis("off")
    summary_lines = [
        title,
        f"auto_pose_call: {row['auto_pose_call']}",
        f"hotspot_contact_5A_count: {row['hotspot_contact_5a_count']}",
        f"min_hotspot_distance_A: {float(row['min_hotspot_distance_a']):.2f}",
        f"ranking_score: {row.get('ranking_score', '')}",
        f"ipTM: {row.get('iptm', '')}",
        "",
        "Visual check guide:",
        "1. Orange hotspots should remain next to magenta AF3 antibody.",
        "2. Magenta AF3 antibody should overlap/neighbor cyan design antibody.",
        "3. Large shift away from orange hotspot = reject/manual concern.",
        "",
        "This PNG is a review aid; formal PyMOL review can use the paired .pml file.",
    ]
    ax.text(0.02, 0.98, "\n".join(summary_lines), va="top", ha="left", fontsize=11, family="monospace")
    handles, labels = fig.axes[0].get_legend_handles_labels()
    fig.legend(handles[:6], labels[:6], loc="lower center", ncol=3, fontsize=10)
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return str(out_path)


def make_montages(png_paths: list[Path], montage_dir: Path, per_sheet: int = 12) -> list[Path]:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return []
    montage_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    thumb_w, thumb_h = 480, 300
    cols = 3
    rows = math.ceil(per_sheet / cols)
    for sheet_idx, start in enumerate(range(0, len(png_paths), per_sheet), start=1):
        chunk = png_paths[start : start + per_sheet]
        canvas = Image.new("RGB", (cols * thumb_w, rows * thumb_h), "white")
        draw = ImageDraw.Draw(canvas)
        for j, path in enumerate(chunk):
            img = Image.open(path).convert("RGB")
            img.thumbnail((thumb_w, thumb_h - 22))
            x = (j % cols) * thumb_w
            y = (j // cols) * thumb_h
            canvas.paste(img, (x + (thumb_w - img.width) // 2, y + 20))
            draw.text((x + 8, y + 4), path.stem[:60], fill=(0, 0, 0))
        out = montage_dir / f"acceptable_montage_{sheet_idx:02d}.png"
        canvas.save(out)
        outputs.append(out)
    return outputs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--input-csv", default="", help="Default: af3_pose_review_acceptable.csv")
    ap.add_argument("--out-subdir", default="pymol_visual_review")
    ap.add_argument("--label", default="acceptable_candidates")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    root = Path(args.root).resolve()
    input_csv = Path(args.input_csv) if args.input_csv else root / "af3_pose_review_acceptable.csv"
    out_root = root / args.out_subdir
    pml_dir = out_root / "pml"
    png_dir = out_root / "png"
    montage_dir = out_root / "montages"
    for d in [out_root, pml_dir, png_dir, montage_dir]:
        d.mkdir(parents=True, exist_ok=True)

    checkpoint_path = out_root / "checkpoint.json"
    partial_path = out_root / "partial_results.csv"
    all_path = out_root / f"pymol_visual_review_{args.label}.csv"
    summary_path = out_root / "pymol_visual_review_summary.md"

    df = pd.read_csv(input_csv)
    df = df.sort_values("optimized_prodigy_delta_g", key=lambda s: s.astype(float)).reset_index(drop=True)
    if args.limit and args.limit > 0:
        df = df.head(args.limit).copy()

    print(f"[{now_text()}] PyMOL visual review prep start total={len(df)} input={input_csv} out={out_root}", flush=True)
    print(f"[{now_text()}] outputs: pml={pml_dir}; png={png_dir}; checkpoint={checkpoint_path}", flush=True)

    parser = PDBParser(QUIET=True)
    rows: list[ReviewArtifact] = []
    png_paths: list[Path] = []
    start = time.time()
    for idx, (_, row) in enumerate(df.iterrows(), start=1):
        design_id = str(row["design_id"])
        rank = int(row["submission_rank"])
        pml_path = pml_dir / f"{rank:03d}_{design_id}_visual_check.pml"
        png_path = png_dir / f"{rank:03d}_{design_id}_visual_check.png"
        try:
            write_pymol_script(row, pml_path)
            make_png(row, png_path, parser)
            png_paths.append(png_path)
            priority = "pymol_confirm_priority"
            notes = "PML and PNG generated. Candidate remains acceptable by automated AF3 pose metrics; formal PyMOL visual confirmation still required."
        except Exception as exc:
            priority = "artifact_error"
            notes = repr(exc)
        artifact = ReviewArtifact(
            submission_rank=rank,
            design_id=design_id,
            optimized_prodigy_delta_g=float(row["optimized_prodigy_delta_g"]),
            hotspot_contact_6a_count=int(float(row["hotspot_contact_6a_count"])),
            centroid_distance=float(row["interface_centroid_distance_to_design_a"]),
            ranking_score=float(row["ranking_score"]) if not pd.isna(row.get("ranking_score")) and str(row.get("ranking_score")) != "" else None,
            iptm=float(row["iptm"]) if not pd.isna(row.get("iptm")) and str(row.get("iptm")) != "" else None,
            pml_path=str(pml_path),
            png_path=str(png_path),
            visual_priority=priority,
            notes=notes,
        )
        rows.append(artifact)
        out_df = pd.DataFrame([r.__dict__ for r in rows])
        out_df.to_csv(partial_path, index=False)
        elapsed = time.time() - start
        rate = idx / elapsed if elapsed else 0
        eta = (len(df) - idx) / rate if rate else 0
        checkpoint = {
            "updated_at": now_text(),
            "stage": "pymol_visual_review_preparation",
            "done": idx,
            "total": len(df),
            "current_design_id": design_id,
            "partial_results": str(partial_path),
            "all_results": str(all_path),
            "pml_dir": str(pml_dir),
            "png_dir": str(png_dir),
        }
        checkpoint_path.write_text(json.dumps(checkpoint, indent=2, ensure_ascii=False), encoding="utf-8")
        print(
            f"[{now_text()}] done={idx}/{len(df)} ({idx/len(df)*100:.1f}%) "
            f"current={design_id} dG={float(row['optimized_prodigy_delta_g']):.3f} "
            f"HS6A={row['hotspot_contact_6a_count']} elapsed={elapsed:.1f}s eta={eta:.1f}s status={priority}",
            flush=True,
        )

    all_df = pd.DataFrame([r.__dict__ for r in rows])
    all_df.to_csv(all_path, index=False)
    montages = make_montages(png_paths, montage_dir)

    top = all_df.sort_values("optimized_prodigy_delta_g").head(12)
    summary = [
        "# Task 4 PyMOL Visual Review Preparation",
        "",
        f"Updated: {now_text()}",
        "",
        f"Input: `{input_csv}`",
        f"Candidates prepared: {len(all_df)}",
        "",
        "Local PyMOL status: no PyMOL executable/module was available in this Windows environment, so this step generated PyMOL-compatible `.pml` scripts plus static review PNGs. Open the `.pml` files in PyMOL for formal visual confirmation.",
        "",
        "Outputs:",
        f"- PML scripts: `{pml_dir}`",
        f"- PNG inspection panels: `{png_dir}`",
        f"- Montage sheets: `{montage_dir}`",
        f"- Review table: `{all_path}`",
        "",
        "Top candidates by PRODIGY Delta G in this input table:",
        "",
    ]
    for _, row in top.iterrows():
        summary.append(
            f"- rank {int(row['submission_rank'])}: `{row['design_id']}`, "
            f"dG={row['optimized_prodigy_delta_g']:.3f}, "
            f"HS6A={int(row['hotspot_contact_6a_count'])}/9, "
            f"centroid={row['centroid_distance']:.2f} A"
        )
    summary_path.write_text("\n".join(summary) + "\n", encoding="utf-8")

    final_checkpoint = {
        "updated_at": now_text(),
        "stage": "pymol_visual_review_preparation_complete",
        "total": len(all_df),
        "all_results": str(all_path),
        "summary": str(summary_path),
        "pml_dir": str(pml_dir),
        "png_dir": str(png_dir),
        "montages": [str(p) for p in montages],
    }
    checkpoint_path.write_text(json.dumps(final_checkpoint, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[{now_text()}] PyMOL visual review prep complete summary={summary_path}", flush=True)
    print(f"[{now_text()}] montage_count={len(montages)} all_results={all_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
