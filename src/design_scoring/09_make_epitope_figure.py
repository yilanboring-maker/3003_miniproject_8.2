#!/usr/bin/env python3
"""Create a reproducible epitope-location figure from Antigen1 coordinates."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HOTSPOTS = {
    "A": {87, 88, 89, 90, 91, 114, 115, 116, 117},
    "C": {87, 88, 89, 90, 91, 114, 115, 116, 117},
}


def read_ca_atoms(path: Path) -> dict[str, list[tuple[int, np.ndarray]]]:
    chains: dict[str, list[tuple[int, np.ndarray]]] = {}
    seen: set[tuple[str, int]] = set()
    with path.open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom = line[12:16].strip()
            if atom != "CA":
                continue
            chain = line[21].strip()
            try:
                resi = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            key = (chain, resi)
            if key in seen:
                continue
            seen.add(key)
            chains.setdefault(chain, []).append((resi, xyz))
    for chain in chains:
        chains[chain].sort(key=lambda item: item[0])
    return chains


def set_equal_axes(ax) -> None:
    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()
    ranges = [abs(xlim[1] - xlim[0]), abs(ylim[1] - ylim[0]), abs(zlim[1] - zlim[0])]
    radius = max(ranges) / 2
    centers = [sum(xlim) / 2, sum(ylim) / 2, sum(zlim) / 2]
    ax.set_xlim3d(centers[0] - radius, centers[0] + radius)
    ax.set_ylim3d(centers[1] - radius, centers[1] + radius)
    ax.set_zlim3d(centers[2] - radius, centers[2] + radius)


def plot_chain(ax, chain_id: str, atoms: list[tuple[int, np.ndarray]], title: str, color: str) -> None:
    coords = np.vstack([xyz for _, xyz in atoms])
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color="#9aa0a6", linewidth=1.0, alpha=0.7)
    hotspot_coords = np.vstack([xyz for resi, xyz in atoms if resi in HOTSPOTS[chain_id]])
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=6, color="#b8bec4", alpha=0.55)
    ax.scatter(
        hotspot_coords[:, 0],
        hotspot_coords[:, 1],
        hotspot_coords[:, 2],
        s=56,
        color=color,
        edgecolor="#222222",
        linewidth=0.45,
        depthshade=False,
        label="Hotspot residues 87-91 and 114-117",
    )
    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_box_aspect((1, 1, 1))
    set_equal_axes(ax)
    ax.view_init(elev=18, azim=-58)
    ax.legend(loc="lower left", fontsize=7, frameon=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-pdb", required=True)
    parser.add_argument("--output-png", required=True)
    args = parser.parse_args()

    chains = read_ca_atoms(Path(args.input_pdb))
    missing = [chain for chain in ["A", "C"] if chain not in chains]
    if missing:
        raise SystemExit(f"Missing chains in PDB: {missing}")

    fig = plt.figure(figsize=(10.5, 4.2))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    plot_chain(ax1, "C", chains["C"], "IgG design epitope on antigen chain C", "#d44f3a")
    plot_chain(ax2, "A", chains["A"], "Nanobody design epitope on antigen chain A", "#2f7f5f")
    fig.suptitle("Selected Antigen1 hotspot patch used for PPIFlow design", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    output = Path(args.output_png)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    print(f"saved {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
