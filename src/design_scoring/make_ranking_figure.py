#!/usr/bin/env python3
"""Create a compact ranking figure from PRODIGY-scored candidates."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output-png", required=True)
    parser.add_argument("--top-n", type=int, default=6)
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)
    df = df[df["prodigy_status"].eq("ok")].copy()
    df["prodigy_delta_g"] = pd.to_numeric(df["prodigy_delta_g"], errors="coerce")
    df["cdr_interface_ratio"] = pd.to_numeric(df["cdr_interface_ratio"], errors="coerce")
    df["hotspot_coverage"] = pd.to_numeric(df["hotspot_coverage"], errors="coerce")
    df = df.dropna(subset=["prodigy_delta_g"])

    panels = []
    for modality in ["antibody", "nanobody"]:
        sub = df[df["modality"].eq(modality)].sort_values("prodigy_delta_g").head(args.top_n)
        panels.append((modality, sub))

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharex=True)
    colors = {"antibody": "#3b7ea1", "nanobody": "#b55d4c"}
    for ax, (modality, sub) in zip(axes, panels):
        y = range(len(sub))
        prefix = "Ab" if modality == "antibody" else "Nb"
        labels = [f"{prefix} {label.removeprefix(modality + '_')}" for label in sub["design_id"]]
        ax.barh(y, sub["prodigy_delta_g"], color=colors[modality], alpha=0.88)
        ax.set_yticks(list(y), labels=labels, fontsize=8)
        ax.invert_yaxis()
        ax.axvline(0, color="#333333", linewidth=0.8)
        ax.grid(axis="x", color="#d8d8d8", linewidth=0.7)
        ax.set_title(f"{modality.capitalize()} top candidates", fontsize=11, weight="bold")
        ax.set_xlabel("PRODIGY predicted delta G (kcal/mol)")
        for i, (_, row) in enumerate(sub.iterrows()):
            ax.text(
                row["prodigy_delta_g"] + 0.08,
                i,
                f"{row['prodigy_delta_g']:.2f}",
                va="center",
                ha="left",
                fontsize=7,
                color="white",
            )

    fig.suptitle("Top PPIFlow designs after interface filtering and PRODIGY scoring", fontsize=12)
    fig.tight_layout(rect=(0.02, 0, 1, 0.93), w_pad=3.0)
    output = Path(args.output_png)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    print(f"saved {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
