#!/usr/bin/env python3
"""Plot GROMACS XVG outputs and derive simple MD validation summaries."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_xvg(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    if not path.exists():
        return rows
    for line in path.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    return rows


def stats(rows: list[tuple[float, float]]) -> dict[str, float | None]:
    vals = [value for _, value in rows if math.isfinite(value)]
    if not vals:
        return {"mean": None, "min": None, "max": None, "last": None}
    return {
        "mean": sum(vals) / len(vals),
        "min": min(vals),
        "max": max(vals),
        "last": vals[-1],
    }


def save_csv(path: Path, rows: list[tuple[float, float]], x_name: str = "time", y_name: str = "value") -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([x_name, y_name])
        writer.writerows(rows)


def plot_series(output: Path, title: str, ylabel: str, series: list[tuple[str, list[tuple[float, float]]]]) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for label, rows in series:
        if rows:
            xs = [x for x, _ in rows]
            ys = [y for _, y in rows]
            ax.plot(xs, ys, label=label, linewidth=1.6)
    ax.set_title(title)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    if len(series) > 1:
        ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    analysis_dir = Path(args.analysis_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "rmsd_complex": analysis_dir / "rmsd_complex.xvg",
        "rmsd_antibody": analysis_dir / "rmsd_antibody.xvg",
        "rmsd_antigen": analysis_dir / "rmsd_antigen.xvg",
        "rmsf_complex_ca": analysis_dir / "rmsf_complex_ca.xvg",
        "rmsf_antibody_ca": analysis_dir / "rmsf_antibody_ca.xvg",
        "rmsf_antigen_ca": analysis_dir / "rmsf_antigen_ca.xvg",
        "rg_complex": analysis_dir / "rg_complex.xvg",
        "rg_antibody": analysis_dir / "rg_antibody.xvg",
        "rg_antigen": analysis_dir / "rg_antigen.xvg",
        "sasa_complex": analysis_dir / "sasa_complex.xvg",
        "sasa_antibody": analysis_dir / "sasa_antibody.xvg",
        "sasa_antigen": analysis_dir / "sasa_antigen.xvg",
        "hbonds": analysis_dir / "hbonds_antibody_antigen.xvg",
        "interface_contacts": analysis_dir / "interface_contacts.xvg",
        "interface_mindist": analysis_dir / "interface_mindist.xvg",
    }
    data = {name: read_xvg(path) for name, path in files.items()}

    for name, rows in data.items():
        save_csv(analysis_dir / f"{name}.csv", rows)

    plot_series(
        output_dir / "rmsd.png",
        "Backbone RMSD",
        "RMSD (nm)",
        [
            ("complex", data["rmsd_complex"]),
            ("antibody", data["rmsd_antibody"]),
            ("antigen", data["rmsd_antigen"]),
        ],
    )
    plot_series(
        output_dir / "rmsf.png",
        "Residue RMSF (C-alpha)",
        "RMSF (nm)",
        [
            ("complex", data["rmsf_complex_ca"]),
            ("antibody", data["rmsf_antibody_ca"]),
            ("antigen", data["rmsf_antigen_ca"]),
        ],
    )
    plot_series(
        output_dir / "rg.png",
        "Radius of Gyration",
        "Rg (nm)",
        [
            ("complex", data["rg_complex"]),
            ("antibody", data["rg_antibody"]),
            ("antigen", data["rg_antigen"]),
        ],
    )
    plot_series(
        output_dir / "sasa.png",
        "Solvent Accessible Surface Area",
        "SASA (nm^2)",
        [
            ("complex", data["sasa_complex"]),
            ("antibody", data["sasa_antibody"]),
            ("antigen", data["sasa_antigen"]),
        ],
    )
    plot_series(output_dir / "hbonds.png", "Antibody-Antigen Hydrogen Bonds", "H-bonds", [("h-bonds", data["hbonds"])])
    plot_series(
        output_dir / "interface_contacts.png",
        "Antibody-Antigen Interface Contacts",
        "Contacts within 0.45 nm",
        [("contacts", data["interface_contacts"])],
    )
    plot_series(
        output_dir / "interface_mindist.png",
        "Antibody-Antigen Minimum Distance",
        "Distance (nm)",
        [("minimum distance", data["interface_mindist"])],
    )

    buried_sasa: list[tuple[float, float]] = []
    for (t1, ab), (t2, ag), (t3, cx) in zip(data["sasa_antibody"], data["sasa_antigen"], data["sasa_complex"]):
        if abs(t1 - t2) < 1e-6 and abs(t1 - t3) < 1e-6:
            buried_sasa.append((t1, ab + ag - cx))
    save_csv(analysis_dir / "buried_sasa.csv", buried_sasa, "time_ns", "buried_sasa_nm2")
    plot_series(output_dir / "buried_sasa.png", "Buried SASA", "Buried SASA (nm^2)", [("buried SASA", buried_sasa)])

    # NIS proxy: solvent-exposed surface not buried in the interface.
    nis: list[tuple[float, float]] = []
    for (t1, cx), (t2, buried) in zip(data["sasa_complex"], buried_sasa):
        if abs(t1 - t2) < 1e-6:
            nis.append((t1, max(cx - buried, 0.0)))
    save_csv(analysis_dir / "nis_proxy.csv", nis, "time_ns", "nis_proxy_nm2")
    plot_series(output_dir / "nis.png", "NIS Proxy from SASA and Interface Burial", "NIS proxy (nm^2)", [("NIS proxy", nis)])

    summary_stats = {name: stats(rows) for name, rows in data.items()}
    summary_stats["buried_sasa"] = stats(buried_sasa)
    summary_stats["nis_proxy"] = stats(nis)
    (analysis_dir / "summary_stats.json").write_text(json.dumps(summary_stats, indent=2), encoding="utf-8")

    lines = [
        "# GROMACS MD summary for antibody_b007_3",
        "",
        "- Candidate: antibody_b007_3",
        "- Antibody chains: A+B",
        "- Antigen chain: C",
        "- Production MD target: 10 ns",
        "- Analyses generated: ICs, NIS proxy, RMSD, hydrogen bonds, Rg, SASA, RMSF",
        "",
        "## Key statistics",
    ]
    for name in [
        "rmsd_complex",
        "rmsd_antibody",
        "rmsd_antigen",
        "hbonds",
        "interface_contacts",
        "interface_mindist",
        "rg_complex",
        "sasa_complex",
        "buried_sasa",
        "nis_proxy",
    ]:
        item = summary_stats.get(name, {})
        lines.append(f"- {name}: mean={item.get('mean')} min={item.get('min')} max={item.get('max')} last={item.get('last')}")
    lines.extend(
        [
            "",
            "## Output figures",
            "- rmsd.png",
            "- rmsf.png",
            "- hbonds.png",
            "- rg.png",
            "- sasa.png",
            "- buried_sasa.png",
            "- interface_contacts.png",
            "- nis.png",
        ]
    )
    Path(args.summary).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary_stats, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
