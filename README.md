# ICA2 Antibody Design Code

This repository contains the core code/technical files for the ICA2 antibody
design project. It is deliberately minimal for submission: it does not include
the report text, submission PDFs, browser automation, Jupyter API helpers,
download/monitoring scripts, package builders, large AlphaFold Server result
archives, GROMACS trajectories or model weights.

## What is included

- `src/design_scoring/`: PPIFlow metric parsing, PRODIGY scoring, epitope checks
  and figure/table generation.
- `src/binder/`: mini-protein binder result collection and scoring.
- `src/optimization/`: corrected AbMPNN + FlowPacker optimisation and PRODIGY
  rescoring scripts.
- `src/alphafold_validation/`: AlphaFold Server input preparation, pose review
  and PyMOL review helper scripts.
- `src/md_analysis/`: GROMACS index generation, MD analysis, plotting and staged
  decision scripts.
- `example_tables/`: small CSV tables that show the main output formats used by
  the report. Full raw AF3/MD outputs are intentionally not included.

## External tools

The full workflow used PPIFlow, PRODIGY, ProteinMPNN/AbMPNN, FlowPacker,
AlphaFold Server, PyMOL and GROMACS. These third-party tools, databases, model
weights and raw simulation outputs are not stored in this repository.

## Reproducibility note

Some scripts preserve the original local or Matpool paths used in the project.
When rerunning, update those paths or pass equivalent command-line arguments.
Candidate IDs such as `antibody_b003_8` are internal workflow IDs kept for
traceability across PPIFlow, PRODIGY, AlphaFold Server/PyMOL and MD results.
