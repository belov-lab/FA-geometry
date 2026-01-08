# FA Geometry — Ensemble Summaries from CSV

**Date:** 2026-01-08

This workflow computes Boltzmann-weighted geometric and volume metrics directly from an input CSV (columns: FA, SMILES).

**Parameters (defaults):** 64 conformers; voxel grid 0.5 Å; SAS probe 1.4 Å; T 298 K; MMFF94(100); ETKDGv3.

**Software versions:** RDKit 2024.09.6; Python 3.x.

**Outputs (when you run locally):** `FA_boltzmann_summaries.csv` and `run_metadata.json` (records RDKit/Python/OS, date, and parameters).

**Notes:** NPR1/NPR2 are shape ratios (I1/I3, I2/I3). L/W/T are principal-axis extents (Å). Volumes use a voxelized union-of-spheres (VDW) and solvent-accessible volume with 1.4 Å probe (SAS).
