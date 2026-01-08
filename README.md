# FA-geometry

Code and methods for comparing fatty-acid (SMILES) 3D geometry from conformer ensembles (RDKit).

## Software & Provenance
- **RDKit:** 2024.09.6
- **Python:** 3.x
- Conformers: ETKDGv3 → MMFF94; ensemble statistics are Boltzmann-weighted at 298 K.

## Reproduce (summary)
- Prepare a CSV with two columns: `FA,SMILES`.
- Run `analysis/fatty_acids/fa_volume_compare.py` to produce ensemble summaries.
- See `docs/FA-geometry-method-csv.md` for details.
