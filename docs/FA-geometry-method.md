# Ensemble-Based Size and Volume Comparison (General)

**Approach (no MD):** RDKit ETKDGv3 to generate 3D conformers; MMFF94 minimization; prune near-duplicates by heavy-atom RMSD; compute metrics per conformer; Boltzmann-average at 298 K.

**Why:** Rapid, reproducible shape/size comparisons across related fatty acids.

**Key metrics:** L/W/T, End-to-End, Rg, NPR1/NPR2; VDW/SAS volumes.

**References (PMIDs):**
- ETKDG: 26575315
- ETKDG updates: 32155061
- VMD: 8744570
- MDAnalysis: 21500218
- FreeSASA: 26973785
