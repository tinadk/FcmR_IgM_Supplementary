# Hierarchical organization of FcµR–IgM recognition across distinct IgM oligomeric states

This repository contains the complete set of Python scripts, core modules,
and configuration files required to reproduce all main‑text and
supplementary figures and tables in the associated manuscript.

## Repository structure
.
├── config.py                          # Global configuration (paths, colours, DPI, matplotlib params)
├── core/                              # Shared computational modules
│   ├── __init__.py
│   ├── compute.py                     # Trajectory trimming, equilibration
│   └── stats.py                       # Bootstrap CI, etc.
├── scripts/                           # Analysis and plotting scripts
│   ├── FoldX_input_config.txt         # FoldX alanine scanning run commands
│   ├── Fig1A_RMSD_heatmap.py
│   ├── Fig2_occurrence_barplot.py
│   ├── Fig3_high_state_dependence.py
│   ├── Fig4_stoichiometry.py
│   ├── Fig7A_stability_RMSD_Rg_SASA_combine.py
│   ├── Fig7B_RMSF_chainC.py
│   ├── Fig7C_COM_distance.py
│   ├── Fig8A_TICA_FELs.py
│   ├── Fig8B_FELs.py
│   ├── Fig9A_equilibrium_macrostate_populations.py
│   ├── Fig9B_TPT_flux_networks.py
│   ├── Fig9C_MFPT.py
│   ├── FigS1_RMSF_both_chains.py
│   ├── FigS2_Hbond_lifetime_salt_survival.py
│   ├── FigS3_DSSP_secondary_structure.py
│   ├── FigS4_PCA_FELs.py
│   ├── FigS5_MSM_validation.py
│   ├── TableS1_global_MD_metrics.py
│   └── TableS2_MSM_validation.py
├── main_figures/                      # High‑resolution main‑text figures (1600 dpi)
│   ├── Fig1A.jpg
│   ├── Fig1B.jpg
│   ├── Fig2.jpg
│   ├── Fig3.jpg
│   ├── Fig4A.jpg
│   ├── Fig4B.jpg
│   ├── Fig5A.jpg
│   ├── Fig5B.jpg
│   ├── Fig6A.jpg
│   ├── Fig6B.jpg
│   ├── Fig6C.jpg
│   ├── Fig7A.jpg
│   ├── Fig7B.jpg
│   ├── Fig7C.jpg
│   ├── Fig8A.jpg
│   ├── Fig8B.jpg
│   ├── Fig9A.jpg
│   ├── Fig9B.jpg
│   ├── Fig9C.jpg
│   ├── Fig10A.jpg
│   ├── Fig10B.jpg
│   ├── Fig10C.jpg
│   ├── Fig10D.jpg
│   ├── Fig10E.jpg
│   ├── Fig10F.jpg
│   └── Fig10G.jpg
├── README.md
├── LICENSE
└── .gitignore                         # Ignore output/, *.xvg, *.npy, *.pdb, etc.


All scripts import shared settings from `config.py` and use functions
defined in the `core` module. Running any script from the project root
creates publication‑quality figures and numerical data tables inside a
local `output/` folder (not tracked by Git).

## What is NOT in this repository

| Item | Where to find it |
|------|------------------|
| MD simulation data (.xvg, .npy, .dat) | Zenodo (DOI: https://doi.org/10.5281/zenodo.21693996) |
| PDB structure files | Zenodo or RCSB PDB (7YTE, 7YSG, 7YTC, 7YTD, 8BPE, 8BPF, 8BPG) |
| FoldX alanine scanning outputs (.fxout) | Supplementary material package (submitted with manuscript) |
| Source data (.csv) for all figures | Supplementary material package |
| Supplementary tables (.xlsx) | Supplementary material package |
| Supplementary figures (.jpg) | Supplementary material package |
| Supplementary information PDF | Supplementary material package |

## Reproducing the analyses

### 1. Clone the repository
git clone https://github.com/tinadk/FcmR_IgM_Supplementary
cd FcmR_IgM_Supplementary

### 2. Install dependencies
Python ≥ 3.10 is required. Install the necessary packages:
pip install numpy pandas matplotlib seaborn mdanalysis networkx biopython scipy


### 3. Obtain the data
Download the Zenodo dataset (DOI: 10.5281/zenodo.21693996) and extract it.
The dataset includes:
- `data/` – MD analysis files (.xvg, .dat)
- `MSM_final/archive/complete/data/` – MSM input files (.npy)
- `source_file/` – PDB structures

### 4. Configure paths
Edit `config.py` and set:
- `PDB_DIR` → where you extracted `source_file/`
- `DATA_DIR` → where you extracted `data/`
- `OUT_DIR`  → where outputs should be saved

All other paths used by the scripts (e.g., MSM data, TICA data) are
defined relative to the project root and normally do not need to be changed.

### 5. Run a script
python scripts/Fig9C_MFPT.py

Generated figures (`.png`, `.jpg`, `.pdf`) and data tables (`.csv`,
`.xlsx`) are saved in the `OUT_DIR` subfolder.

### Script‑to‑figure correspondence
| Figure / Table | Script                                        |
|----------------|-----------------------------------------------|
| Fig.1A         | `Fig1A_RMSD_heatmap.py`                       |
| Fig.2          | `Fig2_occurrence_barplot.py`                  |
| Fig.3          | `Fig3_high_state_dependence.py`               |
| Fig.4A/B       | `Fig4_stoichiometry.py`                       |
| Fig.7A         | `Fig7A_stability_RMSD_Rg_SASA_combine.py`     |
| Fig.7B         | `Fig7B_RMSF_chainC.py`                        |
| Fig.7C         | `Fig7C_COM_distance.py`                       |
| Fig.8A         | `Fig8A_TICA_FELs.py`                          |
| Fig.8B         | `Fig8B_FELs.py`                               |
| Fig.9A         | `Fig9A_equilibrium_macrostate_populations.py` |
| Fig.9B         | `Fig9B_TPT_flux_networks.py`                  |
| Fig.9C         | `Fig9C_MFPT.py`                               |
| Fig.S1         | `FigS1_RMSF_both_chains.py`                   |
| Fig.S2         | `FigS2_Hbond_lifetime_salt_survival.py`       |
| Fig.S3         | `FigS3_DSSP_secondary_structure.py`           |
| Fig.S4         | `FigS4_PCA_FELs.py`                           |
| Fig.S5         | `FigS5_MSM_validation.py`  
| Fig.S6         | `FigS6_Integrated_workflow.mermaid 
| Table S2       | `TableS2_MSM_validation.py`                   |

Tables S3 and S4 are exported automatically by the Fig.9A and Fig.9C
scripts, respectively.

## Data availability
- All processed simulation trajectories, topologies, and feature matrices
  are deposited at Zenodo under DOI:https://doi.org/10.5281/zenodo.21693996).
- PDB structures are available from the RCSB Protein Data Bank (accession
  codes: 7YTE, 7YSG, 7YTC, 7YTD, 8BPE, 8BPF, 8BPG).


## Supplementary material package
The supplementary material submitted to the journal (including the
supplementary information PDF, supplementary figures, supplementary tables,
source data, and FoldX output files) is archived separately and is not part
of this code repository. The scripts included here are sufficient to
regenerate all results from the raw data.

## License
This project is licensed under the MIT License.
