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
├── Scripts/                           # Analysis and plotting scripts
│   ├── Figure_1A_pairwise_Cα_RMSD_heatmap_generate.py
│   ├── Figure_2_occurrence_barplot_generate.py
│   ├── Figure_3_stoichiometry_oligomeric_state_generate.py
│   ├── Figure_5_dynamic_global_structure_stability_generate.py
│   ├── Figure_6_TICA_and_relative_free_energies_generate.py
│   ├── Figure_7_MSM_and_TPT_generate.py
│   ├── Figure_S1_RMSF_both_chains_generate.py
│   ├── Figure_S2_Hbond_lifetime_salt_survival_generate.py
│   ├── Figure_S3_DSSP_secondary_structure_generate.py
│   ├── Figure_S4_PCA_FELs_generate.py
│   ├── Figure_S5_MSM_validation_generate.py
│   ├── Figure_S6_Integrated_workflow_generate.mermaid
│   ├── Table_S1_global_MD_metrics_generate.py
│   └── Table_S2_MSM_validation_generate.py
├── FoldX_Config/                      # FoldX alanine scanning commands and outputs
│   ├── FoldX_input_config.txt         # Run commands for all seven PDB structures
│   └── Outputs/
│       ├── FoldX_alanine_scanning_7YSG.fxout
│       ├── FoldX_alanine_scanning_7YTE.fxout
│       ├── FoldX_alanine_scanning_7YTC.fxout
│       ├── FoldX_alanine_scanning_7YTD.fxout
│       ├── FoldX_alanine_scanning_8BPE.fxout
│       ├── FoldX_alanine_scanning_8BPF.fxout
│       └── FoldX_alanine_scanning_8BPG.fxout
├── Source_Data/                       # Source data (CSV, FASTA) underlying every figure
│   ├── FFigure_1A_RMSD_matrix_data.xlsx
│   ├── Figure_2_occurrence_barplot_data.xlsx
│   ├── Figure_3_stoichiometry_oligomeric_state_data.xlsx
│   ├── Figure_4A_4B_MAFFT_alignment.fasta
│   ├── Figure_5_dynamic_global_structure_stability_data.xlsx
│   ├── Figure_6_TICA_and_relative_free_energies_data.xlsx
│   ├── Figure_7_MSM_and_TPT_data.xlsx
│   ├── Figure_S1_RMSF_both_chains_all_residues_data.xlsx
│   ├── Figure_S2_saltbridge_survival_hbond_lifetime_data.xlsx
│   ├── Figure_S3_DSSP_data.xlsx
│   ├── Figure_S4_PCA_FEL_data.xlsx
│   └── Figure_S5_MSM_validation.xlsx
├── Supplementary_Materials_for_review/           # Supplementary material package (also submitted to journal)
│   ├── Supplementary_Figures/
│   │   ├── Figure_S1.pdf         # RMSF onto both FcµR-D1 chain and IgM-Cµ4 chain
│   │   ├── Figure_S2.pdf         # Hbond and salt survival
│   │   ├── Figure_S3.pdf         # DSSP heatmap
│   │   ├── Figure_S4.pdf         # PCA and FEL heatmap
│   │   ├── Figure_S5.pdf         # MSM validation
│   │   └── Figure_S6.pdf         # Integrated workflow
│   ├── Supplementary_Tables.xlsx           # Include S1/S2/S3/S4 
│   └── Supplementary_Figure_Legends.docx      
├── Main_Figures/                      # High‑resolution main‑text figures (1400 dpi)
│   ├── Figure_1.pdf
│   ├── Figure_2.pdf
│   ├── Figure_3.pdf
│   ├── Figure_4.pdf
│   ├── Figure_5.pdf
│   ├── Figure_6.pdf
│   ├── Figure_7.pdf
│   └── Figure_8.pdf
├── README.md
├── LICENSE
└── .gitignore                         # Ignore output/, *.xvg, *.npy, *.pdb, etc.

All scripts import shared settings from `config.py` and use functions
defined in the `core` module. Running any script from the project root
creates publication‑quality figures and numerical data tables inside a
local `output/` folder (not tracked by Git).

## What is NOT in this repository

| Item                                    | Where to find it                                              |
|-----------------------------------------|---------------------------------------------------------------|
| MD simulation data (.xvg, .npy, .dat)   | Zenodo (DOI: https://doi.org/10.5281/zenodo.22132719)         |
| PDB structure files                     | Zenodo or RCSB PDB (7YTE, 7YSG, 7YTC, 7YTD, 8BPE, 8BPF, 8BPG) | 
| FoldX alanine scanning outputs (.fxout) | FoldX_Config                                                  |
| Source data (.xlsx) for all figures     | Source_Data                                                   |
| Supplementary Tables (.xlsx)            | Supplementary_Materials_for_review                            |
| Supplementary Figures (.pdf)            | Supplementary_Materials_for_review/Supplementary_Figures      |
| Supplementary Figure Legends (.docx)    | Supplementary_Materials_for_review                            |

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
| Figure / Table  | Script                                                    |
|-----------------|-----------------------------------------------------------|
| Figure_1.pdf    | `Figure_1A_pairwise_Cα_RMSD_heatmap_generate.py`          |
| Figure_2.pdf    | `Figure_2_occurrence_barplot_generate.py`                 |
| Figure_3.pdf    | `Figure_3_stoichiometry_oligomeric_state_generate.py`     |
| Figure_4A/B.pdf | `Figure_4A_4B_MAFFT_alignment.fasta`                      |
| Figure_5.pdf    | `Figure_5_dynamic_global_structure_stability_generate.py` |
| Figure_6.pdf    | `Figure_6_TICA_and_relative_free_energies_generate.py     |
| Figure_7.pdf    | `Figure_7_MSM_and_TPT_generate.py`                        |
| Figure_S1.pdf   | `Figure_S1_RMSF_both_chains_generate.py`                  |
| Figure_S2.pdf   | `Figure_S2_Hbond_lifetime_salt_survival_generate.py`      |
| Figure_S3.pdf   | `Figure_S3_DSSP_secondary_structure_generate.py`          |
| Figure_S4.pdf   | `Figure_S4_PCA_FELs_generate.py`                          |
| Figure_S5.pdf   | `Figure_S5_MSM_validation_generate.py`                    |
| Figure_S6.pdf   | `Figure_S6_Integrated_workflow_generate.mermaid`          |
| Table_S1.xlsx   | `Table_S1_global_MD_metrics_generate.py`                  |
| Table_S2.xlsx   | `Table_S2_MSM_validation_generate.py`                     |

Tables S3 and S4 are exported automatically by the Figure 7A and Figure 7F-I
scripts, respectively.

## Data availability
- All processed simulation trajectories, topologies, and feature matrices
  are deposited at Zenodo under DOI:https://doi.org/10.5281/zenodo.22132719).
- PDB structures are available from the RCSB Protein Data Bank (accession
  codes: 7YTE, 7YSG, 7YTC, 7YTD, 8BPE, 8BPF, 8BPG).
- FoldX alanine scanning outputs and run commands, custom analysis scripts,
  configuration files, source data, and the complete supplementary material
  package are all contained within this GitHub repository.

## Supplementary material package
The supplementary material submitted to the journal (including the
supplementary information PDF, supplementary figures, supplementary tables) is archived separately and is not part
of this code repository. The scripts included here are sufficient to regenerate all results from the raw data.

## License
This project is licensed under the MIT License.
