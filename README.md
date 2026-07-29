# Supplementary Material for "Hierarchical organization of FcµR–IgM recognition across distinct IgM oligomeric states"
## Citation
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21247258.svg)](https://doi.org/10.5281/zenodo.21247258)

This directory contains supplementary information, source data, analysis scripts,
and FoldX output files accompanying the manuscript.

================================================================================
FOLDER STRUCTURE
================================================================================

supplementary_information.pdf
    Main supplementary document containing supplementary figures (Figs. S1–S6),
    supplementary tables (Tables S1–S4), and detailed methods.

supplementary_tables/
    Editable Excel tables referenced in the main text and supplementary information.

source_data/
    Raw data underlying the main‑text and supplementary figures.

    - source_data_fig1a_rmsd_matrix.csv
    - source_data_fig2.csv
    - source_data_fig3_three_groups.csv
    - source_data_fig4a_stoichiometry_full.csv
    - source_data_fig4b_pentamer_j_chain_stoichiometry.csv
    - source_data_dimer_vs_pentamer_diff.csv
    - source_data_sIgM_vs_pentamer_diff.csv
    - source_data_fig7_mafft_alignment.fasta
    - source_data_fig7b_RMSF_chainC.csv
    - source_data_fig7c_COM_distance.csv
    - source_data_fig8a_tica_overlay.csv
    - source_data_fig8b_fel_2x2.csv
    - source_data_fig9a_macrostate_populations.csv
    - source_data_fig9b_tpt_flux.csv
    - source_data_fig9c_mfpt_all_systems.csv
    - source_data_figS1_RMSF_both_chains.csv

scripts/
    Python scripts to reproduce the main‑text figures. All scripts import shared
    parameters from the project‑level `config.py`.

    | Figure | Script | Output folder (`output/`) |
    |--------|--------|---------------------------|
    | Fig. 1A | `Fig1A_RMSD_heatmap.py` | `fig1a_rmsd_heatmap/` |
    | Fig. 2  | `fig2-4_contact_analysis.py` | `fig2_output/` |
    | Fig. 7B | `Fig7B_RMSF_chainC_annotated.py` | `fig7b_output/` |
    | Fig. 7C | `Fig7C_COM_distance.py` | `fig7c_output/` |
    | Fig. 8A | `Fig8A_TICA_FELs.py` | `fig8a_tica_fels_output/` |
    | Fig. 8B | `Fig8B_FELs.py` | `fig8b_fels_output/` |
    | Fig. 9A | `Fig9A_Macrostate_populations.py` | `fig9a_macrostate_populations/` |
    | Fig. 9B | `Fig9B_TPT_net_flux_networks.py` | `fig9b_tpt_flux/` |
    | Fig. 9C | `Fig9C_MFPT_matrices.py` | `fig9c_mfpt/` |
    | Fig. S1 | `FigS1_RMSF_both_chains.py` | `figS1_RMSF_both_chains/` |

    - `renumber_chains.py`: PDB chain renumbering (pre‑FoldX).
    - `fig8_workflow.mermaid`: Mermaid source for computational workflow figure.

    Each script saves high‑resolution figures (`.png`, `.jpg`, `.pdf`, 1600 dpi)
    and underlying numerical data (`.csv` / `.xlsx`) in the designated output folder.

foldx_config/
    FoldX alanine scanning raw outputs and execution commands.

    - `foldx_input_config.txt` – Run commands and parameters
    - `foldx_alanine_scanning_7yte.fxout`
    - `foldx_alanine_scanning_7ytc.fxout`
    - `foldx_alanine_scanning_7ytd.fxout`
    - `foldx_alanine_scanning_7ysg.fxout`
    - `foldx_alanine_scanning_8bpe.fxout`
    - `foldx_alanine_scanning_8bpf.fxout`
    - `foldx_alanine_scanning_8bpg.fxout`

================================================================================
DATA & CODE AVAILABILITY
================================================================================

All processed simulation trajectories, topologies, and feature matrices have been
deposited at Zenodo under DOI: 10.5281/zenodo.21247258.

The automated Snakemake workflow and all custom analysis scripts are available on
GitHub at https://github.com/tinadk/FcmR_IgM_Supplementary.git under an MIT license.

================================================================================
USAGE NOTES
================================================================================

- Python scripts require Python 3.10+ with NumPy, Pandas, Matplotlib, Seaborn,
  MDAnalysis, NetworkX, Biopython.
- FoldX analysis requires a license for FoldX 5.0. The provided `.fxout` files
  are the original program output.
- PDB structures are available from the RCSB PDB under accession codes:
  7YTE, 7YTC, 7YTD, 7YSG, 8BPE, 8BPF, 8BPG.
- All supplementary tables are provided as editable `.xlsx` files in the
  `supplementary_tables/` folder.