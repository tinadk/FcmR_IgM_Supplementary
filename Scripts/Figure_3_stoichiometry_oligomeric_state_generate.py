#!/usr/bin/env python3
"""
Figure 3:
A) FcµR-Fcµ stoichiometry heatmap (left, spans both rows)
B) FcµR-J chain stoichiometry heatmap (top right)
C) State-dependent contact proportions (bottom right)
All data exported to Figure_3_data.xlsx (multi-sheet) and individual CSV files.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import OUT_DIR, DPI, MPL_RCPARAMS, PDB_DIR, PDB_INFO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from Bio.PDB import PDBParser, NeighborSearch, is_aa
from collections import defaultdict
import warnings

warnings.filterwarnings('ignore')

OUTPUT_DIR = OUT_DIR / "fig3_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(MPL_RCPARAMS)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

LABEL_FONTSIZE = 9
BASE_DX = -0.12
BASE_DY = 1.06

SUBPLOT_ADJUST = dict(left=0.10, right=0.88, top=0.94, bottom=0.06,
                      wspace=0.65, hspace=0.45)

DIST_CUT = 4.5
REPR = {"Dimer":"7YTE", "Pentamer":"7YTC", "sIgM":"7YSG"}
CHAIN_REPR = {
    "7YTE": {"Fc":["C","D"], "Ig":["A","B"]},
    "7YTC": {"Fc":["R"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
    "7YSG": {"Fc":["U","R","S","V"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
}

def fmt(res):
    return res.resname.capitalize() + str(res.id[1])

def sidechain(res):
    backbone = {'N','CA','C','O'}
    if res.resname == 'GLY':
        return [a for a in res.get_atoms() if a.name == 'CA']
    return [a for a in res.get_atoms() if a.name not in backbone]

def get_pair_proportion(structure, fc_chains, target_chains, cutoff):
    target_atoms = []
    for c in target_chains:
        if c in structure[0]:
            target_atoms.extend(list(structure[0][c].get_atoms()))
    if not target_atoms:
        return {}
    ns = NeighborSearch(target_atoms)
    total_fc = 0
    pair_cnt = defaultdict(int)
    for fc in fc_chains:
        if fc not in structure[0]:
            continue
        total_fc += 1
        seen = defaultdict(set)
        for res in structure[0][fc]:
            if not is_aa(res) or not (18 <= res.id[1] <= 124):
                continue
            sc = sidechain(res)
            if not sc:
                continue
            fc_key = fmt(res)
            for a in sc:
                for nb in ns.search(a.coord, cutoff):
                    nb_res = nb.get_parent()
                    if not is_aa(nb_res):
                        continue
                    tgt_key = fmt(nb_res)
                    seen[fc_key].add(tgt_key)
        for fc_key, tgts in seen.items():
            for tgt_key in tgts:
                pair_cnt[(fc_key, tgt_key)] += 1
    if total_fc == 0:
        return {}
    return {p: cnt / total_fc for p, cnt in pair_cnt.items()}

def get_state_data_high():
    parser = PDBParser(QUIET=True)
    state_data = {}
    for state, pdb in REPR.items():
        path = PDB_DIR / f"{pdb}.pdb"
        if not path.exists():
            print(f"Missing {path}")
            continue
        struct = parser.get_structure(pdb, str(path))
        prop = get_pair_proportion(struct, CHAIN_REPR[pdb]["Fc"], CHAIN_REPR[pdb]["Ig"], DIST_CUT)
        state_data[state] = prop

    all_pairs = set()
    for d in state_data.values():
        all_pairs.update(d.keys())

    pair_mean = {}
    for p in all_pairs:
        dimer_val = state_data["Dimer"].get(p, 0)
        pent_val  = state_data["Pentamer"].get(p, 0)
        sigm_val  = state_data["sIgM"].get(p, 0)
        pair_mean[p] = (dimer_val + pent_val + sigm_val) / 3.0

    grouped = defaultdict(list)
    for p in all_pairs:
        grouped[p[0]].append(p)
    sorted_fc_residues = sorted(grouped.keys(),
                                key=lambda fc: max(pair_mean[p] for p in grouped[fc]),
                                reverse=True)
    pairs_sorted_full = []
    for fc in sorted_fc_residues:
        group_pairs = sorted(grouped[fc], key=lambda p: pair_mean[p], reverse=True)
        pairs_sorted_full.extend(group_pairs)

    high_pairs = [p for p in pairs_sorted_full if pair_mean[p] >= 0.8]
    if len(high_pairs) < 5:
        high_pairs = pairs_sorted_full[:10]

    return state_data, high_pairs

def get_stoich_data():
    required = ["7YTC", "7YTD", "8BPE"]
    parser = PDBParser(QUIET=True)
    data_igm = {}
    for pdb in required:
        path = PDB_DIR / f"{pdb}.pdb"
        if not path.exists():
            print(f"Missing {path}")
            return None
        struct = parser.get_structure(pdb, str(path))
        info = PDB_INFO[pdb]
        prop = get_pair_proportion(struct, info["Fc"], info["Ig"], DIST_CUT)
        data_igm[pdb] = prop

    grp1, grp4, grp8 = data_igm["7YTC"], data_igm["7YTD"], data_igm["8BPE"]
    all_pairs = set(grp1) | set(grp4) | set(grp8)
    score = {p: 0.4 * grp1.get(p, 0) + 0.3 * grp4.get(p, 0) + 0.3 * grp8.get(p, 0) for p in all_pairs}
    keep = [(p, grp1.get(p, 0), grp4.get(p, 0), grp8.get(p, 0)) for p in all_pairs
            if grp1.get(p, 0) > 0 and grp4.get(p, 0) > 0 and grp8.get(p, 0) > 0]
    keep.sort(key=lambda x: score[x[0]], reverse=True)
    if not keep:
        return None
    labels = [f"{p[0][0]}-{p[0][1]}" for p in keep]
    mat = [[p[1], p[2], p[3]] for p in keep]
    df_igm = pd.DataFrame(mat, index=labels, columns=["1 FcμR", "4 FcμR", "8 FcμR"])
    return df_igm

def get_jchain_data():
    j_pdbs = ["7YTC", "7YTD", "8BPE"]
    parser = PDBParser(QUIET=True)
    data_j = {}
    for pdb in j_pdbs:
        path = PDB_DIR / f"{pdb}.pdb"
        if not path.exists():
            print(f"Missing {path}")
            return None
        struct = parser.get_structure(pdb, str(path))
        info = PDB_INFO[pdb]
        if info.get("J") is None:
            print(f"No J chain info for {pdb}")
            return None
        prop = get_pair_proportion(struct, info["Fc"], [info["J"]], DIST_CUT)
        data_j[pdb] = prop

    j1, j4, j8 = data_j["7YTC"], data_j["7YTD"], data_j["8BPE"]
    j_set = set(j1) | set(j4) | set(j8)
    if not j_set:
        return None
    j_score = {p: 0.4 * j1.get(p, 0) + 0.3 * j4.get(p, 0) + 0.3 * j8.get(p, 0) for p in j_set}
    keep_j = [(p, j1.get(p, 0), j4.get(p, 0), j8.get(p, 0)) for p in j_set
              if max(j1.get(p, 0), j4.get(p, 0), j8.get(p, 0)) > 0.3]
    keep_j.sort(key=lambda x: j_score[x[0]], reverse=True)
    if not keep_j:
        return None
    j_labels = [f"{p[0][0]}-{p[0][1]}" for p in keep_j]
    j_mat = [[p[1], p[2], p[3]] for p in keep_j]
    df_j = pd.DataFrame(j_mat, index=j_labels, columns=["1 FcμR", "4 FcμR", "8 FcμR"])
    return df_j

def main():
    state_data, high_pairs = get_state_data_high()
    df_igm = get_stoich_data()
    df_j = get_jchain_data()

    if not state_data or not high_pairs:
        print("No state data or high pairs found.")
        return

    nA = len(df_igm) if df_igm is not None else 0
    nB = len(df_j) if df_j is not None else 0

    if nA > nB and nB > 0:
        height_ratios = [nB, nA - nB]
    else:
        height_ratios = [1, 1]
        if nA <= nB and nA > 0 and nB > 0:
            print("Warning: nA <= nB, cannot make cell sizes identical while keeping positive C height. Using default ratios.")

    fig = plt.figure(figsize=(11, 7))
    fig.patch.set_facecolor('white')
    fig.subplots_adjust(**SUBPLOT_ADJUST)

    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=height_ratios)

    cmap = sns.color_palette("YlGnBu", as_cmap=True)

    # Panel A
    axA = fig.add_subplot(gs[:, 0])
    if df_igm is not None:
        sns.heatmap(df_igm, annot=False, cmap=cmap, vmin=0, vmax=1.0,
                    linewidths=0.5, cbar=False, square=False, ax=axA)
        axA.set_ylabel('Residue pairs (FcμR-D1 – Fcμ Cμ4)')
        axA.set_xlabel('Stoichiometry of FcμR in FcμR-D1–Fcμ–J complex')
        axA.tick_params(axis='y', rotation=0)
        axA.yaxis.tick_right()
        axA.yaxis.set_label_position('right')
    else:
        axA.axis('off')
        axA.text(0.5, 0.5, 'Data missing', ha='center', va='center')

    # Panel B
    axB = fig.add_subplot(gs[0, 1])
    if df_j is not None:
        sns.heatmap(df_j, annot=False, cmap=cmap, vmin=0, vmax=1.0,
                    linewidths=0.5, cbar=False, square=False, ax=axB)
        axB.set_ylabel('Residue pairs (FcμR-D1 – J)')
        axB.set_xlabel('Stoichiometry of FcμR in the FcµR–Fcµ–J complex')
        axB.tick_params(axis='y', rotation=0)
        axB.yaxis.tick_right()
        axB.yaxis.set_label_position('right')
    else:
        axB.axis('off')
        axB.text(0.5, 0.5, 'Data missing', ha='center', va='center')

    # Panel C - width increased to extend beyond panel B's right spine
    axC = fig.add_subplot(gs[1, 1])
    pos = axC.get_position()
    new_width = pos.width * 1.48
    axC.set_position([pos.x0, pos.y0, new_width, pos.height])
    axC.set_clip_on(False)

    # Get axis positions once (used for labels, colorbar, and title)
    posA = axA.get_position()
    posB = axB.get_position()
    posC = axC.get_position()

    # Draw panel C bar plot
    labelsC = [f"{p[0]}-{p[1]}" for p in high_pairs]
    x = np.arange(len(labelsC))
    w = 0.25
    dimer_vals = [state_data["Dimer"].get(p, 0) for p in high_pairs]
    pent_vals  = [state_data["Pentamer"].get(p, 0) for p in high_pairs]
    sigm_vals  = [state_data["sIgM"].get(p, 0) for p in high_pairs]
    axC.bar(x - w, dimer_vals, w, label='Dimer (7YTE)', color='#1f77b4', edgecolor='black')
    axC.bar(x, pent_vals, w, label='Pentamer (7YTC)', color='#ff7f0e', edgecolor='black')
    axC.bar(x + w, sigm_vals, w, label='sIgM (7YSG)', color='#2ca02c', edgecolor='black')
    axC.set_xticks(x)
    axC.set_xticklabels(labelsC, rotation=45, ha='right')
    axC.set_ylabel('Contact proportion')

    # ---- Use figure coordinates for x-axis title ----
    x_title = posC.x0 + posC.width / 2
    y_title = posC.y0 - 0.145   # adjust this value (more negative moves down)
    fig.text(x_title, y_title, 'Residue pairs (FcμR D1 – Fcμ Cμ4)',
             ha='center', va='bottom', fontsize=LABEL_FONTSIZE)

    axC.set_ylim(0, 1.3)
    axC.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axC.tick_params(labelleft=True)
    margin = 0.05 * (len(labelsC) - 1)
    axC.set_xlim(-0.5 - margin, len(labelsC) - 0.5 + margin)
    axC.legend(frameon=False, loc='upper right', title='IgM oligomeric state')
    axC.tick_params(axis='x', pad=8)
    for spine in axC.spines.values():
        spine.set_linewidth(0.8)
    axC.tick_params(width=0.8)

    # Panel labels using absolute coordinates (outside top-left)
    x_margin = -0.02
    y_margin = 0.02

    xA = posA.x0 + x_margin
    yA = posA.y0 + posA.height + y_margin
    xB = posB.x0 + x_margin
    yB = posB.y0 + posB.height + y_margin
    xC = posC.x0 + x_margin
    yC = posC.y0 + posC.height + y_margin

    axA.text(xA, yA, 'A', transform=fig.transFigure, fontweight='bold',
             fontsize=LABEL_FONTSIZE, va='bottom', ha='right', color='black')
    axB.text(xB, yB, 'B', transform=fig.transFigure, fontweight='bold',
             fontsize=LABEL_FONTSIZE, va='bottom', ha='right', color='black')
    axC.text(xC, yC, 'C', transform=fig.transFigure, fontweight='bold',
             fontsize=LABEL_FONTSIZE, va='bottom', ha='right', color='black')

    # Colorbar
    cbar_ax = fig.add_axes([posB.x1 + 0.13, posB.y0, 0.015, posB.height])
    norm = plt.Normalize(vmin=0, vmax=1.0)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Contact proportion')

    # Save figure
    base = OUTPUT_DIR / 'Figure_3'
    for ext in ('png', 'jpg', 'pdf'):
        fig.savefig(f"{base}.{ext}", dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print(f"Figure 3 saved to {OUTPUT_DIR}")

    # Export data
    if df_igm is not None:
        df_igm.to_csv(OUTPUT_DIR / 'Figure_3a_stoichiometry.csv')
    if df_j is not None:
        df_j.to_csv(OUTPUT_DIR / 'Figure_3b_J_chain.csv')

    df_state = pd.DataFrame({
        'pair': [f"{p[0]}-{p[1]}" for p in high_pairs],
        'Dimer': [state_data["Dimer"].get(p, 0) for p in high_pairs],
        'Pentamer': [state_data["Pentamer"].get(p, 0) for p in high_pairs],
        'sIgM': [state_data["sIgM"].get(p, 0) for p in high_pairs]
    })
    df_state.to_csv(OUTPUT_DIR / 'Figure_3c_state_dependence.csv', index=False)

    # ---------- Excel export with custom column names ----------
    excel_path = OUTPUT_DIR / 'Figure_3_data.xlsx'
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        if df_igm is not None:
            # Set index name for sheet A
            df_igm.index.name = 'FcμR-D1–IgM-Cμ4 Residue Pair'
            df_igm.to_excel(writer, sheet_name='A_FcµR-Fcµ_stoichiometry')
        if df_j is not None:
            # Set index name for sheet B
            df_j.index.name = 'FcμR-D1–J-chain Residue Pair'
            df_j.to_excel(writer, sheet_name='B_FcµR-J_chain')
        # For sheet C, rename 'pair' column
        df_state_renamed = df_state.rename(columns={'pair': 'FcμR-D1–IgM-Cμ4 Residue Pair'})
        df_state_renamed.to_excel(writer, sheet_name='C_state_dependence', index=False)

    print(f"Data exported to {OUTPUT_DIR}")
    print(f"Combined Excel workbook saved to {excel_path}")

if __name__ == "__main__":
    main()