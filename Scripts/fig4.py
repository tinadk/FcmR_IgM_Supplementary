#!/usr/bin/env python3
"""
Fig4 : Stoichiometry (A), J‑chain (B), SC contacts (C) + Table1.
Optimizations:
- Heatmap rows sorted by weighted occupancy: score = 0.4*1FcμR + 0.3*4FcμR + 0.3*8FcμR
- Residue roles defined by mean contact proportion: Core anchor (≥0.85), Semi‑core (≥0.5), Peripheral
- Outputs CSV and Excel for supplementary tables
"""

import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from Bio.PDB import PDBParser, NeighborSearch, is_aa
from collections import defaultdict

# ---------- paths & settings ----------
BASE = "/Users/tina/PycharmProjects/PythonProject/source_file"
OUT  = "/Users/tina/PycharmProjects/PythonProject/contact_analysis_manuscript2/fig4_output"
os.makedirs(OUT, exist_ok=True)

PDBS = ["7YTE","7YTC","7YTD","7YSG","8BPE","8BPF","8BPG"]
CHAIN = {
    "7YTE": {"Fc":["C","D"], "Ig":["A","B"]},
    "7YTC": {"Fc":["R"], "Ig":["A","B","C","D","E","F","G","H","K","L"], "J":["J"]},
    "7YTD": {"Fc":["R","S","U","V"], "Ig":["A","B","C","D","E","F","G","H","K","L"], "J":["J"]},
    "7YSG": {"Fc":["U","R","S","V"], "Ig":["A","B","C","D","E","F","G","H","K","L"], "SC":["P"]},
    "8BPE": {"Fc":["I","M","N","O","P","Q","R","S"], "Ig":["A","B","C","D","E","F","G","H","K","L"], "J":["J"]},
    "8BPF": {"Fc":["I"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
    "8BPG": {"Fc":["A","B"], "Ig":["C","D","E","F"]}
}
DIST = 4.5

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 6
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ---------- helper functions ----------
def fmt(res):
    return res.resname.capitalize() + str(res.id[1])

def sidechain(res):
    bb = {'N','CA','C','O'}
    if res.resname == 'GLY':
        return [a for a in res.get_atoms() if a.name == 'CA']
    return [a for a in res.get_atoms() if a.name not in bb]

def get_pair_proportion(structure, fc_chains, target_chains, cutoff):
    tgt_atoms = []
    for c in target_chains:
        if c in structure[0]:
            tgt_atoms.extend(list(structure[0][c].get_atoms()))
    if not tgt_atoms:
        return {}
    ns = NeighborSearch(tgt_atoms)
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
    return {p: cnt/total_fc for p, cnt in pair_cnt.items()}

# ---------- main ----------
def main():
    parser = PDBParser(QUIET=True)
    data_igm, data_j, data_sc = {}, {}, {}
    for pdb in PDBS:
        path = os.path.join(BASE, f"{pdb}.pdb")
        if not os.path.exists(path):
            continue
        struct = parser.get_structure(pdb, path)
        entry = CHAIN[pdb]
        data_igm[pdb] = get_pair_proportion(struct, entry["Fc"], entry["Ig"], DIST)
        if "J" in entry:
            data_j[pdb] = get_pair_proportion(struct, entry["Fc"], entry["J"], DIST)
        if "SC" in entry:
            data_sc[pdb] = get_pair_proportion(struct, entry["Fc"], entry["SC"], DIST)

    # Custom colormap (light to dark blue)
    cmap = LinearSegmentedColormap.from_list('custom_blue', ['#cce5ff', '#08306b'])

    # ========== Fig4A : Stoichiometry (weighted sorting) ==========
    grp1, grp4, grp8 = data_igm["7YTC"], data_igm["7YTD"], data_igm["8BPE"]
    all_pairs = set(grp1) | set(grp4) | set(grp8)
    # Weighted score emphasises overall trend across stoichiometries
    score = {p: 0.4*grp1.get(p,0) + 0.3*grp4.get(p,0) + 0.3*grp8.get(p,0) for p in all_pairs}
    keep_igm = [(p, grp1.get(p,0), grp4.get(p,0), grp8.get(p,0)) for p in all_pairs
                if grp1.get(p,0)>0 and grp4.get(p,0)>0 and grp8.get(p,0)>0]
    keep_igm.sort(key=lambda x: score[x[0]], reverse=True)

    if keep_igm:
        labels = [f"{p[0][0]}-{p[0][1]}" for p in keep_igm]
        mat = [[p[1], p[2], p[3]] for p in keep_igm]
        df_igm = pd.DataFrame(mat, index=labels, columns=["1 FcμR","4 FcμR","8 FcμR"])
        fig, ax = plt.subplots(figsize=(4, max(3, len(keep_igm)*0.2)))
        sns.heatmap(df_igm, annot=False, cmap=cmap, vmin=0, vmax=1.0,
                    linewidths=0.5, cbar_kws={'label':'Contact proportion','aspect': 55},
                    square=False, ax=ax)
        ax.set_ylabel('Residue pair / FcμR-D1 – Fcμ-Cμ4')
        ax.set_xlabel('Number of FcμR-D1 chains')
        ax.tick_params(axis='y', rotation=0)
        plt.tight_layout()
        base = os.path.join(OUT, 'Fig4A_stoichiometry')
        for ext in ['png','jpg','pdf']:
            fig.savefig(f"{base}.{ext}", dpi=1200)
        plt.close('all')
        df_igm.to_csv(os.path.join(OUT, 'supp_stoichiometry.csv'))
        df_igm.to_excel(os.path.join(OUT, 'supp_stoichiometry.xlsx'))

    # ========== Fig4B : J‑chain (weighted sorting) ==========
    j1, j4, j8 = data_j.get("7YTC",{}), data_j.get("7YTD",{}), data_j.get("8BPE",{})
    j_set = set(j1) | set(j4) | set(j8)
    if j_set:
        j_score = {p: 0.4*j1.get(p,0) + 0.3*j4.get(p,0) + 0.3*j8.get(p,0) for p in j_set}
        keep_j = [(p, j1.get(p,0), j4.get(p,0), j8.get(p,0)) for p in j_set
                  if max(j1.get(p,0), j4.get(p,0), j8.get(p,0)) > 0.3]
        keep_j.sort(key=lambda x: j_score[x[0]], reverse=True)
        if keep_j:
            j_labels = [f"{p[0][0]}-{p[0][1]}" for p in keep_j]
            j_mat = [[p[1], p[2], p[3]] for p in keep_j]
            df_j = pd.DataFrame(j_mat, index=j_labels, columns=["1 FcμR","4 FcμR","8 FcμR"])
            fig, ax = plt.subplots(figsize=(4, max(3, len(keep_j)*0.2)*0.5))
            sns.heatmap(df_j, annot=False, cmap=cmap, vmin=0, vmax=1.0,
                        linewidths=0.5, cbar_kws={'label':'Contact proportion','aspect': 15},
                        square=False, ax=ax)
            ax.set_ylabel('Residue pair / FcμR-D1 – J')
            ax.set_xlabel('Number of FcμR-D1 chains')
            ax.tick_params(axis='y', rotation=0)
            plt.tight_layout()
            base = os.path.join(OUT, 'Fig4B_J_chain')
            for ext in ['png','jpg','pdf']:
                fig.savefig(f"{base}.{ext}", dpi=1200)
            plt.close('all')
            df_j.to_csv(os.path.join(OUT, 'supp_j_chain.csv'))
            df_j.to_excel(os.path.join(OUT, 'supp_j_chain.xlsx'))

    # ========== Fig4C : SC contacts ==========
    if "7YSG" in data_sc and data_sc["7YSG"]:
        sc = data_sc["7YSG"]
        sc_labels = [f"{fc}-{sc_}" for fc, sc_ in sc]
        sc_vals = list(sc.values())
        fig, ax = plt.subplots(figsize=(4, 0.6))
        ax.bar(range(len(sc_labels)), sc_vals, color='#2ca02c', edgecolor='black')
        ax.set_xticks(range(len(sc_labels)))
        ax.set_xticklabels(sc_labels, rotation=45, ha='right')
        ax.set_ylabel('Contact proportion')
        ax.set_xlabel('FcμR-D1 – SC pair')
        ax.set_ylim(0, 1.05)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(True)
        plt.tight_layout()
        base = os.path.join(OUT, 'Fig4C_SC_contacts')
        for ext in ['png','jpg','pdf']:
            fig.savefig(f"{base}.{ext}", dpi=1200)
        plt.close('all')
        pd.DataFrame({'pair': sc_labels, 'Contact proportion': sc_vals}).to_csv(
            os.path.join(OUT, 'supp_sc_contacts.csv'), index=False)
        pd.DataFrame({'pair': sc_labels, 'Contact proportion': sc_vals}).to_excel(
            os.path.join(OUT, 'supp_sc_contacts.xlsx'), index=False)

    # ========== Table1 : Residue classification ==========
    all_global = set()
    for d in data_igm.values():
        all_global.update(d)
    occur = defaultdict(int)
    occ_list = defaultdict(list)
    for p in all_global:
        for pdb in PDBS:
            if pdb in data_igm and p in data_igm[pdb]:
                occur[p] += 1
                occ_list[p].append(data_igm[pdb][p])
    fc_res = sorted({p[0] for p in all_global}, key=lambda x: int(re.search(r'\d+', x).group()))
    rows = []
    for fc in fc_res:
        partners = [(p[1], occur[p], np.mean(occ_list[p])) for p in all_global if p[0]==fc]
        if not partners:
            continue
        avg_cnt = np.mean([c for _,c,_ in partners])
        avg_occ = np.mean([o for _,_,o in partners])
        if avg_occ >= 0.85:
            role = 'Core anchor'
        elif avg_occ >= 0.5:
            role = 'Semi-core'
        else:
            role = 'Peripheral'
        rows.append({
            'FcμR-D1 residue': fc,
            'Mean occurrence (max 7)': round(avg_cnt, 1),
            'Mean contact proportion': round(avg_occ, 2),
            'Role': role
        })
    df_table = pd.DataFrame(rows)
    df_table.to_csv(os.path.join(OUT, 'Table1_residue_classification.csv'), index=False)
    df_table.to_excel(os.path.join(OUT, 'Table1_residue_classification.xlsx'), index=False)

    print("Fig4 & Table1 saved →", OUT)

if __name__ == "__main__":
    main()