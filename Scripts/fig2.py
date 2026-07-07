#!/usr/bin/env python3
"""
Fig2 : Occurrence of residue pairs across 7 structures (core/semi-core/peripheral).
Dark blue = 7/7, medium = 5‑6/7, light = ≤4/7.
"""

import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, NeighborSearch, is_aa
from collections import defaultdict

BASE = "/Users/tina/PycharmProjects/PythonProject/source_file"
OUT  = "/Users/tina/PycharmProjects/PythonProject/contact_analysis_manuscript2/fig2_output"
os.makedirs(OUT, exist_ok=True)

PDBS = ["7YTE","7YTC","7YTD","7YSG","8BPE","8BPF","8BPG"]
CHAIN = {
    "7YTE": {"Fc":["C","D"], "Ig":["A","B"]},
    "7YTC": {"Fc":["R"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
    "7YTD": {"Fc":["R","S","U","V"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
    "7YSG": {"Fc":["U","R","S","V"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
    "8BPE": {"Fc":["I","M","N","O","P","Q","R","S"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
    "8BPF": {"Fc":["I"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
    "8BPG": {"Fc":["A","B"], "Ig":["C","D","E","F"]}
}
DIST_CUT = 4.5

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 7
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def fmt(res): return res.resname.capitalize() + str(res.id[1])

def sidechain(res):
    backbone = {'N','CA','C','O'}
    if res.resname == 'GLY':
        return [a for a in res.get_atoms() if a.name == 'CA']
    return [a for a in res.get_atoms() if a.name not in backbone]

def main():
    parser = PDBParser(QUIET=True)
    occurrence = defaultdict(int)
    for pdb in PDBS:
        fname = os.path.join(BASE, f"{pdb}.pdb")
        if not os.path.exists(fname): continue
        struct = parser.get_structure(pdb, fname)
        fc_chains = CHAIN[pdb]["Fc"]; ig_chains = CHAIN[pdb]["Ig"]
        ig_atoms = []
        for c in ig_chains:
            if c in struct[0]: ig_atoms.extend(list(struct[0][c].get_atoms()))
        if not ig_atoms: continue
        ns = NeighborSearch(ig_atoms)
        seen = set()
        for fc in fc_chains:
            if fc not in struct[0]: continue
            for res in struct[0][fc]:
                if not is_aa(res) or not (18 <= res.id[1] <= 124): continue
                sc = sidechain(res)
                if not sc: continue
                fc_key = fmt(res)
                for a in sc:
                    for nb in ns.search(a.coord, DIST_CUT):
                        nb_res = nb.get_parent()
                        if not is_aa(nb_res): continue
                        ig_key = fmt(nb_res)
                        seen.add((fc_key, ig_key))
        for pair in seen: occurrence[pair] += 1

    if not occurrence: print("No contacts."); return

    # Sort: tier first (core=0, semi-core=1, peripheral=2), then count descending, then residue number
    def tier(cnt):
        if cnt == 7: return 0
        elif cnt >= 5: return 1
        else: return 2

    pairs_sorted = sorted(occurrence.items(),
                          key=lambda kv: (tier(kv[1]), -kv[1],
                                          int(re.search(r'\d+', kv[0][0]).group()),
                                          int(re.search(r'\d+', kv[0][1]).group())))
    labels = [f"{p[0]}-{p[1]}" for p,_ in pairs_sorted]
    counts = [cnt for _,cnt in pairs_sorted]

    # Color map
    colors = ['#08306b' if c==7 else '#2171b5' if c>=5 else '#6baed6' for c in counts]

    fig, ax = plt.subplots(figsize=(max(9, len(labels)*0.3), 5))
    fig.patch.set_facecolor('white')
    ax.bar(range(len(labels)), counts, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Number of complexes (out of 7) with ≥ 1 FcμR-D1 chain contact')
    ax.set_xlabel('Residue pair (FcμR D1 – Fcμ Cμ4)')
    ax.set_ylim(0, 7.5)
    ax.set_xlim(-0.6, len(labels)-0.4)

    # Optional horizontal lines
    # ax.axhline(y=7, color='#08306b', linestyle='--', linewidth=0.5, alpha=0.5)
    # ax.axhline(y=5, color='#2171b5', linestyle='--', linewidth=0.5, alpha=0.5)
    # ax.axhline(y=3, color='#6baed6', linestyle='--', linewidth=0.5, alpha=0.5)

    for spine in ax.spines.values(): spine.set_visible(True); spine.set_linewidth(0.8)
    ax.tick_params(width=0.8)
    ax.grid(False)
    plt.tight_layout()

    base = os.path.join(OUT, 'Fig2_occurrence_barplot')
    for ext in ['png','jpg','pdf']:
        fig.savefig(f"{base}.{ext}", dpi=1600, facecolor='white')
    plt.close('all')

    pd.DataFrame({'pair': labels, 'occurrence': counts}).to_csv(
        os.path.join(OUT, 'supp_occurrence_counts.csv'), index=False)
    print("Fig2 saved →", OUT)

if __name__ == "__main__":
    main()