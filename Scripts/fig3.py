#!/usr/bin/env python3
"""
Fig3 : State‑dependent contact proportions (full + high‑occurrence).
Full plot: legend inside upper right.
High plot: legend outside right.
Both sorted by mean proportion descending, grouping by FcμR residue.
Outputs: Fig3_full_state_dependence, Fig3_high_state_dependence
         + corresponding CSV and Excel tables.
"""

import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, NeighborSearch, is_aa
from collections import defaultdict

# ---------- settings ----------
BASE = "/Users/tina/PycharmProjects/PythonProject/source_file"
OUT  = "/Users/tina/PycharmProjects/PythonProject/contact_analysis_manuscript2/fig3_output"
os.makedirs(OUT, exist_ok=True)

REPR = {"Dimer":"7YTE", "Pentamer":"7YTC", "sIgM":"7YSG"}
CHAIN = {
    "7YTE": {"Fc":["C","D"], "Ig":["A","B"]},
    "7YTC": {"Fc":["R"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
    "7YSG": {"Fc":["U","R","S","V"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
}
DIST = 4.5

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 7
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ---------- helper functions ----------
def fmt(res): return res.resname.capitalize() + str(res.id[1])

def sidechain(res):
    bb = {'N','CA','C','O'}
    if res.resname == 'GLY':
        return [a for a in res.get_atoms() if a.name == 'CA']
    return [a for a in res.get_atoms() if a.name not in bb]

def get_pair_proportion(structure, fc_chains, ig_chains, cutoff):
    ig_atoms = []
    for c in ig_chains:
        if c in structure[0]: ig_atoms.extend(list(structure[0][c].get_atoms()))
    if not ig_atoms: return {}
    ns = NeighborSearch(ig_atoms)
    total_fc = 0; pair_cnt = defaultdict(int)
    for fc in fc_chains:
        if fc not in structure[0]: continue
        total_fc += 1
        seen = defaultdict(set)
        for res in structure[0][fc]:
            if not is_aa(res) or not (18 <= res.id[1] <= 124): continue
            sc = sidechain(res)
            if not sc: continue
            fc_key = fmt(res)
            for a in sc:
                for nb in ns.search(a.coord, cutoff):
                    nb_res = nb.get_parent()
                    if not is_aa(nb_res): continue
                    ig_key = fmt(nb_res)
                    seen[fc_key].add(ig_key)
        for fc_key, ig_set in seen.items():
            for ig_key in ig_set:
                pair_cnt[(fc_key, ig_key)] += 1
    if total_fc == 0: return {}
    return {p: cnt/total_fc for p, cnt in pair_cnt.items()}

def main():
    parser = PDBParser(QUIET=True)
    state_data = {}
    for state, pdb in REPR.items():
        path = os.path.join(BASE, f"{pdb}.pdb")
        struct = parser.get_structure(pdb, path)
        prop = get_pair_proportion(struct, CHAIN[pdb]["Fc"], CHAIN[pdb]["Ig"], DIST)
        state_data[state] = prop

    all_pairs = set()
    for d in state_data.values(): all_pairs.update(d.keys())

    # Compute mean proportion and group by FcμR residue
    pair_mean = {}
    for p in all_pairs:
        dimer_val = state_data["Dimer"].get(p, 0)
        pent_val  = state_data["Pentamer"].get(p, 0)
        sigm_val  = state_data["sIgM"].get(p, 0)
        pair_mean[p] = (dimer_val + pent_val + sigm_val) / 3.0

    grouped = defaultdict(list)
    for p in all_pairs:
        grouped[p[0]].append(p)

    # Sort groups by the maximum mean within each group (descending)
    sorted_fc_residues = sorted(grouped.keys(),
                                key=lambda fc: max(pair_mean[p] for p in grouped[fc]),
                                reverse=True)

    # Build full sorted list
    pairs_sorted_full = []
    for fc in sorted_fc_residues:
        # within group, sort by mean descending
        group_pairs = sorted(grouped[fc], key=lambda p: pair_mean[p], reverse=True)
        pairs_sorted_full.extend(group_pairs)

    # High occurrence subset (mean ≥ 0.8)
    high_pairs = [p for p in pairs_sorted_full if pair_mean[p] >= 0.8]
    if len(high_pairs) < 5:   # fallback to top 10
        high_pairs = pairs_sorted_full[:10]

    # Helper function to draw and save
    def draw_and_save(pair_list, suffix, legend_outside=False):
        labels = [f"{p[0]}-{p[1]}" for p in pair_list]
        dimer_vals = [state_data["Dimer"].get(p, 0) for p in pair_list]
        pent_vals  = [state_data["Pentamer"].get(p, 0) for p in pair_list]
        sigm_vals  = [state_data["sIgM"].get(p, 0) for p in pair_list]

        fig, ax = plt.subplots(figsize=(max(9, len(labels)*0.25), 5))
        fig.patch.set_facecolor('white')
        x = np.arange(len(labels)); w = 0.25
        ax.bar(x - w, dimer_vals, w, label='Dimer (7YTE)', color='#1f77b4', edgecolor='black')
        ax.bar(x, pent_vals, w, label='Pentamer (7YTC)', color='#ff7f0e', edgecolor='black')
        ax.bar(x + w, sigm_vals, w, label='sIgM (7YSG)', color='#2ca02c', edgecolor='black')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Contact proportion (per FcμR chain)')
        ax.set_xlabel('Residue pair (FcμR D1 – Fcμ Cμ4)')
        ax.set_ylim(0, 1.05)
        ax.set_xlim(-0.6, len(labels) - 0.4)

        if legend_outside:
            ax.legend(frameon=False, fontsize=8, loc='upper left',
                      bbox_to_anchor=(1.01, 1.0), borderaxespad=0.)
            plt.tight_layout(rect=[0, 0, 0.83, 1])   # leave room for outside legend
        else:
            ax.legend(frameon=False, fontsize=8, loc='upper right')
            plt.tight_layout()

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.8)
        ax.tick_params(width=0.8)
        ax.grid(False)

        base = os.path.join(OUT, f'Fig3_{suffix}_state_dependence')
        for ext in ['png','jpg','pdf']:
            fig.savefig(f"{base}.{ext}", dpi=1200)
        plt.close('all')

        # Save table
        df = pd.DataFrame({'pair': labels,
                           'dimer': dimer_vals,
                           'pentamer': pent_vals,
                           'sIgM': sigm_vals,
                           'mean': [pair_mean.get(p, 0) for p in pair_list]})
        df.to_csv(os.path.join(OUT, f'supp_state_comparison_{suffix}.csv'), index=False)
        df.to_excel(os.path.join(OUT, f'supp_state_comparison_{suffix}.xlsx'), index=False)

    # Generate full figure (legend inside)
    draw_and_save(pairs_sorted_full, 'full', legend_outside=False)

    # Generate high figure (legend outside)
    draw_and_save(high_pairs, 'high', legend_outside=True)

    print("Fig3 (full + high) saved →", OUT)

if __name__ == "__main__":
    main()