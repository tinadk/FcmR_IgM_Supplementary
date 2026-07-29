#!/usr/bin/env python3
"""
Fig2 : Occurrence of residue pairs across 7 structures (core/semi-core/peripheral).
Dark blue = 7/7, medium = 5‑6/7, light = ≤4/7.
Uses shared config.py for output, DPI, plot params, PDB directory, and chain info.
Exports: figure (png/jpg/pdf), CSV, Excel.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import OUT_DIR, DPI, MPL_RCPARAMS, PDB_DIR, PDB_INFO

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, NeighborSearch, is_aa
from collections import defaultdict
import warnings

warnings.filterwarnings('ignore')

# ------------------- Output folder -------------------
OUTPUT_DIR = OUT_DIR / "fig2_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ------------------- Global plot settings -------------------
plt.rcParams.update(MPL_RCPARAMS)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ------------------- Constants from config -------------------
BASE_PDB = PDB_DIR
PDBS = list(PDB_INFO.keys())                     # all 7 structures
# Build Fc/Ig chain dict (ignore J and SC)
CHAIN = {pdb: {"Fc": info["Fc"], "Ig": info["Ig"]}
         for pdb, info in PDB_INFO.items()}
DIST_CUT = 4.5

def fmt(res):
    return res.resname.capitalize() + str(res.id[1])

def sidechain(res):
    backbone = {'N', 'CA', 'C', 'O'}
    if res.resname == 'GLY':
        return [a for a in res.get_atoms() if a.name == 'CA']
    return [a for a in res.get_atoms() if a.name not in backbone]

def main():
    parser = PDBParser(QUIET=True)
    occurrence = defaultdict(int)

    for pdb in PDBS:
        fname = BASE_PDB / f"{pdb}.pdb"
        if not fname.exists():
            print(f"Missing {fname}")
            continue
        struct = parser.get_structure(pdb, str(fname))
        fc_chains = CHAIN[pdb]["Fc"]
        ig_chains = CHAIN[pdb]["Ig"]

        ig_atoms = []
        for c in ig_chains:
            if c in struct[0]:
                ig_atoms.extend(list(struct[0][c].get_atoms()))
        if not ig_atoms:
            continue

        ns = NeighborSearch(ig_atoms)
        seen = set()
        for fc in fc_chains:
            if fc not in struct[0]:
                continue
            for res in struct[0][fc]:
                if not is_aa(res) or not (18 <= res.id[1] <= 124):
                    continue
                sc = sidechain(res)
                if not sc:
                    continue
                fc_key = fmt(res)
                for a in sc:
                    for nb in ns.search(a.coord, DIST_CUT):
                        nb_res = nb.get_parent()
                        if not is_aa(nb_res):
                            continue
                        ig_key = fmt(nb_res)
                        seen.add((fc_key, ig_key))
        for pair in seen:
            occurrence[pair] += 1

    if not occurrence:
        print("No contacts found.")
        return

    def tier(cnt):
        if cnt == 7: return 0
        elif cnt >= 5: return 1
        else: return 2

    pairs_sorted = sorted(occurrence.items(),
                          key=lambda kv: (tier(kv[1]), -kv[1],
                                          int(re.search(r'\d+', kv[0][0]).group()),
                                          int(re.search(r'\d+', kv[0][1]).group())))
    labels = [f"{p[0]}-{p[1]}" for p, _ in pairs_sorted]
    counts = [cnt for _, cnt in pairs_sorted]
    colors = ['#08306b' if c == 7 else '#2171b5' if c >= 5 else '#6baed6' for c in counts]

    fig, ax = plt.subplots(figsize=(max(9, len(labels) * 0.3), 5))
    fig.patch.set_facecolor('white')
    ax.bar(range(len(labels)), counts, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('Number of complexes (out of 7) with ≥ 1 FcμR-D1 chain contact')
    ax.set_xlabel('Residue pairs (FcµR residue – Cµ4 residue)')
    ax.set_ylim(0, 7.5)
    ax.set_xlim(-0.6, len(labels) - 0.4)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
    ax.tick_params(width=0.8)
    ax.grid(False)
    plt.tight_layout()

    base_name = OUTPUT_DIR / 'Fig2_occurrence_barplot'
    for ext in ('png', 'jpg', 'pdf'):
        fig.savefig(f"{base_name}.{ext}", dpi=DPI, facecolor='white')
    plt.close(fig)
    print(f"Figures saved to {OUTPUT_DIR}")

    df_out = pd.DataFrame({'pair': labels, 'occurrence': counts})
    csv_path = OUTPUT_DIR / 'Fig2_occurrence_barplot.csv'
    df_out.to_csv(csv_path, index=False)
    print(f"CSV saved → {csv_path}")

    xlsx_path = OUTPUT_DIR / 'Fig2_occurrence_barplot.xlsx'
    try:
        df_out.to_excel(xlsx_path, index=False)
        print(f"Excel saved → {xlsx_path}")
    except ImportError:
        print("pandas/openpyxl missing; Excel export skipped.")
    except Exception as e:
        print(f"Excel export failed: {e}")

if __name__ == "__main__":
    main()