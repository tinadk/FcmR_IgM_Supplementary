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


# ---------- main ----------

if __name__ == "__main__":
    main()