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

if __name__ == "__main__":
    main()