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


if __name__ == "__main__":
    main()