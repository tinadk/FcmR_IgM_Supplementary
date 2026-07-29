#!/usr/bin/env python3
"""
Fig3 : State‑dependent contact proportions (full + high‑occurrence).
Uses shared config.py for output path, DPI, rcParams and PDB directory.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import OUT_DIR, DPI, MPL_RCPARAMS, PDB_DIR

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, NeighborSearch, is_aa
from collections import defaultdict

# ------------------- Output directory -------------------
OUTPUT_DIR = OUT_DIR / "fig3_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ------------------- Global plot settings -------------------
plt.rcParams.update(MPL_RCPARAMS)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ------------------- Constants -------------------
REPR = {"Dimer":"7YTE", "Pentamer":"7YTC", "sIgM":"7YSG"}
CHAIN = {
    "7YTE": {"Fc":["C","D"], "Ig":["A","B"]},
    "7YTC": {"Fc":["R"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
    "7YSG": {"Fc":["U","R","S","V"], "Ig":["A","B","C","D","E","F","G","H","K","L"]},
}
DIST_CUT = 4.5

# ------------------- Helper functions -------------------


# ------------------- Main -------------------

if __name__ == "__main__":
    main()