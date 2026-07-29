#!/usr/bin/env python3
"""
Fig4 : Stoichiometry (A), J‑chain (B), SC contacts (C) + Table1.
Uses shared config.py for paths, DPI, plot parameters, and PDB info.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import OUT_DIR, DPI, MPL_RCPARAMS, PDB_DIR, PDB_INFO

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from Bio.PDB import PDBParser, NeighborSearch, is_aa
from collections import defaultdict

# ------------------- Local paths -------------------
FIG4_OUT = OUT_DIR / "fig4_output"
FIG4_OUT.mkdir(parents=True, exist_ok=True)

# ------------------- Global plot settings -------------------
plt.rcParams.update(MPL_RCPARAMS)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ------------------- Constants -------------------
DIST_CUT = 4.5        # contact distance threshold (Å)
DEBUG = False         # set to True to print per‑structure contact info

# ------------------- Helper functions -------------------


# ------------------- Main -------------------

if __name__ == "__main__":
    main()