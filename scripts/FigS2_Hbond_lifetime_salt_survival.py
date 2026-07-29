#!/usr/bin/env python3
"""
Figure S2: Hydrogen bond and salt bridge survival for WT and three mutants.
Uses shared config.py.  Raw ACF data from XVG files (no smoothing).
Exports:
  - figure (png, jpg, pdf)
  - lifetimes (τ) as Excel
  - combined raw survival curves (long‑format CSV)
  - separate wide‑format CSV/Excel for each interaction type
All outputs saved to output/figS2_Hbond_salt_survival/
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import DATA_DIR, OUT_DIR, SYSTEMS, COLORS, DPI, MPL_RCPARAMS

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import warnings

warnings.filterwarnings('ignore')
plt.rcParams.update(MPL_RCPARAMS)

# ------------------- Output directory -------------------
OUTPUT_DIR = OUT_DIR / "Supplementary_figS2_Hbond_salt_survival"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ------------------- Figure‑specific settings -------------------
LINE_WIDTH = 1.0
ALPHA = 1.0

YLABEL      = 'Survival probability S(t)'
SUPTITLE    = ''            # no overall title

OUTPUT_BASENAME = 'FigS2_Hbond_salt_survival'
OUTPUT_EXCEL_TAU = OUTPUT_DIR / 'Supplementary_FigS2_Hbond_lifetime_all_systems.xlsx'
OUTPUT_CSV_ALL   = OUTPUT_DIR / 'FigS2_Saltbridge_survival.csv'

OUTPUT_HBOND_CSV   = OUTPUT_DIR / 'FigS2_Hbond_lifetime.csv'
OUTPUT_HBOND_EXCEL = OUTPUT_DIR / 'Supplementary_FigS2_Hbond_lifetime_all_systems.xlsx'
OUTPUT_SALT_CSV    = OUTPUT_DIR / 'FigS2_Saltbridge_survival.csv'
OUTPUT_SALT_EXCEL  = OUTPUT_DIR / 'Supplementary_FigS2_Saltbridge_survival_all_systems.xlsx'

# ------------------- Data loading -------------------


# ------------------- Wide‑format DataFrame builder -------------------

# ------------------- Plotting -------------------

# ------------------- Long‑format CSV -------------------

# ------------------- Wide‑format exports -------------------

# ------------------- Main -------------------

if __name__ == "__main__":
    main()