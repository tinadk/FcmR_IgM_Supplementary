#!/usr/bin/env python3
"""
Global configuration for MD simulation analysis.
All paths are relative to this file's location for portability.
"""

from pathlib import Path

# ---------- Paths ----------
BASE_DIR = Path(__file__).resolve().parent

DATA_DIR = BASE_DIR / "data"
OUT_DIR = BASE_DIR / "output"
OUT_DIR.mkdir(exist_ok=True)


# ---------- Systems & Colors ----------
SYSTEMS = ['WT', 'D111A', 'T110A_D111A', 'T60A_S63A']

COLORS = {
    'WT': '#2C3E50',
    'D111A': '#E74C3C',
    'T110A_D111A': '#27AE60',
    'T60A_S63A': '#2980B9'
}

# ---------- Plotting parameters ----------
DPI = 1400
FONT_FAMILY = 'Arial'
FONT_SIZE = 9

MPL_RCPARAMS = {
    'font.family': FONT_FAMILY,
    'font.size': FONT_SIZE,
    'axes.labelsize': FONT_SIZE,
    'axes.titlesize': FONT_SIZE + 1,
    'xtick.labelsize': FONT_SIZE - 1,
    'ytick.labelsize': FONT_SIZE - 1,
    'legend.fontsize': FONT_SIZE - 1,
    'figure.titlesize': FONT_SIZE + 2,
    'savefig.dpi': DPI,
    'savefig.facecolor': 'white',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'grid.alpha': 0.15,
}

# ---------- Simulation parameters ----------
EQUILIBRATION_NS = 10.0
TIME_MAX_NS = 200.0

# ---------- External data directories ----------
PDB_DIR = Path("/Users/tina/PycharmProjects/PythonProject/source_file")

# ---------- MSM & TICA data directories (simplified Zenodo structure) ----------
MSM_DIR = BASE_DIR / "msm"
TICA_DIR = BASE_DIR / "tica"

# ---------- PDB structure assignments ----------
# Each entry provides:
#   type: "dimer", "pentamer", or "sIgM"
#   Fc: list of chain IDs for Fc region
#   Ig: list of chain IDs for Ig domain (Cμ4)
#   J: J-chain chain ID (None if absent)
# Please verify chain IDs against the actual PDB headers; J-chain often appears as "J" or "K".
PDB_INFO = {
    # Dimers (no J-chain, no SC)
    "7YTE": {"type": "dimer",    "Fc": ["C","D"],                     "Ig": ["A","B"],                          "J": None, "SC": None},
    "8BPG": {"type": "dimer",    "Fc": ["A","B"],                     "Ig": ["C","D","E","F"],                  "J": None, "SC": None},
    # Pentamers (J-chain present, no SC)
    "7YTC": {"type": "pentamer", "Fc": ["R"],                         "Ig": ["A","B","C","D","E","F","G","H","K","L"], "J": "J", "SC": None},
    "7YTD": {"type": "pentamer", "Fc": ["R","S","U","V"],             "Ig": ["A","B","C","D","E","F","G","H","K","L"], "J": "J", "SC": None},
    "8BPE": {"type": "pentamer", "Fc": ["I","M","N","O","P","Q","R","S"], "Ig": ["A","B","C","D","E","F","G","H","K","L"], "J": "J", "SC": None},
    "8BPF": {"type": "pentamer", "Fc": ["I"],                         "Ig": ["A","B","C","D","E","F","G","H","K","L"], "J": "J", "SC": None},
    # sIgM (J-chain + SC present)
    "7YSG": {"type": "sIgM",     "Fc": ["U","R","S","V"],             "Ig": ["A","B","C","D","E","F","G","H","K","L"], "J": "J", "SC": "P"},
}