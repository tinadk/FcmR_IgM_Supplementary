#!/usr/bin/env python3
"""
Figure 7 - Combined MSM analysis: macrostate populations, probability flux networks,
and MFPT matrices for all four systems.
Layout: 3 rows x 3 columns:
Row1: A (macrostate population bar chart), B (WT probability flux), C (D111A probability flux)
Row2: D (T110A_D111A probability flux), E (T60A_S63A probability flux), F (WT MFPT)
Row3: G (D111A MFPT), H (T110A_D111A MFPT), I (T60A_S63A MFPT)
Style parameters match Figure 8 where applicable.
MFPT matrices use plasma colormap; diagonal=0 for valid states, missing macrostates='light gray'.
Exports: figure (png/jpg/pdf), per-system CSV, combined CSV/Excel.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import pandas as pd
import warnings
import networkx as nx

from config import SYSTEMS, COLORS, OUT_DIR, DPI, MPL_RCPARAMS

warnings.filterwarnings('ignore')
plt.rcParams.update(MPL_RCPARAMS)

# ---------- Paths ----------
MSM_DIR = ROOT / 'MSM' / 'archive' / 'complete' / 'data' / 'msm'
TICA_DIR = ROOT / 'MSM' / 'archive' / 'complete' / 'data' / 'tica'
OUTPUT_DIR = OUT_DIR / 'fig7_output'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

LABEL_FONTSIZE = 9
PANEL_LABEL_POS = (-0.12, 1.06)
LAG_TIME = 10
FRAME_TO_NS = 0.1
STATE_COLORS = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00']

# ---------- Load common data ----------
pcca_membership = np.load(MSM_DIR / 'pcca_membership.npy')
dtraj_all = np.load(MSM_DIR / 'dtraj_k150.npy')
system_labels = np.load(TICA_DIR / 'system_labels.npy')

# Handle both 1D and 2D PCCA membership arrays
if pcca_membership.ndim == 2:
    microstate_macro = np.argmax(pcca_membership, axis=1)
elif pcca_membership.ndim == 1:
    microstate_macro = pcca_membership
else:
    raise ValueError('PCCA membership must be 1D or 2D.')

macro_assign = microstate_macro[dtraj_all]
n_macro = microstate_macro.max() + 1
macro_labels = [f'M{i+1}' for i in range(n_macro)]
sys_to_id = {name: i for i, name in enumerate(SYSTEMS)}

# ---------- Helper: build macrostate transition matrix ----------
def build_macro_P(sid):
    mask = system_labels == sid
    traj = macro_assign[mask]
    C = np.zeros((n_macro, n_macro))
    for t in range(len(traj) - 1):
        C[traj[t], traj[t+1]] += 1
    row_sum = C.sum(axis=1)
    P = np.full((n_macro, n_macro), np.nan)
    valid = row_sum > 0
    P[valid] = C[valid] / row_sum[valid, None]
    return P

# ---------- Compute macrostate populations ----------
pop_fractions = {}
for si, sys_name in enumerate(SYSTEMS):
    mask = system_labels == si
    counts = np.bincount(macro_assign[mask], minlength=n_macro)
    pop_fractions[sys_name] = counts / counts.sum()

pop_percent = {sys: pop * 100 for sys, pop in pop_fractions.items()}

# ---------- Compute MFPT matrices ----------
def compute_mfpt(P, tau=1):
    n = P.shape[0]
    mfpt = np.zeros((n, n))
    I = np.eye(n)
    for target in range(n):
        A = I - P.copy()
        A[target, :] = 0.0
        A[target, target] = 1.0
        b = np.ones(n) * tau
        b[target] = 0.0
        mfpt[:, target] = np.linalg.solve(A, b)
    return mfpt

def mfpt_for_system(sid):
    P = build_macro_P(sid)
    valid = ~np.isnan(P).all(axis=1)
    if valid.sum() == 0:
        return np.full((n_macro, n_macro), np.nan)
    P_sub = P[np.ix_(valid, valid)]
    P_sub = np.nan_to_num(P_sub, nan=0.0)
    row_sum = P_sub.sum(axis=1)
    zero_rows = row_sum == 0
    P_sub[zero_rows, zero_rows] = 1.0
    P_sub = P_sub / P_sub.sum(axis=1, keepdims=True)
    mfpt_frames = compute_mfpt(P_sub, tau=LAG_TIME)
    mfpt_ns = mfpt_frames * FRAME_TO_NS
    mfpt_full = np.full((n_macro, n_macro), np.nan)
    vi = np.where(valid)[0]
    for idx_i, si in enumerate(vi):
        for idx_j, sj in enumerate(vi):
            if si != sj:
                mfpt_full[si, sj] = mfpt_ns[idx_i, idx_j]
    # Set diagonal to 0 only for valid states; missing states remain NaN
    for idx in vi:
        mfpt_full[idx, idx] = 0.0
    return mfpt_full

mfpt_data = {}
for sys_name in SYSTEMS:
    mfpt_data[sys_name] = mfpt_for_system(sys_to_id[sys_name])

# ---------- Build flux edges for each system ----------
def get_flux_edges(sys_name):
    pi = pop_fractions[sys_name]
    P = build_macro_P(sys_to_id[sys_name])
    edges = []
    for i in range(n_macro):
        for j in range(n_macro):
            if i != j and not np.isnan(P[i, j]) and pi[i] > 1e-6 and pi[j] > 1e-6:
                flux = pi[i] * P[i, j]
                if flux > 1e-6:
                    edges.append((i, j, flux))
    return edges

flux_edges = {}
for sys_name in SYSTEMS:
    flux_edges[sys_name] = get_flux_edges(sys_name)

# ---------- Export data to a single Excel file with clear sheet names ----------
excel_path = OUTPUT_DIR / 'Figure_7_data.xlsx'
with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
    # Sheet A: macrostate populations (all systems)
    pop_df = pd.DataFrame(pop_percent, index=macro_labels).T
    pop_df.index.name = 'System'
    pop_df.to_excel(writer, sheet_name='A_Pop')

    # Sheets B-E: flux edges for each system
    flux_sheets = {
        'WT': 'B_WT_Flux',
        'D111A': 'C_D111A_Flux',
        'T110A_D111A': 'D_T110A_D111A_Flux',
        'T60A_S63A': 'E_T60A_S63A_Flux'
    }
    for sys_name, sheet_name in flux_sheets.items():
        edges = flux_edges[sys_name]
        if edges:
            df = pd.DataFrame(edges, columns=['Source', 'Target', 'Flux'])
        else:
            df = pd.DataFrame(columns=['Source', 'Target', 'Flux'])
        df.to_excel(writer, sheet_name=sheet_name, index=False)

    # Sheets F-I: MFPT matrices
    mfpt_sheets = {
        'WT': 'F_WT_MFPT',
        'D111A': 'G_D111A_MFPT',
        'T110A_D111A': 'H_T110A_D111A_MFPT',
        'T60A_S63A': 'I_T60A_S63A_MFPT'
    }
    for sys_name, sheet_name in mfpt_sheets.items():
        mat = mfpt_data[sys_name]
        df = pd.DataFrame(mat, index=macro_labels, columns=macro_labels)
        df.index.name = 'Source'
        df.to_excel(writer, sheet_name=sheet_name)

print(f'Excel data saved to {excel_path}')

# ---------- Plot 3x3 ----------
fig, axes = plt.subplots(3, 3, figsize=(18, 15))
fig.subplots_adjust(left=0.06, right=0.88, top=0.94, bottom=0.06,
                    wspace=0.25, hspace=0.30)

def add_label(ax, text):
    ax.text(PANEL_LABEL_POS[0], PANEL_LABEL_POS[1], text,
            transform=ax.transAxes, fontweight='bold', fontsize=LABEL_FONTSIZE,
            va='bottom', ha='center', color='black')

# ----- A: Macrostate population bar chart -----
axA = axes[0, 0]
x = np.arange(n_macro)
bar_width = 0.15
all_bars = []
for si, sys_name in enumerate(SYSTEMS):
    bars = axA.bar(x + si * bar_width, pop_percent[sys_name], bar_width,
                   color=COLORS[sys_name], label=sys_name)
    all_bars.append(bars)

manual_offsets = {
    (0, 0): (0.00, 0.0),
    (0, 1): (-0.02, 6.2),
    (0, 2): (-0.10, 7.0),
    (3, 1): (0.06, 0.0),
    (3, 2): (0.06, 0.0),
    (4, 0): (-0.05, 0.0),
    (4, 1): (-0.05, 0.0),
}
for state_idx in range(n_macro):
    heights = [pop_percent[sys][state_idx] for sys in SYSTEMS]
    max_idx = np.argmax(heights)
    for si, sys_name in enumerate(SYSTEMS):
        bar = all_bars[si][state_idx]
        height = heights[si]
        if height == 0:
            continue
        x_center = bar.get_x() + bar.get_width() / 2.0
        if state_idx == 0:
            dx, dy = manual_offsets.get((0, si), (0.0, 0.0))
            final_x = x_center + dx
            final_y = height + 0.3 + dy
        else:
            shift = -0.08 if si < max_idx else (0.08 if si > max_idx else 0.0)
            dx_manual, dy_manual = manual_offsets.get((state_idx, si), (0.0, 0.0))
            final_x = x_center + shift + dx_manual
            final_y = height + 0.3 + dy_manual
        axA.text(final_x, final_y, f'{height:.1f}%',
                 ha='center', va='bottom', fontsize=7)

axA.set_xticks(x + bar_width * 1.5)
axA.set_xticklabels(macro_labels)
axA.set_xlabel('Macrostate')
axA.set_ylabel('Population (%)')
axA.set_ylim(0, 65)
legA = axA.legend(loc='upper right', fontsize=LABEL_FONTSIZE, title='System', frameon=False,
                  bbox_to_anchor=(1.0, 1.02))
legA._legend_box.align = 'left'
axA.margins(x=0.05)
add_label(axA, 'A')

# ----- Probability flux networks -----
def plot_probability_flux(ax, sys_name):
    pi = pop_fractions[sys_name]
    P = build_macro_P(sys_to_id[sys_name])
    G = nx.DiGraph()
    for i in range(n_macro):
        if pi[i] > 1e-6:
            G.add_node(i, weight=pi[i])
    for i in range(n_macro):
        for j in range(n_macro):
            if i != j and not np.isnan(P[i, j]) and pi[i] > 1e-6 and pi[j] > 1e-6:
                flux = pi[i] * P[i, j]
                if flux > 1e-6:
                    G.add_edge(i, j, flux=flux)
    pos = nx.circular_layout(G)
    node_size = [max(200, G.nodes[i]['weight'] * 5000) for i in G.nodes()]
    node_color_list = [STATE_COLORS[i] for i in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=node_color_list,
                           ax=ax, edgecolors='none')
    fluxes = [G[u][v]['flux'] for u, v in G.edges()]
    if fluxes:
        max_flux = max(fluxes)
        edges = [(u, v) for u, v in G.edges() if G[u][v]['flux'] / max_flux > 0.01]
        widths = [max(0.5, 3 * G[u][v]['flux'] / max_flux) for u, v in edges]
        nx.draw_networkx_edges(G, pos, edgelist=edges, width=widths,
                               edge_color='gray', alpha=0.6, ax=ax)
    labels = {i: f'M{i+1}' for i in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, ax=ax)
    ax.set_title(sys_name, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.axis('off')

# ----- MFPT matrix -----
def plot_mfpt(ax, sys_name):
    data = mfpt_data[sys_name].copy()
    masked = np.ma.array(data, mask=np.isnan(data))
    im = ax.imshow(masked, cmap=cmap, aspect='auto', origin='lower',
                   norm=LogNorm(vmin=vmin, vmax=vmax))
    def fmt_mfpt(val):
        if abs(val) < 10000:
            return f'{val:.1f}'
        else:
            return f'{val:.1e}'
    for i in range(n_macro):
        for j in range(n_macro):
            if i == j:
                if np.isnan(data[i, j]):
                    txt = ''
                    color = 'white'
                else:
                    txt = '0'
                    color = 'white'
            elif np.isnan(data[i, j]):
                txt = ''
                color = 'white'
            else:
                val = data[i, j]
                txt = fmt_mfpt(val)
                color = 'white'
            if txt:
                ax.text(j, i, txt, ha='center', va='center', fontsize=7,
                        color=color, fontweight='bold')
    ax.set_xticks(range(n_macro))
    ax.set_xticklabels(macro_labels)
    ax.set_yticks(range(n_macro))
    ax.set_yticklabels(macro_labels)
    ax.set_xlabel('Target')
    ax.set_ylabel('Source')
    ax.set_title(sys_name, fontweight='bold')
    return im

# Precompute global vmin/vmax for MFPT matrices
all_mfpt_vals = np.concatenate([m[~np.isnan(m)] for m in mfpt_data.values()])
vmin = max(1e-3, np.nanmin(all_mfpt_vals))
vmax = np.nanmax(all_mfpt_vals)
cmap = plt.cm.plasma.copy()
cmap.set_bad('#f0f0f0')  # light gray for NaN (missing macrostates)

# Place panels
axB = axes[0, 1]; plot_probability_flux(axB, 'WT'); add_label(axB, 'B')
axC = axes[0, 2]; plot_probability_flux(axC, 'D111A'); add_label(axC, 'C')
axD = axes[1, 0]; plot_probability_flux(axD, 'T110A_D111A'); add_label(axD, 'D')
axE = axes[1, 1]; plot_probability_flux(axE, 'T60A_S63A'); add_label(axE, 'E')
axF = axes[1, 2]; imF = plot_mfpt(axF, 'WT'); add_label(axF, 'F')
axG = axes[2, 0]; imG = plot_mfpt(axG, 'D111A'); add_label(axG, 'G')
axH = axes[2, 1]; imH = plot_mfpt(axH, 'T110A_D111A'); add_label(axH, 'H')
axI = axes[2, 2]; imI = plot_mfpt(axI, 'T60A_S63A'); add_label(axI, 'I')

fig.align_labels()

# Shared colorbar for MFPT matrices
pos_i = axI.get_position()
cbar_ax = fig.add_axes([pos_i.x1 + 0.02, pos_i.y0, 0.02, pos_i.height])
cbar = fig.colorbar(imI, cax=cbar_ax, label='MFPT (ns)')

# Legend for macrostate colors
macro_handles = [Line2D([0], [0], marker='o', color='w',
                        markerfacecolor=STATE_COLORS[i], markersize=8,
                        label=macro_labels[i]) for i in range(n_macro)]
pos_c = axC.get_position()
legend_ax = fig.add_axes([pos_c.x1 + 0.02, pos_c.y1 - 0.25, 0.06, 0.25])
legend_ax.axis('off')
legend_ax.legend(handles=macro_handles, loc='upper left', frameon=False,
                 title='Macrostate', fontsize=LABEL_FONTSIZE, title_fontsize=LABEL_FONTSIZE)

# Save figure with name "Figure_7"
for ext in ('png', 'jpg', 'pdf'):
    fig.savefig(OUTPUT_DIR / f'Figure_7.{ext}', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Figure 7 saved to {OUTPUT_DIR}')