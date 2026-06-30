"""
Nomogram for Firth logistic regression (publication quality)
Model: N+ + BRAF V600E + TERT
Logic: standard rms::nomogram method
  points_i = (beta_i * x_i - min(beta_i * x_i)) * scale_factor
  risk = 1/(1+exp(-(intercept + min_total_contrib + total_points / scale_factor)))
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
})

# ===================== MODEL COEFFICIENTS =====================
# From Model 5: N+ + BRAF + TERT (Firth)
intercept = -3.442  # ln(0.032)
beta = {
    'N+': 1.486,           # ln(4.418)
    'BRAF V600E': -0.122,  # ln(0.885)
    'TERT': -0.059,        # ln(0.943)
}

var_names = ['N+', 'BRAF V600E', 'TERT']
var_levels = {
    'N+': [0, 1],
    'BRAF V600E': [0, 1],
    'TERT': [0, 1],
}
var_labels = {
    'N+': {0: 'No', 1: 'Yes'},
    'BRAF V600E': {0: 'Wild-type', 1: 'Mutant'},
    'TERT': {0: 'Wild-type', 1: 'Mutant'},
}
var_display = {
    'N+': 'Lymph Node Metastasis (N+)',
    'BRAF V600E': 'BRAF V600E',
    'TERT': 'TERT Promoter',
}

# Compute proper points (rms::nomogram method)
var_min_contrib = {}
var_points = {}
all_contribs = []

for vn in var_names:
    b = beta[vn]
    levels = var_levels[vn]
    contribs = [b * val for val in levels]
    min_c = min(contribs)
    var_min_contrib[vn] = min_c
    all_contribs.extend(contribs)

# Scale: 100 points = max range across all variables
max_range = max(
    abs(b * (max(var_levels[vn]) - min(var_levels[vn])))
    for vn, b in beta.items()
)
scale = 100 / max_range if max_range > 0 else 1

for vn in var_names:
    b = beta[vn]
    levels = var_levels[vn]
    pts = [(b * val - var_min_contrib[vn]) * scale for val in levels]
    var_points[vn] = pts

min_total_contrib = sum(var_min_contrib.values())

# ===================== NOMOGRAM FIGURE =====================
fig, ax = plt.subplots(figsize=(10, 5.5))
ax.set_xlim(-2, 102)
ax.set_ylim(-2.2, len(var_names) * 2.5 + 2.5)
ax.axis('off')

# --- TITLE ---
ax.text(50, len(var_names) * 2.5 + 2.0, 'Nomogram for Recurrence Risk',
        fontsize=14, fontweight='bold', ha='center', va='center')
ax.text(50, len(var_names) * 2.5 + 1.3, 'Model: N+ + BRAF V600E + TERT (Firth Logistic)',
        fontsize=9, ha='center', va='center', style='italic', color='gray')

# --- POINTS RULER (top) ---
y_ruler = len(var_names) * 2.5
ax.text(-1.5, y_ruler, 'Points', fontsize=10, fontweight='bold', va='center', ha='left')
for pt_val in range(0, 101, 10):
    x = pt_val
    ax.plot([x, x], [y_ruler - 0.15, y_ruler + 0.15], color='gray', lw=0.8)
    ax.text(x, y_ruler - 0.5, str(pt_val), fontsize=7.5, ha='center', va='top', color='dimgray')
# Full ruler line
ax.plot([0, 100], [y_ruler, y_ruler], color='black', lw=1.5)

# --- VARIABLE ROWS ---
for idx, vn in enumerate(var_names):
    y_center = (len(var_names) - 1 - idx) * 2.2
    
    # Variable label
    ax.text(-1.5, y_center, var_display[vn], fontsize=10, fontweight='bold',
            va='center', ha='left')
    
    # Horizontal line
    ax.plot([0, 100], [y_center, y_center], color='gray', lw=0.8, alpha=0.4)
    
    # Points at each level
    pts = var_points[vn]
    levels = var_levels[vn]
    labels = var_labels[vn]
    
    for pi, (lv, pt_v) in enumerate(zip(levels, pts)):
        x_pos = pt_v  # 0-100 scale
        label = labels[lv]
        # Dot marker
        ax.plot(x_pos, y_center, 'o', color='#2c3e50', markersize=8, zorder=5)
        # Label text
        ax.text(x_pos, y_center - 0.6, label, fontsize=8.5, ha='center', va='top',
                fontweight='bold')
        # Point value
        ax.text(x_pos + 3, y_center + 0.4, f'{pt_v:.0f}', fontsize=7,
                va='center', color='#7f8c8d')

# --- GRID LINES (vertical alignment guides) ---
for pt_val in range(0, 101, 10):
    x = pt_val
    y_start = len(var_names) * 2.5 - 0.15
    y_end = -1.8
    ax.plot([x, x], [y_start, y_end], color='lightgray', lw=0.3, alpha=0.5, zorder=0)

# --- TOTAL POINTS RULER ---
y_tp = -0.8
ax.text(-1.5, y_tp, 'Total Points', fontsize=10, fontweight='bold', va='center', ha='left')
# Compute total points range
max_total_pts = sum(max(var_points[vn]) for vn in var_names)
# Draw marks at 0, 20, 40, ... max_total_pts
tp_values = np.linspace(0, max_total_pts, 11)
for tp in tp_values:
    x_pct = tp / max_total_pts * 100 if max_total_pts > 0 else 0
    ax.plot([x_pct, x_pct], [y_tp - 0.1, y_tp + 0.1], color='gray', lw=0.8)
    ax.text(x_pct, y_tp - 0.4, f'{tp:.0f}', fontsize=7.5, ha='center', va='top', color='dimgray')
ax.plot([0, 100], [y_tp, y_tp], color='black', lw=1.5)

# --- RISK RULER ---
y_risk = -1.6
ax.text(-1.5, y_risk, 'Risk', fontsize=10, fontweight='bold', va='center', ha='left')
ax.text(-1.5, y_risk - 0.5, 'of Recurrence', fontsize=7, va='center', ha='left', color='gray')

# Risk values
risk_values = []
risk_positions = []
for tp in tp_values:
    lp = intercept + min_total_contrib + tp / scale
    risk = 1 / (1 + np.exp(-lp))
    risk_values.append(risk)
    x_pct = tp / max_total_pts * 100 if max_total_pts > 0 else 0
    risk_positions.append(x_pct)

for i, (tp, risk, x_pct) in enumerate(zip(tp_values, risk_values, risk_positions)):
    risk_pct = risk * 100
    if risk_pct < 1:
        label = f'{risk_pct:.1f}%'
    elif risk_pct < 10:
        label = f'{risk_pct:.1f}%'
    else:
        label = f'{risk_pct:.0f}%'
    ax.text(x_pct, y_risk - 0.3, label, fontsize=7.5, ha='center', va='top', color='#c0392b',
            fontweight='bold')
    ax.plot([x_pct, x_pct], [y_risk + 0.1, y_risk - 0.1], color='#c0392b', lw=0.8)

# Dashed connection from total points to risk
for tp_val in tp_values[::2]:
    x_pct = tp_val / max_total_pts * 100 if max_total_pts > 0 else 0
    ax.plot([x_pct, x_pct], [y_tp + 0.1, y_risk - 0.1], color='gray', lw=0.3, ls=':', alpha=0.5)

# --- HOW TO USE ---
ax.text(50, y_risk - 1.0,
        'Instructions: For each variable, find the point value. Sum all points. Read risk at the bottom.',
        fontsize=7, ha='center', va='top', style='italic', color='gray')

# --- EXAMPLE ANNOTATION ---
# N+, BRAF WT, TERT WT = 100 pts
ex_tp = (var_points['N+'][1] + var_points['BRAF V600E'][0] + var_points['TERT'][0])
ex_x = ex_tp / max_total_pts * 100 if max_total_pts > 0 else 0
ex_risk = 1 / (1 + np.exp(-(intercept + min_total_contrib + ex_tp / scale)))
ax.annotate(f'Example: N+ = 100 pts\nRisk = {ex_risk*100:.1f}%',
            xy=(ex_x, y_risk + 0.3), xytext=(85, y_risk + 1.5),
            fontsize=7.5, ha='center',
            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=0.8),
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#fef9e7', edgecolor='gray', alpha=0.8))

plt.tight_layout()
fig.savefig('H:\\Q3 MNGOC\\nomogram.png', dpi=300)
plt.close(fig)
print(f'  Saved: nomogram.png')
print(f'  Scale factor: {scale:.2f}')
print(f'  Max total points: {max_total_pts:.0f}')
print(f'  Min total contrib: {min_total_contrib:.3f}')
print(f'  Intercept: {intercept:.3f}')

# ===================== VERIFICATION =====================
print('\nVerification:')
test_cases = [
    ('N0, BRAF WT, TERT WT',     [0, 0, 0]),
    ('N+, BRAF WT, TERT WT',      [1, 0, 0]),
    ('N+, BRAF V600E, TERT WT',   [1, 1, 0]),
    ('N0, BRAF V600E, TERT WT',   [0, 1, 0]),
    ('N+, BRAF WT, TERT mutant',  [1, 0, 1]),
    ('N+, BRAF V600E, TERT mut',  [1, 1, 1]),
    ('N0, BRAF WT, TERT mutant',  [0, 0, 1]),
]

for label, vals in test_cases:
    # Actual risk from logistic regression
    lp_actual = intercept + sum(beta[vn] * v for vn, v in zip(var_names, vals))
    risk_actual = 1 / (1 + np.exp(-lp_actual))
    # Nomogram total points
    tp = 0
    for vn, v in zip(var_names, vals):
        pts = (beta[vn] * v - var_min_contrib[vn]) * scale
        tp += pts
    lp_nomo = intercept + min_total_contrib + tp / scale
    risk_nomo = 1 / (1 + np.exp(-lp_nomo))
    match = '[OK]' if abs(risk_actual - risk_nomo) < 0.001 else '[FAIL]'
    print(f'  {label:35s}: actual={risk_actual*100:5.1f}%  nomogram={risk_nomo*100:5.1f}%  pts={tp:5.0f}  {match}')

print('\nDone.')
