import pandas as pd
import numpy as np
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.ticker as mticker

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans'],
    'font.size': 10,
    'figure.dpi': 150,
})

# ===================== LOAD DATA =====================
f1 = 'H:\\Q3 MNGOC\\Bản sao của DỮ LIỆU Q3 GỐC MINH NGỌC MỚI.xlsx'
f2 = 'H:\\Q3 MNGOC\\KẾT QUẢ GEN Q3 MINH NGỌC 2.xlsx'
df1 = pd.read_excel(f1)
df2 = pd.read_excel(f2)
df = df1.merge(df2, on='Mã bệnh nhân', how='left')

us = df['Bất thường siêu âm cụ thể'].fillna('').str.contains('Hạch', na=False)
fna = df['Chọc tế bào chẩn đoán'].fillna('').str.contains('Đã chọc \\+', na=False)
reop = df['Xử trí'].fillna('').str.contains('Phẫu thuật', na=False)
fu = df['Bất thường trong lần khám gần nhất sau xử trí'].fillna('').str.contains('Có bất thường', na=False)
recur = (us | fna | reop | fu).astype(int)

def encode_t(val):
    if pd.isna(val): return np.nan
    v = str(val).strip().lower()
    if v in ['x','0']: return 0
    if v in ['1','1a','1b']: return 1
    if v in ['2']: return 2
    if v in ['3','3a','3b']: return 3
    if v in ['4','4a','4b','4c']: return 4
    return np.nan

def encode_n(val):
    if pd.isna(val): return np.nan
    v = str(val).strip().lower()
    if v in ['0']: return 0
    if v in ['1', '1a']: return 1
    if v in ['1b']: return 2
    if v in ['2']: return 3
    return np.nan

t_stage = df['Phân loại T'].apply(encode_t)
n_stage = df['Phân loại N'].apply(encode_n)
braf = (df['BRAFV600E'] == 2).astype(int)
tert = (df['TERT'] == 2).astype(int)
tp53 = (df['TP53'] == 2).astype(int)
nhan = pd.to_numeric(df['Tổng số nhân trc pt'], errors='coerce').fillna(0)
da_o = (nhan > 1).astype(int)

n_plus = (n_stage > 0).astype(int)
t3_plus = (t_stage >= 3).astype(int)

# ===================== FIRTH LOGIT =====================
class FirthLogit:
    def fit(self, X, y, maxiter=200, tol=1e-8):
        n, p = X.shape
        beta = np.zeros(p)
        for it in range(maxiter):
            eta = X @ beta
            pi = np.clip(1/(1+np.exp(-eta)), 1e-12, 1-1e-12)
            W = np.diag(pi*(1-pi))
            I = X.T @ W @ X
            U = X.T @ (y-pi)
            sqrtW = np.sqrt(pi*(1-pi))
            X_tilde = sqrtW.reshape(-1,1) * X
            try:
                H = X_tilde @ np.linalg.solve(X_tilde.T @ X_tilde, X_tilde.T)
                h = np.diag(H)
            except:
                h = np.zeros(n)
            h = np.clip(h, 0, 1)
            U_star = U + 0.5 * X.T @ (h*(1-2*pi))
            try:
                I_inv = np.linalg.solve(I, np.eye(p))
            except:
                I_inv = np.linalg.pinv(I)
            beta_new = beta + I_inv @ U_star
            if np.max(np.abs(beta_new-beta)) < tol:
                break
            beta = beta_new
        self.beta = beta
        self.se = np.sqrt(np.diag(I_inv))
        self.z = self.beta / self.se
        self.pvalues = 2*(1-stats.norm.cdf(np.abs(self.z)))
        ci = np.column_stack([self.beta-1.96*self.se, self.beta+1.96*self.se])
        self.or_val = np.exp(self.beta[1]) if p>1 else np.exp(self.beta[0])
        self.or_ci = np.exp(ci[1]) if p>1 else np.exp(ci[0])
        return self

def firth_or(var, outcome):
    valid = var.notna() & outcome.notna()
    y = outcome[valid].astype(float).values
    x = var[valid].astype(float).values
    n = len(y); ne = int(y.sum())
    if ne < 2 or (n-ne) < 2:
        return None, None, None, None, n, ne
    X = np.column_stack([np.ones(n), x])
    try:
        m = FirthLogit().fit(X, y)
        return m.or_val, m.or_ci[0], m.or_ci[1], m.pvalues[1], n, ne
    except:
        return None, None, None, None, n, ne

# ===================== FOREST PLOT DATA =====================
subgroups = []

# All
or_v, ci_l, ci_h, pv, n, ne = firth_or(braf, recur)
if or_v is not None:
    subgroups.append(('All patients', f'{n} / {ne}', or_v, ci_l, ci_h, pv))

# N+
mask_np = (n_plus == 1)
or_v, ci_l, ci_h, pv, n, ne = firth_or(braf[mask_np], recur[mask_np])
if or_v is not None:
    subgroups.append(('LN metastasis (N+)', f'{n} / {ne}', or_v, ci_l, ci_h, pv))

# N0
mask_n0 = (n_plus == 0)
or_v, ci_l, ci_h, pv, n, ne = firth_or(braf[mask_n0], recur[mask_n0])
if or_v is not None:
    subgroups.append(('Khong di can hach (N0)', f'{n} / {ne}', or_v, ci_l, ci_h, pv))
else:
    subgroups.append(('No LN metastasis (N0)', f'{n} / {ne}', 1, None, None, 1))

# Multifocal
mask_do = (da_o == 1)
or_v, ci_l, ci_h, pv, n, ne = firth_or(braf[mask_do], recur[mask_do])
if or_v is not None:
    subgroups.append(('Multifocal (>1 nodule)', f'{n} / {ne}', or_v, ci_l, ci_h, pv))

# Unifocal
mask_so = (da_o == 0)
or_v, ci_l, ci_h, pv, n, ne = firth_or(braf[mask_so], recur[mask_so])
if or_v is not None:
    subgroups.append(('Unifocal (1 nodule)', f'{n} / {ne}', or_v, ci_l, ci_h, pv))

# T1-T2
mask_t12 = (t3_plus == 0)
or_v, ci_l, ci_h, pv, n, ne = firth_or(braf[mask_t12], recur[mask_t12])
if or_v is not None:
    subgroups.append(('T stage (T1-T2)', f'{n} / {ne}', or_v, ci_l, ci_h, pv))

# T3+
mask_t3p = (t3_plus == 1)
or_v, ci_l, ci_h, pv, n, ne = firth_or(braf[mask_t3p], recur[mask_t3p])
if or_v is not None:
    subgroups.append(('T stage (T3+)', f'{n} / {ne}', or_v, ci_l, ci_h, pv))

# ===================== PLOT 1: FOREST PLOT =====================
fig, ax = plt.subplots(figsize=(8, 5.5))

y_pos = np.arange(len(subgroups))
or_vals = [s[3] for s in subgroups]
ci_lows = [s[3] for s in subgroups]
ci_highs = [s[4] for s in subgroups]

# Use log scale
ax.set_xscale('log')
ax.set_xlim(0.01, 100)

for i, (label, n_str, or_v, ci_l, ci_h, pv) in enumerate(subgroups):
    y = len(subgroups) - 1 - i
    if ci_l is not None and ci_h is not None:
        ax.plot([ci_l, ci_h], [y, y], color='navy', lw=2, zorder=3)
    ax.plot(or_v, y, marker='s', color='crimson', markersize=8, zorder=4)
    p_text = f'p={pv:.3f}' if pv is not None else ''
    ax.text(100, y, f'  {n_str}  OR={or_v:.2f}  {p_text}', va='center', fontsize=8.5)

ax.axvline(1, color='gray', ls='--', lw=0.8, zorder=1)
ax.set_yticks(range(len(subgroups)))
ax.set_yticklabels([s[0] for s in subgroups][::-1])
ax.set_xlabel('Odds Ratio (log scale) — Firth Logistic', fontsize=10)
ax.set_title('Forest Plot: BRAF V600E and Recurrence\nby Subgroup', fontsize=11, weight='bold')
ax.tick_params(labelsize=9)
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
fig.savefig('H:\\Q3 MNGOC\\braf_forest_plot.png', dpi=200)
plt.close(fig)
print(f'Saved: braf_forest_plot.png ({subgroups.__len__()} subgroups)')

# ===================== PLOT 2: STACKED BAR (BRAF x N) =====================
groups = [
    ('BRAF- / N0', (braf == 0) & (n_plus == 0)),
    ('BRAF+ / N0', (braf == 1) & (n_plus == 0)),
    ('BRAF+ / N+', (braf == 1) & (n_plus == 1)),
    ('BRAF- / N+', (braf == 0) & (n_plus == 1)),
]

bar_data = []
for gname, gmask in groups:
    n_tot = gmask.sum()
    n_rec = recur[gmask].sum()
    pct = n_rec / n_tot * 100 if n_tot > 0 else 0
    bar_data.append((gname, n_tot, n_rec, pct))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

colors = ['#2ecc71', '#3498db', '#e74c3c', '#f39c12']
labels = [d[0] for d in bar_data]
totals = [d[1] for d in bar_data]
recs = [d[2] for d in bar_data]
pcts = [d[3] for d in bar_data]

bars = ax1.bar(range(4), pcts, color=colors, edgecolor='white', width=0.6)
for i, (bar, pct, n_rec, n_tot) in enumerate(zip(bars, pcts, recs, totals)):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
             f'{pct:.1f}%\n({n_rec}/{n_tot})', ha='center', va='bottom', fontsize=8.5, fontweight='bold')
ax1.set_xticks(range(4))
ax1.set_xticklabels(labels, fontsize=8.5)
ax1.set_ylabel('Recurrence rate (%)', fontsize=10)
ax1.set_title('Recurrence rate by BRAF × N', fontsize=11, weight='bold')
ax1.set_ylim(0, max(pcts) * 1.5 + 5)
ax1.grid(axis='y', alpha=0.3)

# Right: N counts by BRAF status
n_braf_neg = (braf == 0).sum()
n_braf_pos = (braf == 1).sum()
n0_neg = ((braf == 0) & (n_plus == 0)).sum()
n_neg = ((braf == 0) & (n_plus == 1)).sum()
n0_pos = ((braf == 1) & (n_plus == 0)).sum()
n_pos = ((braf == 1) & (n_plus == 1)).sum()

x = np.arange(2)
width = 0.3
ax2.bar(x - width/2, [n0_neg, n0_pos], width, label='N0', color='#2ecc71', edgecolor='white')
ax2.bar(x + width/2, [n_neg, n_pos], width, label='N+', color='#e74c3c', edgecolor='white')
for i, (v0, v1) in enumerate(zip([n0_neg, n0_pos], [n_neg, n_pos])):
    ax2.text(i - width/2, v0 + 1, str(v0), ha='center', fontsize=9)
    ax2.text(i + width/2, v1 + 1, str(v1), ha='center', fontsize=9)
ax2.set_xticks(x)
ax2.set_xticklabels(['BRAF-', 'BRAF+'], fontsize=10)
ax2.set_ylabel('Number of patients', fontsize=10)
ax2.set_title('N+ / N0 distribution by BRAF', fontsize=11, weight='bold')
ax2.legend(fontsize=9)
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()
fig.savefig('H:\\Q3 MNGOC\\braf_stacked_bar.png', dpi=200)
plt.close(fig)
print('Saved: braf_stacked_bar.png')

# ===================== PLOT 3: MUTATION HEATMAP =====================
import itertools

genes = ['BRAF V600E', 'TERT', 'TP53']
gene_cols = [braf, tert, tp53]
n_genes = len(genes)

heat_data = np.zeros((n_genes, n_genes))
for i in range(n_genes):
    for j in range(n_genes):
        heat_data[i, j] = ((gene_cols[i] == 1) & (gene_cols[j] == 1)).sum()

fig, ax = plt.subplots(figsize=(5, 4))
im = ax.imshow(heat_data, cmap='YlOrRd', vmin=0, vmax=max(heat_data.max(), 5))
for i in range(n_genes):
    for j in range(n_genes):
        ax.text(j, i, int(heat_data[i, j]),
                ha='center', va='center', fontsize=12, fontweight='bold',
                color='white' if heat_data[i, j] > heat_data.max()/2 else 'black')
ax.set_xticks(range(n_genes))
ax.set_yticks(range(n_genes))
ax.set_xticklabels(genes, fontsize=9, rotation=30)
ax.set_yticklabels(genes, fontsize=9, rotation=0)
ax.set_title('Mutation co-occurrence matrix\n(patients with both mutations)', fontsize=10, weight='bold')

# Annotate totals
for i in range(n_genes):
    total = gene_cols[i].sum()
    ax.text(n_genes + 0.3, i, f'n={total}', va='center', fontsize=9)

cbar = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.08)
cbar.set_label('Number of patients', fontsize=9)

plt.tight_layout()
fig.savefig('H:\\Q3 MNGOC\\mutation_heatmap.png', dpi=200)
plt.close(fig)
print('Saved: mutation_heatmap.png')

# ===================== PLOT 4: RECURRENCE BY BRAF STATUS =====================
fig, ax = plt.subplots(figsize=(6, 4))

braf_neg = (braf == 0)
braf_pos = (braf == 1)
n_bn = braf_neg.sum()
n_bp = braf_pos.sum()
r_bn = recur[braf_neg].sum()
r_bp = recur[braf_pos].sum()

bars = ax.bar(['BRAF-', 'BRAF+'], [r_bn/n_bn*100, r_bp/n_bp*100],
              color=['#3498db', '#e74c3c'], edgecolor='white', width=0.4)
for bar, n, r in zip(bars, [n_bn, n_bp], [r_bn, r_bp]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
            f'{r}/{n}\n({r/n*100:.1f}%)', ha='center', va='bottom', fontsize=11, fontweight='bold')

# Also show OR from Firth
or_v, ci_l, ci_h, pv, _, _ = firth_or(braf, recur)
title = f'Recurrence rate by BRAF V600E\nFirth OR={or_v:.2f}, 95%CI {ci_l:.2f}-{ci_h:.2f}, p={pv:.3f}' if or_v else 'Recurrence rate by BRAF V600E'

ax.set_ylabel('Recurrence rate (%)', fontsize=10)
ax.set_title(title, fontsize=11, weight='bold')
ax.set_ylim(0, max(r_bn/n_bn*100, r_bp/n_bp*100) * 2)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
fig.savefig('H:\\Q3 MNGOC\\braf_recurrence_bar.png', dpi=200)
plt.close(fig)
print('Saved: braf_recurrence_bar.png')

print('\nDone! All plots saved.')
