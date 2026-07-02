"""
Figure 1: Mutation Landscape + Co-occurrence Matrix
Panel A: Grouped bar — mutation frequency by recurrence status
Panel B: Co-occurrence heatmap (BRAF V600E, TERT, TP53)
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
})

# ===================== LOAD DATA =====================
f1 = 'H:\\Q3 MNGOC\\Bản sao của DỮ LIỆU Q3 GỐC MINH NGỌC MỚI.xlsx'
f2 = 'H:\\Q3 MNGOC\\KẾT QUẢ GEN Q3 MINH NGỌC 2.xlsx'
df = pd.read_excel(f1).merge(pd.read_excel(f2), on='Mã bệnh nhân', how='left')

# Outcome
us = df['Bất thường siêu âm cụ thể'].fillna('').str.contains('Hạch', na=False)
fna = df['Chọc tế bào chẩn đoán'].fillna('').str.contains('Đã chọc \\+', na=False)
reop = df['Xử trí'].fillna('').str.contains('Phẫu thuật', na=False)
fu = df['Bất thường trong lần khám gần nhất sau xử trí'].fillna('').str.contains('Có bất thường', na=False)
recur = (us | fna | reop | fu).astype(int).values

# Genes
braf_v600e = (df['BRAFV600E'] == 2).astype(int).values
braf_other = ((df['BRAF'] == 2) & (df['BRAFV600E'] != 2)).astype(int).values
tert = (df['TERT'] == 2).astype(int).values
tp53 = (df['TP53'] == 2).astype(int).values

# ===================== COMPUTE FREQUENCIES =====================
genes = [
    ('BRAF V600E', braf_v600e),
    ('BRAF other', braf_other),
    ('TERT', tert),
    ('TP53', tp53),
]

n_rec = int(recur.sum())
n_norec = len(recur) - n_rec

gene_data = []
for name, arr in genes:
    r_mut = arr[recur == 1]
    nr_mut = arr[recur == 0]
    r_count = int(r_mut.sum())
    nr_count = int(nr_mut.sum())
    r_pct = r_count / n_rec * 100 if n_rec > 0 else 0
    nr_pct = nr_count / n_norec * 100 if n_norec > 0 else 0
    # Fisher exact test
    rec_yes = r_count
    rec_no = n_rec - r_count
    norec_yes = nr_count
    norec_no = n_norec - nr_count
    _, p_val = stats.fisher_exact([[rec_yes, rec_no], [norec_yes, norec_no]])
    gene_data.append({
        'name': name,
        'r_count': r_count, 'nr_count': nr_count,
        'r_pct': r_pct, 'nr_pct': nr_pct,
        'p': p_val,
        'r_total': n_rec, 'nr_total': n_norec,
    })

# Co-occurrence matrix
gene_arrs = {
    'BRAF V600E': braf_v600e,
    'TERT': tert,
    'TP53': tp53,
}
gene_order = ['BRAF V600E', 'TERT', 'TP53']
co_matrix = np.zeros((3, 3), dtype=int)
for i, g1 in enumerate(gene_order):
    for j, g2 in enumerate(gene_order):
        if i == j:
            co_matrix[i][j] = int(gene_arrs[g1].sum())
        else:
            co_matrix[i][j] = int(np.sum((gene_arrs[g1] == 1) & (gene_arrs[g2] == 1)))

# ===================== BUILD FIGURE =====================
fig = plt.figure(figsize=(8, 5.5))
gs = fig.add_gridspec(1, 4, width_ratios=[2.5, 0.15, 0.3, 1.1], hspace=0)

# --- Panel A: Grouped bar ---
ax1 = fig.add_subplot(gs[0, 0])
x = np.arange(len(genes))
width = 0.35

colors_rec = '#c0392b'
colors_norec = '#2980b9'

for i, g in enumerate(gene_data):
    ax1.bar(x[i] - width/2, g['r_pct'], width,
            color=colors_rec, alpha=0.85, edgecolor='white', linewidth=0.5)
    ax1.bar(x[i] + width/2, g['nr_pct'], width,
            color=colors_norec, alpha=0.85, edgecolor='white', linewidth=0.5)
    # Labels
    ax1.text(x[i] - width/2, g['r_pct'] + 1.5, f"{g['r_count']}/{g['r_total']}",
             ha='center', va='bottom', fontsize=6.5, color=colors_rec, fontweight='bold')
    ax1.text(x[i] + width/2, g['nr_pct'] + 1.5, f"{g['nr_count']}/{g['nr_total']}",
             ha='center', va='bottom', fontsize=6.5, color=colors_norec, fontweight='bold')
    # p-value
    p_text = f"p={g['p']:.2f}" if g['p'] >= 0.01 else ("p<0.01" if g['p'] < 0.01 else "p=1.00")
    ax1.text(x[i], max(g['r_pct'], g['nr_pct']) + 6, p_text,
             ha='center', va='bottom', fontsize=7, fontstyle='italic', color='dimgray')

ax1.set_xticks(x)
ax1.set_xticklabels([g['name'] for g in gene_data], fontsize=8.5)
ax1.set_ylabel('Mutation frequency (%)', fontsize=9)
ax1.set_ylim(0, 60)
ax1.set_title('A. Mutation Frequency by Recurrence Status', fontsize=10, fontweight='bold', loc='left')

# Legend
legend_patches = [
    mpatches.Patch(color=colors_rec, alpha=0.85, label=f'Recurrence (n={n_rec})'),
    mpatches.Patch(color=colors_norec, alpha=0.85, label=f'No recurrence (n={n_norec})'),
]
ax1.legend(handles=legend_patches, fontsize=7.5, loc='upper right', framealpha=0.9)

# Grid
ax1.set_axisbelow(True)
ax1.yaxis.grid(True, alpha=0.3)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# --- GAP ---

# --- Panel B: Co-occurrence heatmap ---
ax3 = fig.add_subplot(gs[0, 3])

# Create a nicer heatmap
im = ax3.imshow(co_matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=50)

# Annotate cells
for i in range(3):
    for j in range(3):
        val = co_matrix[i][j]
        color = 'white' if val > 25 else 'black'
        weight = 'bold' if i == j else 'normal'
        ax3.text(j, i, str(val), ha='center', va='center', fontsize=11,
                 color=color, fontweight=weight)

ax3.set_xticks(range(3))
ax3.set_yticks(range(3))
ax3.set_xticklabels(['BRAF\nV600E', 'TERT', 'TP53'], fontsize=8)
ax3.set_yticklabels(['BRAF V600E', 'TERT', 'TP53'], fontsize=8)
ax3.set_title('B. Co-occurrence', fontsize=10, fontweight='bold', loc='left', pad=10)

# Colorbar
cbar = fig.colorbar(im, ax=ax3, fraction=0.046, pad=0.04)
cbar.set_label('Patients (n)', fontsize=7.5)
cbar.ax.tick_params(labelsize=6.5)

# Remove spines
for spine in ax3.spines.values():
    spine.set_visible(False)

# Overall title
fig.suptitle('Figure 1. Mutation Landscape', fontsize=12, fontweight='bold', y=1.01)

plt.tight_layout()
fig.savefig('H:\\Q3 MNGOC\\figure1_mutation.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print('Saved: figure1_mutation.png')
