"""
Publication-quality oncoplot (mutation landscape)
and improved forest plot for thyroid cancer analysis.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ===================== LOAD DATA =====================
f1 = 'H:\\Q3 MNGOC\\Bản sao của DỮ LIỆU Q3 GỐC MINH NGỌC MỚI.xlsx'
f2 = 'H:\\Q3 MNGOC\\KẾT QUẢ GEN Q3 MINH NGỌC 2.xlsx'
df1 = pd.read_excel(f1)
df2 = pd.read_excel(f2)
df = df1.merge(df2, on='Mã bệnh nhân', how='left')
n_total = len(df)

# Recurrence
us = df['Bất thường siêu âm cụ thể'].fillna('').str.contains('Hạch', na=False)
fna = df['Chọc tế bào chẩn đoán'].fillna('').str.contains('Đã chọc \\+', na=False)
reop = df['Xử trí'].fillna('').str.contains('Phẫu thuật', na=False)
fu = df['Bất thường trong lần khám gần nhất sau xử trí'].fillna('').str.contains('Có bất thường', na=False)
recur = (us | fna | reop | fu).astype(int).values

# T, N stage
def encode_t(v):
    if pd.isna(v): return np.nan
    s = str(v).strip().lower()
    if s in ['x','0']: return 0
    if s in ['1','1a','1b']: return 1
    if s in ['2']: return 2
    if s in ['3','3a','3b']: return 3
    if s in ['4','4a','4b','4c']: return 4
    return np.nan

def encode_n(v):
    if pd.isna(v): return np.nan
    s = str(v).strip().lower()
    if s in ['0']: return 0
    if s in ['1','1a']: return 1
    if s in ['1b']: return 2
    if s in ['2']: return 3
    return np.nan

t_stage = np.array([encode_t(v) for v in df['Phân loại T']])
n_stage = np.array([encode_n(v) for v in df['Phân loại N']])
n_plus = (n_stage > 0).astype(int)

# Genes
braf = (df['BRAFV600E'] == 2).astype(int).values
braf_other = ((df['BRAF'] == 2) & (df['BRAFV600E'] != 2)).astype(int).values
tert = (df['TERT'] == 2).astype(int).values
tp53 = (df['TP53'] == 2).astype(int).values

# ===================== ONCOPLOT =====================
print('Generating oncoplot...')

gene_names = ['BRAF V600E', 'BRAF other', 'TERT', 'TP53']
gene_data = np.array([braf, braf_other, tert, tp53])

# Sort patients: by BRAF first, then TERT, then TP53
order = np.lexsort((tp53, tert, braf))
gene_data_sorted = gene_data[:, order]
recur_sorted = recur[order]
n_sorted = n_plus[order]
t_sorted = t_stage[order]

n_genes = len(gene_names)

fig = plt.figure(figsize=(14, 4.5))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 0.25], hspace=0.02,
                       left=0.08, right=0.92, bottom=0.12, top=0.92)

# Main heatmap
ax = fig.add_subplot(gs[0])
cmap = ['#e8e8e8', '#e74c3c', '#f39c12', '#3498db', '#2ecc71']
colors_map = {
    0: '#e8e8e8',  # no mutation
    1: '#e74c3c',  # BRAF V600E
    2: '#f39c12',  # BRAF other
    3: '#3498db',  # TERT
    4: '#2ecc71',  # TP53
}
# Create color matrix
color_matrix = np.zeros((n_genes, n_total, 3))
for i in range(n_genes):
    for j in range(n_total):
        if gene_data_sorted[i, j] == 1:
            if i == 0:
                color_matrix[i, j] = [0.906, 0.298, 0.235]  # BRAF V600E red
            elif i == 1:
                color_matrix[i, j] = [0.953, 0.612, 0.071]  # BRAF other orange
            elif i == 2:
                color_matrix[i, j] = [0.204, 0.596, 0.859]  # TERT blue
            elif i == 3:
                color_matrix[i, j] = [0.180, 0.800, 0.443]  # TP53 green
        else:
            color_matrix[i, j] = [0.91, 0.91, 0.91]  # gray

ax.imshow(color_matrix, aspect='auto', interpolation='nearest')

# Grid lines
for i in range(n_genes + 1):
    ax.axhline(i - 0.5, color='white', lw=1.5)
ax.axvline(-0.5, color='white', lw=1.5)
ax.axvline(n_total - 0.5, color='white', lw=1.5)

ax.set_yticks(range(n_genes))
ax.set_yticklabels(gene_names, fontsize=10, style='italic')
ax.set_xticks([])
ax.set_xlim(-0.5, n_total - 0.5)

# Annotation: recurrence status on top
rec_colors = ['#2c3e50' if r == 1 else '#ecf0f1' for r in recur_sorted]
for j in range(n_total):
    ax.add_patch(Rectangle((j - 0.5, -0.5), 1, 0.3,
                           facecolor=rec_colors[j], edgecolor='gray', lw=0.3))

# Label for recurrence track
ax.text(-2.5, -0.35, 'Recurrence', fontsize=7, ha='center', va='center',
        rotation=0, fontweight='bold')

# Annotation: N+ status
for j in range(n_total):
    if n_sorted[j] == 1:
        ax.add_patch(Rectangle((j - 0.5, -0.2), 1, 0.2,
                               facecolor='#e74c3c', edgecolor='gray', lw=0.2, alpha=0.7))
    else:
        ax.add_patch(Rectangle((j - 0.5, -0.2), 1, 0.2,
                               facecolor='#bdc3c7', edgecolor='gray', lw=0.2, alpha=0.5))

ax.text(-2.5, -0.1, 'N+', fontsize=7, ha='center', va='center', fontweight='bold')

ax.set_ylim(n_genes - 0.5, -0.8)

# Legend
legend_elements = [
    mpatches.Patch(color='#e74c3c', label='BRAF V600E'),
    mpatches.Patch(color='#f39c12', label='BRAF other'),
    mpatches.Patch(color='#3498db', label='TERT'),
    mpatches.Patch(color='#2ecc71', label='TP53'),
    mpatches.Patch(color='#2c3e50', label='Recurrence'),
    mpatches.Patch(color='#ecf0f1', label='No recurrence'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=7,
          ncol=2, framealpha=0.9, edgecolor='gray')

# Bar plot: mutation frequency
ax_bar = fig.add_subplot(gs[1])
freqs = [gene_data[i].sum() / n_total * 100 for i in range(n_genes)]
bar_colors = ['#e74c3c', '#f39c12', '#3498db', '#2ecc71']
bars = ax_bar.barh(range(n_genes), freqs, color=bar_colors, edgecolor='white', height=0.6)
for i, (bar, f) in enumerate(zip(bars, freqs)):
    ax_bar.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                f'{f:.1f}%', va='center', fontsize=8, fontweight='bold')
ax_bar.set_xlim(0, 55)
ax_bar.set_yticks(range(n_genes))
ax_bar.set_yticklabels([])
ax_bar.set_xticks([0, 20, 40])
ax_bar.set_xlabel('Frequency (%)', fontsize=8)
ax_bar.invert_yaxis()
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)

fig.savefig('H:\\Q3 MNGOC\\oncoplot.png', dpi=300)
plt.close(fig)
print(f'  Saved: oncoplot.png ({n_total} patients, {n_genes} genes)')

# ===================== FOREST PLOT (improved) =====================
print('Generating forest plot...')

from firth_logit import firth_univariate

outcome_arr = recur

uni_vars_labels = [
    ('GENDER', 'Gender (Male)'),
    ('AGE_55PLUS', 'Age >=55 years'),
    ('MULTIFOCAL', 'Multifocal'),
    ('BILATERAL', 'Bilateral'),
    ('SIZE_OVER_10MM', 'Size >1cm'),
    ('T3_PLUS', 'T stage T3+'),
    ('N_PLUS', 'LN metastasis (N+)'),
    ('N1B', 'Lateral LN (N1b)'),
    ('CENTRAL_LN_5PLUS', 'Central LN >=5'),
    ('LATERAL_LN_PLUS', 'Lateral LN (+)'),
    ('BRAF_V600E', 'BRAF V600E'),
    ('TERT_MUT', 'TERT'),
    ('TP53_MUT', 'TP53'),
    ('CO_MUT', 'Co-mutation'),
]

# Build D dict from scratch (same as thyroid_pub)
D = {}
D['RECURRENCE'] = recur
nhan_tc = pd.to_numeric(df['Tổng số nhân trc pt'], errors='coerce').fillna(0)
D['MULTIFOCAL'] = (nhan_tc > 1).astype(int).values
vt = df['Vị trí nhân trc pt'].fillna('')
vt_b = vt.apply(lambda x: 1 if any(s in str(x) for s in ['T P', 'TP', 'T,P']) else 0)
vt2 = df['Vị trí'].fillna('').apply(lambda x: 1 if 'Thùy trái , Thùy phải' in str(x) else 0)
D['BILATERAL'] = (vt_b | vt2).astype(int).values
def parse_size_mm(val):
    if pd.isna(val): return np.nan
    v = str(val).strip().replace(',', '.')
    import re
    v = re.sub(r'\([^)]*\)', '', v)
    nums = re.findall(r'[\d.]+', v)
    nums = [n for n in nums if n and n != '.']
    if not nums: return np.nan
    s = max(float(n) for n in nums)
    if 'cm' in v.lower(): s *= 10
    if 0 < s < 1: s *= 10
    return s
kt = df['Kích thước nhân (mm) trc pt'].apply(parse_size_mm)
D['SIZE_OVER_10MM'] = (kt > 10).astype(int).values
D['T_STAGE'] = np.array([encode_t(v) for v in df['Phân loại T']])
D['T3_PLUS'] = (D['T_STAGE'] >= 3).astype(int)
n_stage_arr = np.array([encode_n(v) for v in df['Phân loại N']])
D['N_PLUS'] = (n_stage_arr > 0).astype(int)
D['N1B'] = (n_stage_arr >= 2).astype(int)
met_tt = df['Số hạch di căn (trung tâm)'].fillna(df['Số hạch di căn (trung tâm).1'])
D['CENTRAL_LN_5PLUS'] = (met_tt >= 5).fillna(0).astype(int).values
met_cb = df['Số hạch di căn (cùng bên)'].fillna(df['Số hạch di căn (cùng bên).1'])
met_db = df['Số hạch di căn (đối bên)'].fillna(df['Số hạch di căn (đối bên).1'])
D['LATERAL_LN_PLUS'] = ((met_cb.fillna(0) + met_db.fillna(0)) > 0).astype(int).values
D['BRAF_V600E'] = braf
D['TERT_MUT'] = tert
D['TP53_MUT'] = tp53
D['CO_MUT'] = ((braf.astype(int) + tert.astype(int) + tp53.astype(int)) >= 2).astype(int)

# Gender
gender_map = {'Nữ': 0, 'Nam': 1}
D['GENDER'] = np.array([gender_map.get(v, 0) for v in df['3. Giới']])
age = pd.to_numeric(df['Tuổi tại thời điểm chẩn đoán'], errors='coerce')
D['AGE_55PLUS'] = (age >= 55).astype(int).values

# Run univariate
forest_data = []
for var_code, var_label in uni_vars_labels:
    var_arr = np.array(D[var_code]) if var_code in D else np.array([])
    r = firth_univariate(var_label, var_arr, outcome_arr)
    forest_data.append(r)

# Filter to only those with valid OR
forest_plot_data = [r for r in forest_data if r['OR'] is not None]
n_plot = len(forest_plot_data)

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.set_xscale('log')
ax.set_xlim(0.005, 200)

y_pos = list(range(n_plot))
for i, r in enumerate(forest_plot_data):
    y = n_plot - 1 - i
    or_v = r['OR']
    ci_l = r['CI_low']
    ci_h = r['CI_high']
    if ci_l is not None and ci_h is not None:
        ax.plot([ci_l, ci_h], [y, y], color='#2c3e50', lw=2.5, zorder=3)
    ax.plot(or_v, y, marker='s', color='#e74c3c', markersize=9, zorder=4,
            markeredgecolor='white', markeredgewidth=0.5)
    p_v = r['p']
    p_text = f'p={p_v:.3f}' if p_v is not None else ''
    ci_text = f'({ci_l:.2f}-{ci_h:.2f})' if ci_l is not None else ''
    ax.text(210, y, f'{r["N"]}  {or_v:.2f}  {ci_text}  {p_text}', va='center',
            fontsize=8, fontfamily='monospace')

ax.axvline(1, color='gray', ls='--', lw=0.8, zorder=1, alpha=0.6)
ax.set_yticks(range(n_plot))
ax.set_yticklabels([r['Variable'] for r in forest_plot_data][::-1], fontsize=9)
ax.set_xlabel('Odds Ratio (95% CI) — log scale', fontsize=10)
ax.set_title('Forest Plot: Univariate Firth Logistic Regression\n'
             'Thyroid Cancer Recurrence (ATA 2025)', fontsize=11, weight='bold')
ax.tick_params(labelsize=9)
ax.grid(axis='x', alpha=0.3, which='both')
ax.set_xticks([0.01, 0.1, 1, 10, 100])
ax.set_xticklabels(['0.01', '0.1', '1', '10', '100'])

plt.tight_layout()
fig.savefig('H:\\Q3 MNGOC\\forest_plot_pub.png', dpi=300)
plt.close(fig)
print(f'  Saved: forest_plot_pub.png ({n_plot} variables)')

print('Done.')
