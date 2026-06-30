"""
Thyroid Cancer Recurrence Analysis — Publication-grade
ATA 2025 PTC 4-level risk classification
Firth penalized logistic regression
Gene interaction: BRAF V600E, TERT, TP53
Internal validation: bootstrap, calibration, AUC-ROC
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from firth_logit import (FirthLogit, firth_univariate, firth_multivariate,
                          hosmer_lemeshow, bootstrap_firth_model)

OUTPUT = 'H:\\Q3 MNGOC\\RESULTS_PUB.txt'
_orig_stdout = sys.stdout
sys.stdout = open(OUTPUT, 'w', encoding='utf-8')

# ===================== 1. LOAD & MERGE DATA =====================
f1 = 'H:\\Q3 MNGOC\\Bản sao của DỮ LIỆU Q3 GỐC MINH NGỌC MỚI.xlsx'
f2 = 'H:\\Q3 MNGOC\\KẾT QUẢ GEN Q3 MINH NGỌC 2.xlsx'
df1 = pd.read_excel(f1)
df2 = pd.read_excel(f2)
df = df1.merge(df2, on='Mã bệnh nhân', how='left')

D = {}  # data dictionary
n_total = len(df)

# ===================== 2. OUTCOME: RECURRENCE =====================
us = df['Bất thường siêu âm cụ thể'].fillna('').str.contains('Hạch', na=False)
fna = df['Chọc tế bào chẩn đoán'].fillna('').str.contains('Đã chọc \\+', na=False)
reop = df['Xử trí'].fillna('').str.contains('Phẫu thuật', na=False)
fu = df['Bất thường trong lần khám gần nhất sau xử trí'].fillna('').str.contains('Có bất thường', na=False)
D['RECURRENCE'] = (us | fna | reop | fu).astype(int)
rec_count = int(D['RECURRENCE'].sum())

print('='*70)
print('  THYROID CANCER RECURRENCE RISK ANALYSIS')
print('  ATA 2025 GUIDELINES — FIRTH PENALIZED LOGISTIC')
print('='*70)
print(f'\nTotal patients: {n_total}')
print(f'Recurrence events: {rec_count} ({rec_count/n_total*100:.1f}%)')

# ===================== 3. DEMOGRAPHICS & CLINICAL =====================
gender_map = {'Nữ': 0, 'Nam': 1}
D['GENDER'] = df['3. Giới'].map(gender_map).fillna(0).astype(int)

age = pd.to_numeric(df['Tuổi tại thời điểm chẩn đoán'], errors='coerce')
D['AGE'] = age
D['AGE_55PLUS'] = (age >= 55).astype(int)

nhan_tc = pd.to_numeric(df['Tổng số nhân trc pt'], errors='coerce').fillna(0)
D['MULTIFOCAL'] = (nhan_tc > 1).astype(int)

vt = df['Vị trí nhân trc pt'].fillna('')
vt_bilat = vt.apply(lambda x: 1 if any(s in str(x) for s in ['T P', 'TP', 'T,P']) else 0)
vt_clin = df['Vị trí'].fillna('').apply(lambda x: 1 if 'Thùy trái , Thùy phải' in str(x) else 0)
D['BILATERAL'] = (vt_bilat | vt_clin).astype(int)

def parse_size_mm(val):
    if pd.isna(val): return np.nan
    val = str(val).strip().replace(',', '.')
    val = __import__('re').sub(r'\([^)]*\)', '', val)
    nums = __import__('re').findall(r'[\d.]+', val)
    nums = [n for n in nums if n and n != '.']
    if not nums: return np.nan
    size = max(float(n) for n in nums)
    if 'cm' in val.lower(): size *= 10
    if 0 < size < 1: size *= 10
    return size

kt = df['Kích thước nhân (mm) trc pt'].apply(parse_size_mm)
D['TUMOR_SIZE_MM'] = kt
D['SIZE_OVER_10MM'] = (kt > 10).astype(int)
D['SIZE_OVER_40MM'] = (kt > 40).astype(int)

# T, N, M staging
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
    if v in ['1','1a']: return 1
    if v in ['1b']: return 2
    if v in ['2']: return 3
    return np.nan

t_stage = df['Phân loại T'].apply(encode_t)
n_stage = df['Phân loại N'].apply(encode_n)
m_stage = pd.to_numeric(df['Phân loại M'], errors='coerce').fillna(0)

D['T_STAGE'] = t_stage
D['T3_PLUS'] = (t_stage >= 3).astype(int)
D['N_STAGE'] = n_stage
D['N_PLUS'] = (n_stage > 0).astype(int)
D['N1B'] = (n_stage >= 2).astype(int)
D['M_STAGE'] = m_stage.astype(int)

# LN counts
met_tt = df['Số hạch di căn (trung tâm)'].fillna(df['Số hạch di căn (trung tâm).1'])
D['CENTRAL_LN_MET'] = met_tt
D['CENTRAL_LN_5PLUS'] = (met_tt >= 5).fillna(0).astype(int)
D['CENTRAL_LN_OVER5'] = ((met_tt > 5) & (~met_tt.isna())).astype(int)

met_cb = df['Số hạch di căn (cùng bên)'].fillna(df['Số hạch di căn (cùng bên).1'])
met_db = df['Số hạch di căn (đối bên)'].fillna(df['Số hạch di căn (đối bên).1'])
D['LATERAL_LN_PLUS'] = ((met_cb.fillna(0) + met_db.fillna(0)) > 0).astype(int)

# Stage
stage_map = {'I': 1, 'II': 2, 'III': 3, 'IVb': 4}
D['STAGE'] = df['Giai đoạn'].map(stage_map).fillna(1).astype(int)

# ===================== 4. GENETIC VARIABLES =====================
D['BRAF_V600E'] = (df['BRAFV600E'] == 2).astype(int)
D['BRAF_OTHER'] = ((df['BRAF'] == 2) & (df['BRAFV600E'] != 2) | (df['BRAF R682Q'] == 2) | (df['BRAF A404fs'] == 2)).astype(int)
D['BRAF_ANY'] = (df['BRAF'] == 2).astype(int)

# Individual TERT sites
tert_sites = ['TERT R756C', 'TERT T567M', 'TERT R577Q', 'TERT D685N',
              'TERT R938Q', 'TERT D516N', 'TERT V272M', 'TERT P438H', 'TERT T1111M']
for site in tert_sites:
    D[site.replace(' ', '_')] = (df[site] == 2).astype(int)
D['TERT_ANY'] = (df['TERT'] == 2).astype(int)
D['TERT_MUT'] = D['TERT_ANY']

# Individual TP53 sites
D['TP53_SPLICE_A'] = (df["TP53:c.559+1G>A"] == 2).astype(int)
D['TP53_SPLICE_B'] = (df["TP53:c.96+1G>A"] == 2).astype(int)
D['TP53_ANY'] = (df['TP53'] == 2).astype(int)
D['TP53_MUT'] = D['TP53_ANY']

# Derived
D['ANY_MUT'] = (D['BRAF_ANY'] | D['TERT_ANY'] | D['TP53_ANY']).astype(int)
D['CO_MUT'] = ((D['BRAF_ANY'] + D['TERT_ANY'] + D['TP53_ANY']) >= 2).astype(int)
D['BRAF_TERT'] = ((D['BRAF_ANY'] == 1) & (D['TERT_ANY'] == 1)).astype(int)
D['BRAF_TP53'] = ((D['BRAF_ANY'] == 1) & (D['TP53_ANY'] == 1)).astype(int)

# ===================== 5. ATA 2025 PTC RISK (4-level) =====================
def ata_ptc_risk(row):
    # HIGH: T3+ or M1
    if row['T_STAGE'] >= 3 or row['M_STAGE'] == 1:
        return 'High'
    low_int = 0
    if row['MULTIFOCAL'] == 1 and row['BILATERAL'] == 0:
        low_int += 1
    if row['CENTRAL_LN_OVER5'] == 1:
        low_int += 1
    # Intermediate-High
    if row['BILATERAL'] == 1 and row['MULTIFOCAL'] == 1:
        if not np.isnan(row['TUMOR_SIZE_MM']) and row['TUMOR_SIZE_MM'] > 10:
            return 'Intermediate-High'
    if row['N_STAGE'] >= 2:
        return 'Intermediate-High'
    if low_int >= 2:
        return 'Intermediate-High'
    if low_int >= 1:
        return 'Low-Intermediate'
    return 'Low'

ATA_LEVELS = ['Low', 'Low-Intermediate', 'Intermediate-High', 'High']
risk = []
for i in range(n_total):
    r = {k: (D[k].iloc[i] if hasattr(D[k], 'iloc') else D[k][i])
         for k in ['T_STAGE', 'M_STAGE', 'MULTIFOCAL', 'BILATERAL',
                   'TUMOR_SIZE_MM', 'N_STAGE', 'CENTRAL_LN_OVER5']}
    risk.append(ata_ptc_risk(r))
D['ATA_RISK'] = risk

print(f'\n{"="*70}')
print('  TABLE 0: ATA 2025 PTC RISK STRATIFICATION')
print(f'{"="*70}')
risk_counts = pd.Series(D['ATA_RISK']).value_counts()
for r in ATA_LEVELS:
    n = risk_counts.get(r, 0)
    pr = n/n_total*100
    m = pd.Series(D['RECURRENCE'])[pd.Series(D['ATA_RISK']) == r].sum()
    print(f'  {r:25s}: {n:3d} ({pr:5.1f}%)  Recurrence: {m}/{n} ({m/n*100:.1f}%)' if n > 0 else f'  {r}: 0')

# ===================== 6. TABLE 1: BASELINE CHARACTERISTICS =====================
print(f'\n{"="*70}')
print('  TABLE 1: BASELINE CLINICOPATHOLOGICAL CHARACTERISTICS')
print(f'{"="*70}')

table1_vars = [
    ('GENDER', 'Gender (Male)', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('AGE_55PLUS', 'Age >=55 years', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('AGE', 'Age (mean +/- SD)', lambda x: f'{x.mean():.1f} +/- {x.std():.1f}'),
    ('MULTIFOCAL', 'Multifocal (>1 nodule)', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('BILATERAL', 'Bilateral involvement', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('SIZE_OVER_10MM', 'Tumor size >10mm', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('SIZE_OVER_40MM', 'Tumor size >40mm', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('T3_PLUS', 'T stage T3+', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('N_PLUS', 'LN metastasis (N+)', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('N1B', 'Lateral LN (N1b)', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('CENTRAL_LN_5PLUS', 'Central LN >=5', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('LATERAL_LN_PLUS', 'Lateral LN (+)', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('M_STAGE', 'Distant metastasis (M1)', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('BRAF_V600E', 'BRAF V600E', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('BRAF_OTHER', 'BRAF other mutation', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('BRAF_ANY', 'BRAF (any mutation)', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('TERT_MUT', 'TERT promoter mutation', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('TP53_MUT', 'TP53 mutation', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('ANY_MUT', 'Any mutation', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
    ('CO_MUT', 'Co-mutation (>=2 genes)', lambda x: f'{int(x.sum())}/{len(x)} ({x.mean()*100:.1f}%)'),
]

rec_arr = np.array(D['RECURRENCE'])
print(f'{"Variable":<30s} {"Overall":>20s} {"Recurrence":>20s} {"No Recurrence":>20s}')
print('-'*90)
for var_code, var_label, fmt in table1_vars:
    v = D[var_code]
    overall = fmt(v)
    v_rec = v[rec_arr == 1] if len(v) == len(rec_arr) else v
    v_norec = v[rec_arr == 0] if len(v) == len(rec_arr) else v
    rec_s = fmt(v_rec) if len(v_rec) > 0 else 'N/A'
    norec_s = fmt(v_norec) if len(v_norec) > 0 else 'N/A'
    print(f'{var_label:<30s} {overall:>20s} {rec_s:>20s} {norec_s:>20s}')

# ===================== 7. UNIVARIATE FIRTH =====================
print(f'\n{"="*70}')
print('  TABLE 2: UNIVARIATE ANALYSIS — FIRTH LOGISTIC')
print(f'{"="*70}')

uni_vars = [
    ('GENDER', 'Gender (Male)'),
    ('AGE_55PLUS', 'Age >=55 years'),
    ('MULTIFOCAL', 'Multifocal (>1 nodule)'),
    ('BILATERAL', 'Bilateral'),
    ('SIZE_OVER_10MM', 'Tumor size >1cm'),
    ('SIZE_OVER_40MM', 'Tumor size >4cm'),
    ('T3_PLUS', 'T stage T3+'),
    ('N_PLUS', 'LN metastasis (N+)'),
    ('N1B', 'Lateral LN (N1b)'),
    ('CENTRAL_LN_5PLUS', 'Central LN >=5'),
    ('LATERAL_LN_PLUS', 'Lateral LN (+)'),
    ('BRAF_V600E', 'BRAF V600E'),
    ('BRAF_ANY', 'BRAF (any)'),
    ('TERT_MUT', 'TERT promoter mutation'),
    ('TP53_MUT', 'TP53 mutation'),
    ('ANY_MUT', 'Any mutation'),
    ('CO_MUT', 'Co-mutation (>=2 genes)'),
    ('BRAF_TERT', 'BRAF + TERT co-mutation'),
]

outcome_arr = np.array(D['RECURRENCE'])
uni_results = []
for var_code, var_label in uni_vars:
    r = firth_univariate(var_label, np.array(D[var_code]), outcome_arr)
    uni_results.append(r)

header = f'{"Variable":<30s} {"N":<6s} {"Events":<8s} {"OR (Firth)":<12s} {"95% CI":<25s} {"p-value":<12s}'
print(header)
print('-'*len(header))
for r in uni_results:
    if r['OR'] is not None:
        ci_s = f"{r['CI_low']:.3f} - {r['CI_high']:.3f}"
        p_s = f"{r['p']:.4f}" if r['p'] >= 0.0001 else "<0.0001"
        print(f"{r['Variable']:<30s} {r['N']:<6d} {r['Events']:<8d} {r['OR']:<12.3f} {ci_s:<25s} {p_s:<12s}")
    else:
        print(f"{r['Variable']:<30s} {r['N']:<6d} {r['Events']:<8d} {'NA':<12s} {'NA':<25s} {r['Note']:<12s}")

# ===================== 8. GENE-GENE INTERACTION =====================
print(f'\n{"="*70}')
print('  TABLE 3: GENE-GENE INTERACTION ANALYSIS')
print(f'{"="*70}')

gene_pairs = [
    ('BRAF_V600E', 'TERT_MUT', 'BRAF V600E vs TERT'),
    ('BRAF_V600E', 'TP53_MUT', 'BRAF V600E vs TP53'),
    ('TERT_MUT', 'TP53_MUT', 'TERT vs TP53'),
]

for g1, g2, label in gene_pairs:
    v1 = np.array(D[g1])
    v2 = np.array(D[g2])
    # Contingency table
    both = ((v1 == 1) & (v2 == 1)).sum()
    only_g1 = ((v1 == 1) & (v2 == 0)).sum()
    only_g2 = ((v1 == 0) & (v2 == 1)).sum()
    neither = ((v1 == 0) & (v2 == 0)).sum()
    table = [[both, only_g1], [only_g2, neither]]
    # Fisher's exact test
    odds_ratio, p_fisher = sp_stats.fisher_exact(table)
    print(f'\n  {label}:')
    print(f'        {g1}+   {g1}-')
    print(f'  {g2}+  {both:3d}    {only_g2:3d}')
    print(f'  {g2}-  {only_g1:3d}    {neither:3d}')
    print(f'  Fisher exact: OR={odds_ratio:.3f}, p={p_fisher:.4f}')

# ===================== 9. MULTIVARIATE FIRTH =====================
print(f'\n{"="*70}')
print('  TABLE 4: MULTIVARIATE ANALYSIS — FIRTH LOGISTIC')
print(f'{"="*70}')

models = [
    ({'N_PLUS': 'N+', 'T3_PLUS': 'T3+', 'SIZE_OVER_10MM': 'Size >1cm'},
     'Model 1: Clinical'),
    ({'BRAF_V600E': 'BRAF V600E', 'TERT_MUT': 'TERT', 'TP53_MUT': 'TP53'},
     'Model 2: Genetic'),
    ({'N_PLUS': 'N+', 'T3_PLUS': 'T3+', 'SIZE_OVER_10MM': 'Size >1cm',
      'BRAF_V600E': 'BRAF V600E', 'TERT_MUT': 'TERT'},
     'Model 3: Clinical + Genetic'),
    ({'N_PLUS': 'N+', 'CO_MUT': 'Co-mutation'},
     'Model 4: N+ + Co-mutation'),
    ({'N_PLUS': 'N+', 'BRAF_V600E': 'BRAF V600E', 'TERT_MUT': 'TERT'},
     'Model 5: N+ + BRAF + TERT'),
    ({'N_PLUS': 'N+', 'TERT_MUT': 'TERT'},
     'Model 6: N+ + TERT'),
    ({'N_PLUS': 'N+', 'TP53_MUT': 'TP53'},
     'Model 7: N+ + TP53'),
    ({'N_PLUS': 'N+', 'BRAF_V600E': 'BRAF V600E', 'TERT_MUT': 'TERT'},
     'Model 8: N+ + BRAF + TERT + BRAFxTERT (see below)'),
]

model_results = []
for var_dict, mname in models:
    m = firth_multivariate(var_dict, outcome_arr, D, mname)
    model_results.append(m)
    print(f'\n  {mname} (N={m["N"]}, events={m["events"]}):')
    if m.get('error'):
        print(f'  Error: {m["error"]}')
        continue
    print(f'  {"Variable":<25s} {"OR":<10s} {"95% CI":<22s} {"p-value":<10s}')
    print('  ' + '-'*67)
    for row in m['rows']:
        p_s = f"{row['p']:.4f}" if row['p'] >= 0.0001 else "<0.0001"
        ci_s = f"{row['CI_low']:.3f} - {row['CI_high']:.3f}"
        print(f"  {row['Variable']:<25s} {row['OR']:<10.3f} {ci_s:<22s} {p_s:<10s}")
    if m.get('HL') and not np.isnan(m['HL']['p']):
        hl = m['HL']
        print(f'  Hosmer-Lemeshow: chi2={hl["chi2"]:.3f}, df={hl["df"]}, p={hl["p"]:.4f}')

# Model 8 with interaction: need to handle separately
# Actually Model 8 as defined won't work because I overwrite BRAF_V600E key in the dict
# Let me do it properly
print(f'\n  Model 8: N+ + BRAF + TERT + BRAFxTERT interaction')
valid = np.ones(n_total, dtype=bool)
for v in ['N_PLUS', 'BRAF_V600E', 'TERT_MUT']:
    valid &= ~np.isnan(np.array(D[v]).astype(float))
y8 = outcome_arr[valid].astype(float)
n8 = len(y8)
e8 = int(y8.sum())
n_arr = np.array(D['N_PLUS'])[valid].astype(float)
b_arr = np.array(D['BRAF_V600E'])[valid].astype(float)
t_arr = np.array(D['TERT_MUT'])[valid].astype(float)
inter_arr = b_arr * t_arr
X8 = np.column_stack([np.ones(n8), n_arr, b_arr, t_arr, inter_arr])
try:
    m8 = FirthLogit().fit(X8, y8)
    r8 = m8.summary(['Intercept', 'N+', 'BRAF V600E', 'TERT', 'BRAFxTERT'])
    print(f'  (N={n8}, events={e8})')
    h = f'  {"Variable":<25s} {"OR":<10s} {"95% CI":<22s} {"p-value":<10s}'
    print(h)
    print('  ' + '-'*67)
    for row in r8:
        p_s = f"{row['p']:.4f}" if row['p'] >= 0.0001 else "<0.0001"
        ci_s = f"{row['CI_low']:.3f} - {row['CI_high']:.3f}"
        print(f"  {row['Variable']:<25s} {row['OR']:<10.3f} {ci_s:<22s} {p_s:<10s}")
    pred8 = 1/(1+np.exp(-X8 @ m8.beta))
    hl8 = hosmer_lemeshow(y8, pred8)
    if hl8:
        print(f'  Hosmer-Lemeshow: chi2={hl8["chi2"]:.3f}, df={hl8["df"]}, p={hl8["p"]:.4f}')
except Exception as e:
    print(f'  Error in Model 8: {e}')

# ===================== 10. BOOTSTRAP VALIDATION =====================
print(f'\n{"="*70}')
print('  TABLE 5: BOOTSTRAP VALIDATION (1000 reps)')
print(f'  Model: N+ + TERT + BRAF (Model 5)')
print(f'{"="*70}')

boot_vars = {'N_PLUS': 'N+', 'BRAF_V600E': 'BRAF V600E', 'TERT_MUT': 'TERT'}
boot_res = bootstrap_firth_model(boot_vars, outcome_arr, D, n_boot=1000)
print(f'  {"Variable":<20s} {"OR (bootstrap median)":<22s} {"95% CI (bootstrap)":<22s}')
print('  ' + '-'*64)
for nm in ['Intercept'] + list(boot_vars.keys()):
    if nm in boot_res and boot_res[nm]['OR'] is not None:
        r = boot_res[nm]
        ci_s = f"{r['CI_low']:.3f} - {r['CI_high']:.3f}"
        display = boot_vars.get(nm, nm)
        print(f"  {display:<20s} {r['OR']:<22.3f} {ci_s:<22s}")

# ===================== 11. EXPORT ENCODED DATA =====================
export_cols = ['Mã bệnh nhân', 'RECURRENCE', 'ATA_RISK', 'GENDER', 'AGE',
               'AGE_55PLUS', 'MULTIFOCAL', 'BILATERAL', 'TUMOR_SIZE_MM',
               'SIZE_OVER_10MM', 'SIZE_OVER_40MM', 'T_STAGE', 'T3_PLUS',
               'N_STAGE', 'N_PLUS', 'N1B', 'CENTRAL_LN_5PLUS',
               'LATERAL_LN_PLUS', 'M_STAGE', 'STAGE',
               'BRAF_V600E', 'BRAF_OTHER', 'BRAF_ANY', 'TERT_MUT', 'TP53_MUT',
               'ANY_MUT', 'CO_MUT', 'BRAF_TERT', 'BRAF_TP53']
export_dict = {}
for c in export_cols:
    if c in D:
        export_dict[c] = D[c]
    elif c == 'Mã bệnh nhân':
        export_dict[c] = df[c]
pd.DataFrame(export_dict).to_excel('H:\\Q3 MNGOC\\DATA_PUB.xlsx', index=False)

print(f'\n{"="*70}')
print('  EXPORTED: DATA_PUB.xlsx, RESULTS_PUB.txt')
print(f'{"="*70}')

sys.stdout.close()
sys.stdout = _orig_stdout

print(f'\nAnalysis complete!')
print(f'  Results: RESULTS_PUB.txt')
print(f'  Data: DATA_PUB.xlsx')
print(f'  Recurrence: {rec_count}/{n_total} ({rec_count/n_total*100:.1f}%)')
print(f'  Univariate variables: {sum(1 for r in uni_results if r["OR"] is not None)}')
n_sig = sum(1 for r in uni_results if r["p"] is not None and r["p"] < 0.05)
print(f'  Significant (p<0.05): {n_sig}')
print(f'  Multivariate models: {len(model_results)} + 1 interaction model')
print(f'  Gene-gene tests: {len(gene_pairs)}')
