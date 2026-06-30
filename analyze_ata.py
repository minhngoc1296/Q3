import pandas as pd
import numpy as np
import re
import sys
from scipy import stats
import statsmodels.api as sm

def hosmer_lemeshow_test(y, pred, g=10):
    """Manual Hosmer-Lemeshow test"""
    df = pd.DataFrame({'y': y, 'pred': pred})
    df = df.sort_values('pred').reset_index(drop=True)
    n = len(df)
    # Adjust groups if too small
    g = min(g, n // 3)
    if g < 3:
        return 0, 1.0, 1
    df['group'] = pd.cut(np.arange(n), bins=g, labels=False)
    obs = df.groupby('group')['y'].sum()
    exp = df.groupby('group')['pred'].sum()
    # Remove groups with exp=0 to avoid division by zero
    mask = exp > 0
    if mask.sum() < 2:
        return 0, 1.0, 1
    chi2 = ((obs[mask] - exp[mask])**2 / exp[mask]).sum()
    df_g = mask.sum() - 2
    pval = 1 - stats.chi2.cdf(chi2, df_g)
    return chi2, pval, df_g

# ===================== FIRTH LOGISTIC REGRESSION =====================
class FirthLogit:
    """Firth's penalized logistic regression (Jeffreys prior)
    Handles complete separation by adding 0.5*log|I(β)| penalty.
    """
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

# Redirect output to UTF-8 file
sys.stdout = open('H:\\Q3 MNGOC\\RESULTS_ATA2025_FIRTH.txt', 'w', encoding='utf-8')

# ===================== READ DATA =====================
f1 = 'H:\\Q3 MNGOC\\Bản sao của DỮ LIỆU Q3 GỐC MINH NGỌC MỚI.xlsx'
f2 = 'H:\\Q3 MNGOC\\KẾT QUẢ GEN Q3 MINH NGỌC 2.xlsx'
df1 = pd.read_excel(f1)
df2 = pd.read_excel(f2)
df_orig = df1.merge(df2, on='Mã bệnh nhân', how='left')
cols_keep = list(df1.columns) + [c for c in df2.columns if c not in df1.columns]
data = {col: df_orig[col] for col in cols_keep}

# ===================== DEFINE RECURRENCE OUTCOME =====================
us_abnormal = df_orig['Bất thường siêu âm cụ thể'].fillna('').str.contains('Hạch', na=False)
fna_positive = df_orig['Chọc tế bào chẩn đoán'].fillna('').str.contains('Đã chọc \\+', na=False)
re_operation = df_orig['Xử trí'].fillna('').str.contains('Phẫu thuật', na=False)
last_fu_abnormal = df_orig['Bất thường trong lần khám gần nhất sau xử trí'].fillna('').str.contains('Có bất thường', na=False)

data['RECURRENCE'] = (us_abnormal | fna_positive | re_operation | last_fu_abnormal).astype(int)

print('='*80)
print('THYROID CANCER RECURRENCE RISK ANALYSIS')
print('ACCORDING TO ATA 2025 GUIDELINES')
print('='*80)
n_patients = len(data['RECURRENCE'])
rec_count = sum(data['RECURRENCE'])
print(f'\nRecurrence rate: {rec_count}/{n_patients} ({rec_count/n_patients*100:.1f}%)')
print(f'Total patients: {n_patients}')

# ===================== ATA 2025 ENCODING =====================
def parse_size_mm(val):
    if pd.isna(val):
        return np.nan
    val = str(val).strip().replace(',', '.')
    val = re.sub(r'\([^)]*\)', '', val)
    nums = re.findall(r'[\d.]+', val)
    nums = [n for n in nums if n and n != '.']
    if not nums:
        return np.nan
    nums_float = [float(n) for n in nums]
    size_mm = max(nums_float)
    if 'cm' in val.lower():
        size_mm = size_mm * 10
    if 0 < size_mm < 1:
        size_mm = size_mm * 10
    return size_mm

nhan_tc = pd.to_numeric(df_orig['Tổng số nhân trc pt'], errors='coerce').fillna(0)
data['DA_O'] = (nhan_tc > 1).astype(int)

vt = df_orig['Vị trí nhân trc pt'].fillna('')
data['HAI_BEN'] = vt.apply(lambda x: 1 if any(s in str(x) for s in ['T P', 'TP', 'T,P']) else 0)
vt2 = df_orig['Vị trí'].fillna('')
data['HAI_BEN_CLIN'] = vt2.apply(lambda x: 1 if 'Thùy trái , Thùy phải' in str(x) else 0)
data['HAI_BEN_FINAL'] = (data['HAI_BEN'] | data['HAI_BEN_CLIN']).astype(int)

kt_mm = df_orig['Kích thước nhân (mm) trc pt'].apply(parse_size_mm)
data['KT_MM'] = kt_mm
data['KT_CM'] = kt_mm / 10
data['KT_TREN_1CM'] = (kt_mm > 10).astype(int)
data['KT_TREN_4CM'] = (kt_mm > 40).astype(int)
data['KT_MISSING'] = kt_mm.isna().astype(int)

def encode_t(val):
    if pd.isna(val):
        return np.nan
    val = str(val).strip().lower()
    if val in ['x', '0']:
        return 0
    if val in ['1', '1a', '1b']:
        return 1
    if val in ['2']:
        return 2
    if val in ['3', '3a', '3b']:
        return 3
    if val in ['4', '4a', '4b', '4c']:
        return 4
    return np.nan

def encode_n(val):
    if pd.isna(val):
        return np.nan
    val = str(val).strip().lower()
    if val in ['0']:
        return 0
    if val in ['1', '1a']:
        return 1
    if val in ['1b']:
        return 2
    if val in ['2']:
        return 3
    return np.nan

t_stage = df_orig['Phân loại T'].apply(encode_t)
n_stage = df_orig['Phân loại N'].apply(encode_n)

data['T_STAGE'] = t_stage
data['T3_PLUS'] = (t_stage >= 3).astype(int)
data['N_STAGE'] = n_stage
data['N_PLUS'] = (n_stage > 0).astype(int)
data['N1b'] = (n_stage >= 2).astype(int)

met_tt = df_orig['Số hạch di căn (trung tâm)'].fillna(df_orig['Số hạch di căn (trung tâm).1'])
data['MET_HAI_TRUNGTAM'] = met_tt
data['HACH_TRUNGTAM_5PLUS'] = (met_tt >= 5).fillna(0).astype(int)
data['HACH_TRUNGTAM_OVER5'] = ((met_tt > 5) & (~met_tt.isna())).astype(int)

met_cb = df_orig['Số hạch di căn (cùng bên)'].fillna(df_orig['Số hạch di căn (cùng bên).1'])
met_db = df_orig['Số hạch di căn (đối bên)'].fillna(df_orig['Số hạch di căn (đối bên).1'])
data['CO_HACH_COBEN'] = ((met_cb.fillna(0) + met_db.fillna(0)) > 0).astype(int)

m_stage = df_orig['Phân loại M'].fillna(0).astype(int)
data['M_STAGE'] = m_stage

# ===================== GENETIC VARIABLES =====================
data['BRAF_MUT'] = (df_orig['BRAF'] == 2).astype(int)
data['BRAF_V600E'] = (df_orig['BRAFV600E'] == 2).astype(int)
data['BRAF_OTHER'] = ((df_orig['BRAF'] == 2) & (df_orig['BRAFV600E'] != 2) | (df_orig['BRAF R682Q'] == 2) | (df_orig['BRAF A404fs'] == 2)).astype(int)
data['BRAF_ANY'] = (df_orig['BRAF'] == 2).astype(int)
data['TERT_MUT'] = (df_orig['TERT'] == 2).astype(int)
data['TP53_MUT'] = (df_orig['TP53'] == 2).astype(int)
data['ANY_MUT'] = (data['BRAF_MUT'] | data['TERT_MUT'] | data['TP53_MUT']).astype(int)
data['CO_MUT'] = ((data['BRAF_MUT'] + data['TERT_MUT'] + data['TP53_MUT']) >= 2).astype(int)
data['BRAF_TERT'] = ((data['BRAF_MUT'] == 1) & (data['TERT_MUT'] == 1)).astype(int)
data['BRAF_TP53'] = ((data['BRAF_MUT'] == 1) & (data['TP53_MUT'] == 1)).astype(int)
data['TERT_TP53'] = ((data['TERT_MUT'] == 1) & (data['TP53_MUT'] == 1)).astype(int)

# ===================== ATA 2025 PTC RISK STRATIFICATION =====================
# 4 levels: LOW <10%, LOW-INTERMEDIATE 10-15%, INTERMEDIATE-HIGH 16-30%, HIGH >30%
# T3 → T3b (all HIGH). Skip unavailable criteria: margins, vascular invasion,
# aggressive histology, extranodal extension, LN metastasis size.

def ata_ptc_risk(row):
    # === HIGH > 30% ===
    # T3b (T3→T3b) or T4, or M1
    if row['T_STAGE'] >= 3:
        return 'High'
    if row['M_STAGE'] == 1:
        return 'High'

    # Count low-intermediate factors
    low_int_count = 0
    # Unilateral multifocality → LOW-INTERMEDIATE factor
    if row['DA_O'] == 1 and row['HAI_BEN_FINAL'] == 0:
        low_int_count += 1
    # N1a >5 nodes → LOW-INTERMEDIATE factor
    if row['HACH_TRUNGTAM_OVER5'] == 1:
        low_int_count += 1

    # === INTERMEDIATE-HIGH ≥ 16-30% ===
    # T1-T2 + (bilateral multifocal >1cm OR N1b OR 2+ low-intermediate factors)
    if row['HAI_BEN_FINAL'] == 1 and row['DA_O'] == 1:
        if not pd.isna(row['KT_MM']) and row['KT_MM'] > 10:
            return 'Intermediate-High'
    if row['N_STAGE'] >= 2:  # N1b
        return 'Intermediate-High'
    if low_int_count >= 2:
        return 'Intermediate-High'

    # === LOW-INTERMEDIATE 10-15% ===
    if low_int_count >= 1:
        return 'Low-Intermediate'

    # === LOW < 10% ===
    return 'Low'

ATA_LEVELS = ['Low', 'Low-Intermediate',
              'Intermediate-High', 'High']

risk_labels = []
for i in range(n_patients):
    r = {k: (data[k].iloc[i] if hasattr(data[k], 'iloc') else data[k][i])
         for k in ['T_STAGE', 'M_STAGE', 'DA_O', 'HAI_BEN_FINAL', 'KT_MM',
                   'N_STAGE', 'HACH_TRUNGTAM_OVER5']}
    risk_labels.append(ata_ptc_risk(r))
data['ATA_RISK'] = risk_labels

print(f'\n{"="*80}')
print('1. ATA 2025 PTC RISK STRATIFICATION')
print(f'{"="*80}')
risk_counts = pd.Series(data['ATA_RISK']).value_counts()
for r in ATA_LEVELS:
    n = risk_counts.get(r, 0)
    pct = n/n_patients*100
    print(f'  {r}: {n} ({pct:.1f}%)')

print(f'\nRecurrence rate by risk stratum:')
for r in ATA_LEVELS:
    mask = pd.Series(data['ATA_RISK']) == r
    n_total = mask.sum()
    if n_total > 0:
        n_rec = pd.Series(data['RECURRENCE'])[mask].sum()
        print(f'  {r}: {n_rec}/{n_total} ({n_rec/n_total*100:.1f}%)')

# ===================== DESCRIPTIVE TABLE =====================
print(f'\n{"="*80}')
print('2. DESCRIPTIVE STATISTICS')
print(f'{"="*80}')

desc_vars = [
    ('DA_O', 'Multifocal (>1 nodule)'),
    ('HAI_BEN_FINAL', 'Bilateral'),
    ('KT_TREN_1CM', 'Tumor size >1cm'),
    ('KT_TREN_4CM', 'Tumor size >4cm'),
    ('T3_PLUS', 'T stage 3+'),
    ('N_PLUS', 'Lymph node metastasis (N+)'),
    ('N1b', 'Lateral neck metastasis (N1b)'),
    ('HACH_TRUNGTAM_5PLUS', 'Central LN >=5'),
    ('CO_HACH_COBEN', 'Lateral LN (+)'),
    ('BRAF_V600E', 'BRAF V600E'),
    ('BRAF_OTHER', 'BRAF other'),
    ('TERT_MUT', 'TERT mutation'),
    ('TP53_MUT', 'TP53 mutation'),
    ('ANY_MUT', 'Any mutation'),
    ('CO_MUT', 'Co-mutation (>=2 genes)'),
    ('BRAF_TERT', 'BRAF + TERT'),
    ('BRAF_TP53', 'BRAF + TP53'),
]

for var_code, var_label in desc_vars:
    n = data[var_code].sum() if isinstance(data[var_code], (int, np.integer)) else pd.Series(data[var_code]).sum()
    if isinstance(data[var_code], (int, np.integer)):
        total = 1
        pct = n * 100
    else:
        total = len(data[var_code])
        pct = n / total * 100
    print(f'  {var_label}: {n}/{total} ({pct:.1f}%)')

# ===================== UNIVARIATE ANALYSIS (FIRTH) =====================
results = []

def firth_univariate(var_name, var_data, outcome):
    valid = ~np.isnan(outcome) & ~np.isnan(var_data.astype(float))
    y = outcome[valid].astype(float)
    x = var_data[valid].astype(float)
    n = len(y)
    n_events = int(y.sum())
    n_nonevents = n - n_events
    n_pos = int((x > 0).sum())
    if n < 10 or n_events < 2 or n_nonevents < 2:
        return {'Variable': var_name, 'N': n, 'Events': n_events, 'OR': None, 'CI_lower': None, 'CI_upper': None, 'p': None, 'Note': 'Insufficient events'}
    if n_pos < 2:
        return {'Variable': var_name, 'N': n, 'Events': n_events, 'OR': None, 'CI_lower': None, 'CI_upper': None, 'p': None, 'Note': f'Too few positive cases (n={n_pos})'}
    x_a = np.asarray(x).flatten()
    y_a = np.asarray(y).flatten()
    X = np.column_stack([np.ones_like(x_a), x_a])
    try:
        m = FirthLogit().fit(X, y_a)
        or_v = m.or_val
        ci_low = m.or_ci[0]
        ci_high = m.or_ci[1]
        p_v = m.pvalues[1]
        return {'Variable': var_name, 'N': n, 'Events': n_events, 'OR': or_v, 'CI_lower': ci_low, 'CI_upper': ci_high, 'p': p_v, 'Note': 'Firth'}
    except Exception as e:
        return {'Variable': var_name, 'N': n, 'Events': n_events, 'OR': None, 'CI_lower': None, 'CI_upper': None, 'p': None, 'Note': str(e)[:80]}

uni_vars = [
    ('DA_O', 'Multifocal (>1 nodule)'),
    ('HAI_BEN_FINAL', 'Bilateral'),
    ('KT_TREN_1CM', 'Tumor size >1cm'),
    ('KT_TREN_4CM', 'Tumor size >4cm'),
    ('T3_PLUS', 'T stage 3+'),
    ('N_PLUS', 'LN metastasis (N+)'),
    ('N1b', 'Lateral LN (N1b)'),
    ('HACH_TRUNGTAM_5PLUS', 'Central LN >=5'),
    ('CO_HACH_COBEN', 'Lateral LN (+)'),
    ('BRAF_V600E', 'BRAF V600E'),
    ('BRAF_ANY', 'BRAF (any)'),
    ('TERT_MUT', 'TERT mutation'),
    ('TP53_MUT', 'TP53 mutation'),
    ('ANY_MUT', 'Any mutation'),
    ('CO_MUT', 'Co-mutation (>=2 genes)'),
    ('BRAF_TERT', 'BRAF + TERT'),
    ('BRAF_TP53', 'BRAF + TP53'),
]

outcome_arr = np.array(data['RECURRENCE'])
for var_code, var_label in uni_vars:
    var_arr = np.array(data[var_code])
    res = firth_univariate(var_label, var_arr, outcome_arr)
    results.append(res)

print(f'\n{"="*80}')
print('3. UNIVARIATE ANALYSIS — FIRTH LOGISTIC')
print(f'{"="*80}')
print(f'{"Variable":<30} {"N":<6} {"Events":<8} {"OR (Firth)":<12} {"95% CI (Firth)":<25} {"p-value":<12}')
print('-'*80)
for r in results:
    if r['OR'] is not None:
        ci_str = f"{r['CI_lower']:.3f} - {r['CI_upper']:.3f}"
        p_str = f"{r['p']:.4f}" if r['p'] >= 0.0001 else "<0.0001"
        print(f"{r['Variable']:<30} {r['N']:<6} {r['Events']:<8} {r['OR']:<12.3f} {ci_str:<25} {p_str:<12}")
    else:
        print(f"{r['Variable']:<30} {r['N']:<6} {r['Events']:<8} {'NA':<12} {'NA':<25} {r['Note']:<12}")

# ===================== MULTIVARIATE ANALYSIS (FIRTH) =====================
def firth_multivariate(vars_dict, outcome, data_dict, model_name):
    valid_mask = np.ones(len(outcome), dtype=bool)
    for v in vars_dict:
        valid_mask &= ~np.isnan(np.array(data_dict[v]).astype(float))
    y = outcome[valid_mask].astype(float)
    X_list = [np.array(data_dict[v])[valid_mask].astype(float) for v in vars_dict]
    n = len(y)
    n_e = int(y.sum())
    if n < 20 or n_e < 3:
        print(f'\n{model_name}: Insufficient data (n={n}, events={n_e})')
        return
    X_np = np.column_stack([np.ones(n)] + X_list)
    y_np = np.asarray(y).flatten()
    try:
        m = FirthLogit().fit(X_np, y_np)
        print(f'{model_name} (N={n}, events={n_e}):')
        print(f'{"Variable":<25} {"OR (Firth)":<12} {"95% CI (Firth)":<25} {"p-value":<12}')
        print('-'*70)
        or_all = np.exp(m.beta)
        ci_all = np.exp(np.column_stack([m.beta-1.96*m.se, m.beta+1.96*m.se]))
        p_all = m.pvalues
        names = ['Intercept'] + list(vars_dict.keys())
        for i, nm in enumerate(names):
            p_s = f"{p_all[i]:.4f}" if p_all[i] >= 0.0001 else "<0.0001"
            display_name = 'Intercept' if i == 0 else vars_dict.get(nm, nm)
            print(f"{display_name:<25} {or_all[i]:<12.3f} {ci_all[i,0]:.3f} - {ci_all[i,1]:<10.3f}  {p_s:<12}")
        # HL test using Firth predicted probabilities
        try:
            pred = 1 / (1 + np.exp(-X_np @ m.beta))
            if n_e >= 2 and (n - n_e) >= 2:
                hl_stat, hl_pval, hl_df = hosmer_lemeshow_test(y_np, pred, g=10)
                hl_note = ""
                if not np.isnan(hl_pval) and hl_pval < 0.05:
                    hl_note = " (p<0.05: poor fit)"
                if hl_df > 0:
                    print(f'\n  Hosmer-Lemeshow test: chi2={hl_stat:.3f}, df={hl_df}, p={hl_pval:.4f}{hl_note}')
        except Exception as hl_e:
            print(f'\n  Hosmer-Lemeshow test: Cannot compute ({hl_e})')
    except Exception as e:
        print(f'{model_name}: Error - {e}')

print(f'\n{"="*80}')
print('4. MULTIVARIATE ANALYSIS — FIRTH LOGISTIC')
print(f'{"="*80}')

firth_multivariate(
    {'KT_TREN_1CM': 'Size >1cm', 'T3_PLUS': 'T3+', 'N_PLUS': 'N+'},
    outcome_arr, data, 'Model 1: Clinical')
firth_multivariate(
    {'BRAF_V600E': 'BRAF V600E', 'TERT_MUT': 'TERT'},
    outcome_arr, data, 'Model 2: Genetic')
firth_multivariate(
    {'KT_TREN_1CM': 'Size >1cm', 'T3_PLUS': 'T3+', 'N_PLUS': 'N+', 'BRAF_V600E': 'BRAF V600E', 'TERT_MUT': 'TERT'},
    outcome_arr, data, 'Model 3: Clinical + Genetic')
firth_multivariate(
    {'N_PLUS': 'N+', 'CO_MUT': 'Co-mutation'},
    outcome_arr, data, 'Model 4: N+ + Co-mutation')
firth_multivariate(
    {'TERT_MUT': 'TERT', 'N_PLUS': 'N+'},
    outcome_arr, data, 'Model 5: TERT + N+')
firth_multivariate(
    {'TP53_MUT': 'TP53', 'N_PLUS': 'N+'},
    outcome_arr, data, 'Model 6: TP53 + N+')
firth_multivariate(
    {'TERT_MUT': 'TERT', 'TP53_MUT': 'TP53'},
    outcome_arr, data, 'Model 7: TERT + TP53 (mutual)')

# ===================== EXPORT ENCODED DATA =====================
export_data = {k: data[k] for k in ['Mã bệnh nhân', 'RECURRENCE', 'ATA_RISK',
    'DA_O', 'HAI_BEN_FINAL', 'KT_MM', 'KT_TREN_1CM', 'KT_TREN_4CM',
    'T_STAGE', 'T3_PLUS', 'N_STAGE', 'N_PLUS', 'N1b',
    'HACH_TRUNGTAM_5PLUS', 'CO_HACH_COBEN',
    'BRAF_V600E', 'BRAF_OTHER', 'BRAF_ANY', 'TERT_MUT', 'TP53_MUT',
    'ANY_MUT', 'CO_MUT', 'BRAF_TERT', 'BRAF_TP53']}
pd.DataFrame(export_data).to_excel('H:\\Q3 MNGOC\\DATA_MAHOA_ATA2025.xlsx', index=False)
print(f'\n{"="*80}')
print('EXPORTED FILES: DATA_MAHOA_ATA2025.xlsx, RESULTS_ATA2025_FIRTH.txt')
print(f'{"="*80}')

sys.stdout.close()
