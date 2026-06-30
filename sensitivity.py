"""
Sensitivity Analysis for Thyroid Cancer Recurrence
- Stricter recurrence definition
- Age subgroup (<=55 vs >55)
- Standard logit comparison (for supplementary)
- Missing data report
"""
import sys, os, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
import pandas as pd

from firth_logit import (FirthLogit, firth_univariate, firth_multivariate,
                          hosmer_lemeshow)

OUTPUT = 'H:\\Q3 MNGOC\\RESULTS_SENSITIVITY.txt'
_orig_stdout = sys.stdout
sys.stdout = open(OUTPUT, 'w', encoding='utf-8')

# ===================== LOAD =====================
f1 = 'H:\\Q3 MNGOC\\Bản sao của DỮ LIỆU Q3 GỐC MINH NGỌC MỚI.xlsx'
f2 = 'H:\\Q3 MNGOC\\KẾT QUẢ GEN Q3 MINH NGỌC 2.xlsx'
df1 = pd.read_excel(f1)
df2 = pd.read_excel(f2)
df = df1.merge(df2, on='Mã bệnh nhân', how='left')

D = {}
n_total = len(df)

print('='*70)
print('  SENSITIVITY ANALYSIS')
print('  Thyroid Cancer Recurrence — Firth Logistic')
print('='*70)
print(f'\nTotal patients: {n_total}')

# ===================== BASE OUTCOME =====================
us = df['Bất thường siêu âm cụ thể'].fillna('').str.contains('Hạch', na=False)
fna = df['Chọc tế bào chẩn đoán'].fillna('').str.contains('Đã chọc \\+', na=False)
reop = df['Xử trí'].fillna('').str.contains('Phẫu thuật', na=False)
fu = df['Bất thường trong lần khám gần nhất sau xử trí'].fillna('').str.contains('Có bất thường', na=False)
D['RECURRENCE'] = (us | fna | reop | fu).astype(int)
rec_count = int(D['RECURRENCE'].sum())

# ===================== KEY VARIABLES =====================
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

stage_map = {'I': 1, 'II': 2, 'III': 3, 'IVb': 4}
D['STAGE'] = df['Giai đoạn'].map(stage_map).fillna(1).astype(int)

met_tt = df['Số hạch di căn (trung tâm)'].fillna(df['Số hạch di căn (trung tâm).1'])
D['CENTRAL_LN_5PLUS'] = (met_tt >= 5).fillna(0).astype(int)

met_cb = df['Số hạch di căn (cùng bên)'].fillna(df['Số hạch di căn (cùng bên).1'])
met_db = df['Số hạch di căn (đối bên)'].fillna(df['Số hạch di căn (đối bên).1'])
D['LATERAL_LN_PLUS'] = ((met_cb.fillna(0) + met_db.fillna(0)) > 0).astype(int)

D['BRAF_V600E'] = (df['BRAFV600E'] == 2).astype(int)
D['BRAF_ANY'] = (df['BRAF'] == 2).astype(int)
D['TERT_MUT'] = (df['TERT'] == 2).astype(int)
D['TP53_MUT'] = (df['TP53'] == 2).astype(int)
D['ANY_MUT'] = (D['BRAF_ANY'] | D['TERT_MUT'] | D['TP53_MUT']).astype(int)
D['CO_MUT'] = ((D['BRAF_ANY'] + D['TERT_MUT'] + D['TP53_MUT']) >= 2).astype(int)
D['BRAF_TERT'] = ((D['BRAF_ANY'] == 1) & (D['TERT_MUT'] == 1)).astype(int)
D['BRAF_TP53'] = ((D['BRAF_ANY'] == 1) & (D['TP53_MUT'] == 1)).astype(int)

# ===================== SA1: STRICTER RECURRENCE =====================
print(f'\n{"="*70}')
print('  SA1: STRICTER RECURRENCE DEFINITION (>=2 criteria)')
print(f'{"="*70}')

# Count criteria met per patient
criteria = np.column_stack([us.astype(int), fna.astype(int), reop.astype(int), fu.astype(int)])
strict_recur = (criteria.sum(axis=1) >= 2).astype(int)
print(f'  Strict recurrence (>=2 criteria): {strict_recur.sum()}/{n_total} ({strict_recur.sum()/n_total*100:.1f}%)')
print(f'  Original recurrence: {rec_count}/{n_total} ({rec_count/n_total*100:.1f}%)')

y_strict = strict_recur.astype(float)

sa1_vars = [
    ('GENDER', 'Gender (Male)'),
    ('AGE_55PLUS', 'Age >=55 years'),
    ('MULTIFOCAL', 'Multifocal'),
    ('BILATERAL', 'Bilateral'),
    ('SIZE_OVER_10MM', 'Tumor size >1cm'),
    ('T3_PLUS', 'T stage T3+'),
    ('N_PLUS', 'LN metastasis (N+)'),
    ('CENTRAL_LN_5PLUS', 'Central LN >=5'),
    ('LATERAL_LN_PLUS', 'Lateral LN (+)'),
    ('BRAF_V600E', 'BRAF V600E'),
    ('TERT_MUT', 'TERT promoter'),
    ('TP53_MUT', 'TP53'),
    ('ANY_MUT', 'Any mutation'),
    ('CO_MUT', 'Co-mutation'),
]

print(f'\n{"Variable":<30s} {"N":<6s} {"Events":<8s} {"OR":<12s} {"95% CI":<25s} {"p-value":<12s}')
print('-'*93)
for var_code, var_label in sa1_vars:
    v = np.array(D[var_code])
    valid = ~np.isnan(v.astype(float))
    if valid.sum() == 0:
        continue
    r = firth_univariate(var_label, v[valid], y_strict[valid])
    if r['OR'] is not None:
        ci_s = f"{r['CI_low']:.3f} - {r['CI_high']:.3f}"
        p_s = f"{r['p']:.4f}" if r['p'] >= 0.0001 else "<0.0001"
        print(f"{r['Variable']:<30s} {r['N']:<6d} {r['Events']:<8d} {r['OR']:<12.3f} {ci_s:<25s} {p_s:<12s}")
    else:
        print(f"{r['Variable']:<30s} {r['N']:<6d} {r['Events']:<8d} {'NA':<12s} {'NA':<25s} {r['Note']:<12s}")

# ===================== SA2: AGE SUBGROUP =====================
print(f'\n{"="*70}')
print('  SA2: AGE SUBGROUP ANALYSIS')
print(f'{"="*70}')

age_arr = np.array(D['AGE'])
y = np.array(D['RECURRENCE']).astype(float)

for age_cut, label in [(55, '<=55 vs >55')]:
    young = (age_arr <= 55) & (~np.isnan(age_arr))
    old = (age_arr > 55) & (~np.isnan(age_arr))
    print(f'\n  Age <=55 (N={young.sum()}, events={int(y[young].sum())}):')
    for var_code, var_label in sa1_vars:
        v = np.array(D[var_code])
        valid = young & ~np.isnan(v.astype(float))
        if valid.sum() < 5:
            continue
        r = firth_univariate(var_label, v[valid], y[valid])
        if r['OR'] is not None:
            ci_s = f"{r['CI_low']:.3f} - {r['CI_high']:.3f}"
            p_s = f"{r['p']:.4f}" if r['p'] >= 0.0001 else "<0.0001"
            print(f"    {r['Variable']:<30s} {r['N']:<5d} {r['Events']:<5d} {r['OR']:<10.3f} {ci_s:<22s} {p_s:<10s}")
        else:
            print(f"    {r['Variable']:<30s} {r['N']:<5d} {r['Events']:<5d} {'NA':<10s} {'NA':<22s} {r['Note']:<10s}")

    print(f'\n  Age >55 (N={old.sum()}, events={int(y[old].sum())}):')
    for var_code, var_label in sa1_vars:
        v = np.array(D[var_code])
        valid = old & ~np.isnan(v.astype(float))
        if valid.sum() < 5:
            continue
        r = firth_univariate(var_label, v[valid], y[valid])
        if r['OR'] is not None:
            ci_s = f"{r['CI_low']:.3f} - {r['CI_high']:.3f}"
            p_s = f"{r['p']:.4f}" if r['p'] >= 0.0001 else "<0.0001"
            print(f"    {r['Variable']:<30s} {r['N']:<5d} {r['Events']:<5d} {r['OR']:<10.3f} {ci_s:<22s} {p_s:<10s}")
        else:
            print(f"    {r['Variable']:<30s} {r['N']:<5d} {r['Events']:<5d} {'NA':<10s} {'NA':<22s} {r['Note']:<10s}")

# ===================== SA3: STANDARD LOGIT COMPARISON =====================
print(f'\n{"="*70}')
print('  SA3: FIRTH vs STANDARD MLE LOGISTIC COMPARISON')
print(f'{"="*70}')

try:
    from scipy.optimize import minimize
    from scipy.stats import chi2 as sp_chi2

    def standard_logit(X, y_std):
        """Standard MLE logistic regression via scipy"""
        def neg_log_likelihood(beta):
            linpred = X @ beta
            # Clip to avoid overflow
            linpred = np.clip(linpred, -50, 50)
            p = 1 / (1 + np.exp(-linpred))
            p = np.clip(p, 1e-15, 1 - 1e-15)
            return -np.sum(y_std * np.log(p) + (1 - y_std) * np.log(1 - p))

        def neg_log_likelihood_grad(beta):
            linpred = X @ beta
            linpred = np.clip(linpred, -50, 50)
            p = 1 / (1 + np.exp(-linpred))
            grad = X.T @ (p - y_std)
            return grad

        beta0 = np.zeros(X.shape[1])
        result = minimize(neg_log_likelihood, beta0, method='BFGS',
                          jac=neg_log_likelihood_grad,
                          options={'maxiter': 5000, 'gtol': 1e-8})
        beta_est = result.x
        hess = None
        # Numerical Hessian
        eps2 = 1e-5
        n_params = len(beta_est)
        hess = np.zeros((n_params, n_params))
        for i in range(n_params):
            for j in range(n_params):
                beta_pp, beta_pm, beta_mp, beta_mm = beta_est.copy(), beta_est.copy(), beta_est.copy(), beta_est.copy()
                beta_pp[i] += eps2; beta_pp[j] += eps2
                beta_pm[i] += eps2; beta_pm[j] -= eps2
                beta_mp[i] -= eps2; beta_mp[j] += eps2
                beta_mm[i] -= eps2; beta_mm[j] -= eps2
                ll_pp = -neg_log_likelihood(beta_pp)
                ll_pm = -neg_log_likelihood(beta_pm)
                ll_mp = -neg_log_likelihood(beta_mp)
                ll_mm = -neg_log_likelihood(beta_mm)
                hess[i][j] = (ll_pp - ll_pm - ll_mp + ll_mm) / (4 * eps2 * eps2)
        try:
            cov = np.linalg.inv(-hess)
            se = np.sqrt(np.diag(cov))
        except np.linalg.LinAlgError:
            se = np.full(n_params, np.nan)
        return beta_est, se, result

    # Model 5 for comparison: N+ + BRAF + TERT
    outcome_arr = np.array(D['RECURRENCE'])
    valid = ~np.isnan(np.array(D['N_PLUS']).astype(float))
    y5 = outcome_arr[valid].astype(float)
    n5 = int(y5.sum())
    n_arr = np.array(D['N_PLUS'])[valid].astype(float)
    b_arr = np.array(D['BRAF_V600E'])[valid].astype(float)
    t_arr = np.array(D['TERT_MUT'])[valid].astype(float)

    # Firth
    from firth_logit import FirthLogit
    X = np.column_stack([np.ones(len(y5)), n_arr, b_arr, t_arr])
    m_f = FirthLogit().fit(X, y5)
    firth_sum = m_f.summary(['Intercept', 'N+', 'BRAF', 'TERT'])

    # Standard MLE
    beta_mle, se_mle, opt_result = standard_logit(X, y5)
    names = ['Intercept', 'N+', 'BRAF', 'TERT']

    print(f'\n  {"Parameter":<20s} {"Firth OR":<12s} {"MLE OR":<12s} {"Firth p":<12s} {"MLE p":<12s} {"Ratio":<10s}')
    print('  ' + '-'*78)
    for i, nm in enumerate(names):
        f_or = firth_sum[i]['OR']
        f_se = firth_sum[i].get('SE', np.nan)
        f_p = firth_sum[i]['p']
        mle_or = np.exp(beta_mle[i])
        mle_se = se_mle[i] if not np.isnan(se_mle[i]) else np.nan
        mle_z = beta_mle[i] / mle_se if mle_se and mle_se > 0 else np.nan
        mle_p = 2 * (1 - sp_chi2.cdf(abs(mle_z), 1)) if mle_z and not np.isnan(mle_z) else np.nan
        ratio = f_or / mle_or if mle_or > 0 else np.nan
        f_p_s = f"{f_p:.4f}" if f_p >= 0.0001 else "<0.0001"
        mle_p_s = f"{mle_p:.4f}" if (mle_p and mle_p >= 0.0001) else "<0.0001"
        print(f'  {nm:<20s} {f_or:<12.3f} {mle_or:<12.3f} {f_p_s:<12s} {mle_p_s:<12s} {ratio if not np.isnan(ratio) else 0:<10.3f}')
    print(f'\n  Firth log-likelihood: {-m_f.log_likelihood:.2f}')
    print(f'  MLE converged: {opt_result.success}, iterations: {opt_result.nit}')

except ImportError:
    print('  scipy.optimize not available — skipping standard logit')
except Exception as e:
    print(f'  Error in standard logit: {e}')

# ===================== SA4: MISSING DATA REPORT =====================
print(f'\n{"="*70}')
print('  SA4: MISSING DATA ANALYSIS')
print(f'{"="*70}')

print(f'\n  Total patients: {n_total}')
print(f'\n  {"Variable":<50s} {"Missing":<10s} {"%":<8s} {"Type":<15s} {"Used in":<25s}')
print('  ' + '-'*108)

missing_vars = [
    ('Tuổi tại thời điểm chẩn đoán', D['AGE'], 'continuous'),
    ('3. Giới (Gender)', D['GENDER'], 'binary'),
    ('Phân loại T (T stage)', D['T_STAGE'], 'ordinal'),
    ('Phân loại N (N stage)', D['N_STAGE'], 'ordinal'),
    ('Phân loại M (M stage)', D['M_STAGE'], 'binary'),
    ('Giai đoạn (Stage)', D['STAGE'], 'ordinal'),
    ('Kích thước nhân (Tumor size)', D['TUMOR_SIZE_MM'], 'continuous'),
    ('Bilateral (Vị trí)', D['BILATERAL'], 'binary'),
    ('Multifocal (Tổng số nhân)', D['MULTIFOCAL'], 'binary'),
    ('Số hạch di căn (Central LN)', None, 'count'),
    ('TIRADS trc pt', None, 'ordinal'),
    ('Thể mô bệnh học u (Histology)', None, 'categorical'),
    ('BRAF V600E', D['BRAF_V600E'], 'binary'),
    ('TERT', D['TERT_MUT'], 'binary'),
    ('TP53', D['TP53_MUT'], 'binary'),
]

used_vars = {
    'GENDER': ['Table 1, Table 2'],
    'AGE': ['Table 1'],
    'AGE_55PLUS': ['Table 1, Table 2'],
    'T_STAGE': ['Table 1'],
    'N_PLUS': ['Table 1, Table 2, Model 1-8'],
    'M_STAGE': ['Table 1'],
    'TUMOR_SIZE_MM': ['Table 1'],
    'BILATERAL': ['Table 1, Table 2'],
    'MULTIFOCAL': ['Table 1, Table 2'],
    'BRAF_V600E': ['Table 1, Table 2, Model 2,3,5,8'],
    'TERT_MUT': ['Table 1, Table 2, Model 2,3,5,6,8'],
    'TP53_MUT': ['Table 1, Table 2, Model 2,7'],
}

# Load raw data for missingness
raw_vars = {
    'Tumor size': 'Kích thước nhân (mm) trc pt',
    'TIRADS (pre-op)': 'Phân độ TIRADS trc pt',
    'Histology (tumor)': 'Thể mô bệnh học u',
    'Histology (LN)': 'Mô bệnh học của hạch',
    'Central LN met count': 'Số hạch di căn (trung tâm)',
    'Lateral LN met count': 'Số hạch di căn (cùng bên)',
    'Contralateral LN met': 'Số hạch di căn (đối bên)',
    'FNA cytology': 'Chọc tế bào chẩn đoán',
    'TIRADS (post-op)': 'Phân độ TIRADS sau pt',
}

for raw_name, raw_col in raw_vars.items():
    miss_count = df[raw_col].isna().sum() if raw_col in df.columns else n_total
    miss_pct = miss_count / n_total * 100
    dtype = 'continuous' if 'size' in raw_name.lower() or 'count' in raw_name.lower() else 'categorical'
    used_str = 'TIRADS' if 'TIRADS' in raw_name else ('Histology' if 'Histology' in raw_name else ('FNA' if 'FNA' in raw_name else 'LN counts'))
    print(f'  {raw_name:<50s} {miss_count:<10d} {miss_pct:>5.1f}%   {dtype:<15s} {used_str:<25s}')

for var_code, v_arr, dtype in missing_vars:
    if v_arr is not None:
        count = int(np.isnan(np.array(v_arr).astype(float)).sum()) if hasattr(v_arr, '__len__') else 0
        pct = count / n_total * 100
        used = ', '.join(used_vars.get(var_code, ['N/A']))
        # Truncate long descriptions
        label = var_code[:48]
        print(f'  {label:<50s} {count:<10d} {pct:>5.1f}%   {dtype:<15s} {used:<25s}')

print(f'\n  Note: TIRADS pre-op missingness = 72/114 (63.2%)')
print(f'  Note: FNA cytology missingness = 103/114 (90.4%) — not usable for modeling')
print(f'  Note: Histology missing = 40/114 (35.1%) — not used as covariate due to small N')

# ===================== FINAL SUMMARY =====================
print(f'\n{"="*70}')
print('  SENSITIVITY ANALYSIS COMPLETE')
print(f'{"="*70}')
print(f'\n  SA1: Stricter recurrence (>=2 criteria)')
print(f'  SA2: Age subgroup (<=55 vs >55)')
print(f'  SA3: Firth vs Standard MLE comparison')
print(f'  SA4: Missing data report')

sys.stdout.close()
sys.stdout = _orig_stdout

print(f'\nSensitivity analysis complete!')
print(f'  Results: RESULTS_SENSITIVITY.txt')
