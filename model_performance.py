"""
Model Performance Comparison + Diagnostic Metrics
- AUC, Brier, LogLoss, AIC for all 8 models
- Sensitivity/Specificity/PPV/NPV for M5 at various thresholds
"""
import sys, os, re
import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from firth_logit import FirthLogit

OUTPUT = 'H:\\Q3 MNGOC\\RESULTS_PERFORMANCE.txt'
_orig_stdout = sys.stdout
sys.stdout = open(OUTPUT, 'w', encoding='utf-8')

# ===================== LOAD =====================
f1 = 'H:\\Q3 MNGOC\\Bản sao của DỮ LIỆU Q3 GỐC MINH NGỌC MỚI.xlsx'
f2 = 'H:\\Q3 MNGOC\\KẾT QUẢ GEN Q3 MINH NGỌC 2.xlsx'
df = pd.read_excel(f1).merge(pd.read_excel(f2), on='Mã bệnh nhân', how='left')

D = {}
n_total = len(df)

# Outcome
us = df['Bất thường siêu âm cụ thể'].fillna('').str.contains('Hạch', na=False)
fna = df['Chọc tế bào chẩn đoán'].fillna('').str.contains('Đã chọc \\+', na=False)
reop = df['Xử trí'].fillna('').str.contains('Phẫu thuật', na=False)
fu = df['Bất thường trong lần khám gần nhất sau xử trí'].fillna('').str.contains('Có bất thường', na=False)
D['RECURRENCE'] = (us | fna | reop | fu).astype(int)
y = np.array(D['RECURRENCE']).astype(float)

# Variables
gender_map = {'Nữ': 0, 'Nam': 1}
D['GENDER'] = df['3. Giới'].map(gender_map).fillna(0).astype(int)
age = pd.to_numeric(df['Tuổi tại thời điểm chẩn đoán'], errors='coerce')
D['AGE'] = age
D['AGE_55PLUS'] = (age >= 55).astype(int)
nhan_tc = pd.to_numeric(df['Tổng số nhân trc pt'], errors='coerce').fillna(0)
D['MULTIFOCAL'] = (nhan_tc > 1).astype(int)

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

t_stage = df['Phân loại T'].apply(encode_t)
n_stage = df['Phân loại N'].apply(encode_n)
m_stage = pd.to_numeric(df['Phân loại M'], errors='coerce').fillna(0)
D['T_STAGE'] = t_stage
D['T3_PLUS'] = (t_stage >= 3).astype(int)
D['N_STAGE'] = n_stage
D['N_PLUS'] = (n_stage > 0).astype(int)
D['N1B'] = (n_stage >= 2).astype(int)
D['M_STAGE'] = m_stage.astype(int)

met_tt = df['Số hạch di căn (trung tâm)'].fillna(df['Số hạch di căn (trung tâm).1'])
D['CENTRAL_LN_5PLUS'] = (met_tt >= 5).fillna(0).astype(int)
met_cb = df['Số hạch di căn (cùng bên)'].fillna(df['Số hạch di căn (cùng bên).1'])
met_db = df['Số hạch di căn (đối bên)'].fillna(df['Số hạch di căn (đối bên).1'])
D['LATERAL_LN_PLUS'] = ((met_cb.fillna(0) + met_db.fillna(0)) > 0).astype(int)
stage_map = {'I': 1, 'II': 2, 'III': 3, 'IVb': 4}
D['STAGE'] = df['Giai đoạn'].map(stage_map).fillna(1).astype(int)

D['BRAF_V600E'] = (df['BRAFV600E'] == 2).astype(int)
D['BRAF_ANY'] = (df['BRAF'] == 2).astype(int)
D['TERT_MUT'] = (df['TERT'] == 2).astype(int)
D['TP53_MUT'] = (df['TP53'] == 2).astype(int)
D['ANY_MUT'] = (D['BRAF_ANY'] | D['TERT_MUT'] | D['TP53_MUT']).astype(int)
D['CO_MUT'] = ((D['BRAF_ANY'] + D['TERT_MUT'] + D['TP53_MUT']) >= 2).astype(int)

def parse_size_mm(val):
    if pd.isna(val): return np.nan
    s = str(val).strip().replace(',', '.')
    s = re.sub(r'\([^)]*\)', '', s)
    nums = re.findall(r'[\d.]+', s)
    nums = [n for n in nums if n and n != '.']
    if not nums: return np.nan
    size = max(float(n) for n in nums)
    if 'cm' in s.lower(): size *= 10
    if 0 < size < 1: size *= 10
    return size

KT = df['Kích thước nhân (mm) trc pt'].apply(parse_size_mm)
D['SIZE_OVER_10MM'] = (KT > 10).astype(int)

# ===================== DEFINE MODELS =====================
models = [
    ('M1: Clinical', ['N_PLUS', 'T3_PLUS', 'SIZE_OVER_10MM']),
    ('M2: Genetic', ['BRAF_V600E', 'TERT_MUT', 'TP53_MUT']),
    ('M3: Combined', ['N_PLUS', 'T3_PLUS', 'SIZE_OVER_10MM', 'BRAF_V600E', 'TERT_MUT']),
    ('M4: N+ + Co-mut', ['N_PLUS', 'CO_MUT']),
    ('M5: N+ + BRAF + TERT', ['N_PLUS', 'BRAF_V600E', 'TERT_MUT']),
    ('M6: N+ + TERT', ['N_PLUS', 'TERT_MUT']),
    ('M7: N+ + TP53', ['N_PLUS', 'TP53_MUT']),
    ('M8: N+ + BRAF + TERT + INT', ['N_PLUS', 'BRAF_V600E', 'TERT_MUT', 'interaction']),
]

def bootstrap_auc(y_true, y_pred, n_boot=1000):
    rng = np.random.RandomState(42)
    n = len(y_true)
    aucs = np.zeros(n_boot)
    for i in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        if len(np.unique(y_true[idx])) < 2:
            aucs[i] = np.nan
            continue
        aucs[i] = compute_auc(y_true[idx], y_pred[idx])
    aucs = aucs[~np.isnan(aucs)]
    if len(aucs) < 100:
        return np.mean(aucs), np.nan, np.nan
    return np.mean(aucs), np.percentile(aucs, 2.5), np.percentile(aucs, 97.5)

def compute_auc(y_true, y_pred):
    n = len(y_true)
    idx = np.argsort(y_pred)[::-1]
    y_sorted = y_true[idx]
    tpr = np.zeros(n + 2)
    fpr = np.zeros(n + 2)
    tp = 0; fp = 0
    total_p = int(y_true.sum())
    total_n = n - total_p
    if total_p == 0 or total_n == 0:
        return 0.5
    for i in range(n):
        if y_sorted[i] == 1:
            tp += 1
        else:
            fp += 1
        tpr[i+1] = tp / total_p
        fpr[i+1] = fp / total_n
    tpr[-1] = 1; fpr[-1] = 1
    return abs(np.trapezoid(tpr, fpr))

print('='*70)
print('  MODEL PERFORMANCE COMPARISON')
print('='*70)

results = []
for mname, var_names in models:
    # Build predictor matrix
    valid = np.ones(n_total, dtype=bool)
    for v in var_names:
        if v == 'interaction':
            continue
        arr = np.array(D[v]).astype(float)
        valid &= ~np.isnan(arr)

    y_m = y[valid]
    n_m = len(y_m)
    e_m = int(y_m.sum())

    X_list = [np.ones(n_m)]
    for v in var_names:
        if v == 'interaction':
            arr_b = np.array(D['BRAF_V600E'])[valid].astype(float)
            arr_t = np.array(D['TERT_MUT'])[valid].astype(float)
            X_list.append(arr_b * arr_t)
        else:
            X_list.append(np.array(D[v])[valid].astype(float))
    X = np.column_stack(X_list)

    try:
        print(f'  Fitting {mname}...', file=sys.stderr)
        m = FirthLogit().fit(X, y_m)
        pred = 1.0 / (1.0 + np.exp(-(X @ m.beta)))
        pred_clip = np.clip(pred, 1e-15, 1 - 1e-15)
        ll = np.sum(y_m * np.log(pred_clip) + (1 - y_m) * np.log(1 - pred_clip))
        aic = 2 * len(m.beta) - 2 * ll
        auc = compute_auc(y_m, pred)
        print(f'  AUC={auc:.3f}, bootstrapping...', file=sys.stderr)
        auc_mean, auc_lo, auc_hi = bootstrap_auc(y_m, pred, 200)
        print(f'  Bootstrap done.', file=sys.stderr)
        brier = np.mean((y_m - pred) ** 2)
        log_loss_val = -np.mean(y_m * np.log(pred_clip) + (1 - y_m) * np.log(1 - pred_clip))
        results.append({
            'model': mname, 'N': n_m, 'events': e_m, 'vars': len(var_names),
            'AUC': auc, 'AUC_lo': auc_lo, 'AUC_hi': auc_hi,
            'Brier': brier, 'LogLoss': log_loss_val, 'AIC': aic,
            'preds': pred, 'y': y_m, 'beta': m.beta
        })
    except Exception as e:
        import traceback
        print(f'  {mname}: ERROR - {e}', file=sys.stderr)
        traceback.print_exc()
        results.append({'model': mname, 'N': n_m, 'events': e_m, 'vars': len(var_names),
                        'AUC': np.nan, 'AUC_lo': np.nan, 'AUC_hi': np.nan,
                        'Brier': np.nan, 'LogLoss': np.nan, 'AIC': np.nan,
                        'preds': None, 'y': None, 'beta': None})

# ===================== PRINT COMPARISON TABLE =====================
print(f'\n{"Model":<30s} {"N":<5s} {"Events":<7s} {"k":<4s} {"AUC":<8s} {"95% CI (AUC)":<22s} {"Brier":<8s} {"LogLoss":<8s} {"AIC":<8s}')
print('-'*100)
for r in results:
    if np.isnan(r['AUC']):
        print(f'{r["model"]:<30s} {r["N"]:<5d} {r["events"]:<7d} {r["vars"]:<4d} {"FAILED":<38s}')
    else:
        ci_s = f'{r["AUC_lo"]:.3f} - {r["AUC_hi"]:.3f}' if not np.isnan(r['AUC_lo']) else 'NA'
        aic_s = f'{r["AIC"]:.1f}' if not np.isnan(r['AIC']) else 'NA'
        print(f'{r["model"]:<30s} {r["N"]:<5d} {r["events"]:<7d} {r["vars"]:<4d} '
              f'{r["AUC"]:<8.3f} {ci_s:<22s} {r["Brier"]:<8.4f} {r["LogLoss"]:<8.4f} {aic_s:<8s}')

# ===================== DIAGNOSTIC METRICS FOR M5 =====================
print(f'\n{"="*70}')
print('  DIAGNOSTIC METRICS — Model 5 (N+ + BRAF + TERT)')
print(f'{"="*70}')

m5 = [r for r in results if 'M5' in r['model']][0]
if m5['preds'] is not None:
    y5 = m5['y']
    p5 = m5['preds']

    thresholds = [0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20]
    print(f'\n{"Threshold":<12s} {"TP":<5s} {"FP":<5s} {"FN":<5s} {"TN":<5s} '
          f'{"Sensitivity":<12s} {"Specificity":<12s} {"PPV":<12s} {"NPV":<12s} {"Accuracy":<12s} {"Youden":<10s}')
    print('-'*95)

    best_youden = -1
    best_row = None
    for th in thresholds:
        pred_class = (p5 >= th).astype(int)
        tp = int(np.sum((y5 == 1) & (pred_class == 1)))
        fp = int(np.sum((y5 == 0) & (pred_class == 1)))
        fn = int(np.sum((y5 == 1) & (pred_class == 0)))
        tn = int(np.sum((y5 == 0) & (pred_class == 0)))
        sens = tp / (tp + fn) if (tp + fn) > 0 else 0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0
        acc = (tp + tn) / len(y5)
        youden = sens + spec - 1
        print(f'{th:<12.3f} {tp:<5d} {fp:<5d} {fn:<5d} {tn:<5d} '
              f'{sens:<12.3f} {spec:<12.3f} {ppv:<12.3f} {npv:<12.3f} {acc:<12.3f} {youden:<10.4f}')
        if youden > best_youden:
            best_youden = youden
            best_row = (th, tp, fp, fn, tn, sens, spec, ppv, npv, acc, youden)

    print(f'\n  Best threshold (Youden index = {best_youden:.4f}):')
    th, tp, fp, fn, tn, sens, spec, ppv, npv, acc, youden = best_row
    print(f'    Threshold = {th:.3f}')
    print(f'    TP={tp}, FP={fp}, FN={fn}, TN={tn}')
    print(f'    Sensitivity = {sens:.3f} ({sens*100:.1f}%)')
    print(f'    Specificity = {spec:.3f} ({spec*100:.1f}%)')
    print(f'    PPV = {ppv:.3f} ({ppv*100:.1f}%)')
    print(f'    NPV = {npv:.3f} ({npv*100:.1f}%)')
    print(f'    Accuracy = {acc:.3f} ({acc*100:.1f}%)')

# ===================== BRIEF SUMMARY =====================
print(f'\n{"="*70}')
print('  SUMMARY')
print(f'{"="*70}')
best_auc = max((r for r in results if not np.isnan(r['AUC'])), key=lambda r: r['AUC'])
print(f'  Best AUC: {best_auc["model"]} = {best_auc["AUC"]:.3f} '
      f'(95% CI: {best_auc["AUC_lo"]:.3f} - {best_auc["AUC_hi"]:.3f})')
# Find most parsimonious with AUC close to best
for r in results:
    if r['AUC'] >= best_auc['AUC'] - 0.02 and r['vars'] < best_auc['vars']:
        print(f'  Parsimonious alternative: {r["model"]} (AUC={r["AUC"]:.3f}, k={r["vars"]} vars)')
        break

sys.stdout.close()
sys.stdout = _orig_stdout
print(f'\nPerformance analysis complete: {OUTPUT}')
