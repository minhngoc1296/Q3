"""
Internal validation: Bootstrap calibration, AUC-ROC, LOOCV
For Firth logistic regression models with rare events.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from firth_logit import FirthLogit

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
})

# ===================== LOAD DATA =====================
f1 = 'H:\\Q3 MNGOC\\Bản sao của DỮ LIỆU Q3 GỐC MINH NGỌC MỚI.xlsx'
f2 = 'H:\\Q3 MNGOC\\KẾT QUẢ GEN Q3 MINH NGỌC 2.xlsx'
df1 = pd.read_excel(f1)
df2 = pd.read_excel(f2)
df = df1.merge(df2, on='Mã bệnh nhân', how='left')

# Recurrence
us = df['Bất thường siêu âm cụ thể'].fillna('').str.contains('Hạch', na=False)
fna = df['Chọc tế bào chẩn đoán'].fillna('').str.contains('Đã chọc \\+', na=False)
reop = df['Xử trí'].fillna('').str.contains('Phẫu thuật', na=False)
fu = df['Bất thường trong lần khám gần nhất sau xử trí'].fillna('').str.contains('Có bất thường', na=False)
recur = (us | fna | reop | fu).astype(int).values

# N stage
def encode_n(v):
    if pd.isna(v): return np.nan
    s = str(v).strip().lower()
    if s in ['0']: return 0
    if s in ['1','1a']: return 1
    if s in ['1b']: return 2
    if s in ['2']: return 3
    return np.nan

n_stage = np.array([encode_n(v) for v in df['Phân loại N']])
n_plus = (n_stage > 0).astype(int)

# Genes
braf = (df['BRAFV600E'] == 2).astype(int).values
tert = (df['TERT'] == 2).astype(int).values
tp53 = (df['TP53'] == 2).astype(int).values

# ===================== MODEL: N+ + BRAF + TERT (Model 5) =====================
# This is the best model for validation (clinical + genetic)
y = recur.astype(float)
valid = ~np.isnan(n_plus.astype(float))
y = y[valid]
n = len(y)
n_p = n_plus[valid].astype(float)
b = braf[valid].astype(float)
t = tert[valid].astype(float)
X = np.column_stack([np.ones(n), n_p, b, t])

# Fit final model
m_final = FirthLogit().fit(X, y)
pred_final = 1 / (1 + np.exp(-X @ m_final.beta))

print(f'Validation model: N+ + BRAF + TERT')
print(f'N={n}, events={int(y.sum())}')

# ===================== CALIBRATION CURVE =====================
print('\n--- Calibration ---')
n_bins = 5
bin_edges = np.linspace(0, 1, n_bins + 1)
bin_centers = []
obs_props = []
pred_means = []
for i in range(n_bins):
    if i == n_bins - 1:
        mask = (pred_final >= bin_edges[i]) & (pred_final <= bin_edges[i + 1])
    else:
        mask = (pred_final >= bin_edges[i]) & (pred_final < bin_edges[i + 1])
    if mask.sum() >= 3:
        bin_centers.append(bin_edges[i] + (bin_edges[i+1] - bin_edges[i]) / 2)
        obs_props.append(y[mask].mean())
        pred_means.append(pred_final[mask].mean())

fig, ax = plt.subplots(figsize=(6, 5.5))
ax.plot([0, 1], [0, 1], '--', color='gray', lw=1, alpha=0.6, label='Perfect calibration')
if bin_centers:
    ax.plot(pred_means, obs_props, 'o-', color='#2c3e50', lw=2, markersize=8,
            markerfacecolor='#e74c3c', markeredgecolor='white', markeredgewidth=1)
    # Add error bars (Wilson CI for proportions)
    for pm, op in zip(pred_means, obs_props):
        if op > 0 and op < 1:
            se = np.sqrt(op * (1 - op) / max(sum(y == 1) * 0.2, 3))
            ax.plot([pm, pm], [max(0, op - se), min(1, op + se)], color='#2c3e50', lw=1)
    # Calibration intercept and slope
    if len(bin_centers) >= 3:
        logit_obs = np.log(np.clip(np.array(obs_props), 0.01, 0.99) / (1 - np.clip(np.array(obs_props), 0.01, 0.99)))
        logit_pred = np.log(np.clip(np.array(pred_means), 0.01, 0.99) / (1 - np.clip(np.array(pred_means), 0.01, 0.99)))
        slope, intercept = np.polyfit(logit_pred, logit_obs, 1)
        ax.set_title(f'Calibration Curve\nN+ + BRAF + TERT model\n'
                     f'Intercept={intercept:.3f}, Slope={slope:.3f}', fontsize=11, weight='bold')
    else:
        ax.set_title('Calibration Curve (N+ + BRAF + TERT)', fontsize=11, weight='bold')
else:
    ax.set_title('Calibration Curve (insufficient bins)', fontsize=11, weight='bold')

ax.set_xlim(0, 0.35)
ax.set_ylim(0, 0.35)
ax.set_xlabel('Predicted risk', fontsize=10)
ax.set_ylabel('Observed proportion', fontsize=10)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig('H:\\Q3 MNGOC\\calibration_plot.png', dpi=300)
plt.close(fig)
print('  Saved: calibration_plot.png')

# ===================== AUC-ROC =====================
print('\n--- AUC-ROC ---')
# Manual ROC curve implementation (no sklearn dependency)
def roc_curve_manual(y_true, y_score, n_thresh=100):
    thresholds = np.linspace(0, 1, n_thresh + 1)
    tpr_list, fpr_list = [], []
    for thresh in thresholds:
        pred = (y_score >= thresh).astype(float)
        tp = ((y_true == 1) & (pred == 1)).sum()
        fn = ((y_true == 1) & (pred == 0)).sum()
        fp = ((y_true == 0) & (pred == 1)).sum()
        tn = ((y_true == 0) & (pred == 0)).sum()
        tpr_list.append(tp / (tp + fn) if (tp + fn) > 0 else 0)
        fpr_list.append(fp / (fp + tn) if (fp + tn) > 0 else 0)
    return np.array(fpr_list), np.array(tpr_list), thresholds

def auc_manual(fpr, tpr):
    return -np.trapezoid(tpr, fpr)

fpr, tpr, _ = roc_curve_manual(y, pred_final)
auc_val = auc_manual(fpr, tpr)

# Bootstrap AUC CI
rng = np.random.default_rng(42)
boot_aucs = []
for _ in range(1000):
    idx = rng.integers(0, n, size=n)
    yb, pb = y[idx], pred_final[idx]
    if yb.sum() >= 2 and (len(yb) - yb.sum()) >= 2:
        try:
            fpr_b, tpr_b, _ = roc_curve_manual(yb, pb)
            boot_aucs.append(auc_manual(fpr_b, tpr_b))
        except:
            continue
boot_aucs = np.array(boot_aucs)
if len(boot_aucs) >= 100:
    auc_ci = np.percentile(boot_aucs, [2.5, 97.5])
else:
    auc_ci = [np.nan, np.nan]

print(f'  AUC = {auc_val:.3f} (95% CI: {auc_ci[0]:.3f} - {auc_ci[1]:.3f})')

fig, ax = plt.subplots(figsize=(6, 5.5))
ax.plot(fpr, tpr, 'b-', lw=2.5, label=f'AUC = {auc_val:.3f} (95% CI: {auc_ci[0]:.3f}-{auc_ci[1]:.3f})')
ax.plot([0, 1], [0, 1], '--', color='gray', lw=1, alpha=0.6)
ax.fill_between(fpr, tpr, alpha=0.15, color='#3498db')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('1 — Specificity', fontsize=10)
ax.set_ylabel('Sensitivity', fontsize=10)
ax.set_title('ROC Curve: N+ + BRAF + TERT', fontsize=11, weight='bold')
ax.legend(fontsize=9, loc='lower right')
ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig('H:\\Q3 MNGOC\\roc_curve.png', dpi=300)
plt.close(fig)
print('  Saved: roc_curve.png')

# ===================== LOOCV =====================
print('\n--- Leave-One-Out Cross-Validation ---')
pred_loocv = np.zeros(n)
for i in range(n):
    X_train = np.delete(X, i, axis=0)
    y_train = np.delete(y, i)
    X_test = X[i:i+1]
    try:
        m_loocv = FirthLogit().fit(X_train, y_train, maxiter=100)
        pred_loocv[i] = 1 / (1 + np.exp(-X_test @ m_loocv.beta))[0]
    except:
        pred_loocv[i] = np.nan

valid_loocv = ~np.isnan(pred_loocv)
if valid_loocv.sum() > 10:
    y_loocv = y[valid_loocv]
    p_loocv = pred_loocv[valid_loocv]
    fpr_l, tpr_l, _ = roc_curve_manual(y_loocv, p_loocv)
    auc_loocv = auc_manual(fpr_l, tpr_l)
    print(f'  LOOCV AUC = {auc_loocv:.3f}')

    # Brier score
    brier = ((y_loocv - p_loocv) ** 2).mean()
    print(f'  Brier score = {brier:.4f}')

    # Log-loss
    eps = 1e-12
    logloss = -(y_loocv * np.log(np.clip(p_loocv, eps, 1)) +
                (1 - y_loocv) * np.log(np.clip(1 - p_loocv, eps, 1))).mean()
    print(f'  Log-loss = {logloss:.4f}')
else:
    print('  LOOCV: insufficient valid predictions')

print('\nValidation complete.')
