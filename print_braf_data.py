import pandas as pd, numpy as np
from scipy import stats

f1 = 'H:\\Q3 MNGOC\\Bản sao của DỮ LIỆU Q3 GỐC MINH NGỌC MỚI.xlsx'
f2 = 'H:\\Q3 MNGOC\\KẾT QUẢ GEN Q3 MINH NGỌC 2.xlsx'
df1 = pd.read_excel(f1); df2 = pd.read_excel(f2)
df = df1.merge(df2, on='Mã bệnh nhân', how='left')

us = df['Bất thường siêu âm cụ thể'].fillna('').str.contains('Hạch', na=False)
fna = df['Chọc tế bào chẩn đoán'].fillna('').str.contains('Đã chọc \\+', na=False)
reop = df['Xử trí'].fillna('').str.contains('Phẫu thuật', na=False)
fu = df['Bất thường trong lần khám gần nhất sau xử trí'].fillna('').str.contains('Có bất thường', na=False)
rec = (us|fna|reop|fu).astype(int)

braf = (df['BRAFV600E'] == 2).astype(int)
tert = (df['TERT'] == 2).astype(int)
tp53 = (df['TP53'] == 2).astype(int)
n_stage = df['Phân loại N'].apply(lambda v: 0 if pd.isna(v) else (1 if str(v).strip().lower() not in ['0','x'] else 0))
n_plus = n_stage.astype(int)

def encode_t(v):
    if pd.isna(v): return np.nan
    v = str(v).strip().lower()
    if v in ['x','0']: return 0
    if v in ['1','1a','1b']: return 1
    if v in ['2']: return 2
    if v in ['3','3a','3b']: return 3
    if v in ['4','4a','4b','4c']: return 4
    return np.nan
t3 = df['Phân loại T'].apply(lambda v: 1 if pd.notna(v) and encode_t(v) >= 3 else 0)
nhan = pd.to_numeric(df['Tổng số nhân trc pt'], errors='coerce').fillna(0)
da_o = (nhan > 1).astype(int)

class FirthLogit:
    def fit(self, X, y, maxiter=200, tol=1e-8):
        n, p = X.shape
        beta = np.zeros(p)
        for it in range(maxiter):
            eta = X @ beta
            pi = np.clip(1/(1+np.exp(-eta)), 1e-12, 1-1e-12)
            W = np.diag(pi*(1-pi))
            I = X.T @ W @ X; U = X.T @ (y-pi)
            sqrtW = np.sqrt(pi*(1-pi))
            X_tilde = sqrtW.reshape(-1,1) * X
            try: H = X_tilde @ np.linalg.solve(X_tilde.T @ X_tilde, X_tilde.T); h = np.diag(H)
            except: h = np.zeros(n)
            h = np.clip(h, 0, 1)
            U_star = U + 0.5 * X.T @ (h*(1-2*pi))
            try: I_inv = np.linalg.solve(I, np.eye(p))
            except: I_inv = np.linalg.pinv(I)
            beta_new = beta + I_inv @ U_star
            if np.max(np.abs(beta_new-beta)) < tol: break
            beta = beta_new
        self.beta = beta; self.se = np.sqrt(np.diag(I_inv))
        self.z = self.beta / self.se
        self.pvalues = 2*(1-stats.norm.cdf(np.abs(self.z)))
        ci = np.column_stack([self.beta-1.96*self.se, self.beta+1.96*self.se])
        self.or_val = np.exp(self.beta[1]) if p>1 else np.exp(self.beta[0])
        self.or_ci = np.exp(ci[1]) if p>1 else np.exp(ci[0])
        return self

print('=== FOREST PLOT: BRAF V600E OR theo phan nhom ===')
print(f'{"Phan nhom":<30s} {"N/evt":<8s} {"OR":<8s} {"95% CI":<22s} {"p":<8s}')
print('-'*76)

for label, mask in [
    ('1. Toan bo', slice(None)),
    ('2. N+', n_plus==1),
    ('3. N0', n_plus==0),
    ('4. Da o (>1 nhan)', da_o==1),
    ('5. Don o (1 nhan)', da_o==0),
    ('6. T1-T2', t3==0),
    ('7. T3+', t3==1),
]:
    if isinstance(mask, slice):
        y, x = rec.values, braf.values
    else:
        y, x = rec[mask].values, braf[mask].values
    valid = ~np.isnan(x.astype(float))
    y, x = y[valid], x[valid]
    n = len(y); ne = int(y.sum())
    if ne < 2 or (n-ne) < 2:
        print(f'{label:<30s} {n:>3d}/{ne:<3d} {"--":<8s} {"--":<22s} {"--":<8s}')
        continue
    X = np.column_stack([np.ones(n), x.astype(float)])
    try:
        m = FirthLogit().fit(X, y.astype(float))
        or_v = m.or_val; ci_l = m.or_ci[0]; ci_h = m.or_ci[1]; pv = m.pvalues[1]
        print(f'{label:<30s} {n:>3d}/{ne:<3d} {or_v:<8.3f} {ci_l:.3f}-{ci_h:<13.3f} {pv:<8.4f}')
    except Exception as e:
        print(f'{label:<30s} {n:>3d}/{ne:<3d} {"ERR":<8s} {str(e)[:30]}')

print()
print('=== STACKED BAR: Tai phat theo BRAF x N ===')
print(f'{"Nhom":<20s} {"Tong":<6s} {"Tai phat":<10s} {"Ty le":<8s}')
print('-'*44)
for label, m in [
    ('BRAF- / N0', (braf==0)&(n_plus==0)),
    ('BRAF+ / N0', (braf==1)&(n_plus==0)),
    ('BRAF+ / N+', (braf==1)&(n_plus==1)),
    ('BRAF- / N+', (braf==0)&(n_plus==1)),
]:
    t = int(m.sum()); r = int(rec[m].sum())
    print(f'{label:<20s} {t:<6d} {r:<10d} {r/t*100 if t>0 else 0:<6.1f}%')

print()
print('=== DONG MAC DOT BIEN (so BN co ca 2 dot bien) ===')
print(f'{"":>10s} {"BRAF":>8s} {"TERT":>8s} {"TP53":>8s}')
print('-'*34)
for gname, col in [('BRAF', braf), ('TERT', tert), ('TP53', tp53)]:
    b = ((col==1) & (braf==1)).sum()
    t = ((col==1) & (tert==1)).sum()
    p = ((col==1) & (tp53==1)).sum()
    print(f'{gname:>10s} {b:>8d} {t:>8d} {p:>8d}')

print()
print('=== PHAN BO BRAF THEO N ===')
print(f'{"":>10s} {"BRAF-":>8s} {"BRAF+":>8s} {"Tong":>8s}')
b_n0 = ((braf==0)&(n_plus==0)).sum()
b_np = ((braf==0)&(n_plus==1)).sum()
v_n0 = ((braf==1)&(n_plus==0)).sum()
v_np = ((braf==1)&(n_plus==1)).sum()
print(f'{"N0":>10s} {b_n0:>8d} {v_n0:>8d} {b_n0+v_n0:>8d}')
print(f'{"N+":>10s} {b_np:>8d} {v_np:>8d} {b_np+v_np:>8d}')
print(f'{"Tong":>10s} {b_n0+b_np:>8d} {v_n0+v_np:>8d} {114:>8d}')
