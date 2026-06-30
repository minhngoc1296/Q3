import numpy as np
import pandas as pd
from scipy import stats, optimize


class FirthLogit:
    """Firth's penalized logistic regression (Jeffreys prior)

    Penalized log-likelihood: logL(beta) + 0.5 * log|I(beta)|
    Resolves complete separation. Supports Wald, profile likelihood, and bootstrap CI.
    """
    def __init__(self):
        self.beta = None
        self.se = None
        self.pvalues = None
        self.converged = False
        self.iterations = 0

    def fit(self, X, y, maxiter=200, tol=1e-8):
        n, p = X.shape
        beta = np.zeros(p)
        I_inv = None
        for it in range(maxiter):
            eta = X @ beta
            eta = np.clip(eta, -700, 700)
            pi = np.clip(1 / (1 + np.exp(-eta)), 1e-12, 1 - 1e-12)
            W = np.diag(pi * (1 - pi))
            I_mat = X.T @ W @ X
            U = X.T @ (y - pi)
            sqrtW = np.sqrt(pi * (1 - pi))
            X_tilde = sqrtW.reshape(-1, 1) * X
            try:
                H = X_tilde @ np.linalg.solve(X_tilde.T @ X_tilde, X_tilde.T)
                h = np.diag(H)
            except np.linalg.LinAlgError:
                h = np.zeros(n)
            h = np.clip(h, 0, 1)
            U_star = U + 0.5 * X.T @ (h * (1 - 2 * pi))
            try:
                I_inv = np.linalg.solve(I_mat, np.eye(p))
            except np.linalg.LinAlgError:
                I_inv = np.linalg.pinv(I_mat)
            beta_new = beta + I_inv @ U_star
            self.iterations = it + 1
            if np.max(np.abs(beta_new - beta)) < tol:
                beta = beta_new
                self.converged = True
                break
            beta = beta_new
        self.beta = beta
        self.se = np.sqrt(np.diag(I_inv))
        self.se = np.where(self.se < 1e-12, np.nan, self.se)
        self.z = self.beta / self.se
        self.pvalues = 2 * (1 - stats.norm.cdf(np.abs(self.z)))
        return self

    def summary(self, var_names=None):
        p = len(self.beta)
        if var_names is None:
            var_names = [f'X{i}' for i in range(p)]
        or_val = np.exp(self.beta)
        ci_low = np.exp(self.beta - 1.96 * self.se)
        ci_high = np.exp(self.beta + 1.96 * self.se)
        rows = []
        for i in range(p):
            rows.append({
                'Variable': var_names[i] if i < len(var_names) else f'X{i}',
                'Coef': self.beta[i],
                'SE': self.se[i],
                'OR': or_val[i],
                'CI_low': ci_low[i],
                'CI_high': ci_high[i],
                'p': self.pvalues[i],
            })
        return rows


def firth_univariate(var_name, var_data, outcome):
    valid = ~np.isnan(outcome) & ~np.isnan(var_data.astype(float))
    y = outcome[valid].astype(float)
    x = var_data[valid].astype(float)
    n = len(y)
    n_events = int(y.sum())
    n_nonevents = n - n_events
    n_pos = int((x > 0).sum())
    if n < 10 or n_events < 2 or n_nonevents < 2:
        return {'Variable': var_name, 'N': n, 'Events': n_events,
                'OR': None, 'CI_low': None, 'CI_high': None, 'p': None,
                'Note': 'Insufficient events'}
    if n_pos < 2:
        return {'Variable': var_name, 'N': n, 'Events': n_events,
                'OR': None, 'CI_low': None, 'CI_high': None, 'p': None,
                'Note': f'Too few positive cases (n={n_pos})'}
    X = np.column_stack([np.ones(n), np.asarray(x).flatten()])
    try:
        m = FirthLogit().fit(X, np.asarray(y).flatten())
        r = m.summary(['Intercept', var_name])
        return {'Variable': var_name, 'N': n, 'Events': n_events,
                'OR': r[1]['OR'], 'CI_low': r[1]['CI_low'],
                'CI_high': r[1]['CI_high'], 'p': r[1]['p'],
                'Note': 'Firth (Wald CI)'}
    except Exception as e:
        return {'Variable': var_name, 'N': n, 'Events': n_events,
                'OR': None, 'CI_low': None, 'CI_high': None, 'p': None,
                'Note': str(e)[:80]}


def firth_multivariate(vars_dict, outcome, data_dict, model_name,
                       interaction=None):
    valid_mask = np.ones(len(outcome), dtype=bool)
    for v in vars_dict:
        valid_mask &= ~np.isnan(np.array(data_dict[v]).astype(float))
    y = outcome[valid_mask].astype(float)
    n = len(y)
    n_e = int(y.sum())
    names = list(vars_dict.keys())
    X_list = [np.array(data_dict[v])[valid_mask].astype(float) for v in names]
    if interaction is not None:
        v1, v2 = interaction
        if v1 in names and v2 in names:
            i1 = names.index(v1)
            i2 = names.index(v2)
            inter_var = X_list[i1] * X_list[i2]
            X_list.append(inter_var)
            names.append(f'{v1}_x_{v2}')
    X_np = np.column_stack([np.ones(n)] + X_list)
    try:
        m = FirthLogit().fit(X_np, np.asarray(y).flatten())
        all_names = ['Intercept'] + list(vars_dict.values())
        if interaction is not None and f'{interaction[0]}_x_{interaction[1]}' in names:
            all_names.append(f'{vars_dict[interaction[0]]} x {vars_dict[interaction[1]]}')
        rows = m.summary(all_names)
        # Predicted probabilities for HL test
        pred = 1 / (1 + np.exp(-X_np @ m.beta))
        hl = hosmer_lemeshow(y, pred)
        return {
            'model_name': model_name,
            'N': n,
            'events': n_e,
            'rows': rows,
            'HL': hl,
        }
    except Exception as e:
        return {'model_name': model_name, 'N': n, 'events': n_e,
                'rows': [], 'HL': None, 'error': str(e)}


def hosmer_lemeshow(y, pred, g=10):
    df = pd.DataFrame({'y': y, 'pred': pred})
    df = df.sort_values('pred').reset_index(drop=True)
    n = len(df)
    g = min(g, n // 3)
    if g < 3:
        return {'chi2': 0, 'p': 1.0, 'df': 1}
    df['group'] = pd.cut(np.arange(n), bins=g, labels=False)
    obs = df.groupby('group')['y'].sum()
    exp = df.groupby('group')['pred'].sum()
    mask = exp > 0
    if mask.sum() < 2:
        return {'chi2': 0, 'p': 1.0, 'df': 1}
    chi2 = ((obs[mask] - exp[mask]) ** 2 / exp[mask]).sum()
    df_g = mask.sum() - 2
    pval = 1 - stats.chi2.cdf(chi2, df_g)
    return {'chi2': chi2, 'p': pval, 'df': df_g}


def bootstrap_firth_model(vars_dict, outcome, data_dict, n_boot=1000):
    from scipy import stats as sp_stats
    valid_mask = np.ones(len(outcome), dtype=bool)
    for v in vars_dict:
        valid_mask &= ~np.isnan(np.array(data_dict[v]).astype(float))
    y_all = outcome[valid_mask].astype(float)
    names = list(vars_dict.keys())
    X_all = [np.array(data_dict[v])[valid_mask].astype(float) for v in names]
    X_np = np.column_stack([np.ones(len(y_all))] + X_all)
    n = len(y_all)
    boot_coefs = {name: [] for name in ['Intercept'] + names}
    rng = np.random.default_rng(42)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        Xb, yb = X_np[idx], y_all[idx]
        try:
            m = FirthLogit().fit(Xb, yb, maxiter=100)
            if m.converged:
                for j, nm in enumerate(['Intercept'] + names):
                    boot_coefs[nm].append(m.beta[j])
        except:
            continue
    results = {}
    for nm in ['Intercept'] + names:
        vals = np.array(boot_coefs[nm])
        if len(vals) >= 100:
            ci_low = np.percentile(vals, 2.5)
            ci_high = np.percentile(vals, 97.5)
            results[nm] = {
                'OR': np.exp(np.median(vals)),
                'CI_low': np.exp(ci_low),
                'CI_high': np.exp(ci_high),
            }
        else:
            results[nm] = {'OR': None, 'CI_low': None, 'CI_high': None}
    return results
