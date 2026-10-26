"""
Motori PK (1 e 2 compartimenti, infusione a velocità costante, sovrapposizione
lineare delle dosi) e stimatore bayesiano MAP generico.

Una dose è un dict {'ore': t di inizio infusione (h da t0), 'dose': mg,
'tinf': durata infusione in h (facoltativa)}. Un livello è {'ore', 'conc'}.
"""

import numpy as np
from scipy.optimize import minimize


def _tinf(s, t_inf_default):
    v = s.get('tinf')
    return float(v) if v and v > 0 else float(t_inf_default)


def conc_1cmt(CL, V, storia, t, t_inf=0.5):
    """Concentrazione (mg/L) a uno o più tempi t."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    C = np.zeros_like(t)
    if CL <= 0 or V <= 0:
        return C
    ke = CL / V
    for s in storia:
        ti = _tinf(s, t_inf)
        dt = t - s['ore']
        k0 = s['dose'] / ti
        dur = np.clip(dt, 0, ti)
        fase = (k0 / CL) * (1.0 - np.exp(-ke * dur))
        post = np.where(dt > ti, np.exp(-ke * (dt - ti)), 1.0)
        C += np.where(dt > 0, fase * post, 0.0)
    return np.maximum(C, 0.0)


def conc_2cmt(CL, V1, Q, V2, storia, t, t_inf=0.5):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    C = np.zeros_like(t)
    if min(CL, V1, Q, V2) <= 0:
        return conc_1cmt(CL, V1 + max(V2, 0), storia, t, t_inf)
    k10, k12, k21 = CL / V1, Q / V1, Q / V2
    s_ = k10 + k12 + k21
    d_ = max(s_ ** 2 - 4 * k10 * k21, 0.0)
    a = (s_ + np.sqrt(d_)) / 2
    b = (s_ - np.sqrt(d_)) / 2
    if abs(a - b) < 1e-12:
        return conc_1cmt(CL, V1 + V2, storia, t, t_inf)
    A = (k21 - a) / ((b - a) * V1)
    B = (k21 - b) / ((a - b) * V1)
    for s in storia:
        ti = _tinf(s, t_inf)
        dt = t - s['ore']
        k0 = s['dose'] / ti
        dur = np.clip(dt, 0, ti)
        rest = np.clip(dt - ti, 0, None)
        c = k0 * ((A / a) * (1 - np.exp(-a * dur)) * np.exp(-a * rest)
                  + (B / b) * (1 - np.exp(-b * dur)) * np.exp(-b * rest))
        C += np.where(dt > 0, c, 0.0)
    return np.maximum(C, 0.0)


def conc(par, storia, t, t_inf=0.5):
    """par: dict con CL, V1 e opzionalmente Q, V2 (Q=0 → 1 compartimento)."""
    if par.get('Q', 0) and par.get('V2', 0):
        return conc_2cmt(par['CL'], par['V1'], par['Q'], par['V2'], storia, t, t_inf)
    return conc_1cmt(par['CL'], par['V1'], storia, t, t_inf)


def metriche_finestra(par, storia, t_inizio, durata=24.0, t_inf=0.5, n=400):
    """AUC (trapezi), Cmax, Cmin nella finestra [t_inizio, t_inizio+durata]."""
    tg = np.linspace(t_inizio, t_inizio + durata, n)
    cg = conc(par, storia, tg, t_inf)
    auc = float(np.sum((cg[1:] + cg[:-1]) / 2.0 * np.diff(tg)))
    return auc, float(cg.max()), float(cg.min())


def tempo_sotto_soglia(par, storia, soglia, t_partenza, t_max, passo=0.25, t_inf=0.5):
    """Primo tempo ≥ t_partenza in cui C < soglia (dopo eventuali picchi). None se mai."""
    tg = np.arange(t_partenza, t_max + passo, passo)
    cg = conc(par, storia, tg, t_inf)
    idx = np.where(cg < soglia)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    if i == 0:
        return float(tg[0])
    lo, hi = tg[i - 1], tg[i]
    for _ in range(40):
        mid = (lo + hi) / 2
        if conc(par, storia, mid, t_inf)[0] >= soglia:
            lo = mid
        else:
            hi = mid
    return float(hi)


def mappa_bayesiana(prior, omega, storia, livelli, sigma_prop, sigma_add=0.0, t_inf=0.5,
                    limiti=(0.05, 10.0)):
    """
    MAP con prior log-normali.

    prior: dict tipici {'CL','V1','Q','V2'} (Q/V2 assenti o 0 per 1 cmt)
    omega: dict SD su scala log per ciascun parametro; omega 0 = parametro fisso
    sigma_prop: errore residuo proporzionale (SD); sigma_add: additivo (mg/L)
    Con sigma_add = 0 il residuo è log-normale (come nelle versioni precedenti);
    con sigma_add > 0 si usa l'errore combinato su scala lineare.
    """
    stimati = [k for k in ('CL', 'V1', 'Q', 'V2') if omega.get(k, 0) > 0 and prior.get(k, 0) > 0]
    if not livelli or not stimati:
        return dict(prior), "Prior di popolazione (nessun livello inserito)"

    t_obs = np.array([l['ore'] for l in livelli], dtype=float)
    y = np.array([max(l['conc'], 1e-6) for l in livelli], dtype=float)
    lp0 = np.array([np.log(prior[k]) for k in stimati])
    om = np.array([omega[k] for k in stimati])

    def unpack(x):
        p = dict(prior)
        for k, v in zip(stimati, x):
            p[k] = float(np.exp(v))
        return p

    def obj(x):
        p = unpack(x)
        pred = np.maximum(conc(p, storia, t_obs, t_inf), 1e-6)
        pen = np.sum(((x - lp0) / om) ** 2)
        if sigma_add > 0:
            var = (sigma_prop * pred) ** 2 + sigma_add ** 2
            lik = np.sum((y - pred) ** 2 / var + np.log(var))
        else:
            lik = np.sum(((np.log(y) - np.log(pred)) / sigma_prop) ** 2)
        return 0.5 * (pen + lik)

    bnd = [(v + np.log(limiti[0]), v + np.log(limiti[1])) for v in lp0]
    best = None
    for start in (lp0, lp0 + 0.3, lp0 - 0.3):
        r = minimize(obj, start, method='L-BFGS-B', bounds=bnd,
                     options={'maxiter': 3000, 'ftol': 1e-12})
        if best is None or r.fun < best.fun:
            best = r
    p = unpack(best.x)
    n = len(livelli)
    descr = (f"MAP bayesiano su {n} livell{'o' if n == 1 else 'i'}; "
             f"stimati: {', '.join(stimati)}"
             + (f"; fissi al prior: {', '.join(k for k in ('Q', 'V2', 'V1') if k in prior and prior.get(k) and k not in stimati)}"
                if any(k not in stimati and prior.get(k) for k in ('V1', 'Q', 'V2')) else ""))
    return {k: float(v) for k, v in p.items()}, descr


def omega_da_cv(cv_percent):
    """CV% di una distribuzione log-normale → SD su scala log."""
    c = cv_percent / 100.0
    return float(np.sqrt(np.log(1 + c * c)))
