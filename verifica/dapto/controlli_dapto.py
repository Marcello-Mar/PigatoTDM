"""Controlli analitici e casi limite (daptomicina). Uso: python controlli_dapto.py <repo>"""
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
from pigato_bayes import models
from pigato_bayes.pk_engine import conc_1cmt, metriche_finestra, mappa_bayesiana
import scipy.optimize as so

CL, V, D = 0.5, 5.0, 350.0; k = CL / V
# 1. bolo (infusione di 1e-6 h) vs D/V·e^(-kt)
t = np.array([0.5, 2, 8, 23.9])
c = conc_1cmt(CL, V, [{'ore': 0, 'dose': D, 'tinf': 1e-6}], t)
print('1. bolo, scarto max %:', np.max(np.abs(c / (D / V * np.exp(-k * t)) - 1)) * 100)
# 2. infusione singola vs formula chiusa scritta a mano (durante e dopo)
ti = 0.5; k0 = D / ti
ref = np.where(t <= ti, k0 / CL * (1 - np.exp(-k * t)), k0 / CL * (1 - np.exp(-k * ti)) * np.exp(-k * (t - ti)))
c = conc_1cmt(CL, V, [{'ore': 0, 'dose': D, 'tinf': ti}], t)
print('2. infusione, scarto max %:', np.max(np.abs(c / ref - 1)) * 100)
# 3. stato stazionario: AUC24 = dose/CL
st = [{'ore': 24 * i, 'dose': D, 'tinf': 0.5} for i in range(40)]
auc, _, cmin = metriche_finestra({'CL': CL, 'V1': V}, st, 24 * 39, 24, 0.5)
print('3. AUC24 ss (trapezi 400 pt) vs dose/CL: %.3f vs %.3f, scarto %%: %.4f' % (auc, D / CL, (auc / (D / CL) - 1) * 100))
# 3b. Cmin ss vs formula: k0/CL(1-e^-k ti) e^-k(tau-ti)/(1-e^-k tau)
cmin_ref = k0 / CL * (1 - np.exp(-k * ti)) * np.exp(-k * (24 - ti)) / (1 - np.exp(-k * 24))
print('3b. Cmin ss scarto %%: %.5f' % ((cmin / cmin_ref - 1) * 100))

# 4. Casi limite
pr = models.prior_dapto('Maschio', 3300, 30.0, 0.5, False, 0.2, 2.0, 'Nessuna')
st = [{'ore': 24 * i, 'dose': 270, 'tinf': 0.5} for i in range(3)]
for lab, liv in [('livello = 0 (sotto LOQ inserito come 0)', [{'ore': 71.8, 'conc': 0.0}]),
                 ('livello 0,05 mg/L', [{'ore': 71.8, 'conc': 0.05}]),
                 ('livello prima della prima dose (ore -1)', [{'ore': -1.0, 'conc': 5.0}]),
                 ('livello durante l\'infusione (0,25 h)', [{'ore': 48.25, 'conc': 40.0}])]:
    p, d = mappa_bayesiana(pr['par'], pr['omega'], st, liv, pr['sigma_prop'], 0, pr['t_inf'])
    print('4.', lab, '-> CL %.4f (x%.2f prior), V %.3f (x%.2f)' % (p['CL'], p['CL'] / pr['par']['CL'], p['V1'], p['V1'] / pr['par']['V1']))

# 5. convergenza: l'esito di minimize (success) viene scartato?
orig = so.minimize; esiti = []
def spia(*a, **k):
    r = orig(*a, **k); esiti.append((r.success, r.message)); return r
import pigato_bayes.pk_engine as pe; pe.minimize = spia
pe.mappa_bayesiana(pr['par'], pr['omega'], st, [{'ore': 71.8, 'conc': 0.0}], pr['sigma_prop'], 0, pr['t_inf'])
print('5. esiti ottimizzatore (caso livello 0):', esiti)
