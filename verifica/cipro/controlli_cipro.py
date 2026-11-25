"""Controlli analitici 2 cmt e peso della riduzione degli omega Q/V2. Uso: python controlli_cipro.py <repo> casi_cipro.json"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
from pigato_bayes.pk_engine import conc_2cmt, conc_1cmt, metriche_finestra, mappa_bayesiana
CL, V1, Q, V2, D = 10.0, 30.0, 15.0, 40.0, 400.0
# 1. bolo 2 cmt vs A·e^-at + B·e^-bt scritto a mano (Gibaldi)
k10, k12, k21 = CL/V1, Q/V1, Q/V2; s = k10+k12+k21; r = np.sqrt(s*s-4*k10*k21); a, b = (s+r)/2, (s-r)/2
t = np.array([0.25, 1, 4, 12, 23.5])
ref = D/V1*((a-k21)/(a-b)*np.exp(-a*t) + (k21-b)/(a-b)*np.exp(-b*t))
print('1. bolo 2 cmt, scarto max %:', np.max(np.abs(conc_2cmt(CL, V1, Q, V2, [{'ore': 0, 'dose': D, 'tinf': 1e-7}], t)/ref-1))*100)
# 2. stato stazionario: AUCtau = dose/CL
st = [{'ore': 12*i, 'dose': D, 'tinf': 1.0} for i in range(60)]
auc, _, _ = metriche_finestra({'CL': CL, 'V1': V1, 'Q': Q, 'V2': V2}, st, 12*59, 12, 1.0)
print('2. AUC12 ss vs dose/CL, scarto %%: %.4f' % ((auc/(D/CL)-1)*100))
# 3. Q o V2 = 0 -> 1 cmt con V1+V2
print('3. fallback Q=0:', np.allclose(conc_2cmt(CL, V1, 0, V2, st[:2], t), conc_1cmt(CL, V1+V2, st[:2], t)))
# 4. effetto della riduzione degli omega Q/V2 (bayes_page) rispetto agli omega pubblicati
for c in json.load(open(sys.argv[2])):
    if c['omega'] == c['omega_pubblicati']: continue
    p2, _ = mappa_bayesiana(c['prior'], c['omega_pubblicati'], c['storia'], c['livelli'], c['sigma_prop'], c['sigma_add'], c['t_inf'])
    p1 = c['map']
    print('4.', c['id'], len(c['livelli']), 'liv.', ' '.join(f"{k} {100*(p1[k]/p2[k]-1):+.1f}%" for k in ('CL', 'V1', 'Q', 'V2')))
