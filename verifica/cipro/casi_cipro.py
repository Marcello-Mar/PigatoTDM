"""
Verifica PIGATO TDM, ciprofloxacina: casi e risultati PIGATO.
Riproduce _calcola di bayes_page (omega pubblicati).
L'eGFR entra come covariata già calcolata (le formule renali si verificano a parte).
Uso: python casi_cipro.py <repo> <cartella_uscita>
"""
import sys, json, csv
import numpy as np

sys.path.insert(0, sys.argv[1])
from pigato_bayes import models
from pigato_bayes.pk_engine import conc, mappa_bayesiana, metriche_finestra

OUT = sys.argv[2]
H, R_, F = models.CIPRO_MODELLI


def regime(n, tau, dose, tinf, t0=0.0):
    return [{'ore': t0 + i * tau, 'dose': dose, 'tinf': tinf} for i in range(n)]


C = []
def caso(id, descr, modello, sesso, gg, peso, storia, livelli, egfr=None, crea=None, altezza=None,
         crrt=False, Sc=0.2, Qef=2.0):
    C.append(dict(id=id, descr=descr, modello=modello, sesso=sesso, gg=gg, peso=peso, egfr=egfr, crea=crea,
                  altezza=altezza, crrt=crrt, Sc=Sc, Qef=Qef, storia=storia, livelli=livelli))

# --- Hirt 2021 (predefinito < 18 anni): IIV solo su CL, errore proporzionale log-normale
caso('C01', 'Hirt, 3 mesi 5,5 kg, eGFR 80, 10 mg/kg/8h, 1 valle', H, 'Maschio', 91, 5.5,
     regime(4, 8, 55, 1.0), [{'ore': 31.9, 'conc': 0.35}], egfr=80)
caso('C02', 'Hirt, 2 anni 12 kg, eGFR 150, picco + valle', H, 'Femmina', 730, 12.0,
     regime(4, 12, 120, 1.0), [{'ore': 37.0, 'conc': 3.1}, {'ore': 47.8, 'conc': 0.25}], egfr=150)
caso('C03', 'Hirt, 8 anni 25 kg, eGFR mancante (199), 15 mg/kg/12h', H, 'Maschio', 2922, 25.0,
     regime(3, 12, 375, 1.0), [{'ore': 25.5, 'conc': 4.0}, {'ore': 35.5, 'conc': 0.6}])
caso('C04', 'Hirt, 15 anni 55 kg, eGFR 60, 3 livelli', H, 'Femmina', 5479, 55.0,
     regime(4, 12, 400, 1.0), [{'ore': 37.0, 'conc': 4.2}, {'ore': 40.0, 'conc': 2.4}, {'ore': 47.9, 'conc': 0.9}],
     egfr=60)
caso('C05', 'Hirt, neonato 10 gg 3,2 kg, eGFR 30, q12h', H, 'Maschio', 10, 3.2,
     regime(3, 12, 32, 1.0), [{'ore': 35.0, 'conc': 1.2}], egfr=30)
caso('C06', 'Hirt + CRRT (Sc 0,5 × Qef 1,0), 6 anni 20 kg, eGFR 20', H, 'Femmina', 2191, 20.0,
     regime(3, 12, 200, 1.0), [{'ore': 25.2, 'conc': 3.0}, {'ore': 35.0, 'conc': 1.0}], egfr=20,
     crrt=True, Sc=0.5, Qef=1.0)
# --- Rajagopalan 2003: IIV su 4 parametri, errore combinato (prop 0,27 + add 0,04)
caso('C07', 'Rajagopalan, 6 mesi 7 kg, eGFR 110, 1 livello', R_, 'Maschio', 183, 7.0,
     regime(3, 12, 70, 1.0), [{'ore': 33.0, 'conc': 0.4}], egfr=110)
caso('C08', 'Rajagopalan, 5 anni 18 kg, eGFR 120, 2 livelli', R_, 'Femmina', 1826, 18.0,
     regime(4, 12, 180, 1.0), [{'ore': 37.0, 'conc': 2.9}, {'ore': 47.5, 'conc': 0.3}], egfr=120)
caso('C09', 'Rajagopalan, 12 anni 40 kg, eGFR 50 (correzione renale), 3 livelli', R_, 'Maschio', 4383, 40.0,
     regime(4, 12, 400, 1.0), [{'ore': 37.0, 'conc': 4.0}, {'ore': 42.0, 'conc': 1.9}, {'ore': 47.9, 'conc': 1.0}],
     egfr=50)
caso('C10', 'Rajagopalan, 2 mesi 4,5 kg (fuori range < 100 gg)', R_, 'Maschio', 60, 4.5,
     regime(2, 12, 45, 1.0), [{'ore': 23.5, 'conc': 0.5}], egfr=60)
caso('C11', 'Rajagopalan, valle 0,08 mg/L (domina l\'errore additivo)', R_, 'Femmina', 3000, 26.0,
     regime(3, 12, 260, 1.0), [{'ore': 35.9, 'conc': 0.08}], egfr=130)
caso('C12', 'Rajagopalan, 6 livelli (omega Q e V2 non ridotti)', R_, 'Maschio', 2500, 22.0,
     regime(3, 12, 220, 1.0), [{'ore': 25.0, 'conc': 3.5}, {'ore': 25.5, 'conc': 2.6}, {'ore': 26.5, 'conc': 1.8},
                                {'ore': 28.0, 'conc': 1.2}, {'ore': 31.0, 'conc': 0.7}, {'ore': 35.8, 'conc': 0.35}],
     egfr=115)
# --- Forrest 1993 (predefinito >= 18 anni)
caso('C13', 'Forrest, 30 anni M 75 kg 178 cm, crea 0,9 (CKD-EPI calcolata nel modello)', F, 'Maschio', 10958, 75.0,
     regime(4, 8, 400, 1.0), [{'ore': 25.0, 'conc': 3.6}, {'ore': 31.9, 'conc': 0.45}], crea=0.9, altezza=178)
caso('C14', 'Forrest, 65 anni F 60 kg, crea 1,8, senza altezza (BSA 1,73)', F, 'Femmina', 23741, 60.0,
     regime(3, 12, 400, 1.0), [{'ore': 35.5, 'conc': 1.6}], crea=1.8)
caso('C15', 'Forrest, eGFR 25 + CRRT (Sc 0,5 × Qef 2,0)', F, 'Maschio', 18262, 80.0,
     regime(3, 12, 400, 1.0), [{'ore': 25.0, 'conc': 4.5}, {'ore': 35.9, 'conc': 1.5}], egfr=25, altezza=175,
     crrt=True, Sc=0.5, Qef=2.0)
caso('C16', 'Forrest, campionamento ricco (7 livelli)', F, 'Femmina', 14610, 68.0,
     regime(2, 12, 400, 1.0), [{'ore': 13.0, 'conc': 4.8}, {'ore': 13.5, 'conc': 3.6}, {'ore': 14.0, 'conc': 3.0},
                                {'ore': 15.0, 'conc': 2.3}, {'ore': 17.0, 'conc': 1.6}, {'ore': 20.0, 'conc': 1.0},
                                {'ore': 23.9, 'conc': 0.6}], egfr=95, altezza=165)
caso('C17', 'Forrest applicato a 16 anni (fuori range)', F, 'Maschio', 5844, 62.0,
     regime(3, 12, 400, 1.0), [{'ore': 25.0, 'conc': 3.8}], egfr=110, altezza=172)
# --- Regime irregolare e livelli lontani dal prior
caso('C18', 'Hirt, regime irregolare, infusioni 0,5/2/1 h, 3 livelli, 10 anni 32 kg', H, 'Femmina', 3652, 32.0,
     [{'ore': 0.0, 'dose': 320, 'tinf': 0.5}, {'ore': 7.0, 'dose': 320, 'tinf': 2.0},
      {'ore': 17.0, 'dose': 400, 'tinf': 1.0}, {'ore': 24.5, 'dose': 320, 'tinf': 1.0}],
     [{'ore': 9.5, 'conc': 3.0}, {'ore': 18.5, 'conc': 4.5}, {'ore': 32.0, 'conc': 0.7}], egfr=140)
caso('C19', 'Hirt, valle molto alta (8 mg/L)', H, 'Maschio', 4000, 35.0,
     regime(3, 12, 350, 1.0), [{'ore': 35.8, 'conc': 8.0}], egfr=170)
caso('C20', 'Rajagopalan, picco molto basso (0,5 mg/L)', R_, 'Femmina', 3500, 28.0,
     regime(3, 12, 280, 1.0), [{'ore': 25.2, 'conc': 0.5}], egfr=120)
caso('C21', 'Forrest, valle molto alta (6 mg/L)', F, 'Maschio', 25567, 70.0,
     regime(3, 12, 400, 1.0), [{'ore': 35.9, 'conc': 6.0}], egfr=70, altezza=170)


def omega_app(om, n):
    """Come bayes_page._calcola: omega pubblicati senza modifiche (riduzione di Q e V2 tolta il 24/11/2026)."""
    return dict(om)


def metriche_app(par, storia, t_inf):
    t_last = storia[-1]['ore']
    taus = np.diff([s['ore'] for s in storia])
    tau_med = float(np.median(taus)) if len(taus) else None
    proiez = list(storia)
    if tau_med and tau_med < 24:
        k = 1
        while t_last + k * tau_med < t_last + 24 - 1e-6:
            proiez.append({'ore': t_last + k * tau_med, 'dose': storia[-1]['dose'], 'tinf': storia[-1]['tinf']})
            k += 1
    if tau_med and tau_med > 24:
        auc, cmax, cmin = metriche_finestra(par, storia, t_last, tau_med, t_inf)
        return dict(auc24=auc * 24 / tau_med, cmin=cmin, w_start=t_last, w_len=tau_med, dosi=storia, fattore=24 / tau_med)
    auc, cmax, cmin = metriche_finestra(par, proiez, t_last, 24.0, t_inf)
    return dict(auc24=auc, cmin=cmin, w_start=t_last, w_len=24.0, dosi=proiez, fattore=1.0)


PAR = ('CL', 'V1', 'Q', 'V2')
out = []
for c in C:
    pr = models.prior_cipro(c['sesso'], c['gg'], c['peso'], c['crea'], c['crrt'], c['Sc'], c['Qef'],
                            altezza=c['altezza'], egfr=c['egfr'], modello=c['modello'])
    t_inf = pr['t_inf']
    storia = [{**s, 'tinf': s['tinf'] if s['tinf'] > 0 else t_inf} for s in c['storia']]
    om = omega_app(pr['omega'], len(c['livelli']))
    par, descr = mappa_bayesiana(pr['par'], om, storia, c['livelli'], pr['sigma_prop'], pr['sigma_add'], t_inf)
    par_lib, _ = mappa_bayesiana(pr['par'], om, storia, c['livelli'], pr['sigma_prop'], pr['sigma_add'], t_inf,
                                 limiti=(1e-6, 1e6))
    stimati = [k for k in PAR if om.get(k, 0) > 0]
    lim = [k for k in stimati if not (0.0501 < par[k] / pr['par'][k] < 9.99)]
    usa = par_lib if lim else par
    m = metriche_app(usa, storia, t_inf)
    tgrid = np.round(np.linspace(0, storia[-1]['ore'] + 24, 97), 4).tolist()
    cov = {k: c.get(k) for k in ('sesso', 'gg', 'peso', 'crea', 'egfr', 'altezza', 'dialisi', 'crrt', 'Sc', 'Qef')}
    out.append(dict(id=c['id'], cov=cov, descr=c['descr'], modello=pr['nome'], fuori_range=pr['fuori_range'], note=pr['note'],
                    prior=pr['par'], omega_pubblicati=pr['omega'], omega=om, stimati=stimati,
                    sigma_prop=pr['sigma_prop'], sigma_add=pr['sigma_add'], t_inf=t_inf,
                    storia=storia, livelli=c['livelli'], map=par, map_senza_limiti=par_lib, limiti_attivi=lim,
                    auc24=m['auc24'], cmin=m['cmin'],
                    finestra=dict(start=m['w_start'], len=m['w_len'], fattore=m['fattore'], dosi=m['dosi']),
                    tgrid=tgrid, pop_curve=conc(pr['par'], storia, tgrid, t_inf).tolist()))

json.dump(out, open(f'{OUT}/casi_cipro.json', 'w'), indent=1, ensure_ascii=False)
with open(f'{OUT}/pigato_cipro.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['id', 'stimati', 'CL', 'V1', 'Q', 'V2', 'AUC24', 'Cmin', 'limiti_attivi'])
    for r in out:
        p = r['map_senza_limiti'] if r['limiti_attivi'] else r['map']
        w.writerow([r['id'], ';'.join(r['stimati'])] + [p[k] for k in PAR] + [r['auc24'], r['cmin'], ';'.join(r['limiti_attivi'])])
print(f"{len(out)} casi scritti")
