"""
Verifica PIGATO TDM, daptomicina: generazione dei casi e risultati PIGATO.
Chiama direttamente pigato_bayes.models e pigato_bayes.pk_engine (nessuna interfaccia).
Uscita: casi_dapto.json (input per R/mapbayr) e pigato_dapto.csv.

Uso: python casi_dapto.py <percorso_repo_PigatoTDM> <cartella_uscita>
"""
import sys, json, csv
import numpy as np

sys.path.insert(0, sys.argv[1])
from pigato_bayes import models
from pigato_bayes.pk_engine import conc, mappa_bayesiana, metriche_finestra

OUT = sys.argv[2]


def regime(n, tau, dose, tinf, t0=0.0):
    return [{'ore': t0 + i * tau, 'dose': dose, 'tinf': tinf} for i in range(n)]


# Ogni caso: covariate PIGATO, storia delle dosi, livelli (ore dalla prima dose, mg/L)
C = []
def caso(id, descr, sesso, gg, peso, crea, storia, livelli, dialisi="Nessuna", crrt=False, Sc=0.2, Qef=2.0):
    C.append(dict(id=id, descr=descr, sesso=sesso, gg=gg, peso=peso, crea=crea, dialisi=dialisi,
                  crrt=crrt, Sc=Sc, Qef=Qef, storia=storia, livelli=livelli))

# --- Wei 2020, bande d'età, funzione renale normale
caso('D01', 'Wei 3-6m, 4 mesi 6 kg, 12 mg/kg/24h, 1 valle', 'Maschio', 120, 6.0, 0.25,
     regime(3, 24, 72, 1.0), [{'ore': 71.9, 'conc': 4.2}])
caso('D02', 'Wei 7-12m, 9 mesi 8,5 kg, picco + valle', 'Femmina', 275, 8.5, 0.28,
     regime(4, 24, 102, 1.0), [{'ore': 73.5, 'conc': 68.0}, {'ore': 95.8, 'conc': 5.1}])
caso('D03', 'Wei 13-24m, 18 mesi 11 kg, 2 livelli', 'Maschio', 548, 11.0, 0.30,
     regime(3, 24, 132, 1.0), [{'ore': 50.0, 'conc': 75.0}, {'ore': 60.0, 'conc': 22.0}])
caso('D04', 'Wei 2-6y, 4 anni 16 kg, 1 picco (default t_inf da modello)', 'Femmina', 1461, 16.0, 0.35,
     regime(5, 24, 192, 0.0), [{'ore': 97.5, 'conc': 95.0}])
caso('D05', 'Wei 7-11y, 9 anni 28 kg, 9 mg/kg, 2 livelli', 'Maschio', 3287, 28.0, 0.45,
     regime(4, 24, 252, 0.5), [{'ore': 73.0, 'conc': 85.0}, {'ore': 95.5, 'conc': 9.5}])
caso('D06', 'Wei 12-17y, 15 anni 55 kg, 7 mg/kg, 2 livelli', 'Femmina', 5479, 55.0, 0.70,
     regime(4, 24, 385, 0.5), [{'ore': 72.75, 'conc': 92.0}, {'ore': 95.9, 'conc': 11.0}])
caso('D07', 'Età < 91 gg (2 mesi, 4,5 kg): estrapolazione 3-6m, fuori range', 'Maschio', 60, 4.5, 0.30,
     regime(2, 24, 54, 1.0), [{'ore': 47.5, 'conc': 3.0}])
caso('D08', 'Bordo di banda: 2190 gg esatti (t_inf passa a 0,5 h)', 'Maschio', 2190, 19.0, 0.40,
     regime(3, 24, 228, 0.0), [{'ore': 49.0, 'conc': 80.0}])
# --- Correzione renale Ye 2022
caso('D09', 'Ye moderata, 5 anni 18 kg', 'Femmina', 1826, 18.0, 0.9,
     regime(3, 24, 216, 0.5), [{'ore': 71.5, 'conc': 14.0}], dialisi="Moderata (GFR 30-60)")
caso('D10', 'Ye ESRD, 10 anni 30 kg, q48h, 2 livelli', 'Maschio', 3652, 30.0, 5.0,
     regime(3, 48, 210, 0.5), [{'ore': 97.0, 'conc': 75.0}, {'ore': 143.5, 'conc': 16.0}],
     dialisi="ESRD (GFR <15)")
caso('D11', 'Ye grave sotto l\'anno (8 mesi 8 kg, fascia 1-2y estrapolata)', 'Maschio', 243, 8.0, 0.8,
     regime(2, 24, 96, 1.0), [{'ore': 47.0, 'conc': 9.0}], dialisi="Grave (GFR 15-30)")
caso('D12', 'Ye lieve, 14 anni 50 kg', 'Femmina', 5113, 50.0, 0.9,
     regime(3, 24, 350, 0.5), [{'ore': 50.0, 'conc': 60.0}], dialisi="Lieve (GFR 60-90)")
# --- Lim 2021 (adolescenti HD/CCPD)
caso('D13', 'Lim HD, 14 anni 45 kg, q48h, 2 livelli', 'Maschio', 5113, 45.0, 6.0,
     regime(3, 48, 315, 0.5), [{'ore': 96.75, 'conc': 95.0}, {'ore': 143.0, 'conc': 30.0}],
     dialisi="Emodialisi (HD)")
caso('D14', 'Lim CCPD, 16 anni 60 kg, 1 livello', 'Femmina', 5844, 60.0, 7.0,
     regime(2, 48, 360, 0.5), [{'ore': 95.0, 'conc': 25.0}], dialisi="Dialisi Peritoneale (CCPD)")
# --- Soraluce 2018 (>= 18 anni)
caso('D15', 'Soraluce, 19 anni M 70 kg, crea 0,8, 8 mg/kg, 2 livelli', 'Maschio', 6940, 70.0, 0.8,
     regime(4, 24, 560, 0.5), [{'ore': 72.5, 'conc': 85.0}, {'ore': 95.8, 'conc': 12.0}])
caso('D16', 'Soraluce, 25 anni F 60 kg, crea 2,0, q48h', 'Femmina', 9131, 60.0, 2.0,
     regime(3, 48, 480, 0.5), [{'ore': 143.0, 'conc': 22.0}])
caso('D17', 'Soraluce, creatinina mancante (ClCr 60 di default)', 'Maschio', 7300, 80.0, None,
     regime(3, 24, 640, 0.5), [{'ore': 71.0, 'conc': 10.0}])
# --- CRRT
caso('D18', 'CRRT, 6 anni 20 kg (Wei + Sc·Qef 0,2×1,0)', 'Maschio', 2191, 20.0, 1.5,
     regime(3, 24, 240, 0.5), [{'ore': 72.5, 'conc': 70.0}, {'ore': 84.0, 'conc': 20.0}],
     crrt=True, Sc=0.2, Qef=1.0)
caso('D19', 'CRRT adulto, 30 anni 75 kg crea 3,5 (Soraluce + 0,15×2,5)', 'Femmina', 10958, 75.0, 3.5,
     regime(3, 24, 600, 0.5), [{'ore': 49.0, 'conc': 60.0}], crrt=True, Sc=0.15, Qef=2.5)
# --- Regimi irregolari, infusioni di durata diversa, molti livelli
caso('D20', 'Regime irregolare, infusioni 0,5/2/1 h, 3 livelli, 7 anni 22 kg', 'Femmina', 2557, 22.0, 0.40,
     [{'ore': 0.0, 'dose': 220, 'tinf': 0.5}, {'ore': 20.0, 'dose': 260, 'tinf': 2.0},
      {'ore': 47.5, 'dose': 200, 'tinf': 1.0}, {'ore': 70.0, 'dose': 264, 'tinf': 0.5}],
     [{'ore': 21.0, 'conc': 40.0}, {'ore': 48.0, 'conc': 30.0}, {'ore': 93.0, 'conc': 8.0}])
caso('D21', 'Molti livelli (5), adulto Soraluce', 'Maschio', 18262, 85.0, 1.1,
     regime(5, 24, 680, 0.5), [{'ore': 0.75, 'conc': 70.0}, {'ore': 6.0, 'conc': 45.0},
                                {'ore': 23.5, 'conc': 15.0}, {'ore': 72.6, 'conc': 95.0},
                                {'ore': 95.9, 'conc': 18.0}])
# --- Livelli molto lontani dal prior
caso('D22', 'Valle molto alta (40 mg/L), 9 anni 30 kg', 'Maschio', 3300, 30.0, 0.5,
     regime(3, 24, 270, 0.5), [{'ore': 71.8, 'conc': 40.0}])
caso('D23', 'Picco molto basso (5 mg/L), 12 anni 40 kg', 'Femmina', 4400, 40.0, 0.6,
     regime(3, 24, 280, 0.5), [{'ore': 49.0, 'conc': 5.0}])
caso('D24', 'Lim, valle molto alta (bassa omega: limiti attivi?)', 'Maschio', 5000, 40.0, 8.0,
     regime(2, 48, 280, 0.5), [{'ore': 95.5, 'conc': 60.0}], dialisi="Emodialisi (HD)")


def metriche_app(par, storia, t_inf):
    """Stessa logica di bayes_page._calcola: AUC24 e Cmin nelle 24 h dopo l'ultima dose."""
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
        auc_tau, cmax, cmin = metriche_finestra(par, storia, t_last, tau_med, t_inf)
        return dict(auc24=auc_tau * 24 / tau_med, cmax=cmax, cmin=cmin, w_start=t_last, w_len=tau_med,
                    dosi_finestra=storia, fattore=24 / tau_med)
    auc, cmax, cmin = metriche_finestra(par, proiez, t_last, 24.0, t_inf)
    return dict(auc24=auc, cmax=cmax, cmin=cmin, w_start=t_last, w_len=24.0, dosi_finestra=proiez, fattore=1.0)


out, righe = [], []
for c in C:
    pr = models.prior_dapto(c['sesso'], c['gg'], c['peso'], c['crea'], c['crrt'], c['Sc'], c['Qef'], c['dialisi'])
    t_inf = pr['t_inf']
    # tinf effettivo per dose (0 → predefinito del modello), esplicitato per R
    storia = [{**s, 'tinf': s['tinf'] if s['tinf'] > 0 else t_inf} for s in c['storia']]
    par, descr = mappa_bayesiana(pr['par'], pr['omega'], c['storia'], c['livelli'],
                                 pr['sigma_prop'], pr['sigma_add'], t_inf)
    par_lib, _ = mappa_bayesiana(pr['par'], pr['omega'], c['storia'], c['livelli'],
                                 pr['sigma_prop'], pr['sigma_add'], t_inf, limiti=(1e-6, 1e6))
    m_post = metriche_app(par, storia, t_inf)
    m_pop = metriche_app(pr['par'], storia, t_inf)
    m_lib = metriche_app(par_lib, storia, t_inf)
    pred_liv = conc(par, storia, [l['ore'] for l in c['livelli']], t_inf).tolist()
    tgrid = np.round(np.linspace(0, storia[-1]['ore'] + 24, 97), 4).tolist()
    pop_curve = conc(pr['par'], storia, tgrid, t_inf).tolist()
    # limiti di PIGATO attivi?
    lim = []
    for k in ('CL', 'V1'):
        r = par[k] / pr['par'][k]
        if r < 0.0501 or r > 9.99:
            lim.append(k)
    cov = {k: c.get(k) for k in ('sesso', 'gg', 'peso', 'crea', 'egfr', 'altezza', 'dialisi', 'crrt', 'Sc', 'Qef')}
    rec = dict(id=c['id'], cov=cov, descr=c['descr'], modello=pr['nome'], fuori_range=pr['fuori_range'],
               note=pr['note'], prior=pr['par'], omega=pr['omega'], sigma_prop=pr['sigma_prop'],
               sigma_add=pr['sigma_add'], t_inf=t_inf, storia=storia, livelli=c['livelli'],
               map=par, map_senza_limiti=par_lib, pred_livelli=pred_liv, limiti_attivi=lim,
               auc24=m_post['auc24'], cmin=m_post['cmin'], cmax=m_post['cmax'],
               auc24_senza_limiti=m_lib['auc24'], cmin_senza_limiti=m_lib['cmin'], auc24_pop=m_pop['auc24'], cmin_pop=m_pop['cmin'],
               finestra=dict(start=m_post['w_start'], len=m_post['w_len'], fattore=m_post['fattore'],
                             dosi=m_post['dosi_finestra']),
               tgrid=tgrid, pop_curve=pop_curve)
    out.append(rec)
    righe.append([c['id'], par['CL'], par['V1'], m_post['auc24'], m_post['cmin'], ';'.join(lim)])

json.dump(out, open(f'{OUT}/casi_dapto.json', 'w'), indent=1, ensure_ascii=False)
with open(f'{OUT}/pigato_dapto.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['id', 'CL', 'V', 'AUC24', 'Cmin', 'limiti_attivi'])
    w.writerows(righe)
print(f"{len(out)} casi scritti")
