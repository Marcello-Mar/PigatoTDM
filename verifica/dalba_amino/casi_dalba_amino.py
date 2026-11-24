"""
Verifica PIGATO TDM, dalbavancina e aminoglicosidi: casi e risultati PIGATO.
Riproduce _calcola di bayes_page (omega pubblicati; AUC24/Cmin
nella finestra dell'app) e, per la dalbavancina, il tempo sotto soglia del pannello forecast.
Uso: python casi_dalba_amino.py <repo> <cartella_uscita>
"""
import sys, json, csv
import numpy as np

sys.path.insert(0, sys.argv[1])
from pigato_bayes import models
from pigato_bayes.pk_engine import conc, mappa_bayesiana, metriche_finestra, tempo_sotto_soglia

OUT = sys.argv[2]
PAR = ('CL', 'V1', 'Q', 'V2')
MA, CO, BA = list(models.DALBA_MODELLI)
SOGLIA = models.dalba_soglia_auc(111.1, 0.25, models.DALBA_FU)   # 16,53 mg/L
G = 24.0


def regime(n, tau, dose, tinf, t0=0.0):
    return [{'ore': t0 + i * tau, 'dose': dose, 'tinf': tinf} for i in range(n)]


def liv(*coppie):
    return [{'ore': t, 'conc': c} for t, c in coppie]


DALBA = [
    ('B01', 'Mariani, 10 anni 30 kg, 18 mg/kg dose singola, livelli g7 e g14', MA, 30, regime(1, 0, 540, 0.5),
     liv((7 * G, 35.0), (14 * G, 18.0))),
    ('B02', 'Mariani, 3 anni 14 kg, 22,5 mg/kg, 1 livello g7', MA, 14, regime(1, 0, 315, 0.5), liv((7 * G, 30.0))),
    ('B03', 'Mariani, 15 anni 60 kg, 1500 mg g1+g8, 5 livelli fino a g56', MA, 60, regime(2, 7 * G, 1500, 0.5),
     liv((7 * G, 45.0), (14 * G, 60.0), (28 * G, 30.0), (42 * G, 17.0), (56 * G, 9.0))),
    ('B04', 'Mariani, 6 mesi 7 kg, 22,5 mg/kg', MA, 7, regime(1, 0, 158, 0.5), liv((168.0, 28.0))),
    ('B05', 'Mariani, 17 anni 80 kg, 1000 mg ogni 14 gg ×4, 6 livelli (terapia lunga)', MA, 80,
     regime(4, 14 * G, 1000, 0.5), liv((7 * G, 30.0), (14 * G, 20.0), (28 * G, 28.0), (42 * G, 30.0),
                                       (56 * G, 25.0), (70 * G, 14.0))),
    ('B06', 'Cojutti, adulto, 1500 mg g1+g8, livelli g14 e g28', CO, 75, regime(2, 7 * G, 1500, 0.5),
     liv((14 * G, 55.0), (28 * G, 25.0))),
    ('B07', 'Cojutti, 1 livello g7', CO, 70, regime(1, 0, 1500, 0.5), liv((7 * G, 40.0))),
    ('B08', 'Cojutti, regime irregolare 1500/1000/1500 a g0, g8, g22, 4 livelli', CO, 82,
     [{'ore': 0, 'dose': 1500, 'tinf': 0.5}, {'ore': 8 * G, 'dose': 1000, 'tinf': 0.5},
      {'ore': 22 * G + 5, 'dose': 1500, 'tinf': 1.0}], liv((7 * G, 42.0), (15 * G, 50.0), (22 * G, 22.0), (36 * G, 30.0))),
    ('B09', 'Baiardi, 70 kg, 1500 g1+g8 poi ogni 14 gg (6 dosi), 8 livelli', BA, 70,
     regime(2, 7 * G, 1500, 0.5) + regime(4, 14 * G, 1500, 0.5, t0=22 * G),
     liv((7 * G, 50.0), (14 * G, 65.0), (21 * G, 40.0), (35 * G, 45.0), (49 * G, 48.0), (63 * G, 50.0),
         (70 * G, 70.0), (84 * G, 30.0))),
    ('B10', 'Baiardi, 50 kg, 1 livello g14', BA, 50, regime(2, 7 * G, 1500, 0.5), liv((14 * G, 70.0))),
    ('B11', 'Mariani, livello molto alto (120 mg/L a g14)', MA, 35, regime(1, 0, 630, 0.5), liv((14 * G, 120.0))),
    ('B12', 'Baiardi, livello molto basso (5 mg/L a g7)', BA, 72, regime(1, 0, 1500, 0.5), liv((7 * G, 5.0))),
]

GE, AM, TO = 'Gentamicina', 'Amikacina', 'Tobramicina'
LP, AL, HE = models.AMINO_MODELLI[GE][0], models.AMINO_MODELLI[AM][0], models.AMINO_MODELLI[TO][1]
LPA = models.AMINO_MODELLI[AM][1]
AMINO = [  # id, descr, farmaco, modello, sesso, gg, peso, altezza, crea, eg_sett, crrt(Sc,Qef), storia, livelli
    ('A01', 'Llanos-Paez genta, 5 anni M 18 kg 108 cm, 1 livello a 10 h', GE, LP, 'Maschio', 1826, 18, 108, 0.40, None,
     None, regime(1, 0, 126, 0.5), liv((10.0, 0.9))),
    ('A02', 'Llanos-Paez genta, 1 anno F 9,5 kg 75 cm, picco + 10 h', GE, LP, 'Femmina', 365, 9.5, 75, 0.25, None,
     None, regime(1, 0, 67, 0.5), liv((1.0, 14.0), (10.0, 0.6))),
    ('A03', 'Llanos-Paez genta, 14 anni M 50 kg 165 cm, 1 livello a 10 h', GE, LP, 'Maschio', 5113, 50, 165, 0.70,
     None, None, regime(1, 0, 350, 0.5), liv((10.0, 1.4))),
    ('A04', 'Llanos-Paez genta, pretermine 20 gg EG 30 sett 1,3 kg 38 cm (fuori range)', GE, LP, 'Maschio', 20, 1.3,
     38, 0.80, 30, None, regime(2, 36, 6.5, 0.5), liv((46.0, 3.2))),
    ('A05', 'Llanos-Paez genta, creatinina 1,5 mg/dL, 8 anni 25 kg 128 cm, valle a 24 h', GE, LP, 'Femmina', 2922, 25,
     128, 1.5, None, None, regime(1, 0, 175, 0.5), liv((1.0, 16.0), (24.0, 2.5))),
    ('A06', 'Llanos-Paez per amikacina, 10 anni F 32 kg 140 cm, 10 h', AM, LPA, 'Femmina', 3652, 32, 140, 0.5, None,
     None, regime(1, 0, 480, 0.5), liv((10.0, 4.0))),
    ('A07', 'Llanos-Paez per tobramicina + CRRT (0,85 × 1,5), 6 anni 20 kg 115 cm', TO, models.AMINO_MODELLI[TO][0],
     'Maschio', 2191, 20, 115, 1.2, None, (0.85, 1.5), regime(1, 0, 140, 0.5), liv((1.0, 12.0), (10.0, 3.0))),
    ('A08', 'Llanos-Paez genta, 19 anni 70 kg 175 cm (Janmahasatian, fuori range)', GE, LP, 'Maschio', 6940, 70, 175,
     0.9, None, None, regime(1, 0, 490, 0.5), liv((10.0, 1.1))),
    ('A09', 'Alhadab amika, 8 anni 25 kg, picco + 10 h', AM, AL, 'Maschio', 2922, 25, None, 0.4, None, None,
     regime(1, 0, 500, 0.5), liv((1.0, 40.0), (10.0, 3.5))),
    ('A10', 'Alhadab amika, 16 anni 55 kg, 1 livello a 10 h', AM, AL, 'Femmina', 5844, 55, None, 0.7, None, None,
     regime(1, 0, 1100, 0.5), liv((10.0, 6.0))),
    ('A11', 'Alhadab amika, livello 0,3 mg/L (sotto l\'errore additivo 0,5)', AM, AL, 'Maschio', 4000, 35, None, 0.5,
     None, None, regime(1, 0, 700, 0.5), liv((18.0, 0.3))),
    ('A12', 'Hennig tobra FC, 12 anni 35 kg, picco + 10 h', TO, HE, 'Femmina', 4383, 35, None, 0.5, None, None,
     regime(1, 0, 350, 0.5), liv((1.0, 25.0), (10.0, 1.5))),
    ('A13', 'Hennig tobra FC, 4 anni 15 kg, 1 livello', TO, HE, 'Maschio', 1461, 15, None, 0.3, None, None,
     regime(1, 0, 150, 0.5), liv((10.0, 0.8))),
    ('A14', 'Hennig, valle molto alta (6 mg/L a 10 h)', TO, HE, 'Maschio', 3000, 26, None, 0.5, None, None,
     regime(1, 0, 260, 0.5), liv((10.0, 6.0))),
    ('A15', 'Llanos-Paez genta q8h, 3 dosi, 3 livelli', GE, LP, 'Femmina', 2557, 22, 120, 0.4, None, None,
     regime(3, 8, 55, 0.5), liv((1.0, 5.5), (16.5, 6.0), (23.9, 0.6))),
]


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
        auc, _, cmin = metriche_finestra(par, storia, t_last, tau_med, t_inf)
        return dict(auc24=auc * 24 / tau_med, cmin=cmin, w_start=t_last, w_len=tau_med, dosi=storia, fattore=24 / tau_med)
    auc, _, cmin = metriche_finestra(par, proiez, t_last, 24.0, t_inf)
    return dict(auc24=auc, cmin=cmin, w_start=t_last, w_len=24.0, dosi=proiez, fattore=1.0)


def esegui(pr, cov, id, descr, storia, livelli, dalba):
    t_inf = pr['t_inf']
    storia = [{**s, 'tinf': s['tinf'] if s['tinf'] > 0 else t_inf} for s in storia]
    om = omega_app(pr['omega'], len(livelli))
    par, _ = mappa_bayesiana(pr['par'], om, storia, livelli, pr['sigma_prop'], pr['sigma_add'], t_inf)
    par_lib, _ = mappa_bayesiana(pr['par'], om, storia, livelli, pr['sigma_prop'], pr['sigma_add'], t_inf,
                                 limiti=(1e-6, 1e6))
    stimati = [k for k in PAR if om.get(k, 0) > 0]
    lim = [k for k in stimati if not (0.0501 < par[k] / pr['par'][k] < 9.99)]
    usa = par_lib if lim else par
    m = metriche_app(usa, storia, t_inf)
    tgrid = np.round(np.linspace(0, storia[-1]['ore'] + (40 * G if dalba else 24), 97), 4).tolist()
    rec = dict(id=id, cov=cov, descr=descr, modello=pr['nome'], fuori_range=pr['fuori_range'], note=pr['note'],
               prior=pr['par'], omega_pubblicati=pr['omega'], omega=om, stimati=stimati,
               sigma_prop=pr['sigma_prop'], sigma_add=pr['sigma_add'], t_inf=t_inf, storia=storia, livelli=livelli,
               map=par, map_senza_limiti=par_lib, limiti_attivi=lim, auc24=m['auc24'], cmin=m['cmin'],
               finestra=dict(start=m['w_start'], len=m['w_len'], fattore=m['fattore'], dosi=m['dosi']),
               tgrid=tgrid, pop_curve=conc(pr['par'], storia, tgrid, t_inf).tolist())
    if dalba:   # come _pannello_dalba: da max(ultimo livello, ultima dose), passo 2 h, orizzonte 100 gg
        t_val = max(livelli[-1]['ore'], storia[-1]['ore'])
        rec['soglia'] = SOGLIA
        rec['t_val'] = t_val
        rec['t_sotto_soglia'] = tempo_sotto_soglia(usa, storia, SOGLIA, t_val, t_val + 2400, passo=2.0, t_inf=t_inf)
    return rec


out_d = []
for id, descr, mod, peso, storia, livelli in DALBA:
    pr = models.prior_dalba(peso, mod)
    out_d.append(esegui(pr, dict(peso=peso), id, descr, storia, livelli, True))
out_a = []
for id, descr, farm, mod, sesso, gg, peso, alt, crea, eg, crrt, storia, livelli in AMINO:
    Sc, Qef = crrt or (0.85, 2.0)
    pr = models.prior_amino(farm, mod, sesso, gg, peso, alt, crea, bool(crrt), Sc, Qef, eg_sett=eg)
    cov = dict(farmaco=farm, sesso=sesso, gg=gg, peso=peso, altezza=alt, crea=crea, eg_sett=eg,
               crrt=bool(crrt), Sc=Sc, Qef=Qef)
    out_a.append(esegui(pr, cov, id, descr, storia, livelli, False))

for nome, out in (('dalba', out_d), ('amino', out_a)):
    json.dump(out, open(f'{OUT}/casi_{nome}.json', 'w'), indent=1, ensure_ascii=False)
    with open(f'{OUT}/pigato_{nome}.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['id', 'stimati', 'CL', 'V1', 'Q', 'V2', 'AUC24', 'Cmin', 't_sotto_soglia_h', 'limiti_attivi'])
        for r in out:
            p = r['map_senza_limiti'] if r['limiti_attivi'] else r['map']
            w.writerow([r['id'], ';'.join(r['stimati'])] + [p[k] for k in PAR]
                       + [r['auc24'], r['cmin'], r.get('t_sotto_soglia', ''), ';'.join(r['limiti_attivi'])])
    print(nome, len(out), 'casi;', 'limiti attivi:', [r['id'] for r in out if r['limiti_attivi']])
