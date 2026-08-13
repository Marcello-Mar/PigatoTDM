"""
TDM Buddy v5.3
==============

FUNZIONALITÀ PRINCIPALI
  - MAP Bayesiano 1-compartimento (Daptomicina) e 2-compartimento (Ciprofloxacina/Dalbavancin)
  - Integrazione storia dosi: la TDM successiva recupera automaticamente tutte le dosi precedenti 
    dallo storico di sessione e le fonde nella timeline.
  - Archivio clinico SQLite con storico per paziente
  - Export PDF referto strutturato
  - Grafici PK Interattivi via Plotly (zoom, pan, toggle tracce)

MODELLI INTEGRATI — DAPTOMICINA
  [1] Soraluce A, et al. Int J Antimicrob Agents. 2018 Aug;52(2):158-165.
  [2] Wei XC, et al. J Clin Pharmacol. 2020 Jun;60(6):768-774.
  [3] Ye L, et al. Front Pharmacol. 2022 Aug 16;13:838599.
  [4] Ye L, et al. Eur J Pharm Sci. 2025 Mar;208:107058.
  [5] Lim SY, et al. J Pediatr Pharmacol Ther. 2021 Feb 15;26(2):123-132.
  [6] Bhavnani SM, et al. Clin Infect Dis. 2010 Jun 15;50(12):1568-1574.

MODELLI INTEGRATI — CIPROFLOXACINA
  [7] Schlender JF, et al. Clin Pharmacokinet. 2018 Dec;57(12):1613-1634.
  [8] Balbas-Martinez V, et al. Eur J Pharm Sci. 2019 Feb 1;128:171-179.
  [9] Forrest A, et al. Antimicrob Agents Chemother. 1993 May;37(5):1065-1072.
  [10] Ambrose PG, et al. Clin Infect Dis. 2007;44(Suppl 2):S79-S86.
"""

import streamlit as st
import numpy as np
from scipy.optimize import minimize
import plotly.graph_objects as go
import pandas as pd
from datetime import datetime, date, timedelta
import json
from fpdf import FPDF

# ============================================================
# STORAGE
# ============================================================

def storico_vuoto():
    return {
        'codice_paziente': '',
        'sesso': 'Maschio',
        'monitoraggi': [],
    }

if 'storico' not in st.session_state:
    st.session_state['storico'] = storico_vuoto()


# ============================================================
# PARAMETRI PK — DAPTOMICINA
# ============================================================

SORALUCE = {
    'CLNR': 0.160, 'CLR_theta': 0.366, 'Clcr_ref': 49.0,
    'V': 12.30, 'omega_CL': 0.367, 'omega_V': 0.278, 'sigma': 0.123,
    'Sc_default': 0.20,
    'nome': 'Soraluce et al., Int J Antimicrob Agents 2018;52(2):158-165 [Adulto ICU]',
}

WEI = {
    '3-6m'  : {'CL':19.7,'SD':5.5, 'V':0.17,'nome':'Wei et al., J Clin Pharmacol 2020;60(6):768-774 [3-6 mesi]'},
    '7-12m' : {'CL':19.6,'SD':5.8, 'V':0.17,'nome':'Wei et al., J Clin Pharmacol 2020;60(6):768-774 [7-12 mesi]'},
    '13-24m': {'CL':21.8,'SD':3.0, 'V':0.11,'nome':'Wei et al., J Clin Pharmacol 2020;60(6):768-774 [13-24 mesi]'},
    '2-6y'  : {'CL':20.0,'SD':6.1, 'V':0.11,'nome':'Wei et al., J Clin Pharmacol 2020;60(6):768-774 [2-6 anni]'},
    '7-11y' : {'CL':17.0,'SD':5.7, 'V':0.11,'nome':'Wei et al., J Clin Pharmacol 2020;60(6):768-774 [7-11 anni]'},
    '12-17y': {'CL':11.0,'SD':2.1, 'V':0.11,'nome':'Wei et al., J Clin Pharmacol 2020;60(6):768-774 [12-17 anni]'},
}

YE2025 = {
    'CL':19.6,'SD':5.5,'V':0.17,
    'nome': 'Ye L, Zhou H, et al., Eur J Pharm Sci 2025;208:107058 [Lattante <1 anno, PBPK]',
}

YE2022 = {
    '12-17y': {'mild':0.98,'moderate':1.08,'severe':1.24,'esrd':1.33},
    '7-11y' : {'mild':1.00,'moderate':1.19,'severe':1.49,'esrd':1.73},
    '2-6y'  : {'mild':1.02,'moderate':1.28,'severe':1.74,'esrd':2.20},
    '1-2y'  : {'mild':1.04,'moderate':1.29,'severe':1.72,'esrd':2.13},
}

LIM = {
    'Vss':5.90,'CLren':0.023,'CLnr':0.158,'CL_HD':0.619,
    'omega_V1':0.112,'omega_knr':0.088,'sigma':0.085,'peso_ref':40.0,
    'nome': 'Lim SY, et al., J Pediatr Pharmacol Ther 2021;26(2):123-132 [Adolescente HD/CCPD]',
    'nota': ('⚠️ Modello 2-cmt (Lim 2021) approssimato a 1-cmt per MAP. '
             'Q e V2 periferico non stimati. Per analisi dettagliata usare Phoenix NLME.'),
}

PD_DAPTO = {'AUC_MIC_staph':666,'AUC_MIC_entero':143,'Cmin_tox':24.3}

# ============================================================
# PARAMETRI PK — CIPROFLOXACINA
# ============================================================

CIPRO = {
    'CL_kg' :0.50,'V1_kg':0.55,'Q_kg':0.30,'V2_kg':1.90,
    'omega_CL':0.37,'omega_V1':0.32,'omega_Q':0.25,'omega_V2':0.30,
    'sigma':0.15,'fu':0.67,
    'nome':'Schlender JF, et al., Clin Pharmacokinet 2018;57(12):1613-1634 '
            '/ Balbas-Martinez V, et al., Eur J Pharm Sci 2019;128:171-179',
}

CIPRO_PED = {
    (91,  365): 1.40,
    (365, 730): 1.30,
    (730,2190): 1.20,
    (2190,4380):1.10,
    (4380,6570):1.00,
}
CIPRO_FRAZ_RENAL = 0.70

PD_CIPRO = {'AUC_MIC_target':125,'AUC_alta':400}

IRC_MAP = {
    "Insufficienza Renale Lieve (GFR 60-90)"   : 'mild',
    "Insufficienza Renale Moderata (GFR 30-60)" : 'moderate',
    "Insufficienza Renale Grave (GFR 15-30)"    : 'severe',
    "ESRD (GFR <15)"                            : 'esrd',
}

# ============================================================
# PARAMETRI PK — DALBAVANCINA (modelli Mariani e Cojutti 2021)
# ============================================================

DALBA_PARAMS = {
    'CL_typ70': 0.0305, 'V1_typ70': 4.81, 'Q_typ70': 0.0081, 'V2_typ70': 4.81,
    'omega_CL': 0.198, 'sigma': 0.375, 't_inf': 0.5, 'fu': 0.07,
    'nome': ("Modello Mariani (dati IRCCS Gaslini, non ancora pubblicato, "
             "modello interno fornito dall'utente)"),
}

COJUTTI_PARAMS = {
    'CL_typ70': 0.106, 'V1_typ70': 17.4, 'Q_typ70': 0.103, 'V2_typ70': 15.1,
    'omega_CL': 0.362, 'sigma': 0.18, 't_inf': 0.5, 'fu': 0.07,
    'nome': ("Cojutti PG, et al. Antimicrob Agents Chemother. 2021;65(9):"
             "e00879-21 (adulti, infezioni osteoarticolari, n=15; nessuno "
             "scaling per peso fornito — valori di popolazione fissi)"),
}

DALBA_MODELLI = {
    'Mariani (pediatrico, allometrico — Gaslini)': DALBA_PARAMS,
    'Cojutti 2021 (adulti, osteoarticolare — AAC)': COJUTTI_PARAMS,
}

DALBA_MIC_EUCAST = {
    'MIC90 osteoarticolare (0.0625 mg/L)': 0.0625,
    'Breakpoint EUCAST precedente (0.125 mg/L)': 0.125,
    'Breakpoint EUCAST attuale S.aureus (0.25 mg/L)': 0.25,
}

DALBA_TARGETS = {
    'ABSSSI semplice (fAUC24h/MIC >= 27.1)'      : 27.1,
    'Infezione complessa (fAUC24h/MIC >= 111.1)' : 111.1,
}

def dalba_soglia_concentrazione_totale(target_auc_mic, mic, f):
    return (target_auc_mic * mic) / (24.0 * f)

def dalba_soglia_trough(mic, f):
    return (4.0 * mic) / f

# ============================================================
# MOTORI PK
# ============================================================

def pk_1cmt(CL, V, storia, t_obs, t_inf=0.5):
    if CL <= 0 or V <= 0: return 0.0
    ke = CL / V
    C = 0.0
    for s in storia:
        dt = t_obs - s['ore']
        if dt <= 0: continue
        D = s['dose']
        if dt <= t_inf:
            C += (D/(V*ke*t_inf))*(1.0-np.exp(-ke*dt))
        else:
            Cf = (D/(V*ke*t_inf))*(1.0-np.exp(-ke*t_inf))
            C += Cf*np.exp(-ke*(dt-t_inf))
    return max(C, 0.0)


def pk_2cmt(CL, V1, Q, V2, storia, t_obs, t_inf=0.5):
    if CL<=0 or V1<=0 or Q<=0 or V2<=0: return 0.0
    k10=CL/V1; k12=Q/V1; k21=Q/V2
    s_=k10+k12+k21; d_=max(s_**2-4*k10*k21,0)
    alp=(s_+np.sqrt(d_))/2; bet=(s_-np.sqrt(d_))/2
    if abs(alp-bet)<1e-10: return pk_1cmt(CL,V1+V2,storia,t_obs,t_inf)
    A=(k21-alp)/((bet-alp)*V1); B=(k21-bet)/((alp-bet)*V1)
    C=0.0
    for s in storia:
        dt=t_obs-s['ore']
        if dt<=0: continue
        D=s['dose']; k0=D/t_inf
        if dt<=t_inf:
            ca=(A/alp)*(1-np.exp(-alp*dt))
            cb=(B/bet)*(1-np.exp(-bet*dt))
        else:
            ca=(A/alp)*(1-np.exp(-alp*t_inf))*np.exp(-alp*(dt-t_inf))
            cb=(B/bet)*(1-np.exp(-bet*t_inf))*np.exp(-bet*(dt-t_inf))
        C+=k0*(ca+cb)
    return max(C,0.0)


def metriche_esposizione_24h(pk_func, params, storia, t_inf, n_punti=300):
    t_ultima_dose = storia[-1]['ore']
    t_grid = np.linspace(t_ultima_dose, t_ultima_dose + 24, n_punti)
    c_grid = np.array([pk_func(*params, storia, t, t_inf) for t in t_grid])
    
    if hasattr(np, 'trapezoid'):
        auc = float(np.trapezoid(c_grid, x=t_grid))
    else:
        auc = float(np.trapz(c_grid, x=t_grid))
        
    return auc, float(c_grid.max()), float(c_grid.min())


def tau_mediano(storia):
    if len(storia) < 2:
        return None
    ore = sorted(s['ore'] for s in storia)
    diffs = [ore[i+1]-ore[i] for i in range(len(ore)-1)]
    diffs.sort()
    n = len(diffs)
    return diffs[n//2] if n % 2 == 1 else (diffs[n//2-1]+diffs[n//2])/2.0

# ============================================================
# PRIOR DI POPOLAZIONE
# ============================================================

def clcr_cg(eta_a, peso, crea, sesso):
    if crea<=0: return 60.0
    clcr=((140-eta_a)*peso)/(72*crea)
    if sesso=="Femmina": clcr*=0.85
    return max(clcr,1.0)


def prior_dapto(sesso, gg, peso, crea, crrt, Sc, Qef, dialisi):
    eta=gg/365.25; note=""
    if gg<365:
        p=YE2025; CL=(p['CL']*peso)/1000; V=p['V']*peso
        oCL=p['SD']/p['CL']; oV=0.30; sig=0.12; t_inf=1.0; nm=p['nome']
        if crrt: CLEC=Sc*Qef; CL+=CLEC; note+=f"CRRT: CLEC={CLEC:.3f} L/h (Sc={Sc}, Qef={Qef} L/h). "
        return CL,V,oCL,oV,sig,t_inf,nm,note
    if 4380<=gg<6570 and dialisi in ("Emodialisi (HD)","Dialisi Peritoneale (CCPD)"):
        sc=peso/LIM['peso_ref']
        CL=(LIM['CLren']+LIM['CLnr'])*sc; V=LIM['Vss']*sc
        nm=LIM['nome']; note=LIM['nota']
        if dialisi=="Emodialisi (HD)": note+=f" CLHD={LIM['CL_HD']*sc:.3f} L/h (somministrare DOPO HD)."
        return CL,V,LIM['omega_knr'],LIM['omega_V1'],LIM['sigma'],0.5,nm,note
    if 91<=gg<6570:
        if   gg<213:  ag='3-6m'
        elif gg<365:  ag='7-12m'
        elif gg<730:  ag='13-24m'
        elif gg<2190:  ag='2-6y'
        elif gg<4380: ag='7-11y'
        else:         ag='12-17y'
        p=WEI[ag]; CL=(p['CL']*peso)/1000; V=p['V']*peso
        oCL=p['SD']/p['CL']; oV=0.25; sig=0.12
        t_inf=1.0 if gg<730 else 0.5; nm=p['nome']
        age_ye = ag if ag in YE2022 else ('2-6y' if ag in ('13-24m','2-6y') else '1-2y')
        if dialisi in IRC_MAP and age_ye in YE2022:
            f=YE2022[age_ye][IRC_MAP[dialisi]]; CL/=f
            note+=f"Correzione IRC (Ye L, et al., Front Pharmacol 2022;13:838599): ×{f:.2f} AUC → CL={CL:.3f} L/h. "
            if IRC_MAP[dialisi] in ('severe','esrd'): note+="⚠️ Considerare intervallo q48h (Tabella 6, Ye 2022). "
        if crrt: CLEC=Sc*Qef; CL+=CLEC; note+=f"CRRT: CLEC={CLEC:.3f} L/h. "
        return CL,V,oCL,oV,sig,t_inf,nm,note
    s=SORALUCE; clcr=clcr_cg(eta,peso,crea,sesso)
    CL=s['CLNR']+s['CLR_theta']*(clcr/s['Clcr_ref']); V=s['V']
    nm=f"{s['nome']} (ClCr={clcr:.0f} mL/min)"
    if crrt: CLEC=Sc*Qef; CL+=CLEC; note+=f"CRRT: CLEC={CLEC:.3f} L/h (Sc={Sc}, Qef={Qef} L/h). CL totale={CL:.3f} L/h."
    return CL,V,s['omega_CL'],s['omega_V'],s['sigma'],0.5,nm,note


def prior_cipro(sesso, gg, peso, crea, crrt, dialisi):
    eta=gg/365.25; note=""; c=CIPRO
    fc=1.0
    for (lo,hi),f in CIPRO_PED.items():
        if lo<=gg<hi: fc=f; break
    if gg<6570 and fc!=1.0:
        note+=f"Scaling pediatrico (Schlender 2018, Fig.7): CL×{fc:.2f} per età {eta:.1f}a. "
    CL=c['CL_kg']*peso*fc; V1=c['V1_kg']*peso
    Q=c['Q_kg']*peso;       V2=c['V2_kg']*peso
    clcr=clcr_cg(eta,peso,crea,sesso)
    if clcr>=30:
        fr=min(clcr/120.0,1.0)
        CL=(c['CL_kg']*peso*fc*CIPRO_FRAZ_RENAL*fr)+(c['CL_kg']*peso*fc*(1-CIPRO_FRAZ_RENAL))
        if clcr<90: note+=f"IRC lieve-moderata: ClCr={clcr:.0f} mL/min → CL ridotta a {CL:.1f} L/h (Balbas-Martinez, Eur J Pharm Sci 2019;128:171-179). "
    else:
        CL_r=c['CL_kg']*peso*fc*CIPRO_FRAZ_RENAL*0.20
        CL_nr=c['CL_kg']*peso*fc*(1-CIPRO_FRAZ_RENAL)*0.63
        CL=CL_r+CL_nr
        note+=(f"⚠️ IRC grave (ClCr={clcr:.0f} mL/min): CL ridotta ({CL:.1f} L/h). "
               "CYP1A2 -37% (Rowland Yeo et al., Expert Rev Clin Pharmacol 2011;4:261-274). "
               "Considerare intervallo q48h e monitoraggio intensivo.")
    if crrt: note+="⚠️ CRRT: rimozione cipro variabile, non modellata esplicitamente. Monitorare livelli."
    nm=f"{c['nome']} (ClCr={clcr:.0f} mL/min)"
    return CL,V1,Q,V2,c['omega_CL'],c['omega_V1'],c['omega_Q'],c['omega_V2'],c['sigma'],1.0,nm,note


def prior_dalba(peso, modello_scelto):
    p = DALBA_MODELLI[modello_scelto]

    if modello_scelto.startswith('Mariani'):
        fw = peso / 70.0
        CL_pop = p['CL_typ70'] * (fw ** 0.75)
        V1_pop = p['V1_typ70'] * fw
        Q_pop  = p['Q_typ70']  * (fw ** 0.75)
        V2_pop = p['V2_typ70'] * fw
        note = ("Scaling allometrico applicato (CL,Q: peso/70^0.75; V1,V2: peso/70). "
                "Per insufficienza renale grave (ClCr<30 mL/min) fare riferimento "
                "alla scheda tecnica DALVANCE (FDA 2025) per l'aggiustamento posologico.")
    else:
        CL_pop, V1_pop, Q_pop, V2_pop = p['CL_typ70'], p['V1_typ70'], p['Q_typ70'], p['V2_typ70']
        note = ("⚠️ Modello Cojutti 2021: valori di popolazione applicati come fissi, "
                "indipendentemente dal peso inserito. Usare con cautela.")

    return CL_pop, V1_pop, Q_pop, V2_pop, p['omega_CL'], p['sigma'], p['t_inf'], p['nome'], note

# ============================================================
# MAP BAYESIANO
# ============================================================

def map_1cmt(CL_p, V_p, oCL, oV, sig, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0:
        return CL_p, V_p, "Prior di popolazione (nessun livello inserito)"

    def f(p):
        CL_, V_ = np.exp(p[0]), np.exp(p[1])
        pr = (((p[0]-np.log(CL_p))**2)/oCL**2 + ((p[1]-np.log(V_p))**2)/oV**2)
        lk = 0.0
        for liv in livelli:
            cp = pk_1cmt(CL_, V_, storia, liv['ore'], t_inf)
            lk += ((np.log(max(liv['conc'],1e-6)) - np.log(max(cp,1e-6)))**2) / sig**2
        return 0.5*(pr+lk)

    x0 = [np.log(CL_p), np.log(V_p)]
    bnd = [(np.log(CL_p*0.05), np.log(CL_p*8)), (np.log(V_p*0.15), np.log(V_p*5))]
    r = minimize(f, x0, method='L-BFGS-B', bounds=bnd, options={'maxiter':2000,'ftol':1e-12})
    lbl = f"MAP Bayesiano 1-cmt su {n} livell{'o' if n==1 else 'i'} plasmatic{'o' if n==1 else 'i'}"
    return np.exp(r.x[0]), np.exp(r.x[1]), lbl


def map_2cmt(CL_p, V1_p, Q_p, V2_p, oCL, oV1, oQ, oV2, sig, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0:
        return CL_p, V1_p, Q_p, V2_p, "Prior di popolazione (nessun livello inserito)"

    fattore_allentamento = min(1.0, n / 6.0)
    OQ_eff  = 0.10 + (oQ  - 0.10) * fattore_allentamento
    OV2_eff = 0.10 + (oV2 - 0.10) * fattore_allentamento

    def f(p):
        CL_, V1_, Q_, V2_ = np.exp(p[0]), np.exp(p[1]), np.exp(p[2]), np.exp(p[3])
        pr = (((p[0]-np.log(CL_p))**2)/oCL**2 + ((p[1]-np.log(V1_p))**2)/oV1**2
             +((p[2]-np.log(Q_p)) **2)/OQ_eff**2 + ((p[3]-np.log(V2_p))**2)/OV2_eff**2)
        lk = 0.0
        for liv in livelli:
            cp = pk_2cmt(CL_, V1_, Q_, V2_, storia, liv['ore'], t_inf)
            lk += ((np.log(max(liv['conc'],1e-6)) - np.log(max(cp,1e-6)))**2) / sig**2
        return 0.5*(pr+lk)

    x0 = [np.log(CL_p), np.log(V1_p), np.log(Q_p), np.log(V2_p)]
    bnd = [(np.log(CL_p*0.05), np.log(CL_p*10)), (np.log(V1_p*0.10), np.log(V1_p*6)),
           (np.log(Q_p*0.30),  np.log(Q_p*3.0)), (np.log(V2_p*0.30), np.log(V2_p*3.0))]
    r = minimize(f, x0, method='L-BFGS-B', bounds=bnd, options={'maxiter':2000,'ftol':1e-12})
    nota_id = ("Q e V2 ancorati stretti al prior (pochi campioni)" if n <= 2
              else f"Q e V2 stimati con vincolo allentato ({n} campioni)")
    lbl = f"MAP Bayesiano 2-cmt su {n} livell{'o' if n==1 else 'i'} — {nota_id}"
    return np.exp(r.x[0]), np.exp(r.x[1]), np.exp(r.x[2]), np.exp(r.x[3]), lbl


def map_dalba(CL_pop, V1_pop, Q_pop, V2_pop, omega_CL, sigma, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0:
        return CL_pop, "Prior di popolazione (nessun livello inserito)"

    def f(p):
        CL_ = np.exp(p[0])
        prior = ((p[0] - np.log(CL_pop)) ** 2) / omega_CL ** 2
        lik = 0.0
        for liv in livelli:
            cp = pk_2cmt(CL_, V1_pop, Q_pop, V2_pop, storia, liv['ore'], t_inf)
            lik += ((np.log(max(liv['conc'],1e-6)) - np.log(max(cp,1e-6))) ** 2) / sigma ** 2
        return 0.5 * (prior + lik)

    x0 = [np.log(CL_pop)]
    bnd = [(np.log(CL_pop*0.05), np.log(CL_pop*10.0))]
    r = minimize(f, x0, method='L-BFGS-B', bounds=bnd, options={'maxiter': 2000, 'ftol': 1e-12})
    lbl = f"MAP Bayesiano su CL, {n} livell{'o' if n==1 else 'i'} plasmatic{'o' if n==1 else 'i'}; V1/Q/V2 fissi"
    return np.exp(r.x[0]), lbl


def trova_tempo_soglia(CL, V1, Q, V2, storia, soglia, t_inf, t_partenza=0.0, t_max=2000.0, passo=1.0):
    def C(t):
        return pk_2cmt(CL, V1, Q, V2, storia, t, t_inf)

    c0 = C(t_partenza)
    if c0 < soglia:
        return t_partenza

    t_prev = t_partenza
    t = t_partenza + passo
    while t <= t_max:
        c = C(t)
        if c < soglia:
            lo, hi = t_prev, t
            for _ in range(60):
                mid = (lo + hi) / 2.0
                if C(mid) >= soglia:
                    lo = mid
                else:
                    hi = mid
            return hi
        t_prev = t
        t += passo
    return None

# ============================================================
# STORIA DOSI: COSTRUZIONE E FUSIONE
# ============================================================

def carica_storia_precedente(farmaco):
    monitoraggi = st.session_state['storico'].get('monitoraggi', [])
    candidati = [m for m in monitoraggi if m.get('farmaco') == farmaco]
    if not candidati:
        return [], [], None
    ultimo = candidati[-1]
    try:
        storia = ultimo['storia_dosi']
        livelli = ultimo.get('livelli', [])
        t0 = datetime.strptime(ultimo['data_prima_dose'], '%Y-%m-%d %H:%M')
        return storia, livelli, t0
    except Exception:
        return [], [], None

def costruisci_storia(nuove_somm, t0_assoluto, storia_precedente):
    storia = list(storia_precedente)
    for s in nuove_somm:
        ore = (s['dt'] - t0_assoluto).total_seconds() / 3600.0
        storia.append({'ore': ore, 'dose': s['dose']})
    storia.sort(key=lambda x: x['ore'])
    return storia

def costruisci_livelli(nuovi_liv, t0_assoluto, livelli_precedenti):
    livelli = list(livelli_precedenti)
    for l in nuovi_liv:
        ore = (l['dt'] - t0_assoluto).total_seconds() / 3600.0
        livelli.append({'ore': ore, 'conc': l['conc'], 'nota': l.get('nota', '')})
    livelli.sort(key=lambda x: x['ore'])
    return livelli

# ============================================================
# PDF
# ============================================================

def genera_pdf(testo):
    safe = (testo.replace('—','-').replace('–','-')
                 .replace('\u2019',"'").replace('\u201c','"').replace('\u201d','"')
                 .replace('\u03bc','u').replace('\u2265','>='))
    safe = safe.encode('latin-1','replace').decode('latin-1')
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial",'B',13)
    pdf.cell(0, 9, "TDM Buddy v5.3 - Referto di Monitoraggio Terapeutico Farmacologico", ln=1, align='C')
    pdf.set_font("Arial", size=9)
    pdf.ln(2)
    for riga in safe.split('\n'):
        pdf.multi_cell(0, 5, txt=riga)
    return pdf.output(dest='S').encode('latin-1')

# ============================================================
# GRAFICI PLOTLY INTERATTIVI
# ============================================================

def disegna_grafico_plotly(farmaco, storia, livelli,
                           CL_ind, V_ind, Q_ind, V2_ind, CL_pop, V_pop, Q_pop, V2_pop,
                           t_inf, storia_precedente, t_fine_extra=12):
    t_ultimo_livello = max((l['ore'] for l in livelli), default=0)
    t_fine = max(t_ultimo_livello, max(s['ore'] for s in storia) + 12) + t_fine_extra
    t_arr = np.linspace(0, t_fine, 800)

    if farmaco == "Daptomicina":
        C_ind = [pk_1cmt(CL_ind, V_ind, storia, t, t_inf) for t in t_arr]
        C_pop = [pk_1cmt(CL_pop, V_pop, storia, t, t_inf) for t in t_arr]
    else:
        C_ind = [pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia, t, t_inf) for t in t_arr]
        C_pop = [pk_2cmt(CL_pop, V_pop, Q_pop, V2_pop, storia, t, t_inf) for t in t_arr]

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=t_arr, y=C_pop, mode='lines',
        line=dict(dash='dash', color='#94A3B8', width=2),
        name='Profilo Popolazione (Prior)'
    ))
    
    fig.add_trace(go.Scatter(
        x=t_arr, y=C_ind, mode='lines',
        line=dict(color='#0F172A', width=3),
        name='Profilo Individuale (MAP)'
    ))

    max_c = max(max(C_ind), max(C_pop)) if len(C_ind) > 0 else 100

    t_cambio_val = max((s['ore'] for s in storia_precedente), default=-1) if storia_precedente else -1

    if storia_precedente:
        fig.add_vrect(x0=0, x1=t_cambio_val, fillcolor="#6EA8D8", opacity=0.1, 
                      layer="below", line_width=0, name="Terapia precedente")

    x_dosi_vecchie, y_dosi_vecchie = [], []
    x_dosi_nuove, y_dosi_nuove = [], []

    for s in storia:
        if s['ore'] <= t_cambio_val:
            x_dosi_vecchie.extend([s['ore'], s['ore'], None])
            y_dosi_vecchie.extend([0, max_c, None])
        else:
            x_dosi_nuove.extend([s['ore'], s['ore'], None])
            y_dosi_nuove.extend([0, max_c, None])

    if x_dosi_vecchie:
        fig.add_trace(go.Scatter(x=x_dosi_vecchie, y=y_dosi_vecchie, mode='lines',
                                 line=dict(color='#6EA8D8', dash='dot', width=1),
                                 name='Dosi (Storico)', hoverinfo='none'))
    if x_dosi_nuove:
        fig.add_trace(go.Scatter(x=x_dosi_nuove, y=y_dosi_nuove, mode='lines',
                                 line=dict(color='#475569', dash='dot', width=1),
                                 name='Dosi (Nuove)', hoverinfo='none'))

    if livelli:
        x_liv = [l['ore'] for l in livelli]
        y_liv = [l['conc'] for l in livelli]
        testi = [f"Misurazione<br>t={l['ore']:.1f}h<br>C={l['conc']:.2f} mg/L<br>{l.get('nota','')}" for l in livelli]
        fig.add_trace(go.Scatter(x=x_liv, y=y_liv, mode='markers',
                                 marker=dict(symbol='star', size=14, color='#16A34A', line=dict(width=1, color='white')),
                                 name='Livelli misurati', text=testi, hoverinfo='text'))

    if farmaco == "Daptomicina":
        fig.add_hline(y=24.3, line_dash="dashdot", line_color="#DC2626", 
                      annotation_text="Soglia tossicità CPK (24.3 mg/L)", annotation_position="top right")

    fig.update_layout(
        xaxis_title="Ore dall'inizio della terapia (t=0 prima dose)",
        yaxis_title=f"Concentrazione plasmatica {farmaco} (mg/L)",
        plot_bgcolor='#FAFAFA',
        paper_bgcolor='#FAFAFA',
        hovermode="x unified",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=40, r=40, t=40, b=40)
    )
    return fig


def disegna_forecast_plotly(t_arr_sim, c_sim, soglia_target, t_new_ore, dose_new, storia_completa):
    fig2 = go.Figure()
    
    fig2.add_trace(go.Scatter(
        x=t_arr_sim, y=c_sim, mode='lines', 
        line=dict(color='#0F172A', width=2.5), 
        name='Profilo proiettato'
    ))
    
    fig2.add_hline(y=soglia_target, line_dash="dashdot", line_color="#DC2626", 
                   annotation_text=f"Soglia target ({soglia_target:.2f} mg/L)", annotation_position="top right")
                   
    fig2.add_vline(x=t_new_ore, line_dash="dash", line_color="#16A34A", 
                   annotation_text=f"Nuova dose ({dose_new:.0f} mg)")

    x_dosi, y_dosi = [], []
    max_c_sim = max(c_sim) if len(c_sim) > 0 else 100
    for s in storia_completa:
        x_dosi.extend([s['ore'], s['ore'], None])
        y_dosi.extend([0, max_c_sim, None])
    
    if x_dosi:
        fig2.add_trace(go.Scatter(x=x_dosi, y=y_dosi, mode='lines',
                                  line=dict(color='#94A3B8', dash='dot', width=1),
                                  name='Dosi passate', hoverinfo='none'))

    fig2.update_layout(
        xaxis_title="Ore dall'inizio terapia",
        yaxis_title="Dalbavancina (mg/L)",
        plot_bgcolor='#FAFAFA', paper_bgcolor='#FAFAFA',
        hovermode="x unified",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=40, r=40, t=40, b=40)
    )
    return fig2

# ============================================================
# ANALISI PRINCIPALE
# ============================================================

def esegui_analisi(farmaco, id_paz, sesso, gg, peso, dt_val,
                   crea, crrt, Sc, Qef, dialisi,
                   mic, livelli_completi,
                   storia_completa, storia_precedente, livelli_precedenti, t0_assoluto,
                   salva=False, modello_dalba=None):

    if not storia_completa:
        st.error("Nessuna dose inserita."); return

    dose_mg       = storia_completa[-1]['dose']
    dose_mg_kg    = dose_mg / peso
    eta           = gg / 365.25
    tau_h_display = tau_mediano(storia_completa)  

    if farmaco == "Daptomicina":
        CL_p,V_p,oCL,oV,sig,t_inf,nm,note = prior_dapto(sesso,gg,peso,crea,crrt,Sc,Qef,dialisi)
        CL_ind,V_ind,tipo_fit = map_1cmt(CL_p,V_p,oCL,oV,sig,storia_completa,livelli_completi,t_inf)
        Q_ind=V2_ind=Q_p=V2_p=0.0
        AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_1cmt, (CL_ind,V_ind), storia_completa, t_inf)
        auc_mic  = AUC24_v / mic
        emivita  = np.log(2)*V_ind/CL_ind
        Vss_d    = V_ind
        ok_staph = auc_mic>=PD_DAPTO['AUC_MIC_staph']
        ok_ent   = auc_mic>=PD_DAPTO['AUC_MIC_entero']
        tox      = Cmin_ss>=PD_DAPTO['Cmin_tox']
    elif farmaco == "Dalbavancin":
        if modello_dalba is None:
            modello_dalba = list(DALBA_MODELLI.keys())[0]
        CL_p,V1_p,Qp,V2_p,oCL,sig,t_inf,nm,note = prior_dalba(peso, modello_dalba)
        CL_ind, tipo_fit = map_dalba(CL_p,V1_p,Qp,V2_p,oCL,sig,storia_completa,livelli_completi,t_inf)
        V_ind, Q_ind, V2_ind = V1_p, Qp, V2_p
        Q_p = Qp; Vss_d = V_ind + V2_ind
        AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_2cmt, (CL_ind,V_ind,Q_ind,V2_ind), storia_completa, t_inf)
        auc_mic  = AUC24_v / mic
        emivita  = np.log(2)*Vss_d/CL_ind
        ok_target= auc_mic >= DALBA_TARGETS['ABSSSI semplice (fAUC24h/MIC >= 27.1)']
        tox      = False  
    else:
        CL_p,V1_p,Qp,V2_p,oCL,oV1,oQ,oV2,sig,t_inf,nm,note = prior_cipro(sesso,gg,peso,crea,crrt,dialisi)
        CL_ind,V_ind,Q_ind,V2_ind,tipo_fit = map_2cmt(CL_p,V1_p,Qp,V2_p,oCL,oV1,oQ,oV2,sig,storia_completa,livelli_completi,t_inf)
        Q_p=Qp; Vss_d=V_ind+V2_ind
        AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_2cmt, (CL_ind,V_ind,Q_ind,V2_ind), storia_completa, t_inf)
        auc_mic  = AUC24_v/mic
        emivita  = np.log(2)*Vss_d/CL_ind
        ok_target= auc_mic>=PD_CIPRO['AUC_MIC_target']
        tox      = AUC24_v>=PD_CIPRO['AUC_alta']

    st.markdown("---")
    if storia_precedente or livelli_precedenti:
        n_nuove_dosi = len(storia_completa) - len(storia_precedente)
        n_nuovi_liv  = len(livelli_completi) - len(livelli_precedenti)
        st.info(f"🔗 **Analisi integrata con storico:** {len(storia_precedente)} dosi + "
                f"{len(livelli_precedenti)} livelli precedenti, più {n_nuove_dosi} nuove dosi "
                f"e {n_nuovi_liv} nuovi livelli in questa valutazione. "
                f"Timeline totale dal t=0 originale ({t0_assoluto.strftime('%d/%m/%Y %H:%M')}).")

    st.subheader(f"🔬 Risultati PK/PD — {farmaco}")
    st.caption(f"**Modello:** {nm}")
    if note: st.info(note)

    c1,c2,c3 = st.columns(3)
    c1.metric("CL Individuale",f"{CL_ind:.3f} L/h",
              delta=f"{(CL_ind/peso)*1000:.1f} mL/h/kg",
              help=f"Prior di popolazione: {CL_p:.3f} L/h")
    c2.metric("Vss (distribuzione)",f"{Vss_d:.1f} L",
              delta=f"{Vss_d/peso:.2f} L/kg",
              help="1-cmt (dapto): V; 2-cmt (cipro/dalba): V1+V2")
    c3.metric("Emivita terminale (T½)",f"{emivita:.1f} h",
              help="Approssimazione: T½ = ln2 × Vss / CL")
    st.markdown("")
    c4,c5,c6 = st.columns(3)
    c4.metric("AUC₂₄ stimata",f"{AUC24_v:.0f} mg·h/L")

    if farmaco=="Daptomicina":
        c5.metric("AUC₂₄/MIC  (target S.aureus ≥ 666)",f"{auc_mic:.0f}",
                  delta="✅ TARGET" if ok_staph else "⚠️ SOTTO TARGET",
                  delta_color="normal" if ok_staph else "inverse")
        c6.metric(f"Cmin SS prevista  (soglia CPK ≥ 24.3 mg/L)",f"{Cmin_ss:.1f} mg/L",
                  delta="⚠️ RISCHIO CPK" if tox else "✅ Sicuro",
                  delta_color="inverse" if tox else "normal")
    elif farmaco=="Dalbavancin":
        c5.metric("fAUC₂₄h/MIC  (target ABSSSI ≥ 27.1)",f"{auc_mic:.0f}",
                  delta="✅ TARGET" if ok_target else "⚠️ SOTTO TARGET",
                  delta_color="normal" if ok_target else "inverse",
                  help="fAUC calcolata numericamente sulle 24h successive all'ultima dose")
        c6.metric("Cmax / Cmin previste",f"{Cmax_ss:.1f} / {Cmin_ss:.1f} mg/L",
                  help="Vedi sezione Forecast per pianificazione redose")
    else:
        c5.metric("AUC₂₄/MIC  (target gram-neg ≥ 125)",f"{auc_mic:.0f}",
                  delta="✅ TARGET" if ok_target else "⚠️ SOTTO TARGET",
                  delta_color="normal" if ok_target else "inverse")
        c6.metric(f"AUC₂₄  (soglia neurotox ≥ 400 mg·h/L)",f"{AUC24_v:.0f} mg·h/L",
                  delta="⚠️ AUC ELEVATA" if tox else "✅ Esposizione sicura",
                  delta_color="inverse" if tox else "normal")

    st.subheader("📈 Profilo Farmacocinetico Interattivo")
    fig = disegna_grafico_plotly(farmaco, storia_completa, livelli_completi,
                                 CL_ind, V_ind, Q_ind, V2_ind, CL_p, V_p if farmaco=="Daptomicina" else V1_p,
                                 Q_p, V2_p, t_inf, storia_precedente)
    st.plotly_chart(fig, use_container_width=True)

    def E(ok): return "TARGET RAGGIUNTO" if ok else "INSUFFICIENTE — RIVALUTARE POSOLOGIA"
    storico_note = ""
    if storia_precedente:
        storico_note = (f"\n[INTEGRAZIONE STORICO DOSI]\n"
                        f"  Dosi precedenti caricate da storico sessione: {len(storia_precedente)}\n"
                        f"  Nuove dosi correnti:            {len(storia_completa)-len(storia_precedente)}\n"
                        f"  t=0 assoluto terapia:           {t0_assoluto.strftime('%d/%m/%Y %H:%M')}\n"
                        f"  Timeline totale:                {storia_completa[-1]['ore']:.1f} h dal t=0\n")

    if farmaco=="Dalbavancin":
        pd_block=(
            f"  Efficacia (fAUC24h/MIC >= 27.1, target ABSSSI semplice):\n"
            f"    fAUC24h/MIC = {auc_mic:.0f}  ->  {E(ok_target)}\n"
            f"  Cmax prevista (24h dopo l'ultima dose): {Cmax_ss:.1f} mg/L\n"
            f"  Cmin prevista (24h dopo l'ultima dose): {Cmin_ss:.1f} mg/L\n"
        )
    elif farmaco=="Daptomicina":
        pd_block=(
            f"  Efficacia vs S. aureus / MRSA     (AUC24/MIC >= 666):\n"
            f"    AUC24/MIC = {auc_mic:.0f}  ->  {E(ok_staph)}\n"
            f"  Efficacia vs E. faecium           (AUC24/MIC >= 143):\n"
            f"    AUC24/MIC = {auc_mic:.0f}  ->  {E(ok_ent)}\n"
            f"  Sicurezza muscolare  (Cmin_ss >= 24.3 mg/L):\n"
            f"    Cmin_ss = {Cmin_ss:.1f} mg/L  ->  {'ATTENZIONE: MONITORARE CPK' if tox else 'SICURO'}\n"
        )
    else:
        pd_block=(
            f"  Efficacia gram-negativi           (AUC24/MIC >= 125):\n"
            f"    AUC24/MIC = {auc_mic:.0f}  ->  {E(ok_target)}\n"
            f"  Esposizione sistemica: AUC24 = {AUC24_v:.0f} mg*h/L\n"
            f"  Rischio neurotossicita' (AUC24 >= 400 mg*h/L):\n"
            f"    {'MONITORARE: AUC ELEVATA' if tox else 'SICURO'}\n"
        )

    referto = f"""\
=== TDM BUDDY v5.3 — REFERTO {farmaco.upper()} ===
Data valutazione: {dt_val}   ID paziente: {id_paz}
Sesso: {sesso}   Eta': {eta:.1f} anni   Peso: {peso} kg

[MODELLO FARMACOGENETICO APPLICATO]
  {nm}
  Algoritmo: {tipo_fit}
{storico_note}
[DATI CLINICI]
  Creatinina sierica:      {crea} mg/dL
  Dialisi / CRRT:          {dialisi}
  Dose più recente:        {dose_mg:.0f} mg  ({dose_mg_kg:.1f} mg/kg)
  Numero dosi totali:      {len(storia_completa)}
  Intervallo mediano dosi: {(f'{tau_h_display:.1f} h' if tau_h_display else "N/D")}
  MIC patogeno:            {mic} mg/L
  Livelli plasmatici usati per il fit: {len(livelli_completi)}

[PARAMETRI PK INDIVIDUALI STIMATI]
  Clearance (CL):               {CL_ind:.3f} L/h  ({(CL_ind/peso)*1000:.1f} mL/h/kg)
  Volume di distribuzione (Vss):{Vss_d:.2f} L     ({Vss_d/peso:.2f} L/kg)
  Emivita terminale (T1/2):     {emivita:.1f} h

[VALUTAZIONE PK/PD]
  AUC24 stimata (finestra 24h dopo l'ultima dose): {AUC24_v:.0f} mg*h/L

{pd_block}
"""
    st.subheader("📋 Referto per Cartella Clinica Elettronica")
    col_r, col_b = st.columns([4,1])
    with col_r:
        st.text_area("Testo strutturato:", value=referto, height=300)
    with col_b:
        st.write(""); st.write("")
        pdf_bytes = genera_pdf(referto)
        st.download_button("📄 Scarica PDF", data=pdf_bytes,
                           file_name=f"TDMBuddy_{farmaco}_{id_paz}_{dt_val[:10]}.pdf",
                           mime="application/pdf", use_container_width=True)

    if farmaco == "Dalbavancin":
        st.markdown("---")
        st.subheader("🔮 Forecast — Pianificazione Prossima Somministrazione")

        colF1, colF2, colF3 = st.columns(3)
        with colF1:
            target_scelto = st.selectbox(
                "Target PK/PD di riferimento:",
                list(DALBA_TARGETS.keys()) + ["Trough (100% fT > 4×MIC)"]
            )
        with colF2:
            mic_scelta_label = st.selectbox("MIC / breakpoint:", list(DALBA_MIC_EUCAST.keys()) + ["Personalizzato"])
            if mic_scelta_label == "Personalizzato":
                mic_forecast = st.number_input("MIC (mg/L)", 0.001, 4.0, 0.125, step=0.001, format="%.4f")
            else:
                mic_forecast = DALBA_MIC_EUCAST[mic_scelta_label]
        with colF3:
            f_unbound = st.number_input("Frazione libera (f)", 0.01, 0.10, DALBA_PARAMS['fu'], 0.01)

        if target_scelto == "Trough (100% fT > 4×MIC)":
            soglia_target = dalba_soglia_trough(mic_forecast, f_unbound)
            desc_target = f"trough 4×MIC libera → soglia totale = {soglia_target:.2f} mg/L"
        else:
            target_val = DALBA_TARGETS[target_scelto]
            soglia_target = dalba_soglia_concentrazione_totale(target_val, mic_forecast, f_unbound)
            desc_target = f"{target_scelto} → soglia totale equivalente = {soglia_target:.2f} mg/L"

        st.info(f"📐 Soglia calcolata: **{desc_target}**")

        t_ultimo_liv = max((l['ore'] for l in livelli_completi), default=0)
        t_valutazione = max(t_ultimo_liv, storia_completa[-1]['ore'])
        c_attuale = pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia_completa, t_valutazione, t_inf)

        t_crossing = trova_tempo_soglia(CL_ind, V_ind, Q_ind, V2_ind, storia_completa,
                                        soglia_target, t_inf, t_partenza=t_valutazione,
                                        t_max=t_valutazione + 2000, passo=2.0)

        if t_crossing is None:
            st.warning("⚠️ Soglia non raggiunta entro 83 giorni.")
        elif t_crossing <= t_valutazione + 0.5:
            st.error(f"🔴 Concentrazione già sotto soglia ({c_attuale:.2f} mg/L).")
        else:
            data_crossing = t0_assoluto + timedelta(hours=t_crossing)
            giorni_rimanenti = (t_crossing - t_valutazione) / 24.0
            st.success(f"🗓️ La concentrazione scenderà sotto soglia il **{data_crossing.strftime('%d/%m/%Y %H:%M')}** ({giorni_rimanenti:.1f} giorni).")

        st.markdown("#### 💉 Simula nuova somministrazione ipotetica")
        col_sim1, col_sim2, col_sim3 = st.columns(3)
        with col_sim1: simula = st.checkbox("Aggiungi dose al forecast")
        with col_sim2:
            if simula:
                data_new = st.date_input("Data nuova dose", value=(t0_assoluto + timedelta(hours=t_valutazione+24)).date(), format="DD/MM/YYYY")
                ora_new = st.time_input("Ora", value=datetime.strptime("10:00","%H:%M").time())
        with col_sim3:
            if simula: dose_new = st.number_input("Dose ipotetica (mg)", 100.0, 1500.0, 1500.0, 50.0)

        if simula:
            dt_new = datetime.combine(data_new, ora_new)
            t_new_ore = (dt_new - t0_assoluto).total_seconds() / 3600.0
            storia_sim = storia_completa + [{'ore': t_new_ore, 'dose': dose_new}]
            storia_sim.sort(key=lambda x: x['ore'])

            t_crossing_sim = trova_tempo_soglia(CL_ind, V_ind, Q_ind, V2_ind, storia_sim,
                                                soglia_target, t_inf, t_partenza=t_new_ore, t_max=t_new_ore+2000, passo=2.0)
            if t_crossing_sim:
                data_crossing_sim = t0_assoluto + timedelta(hours=t_crossing_sim)
                st.success(f"✅ Con dose di {dose_new:.0f} mg, concentrazione sopra soglia fino al **{data_crossing_sim.strftime('%d/%m/%Y')}**.")

            t_fine_sim = t_crossing_sim + 48 if t_crossing_sim else t_new_ore + 500
            t_arr_sim = np.linspace(0, t_fine_sim, 1000)
            c_sim = [pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia_sim, t, t_inf) for t in t_arr_sim]
            
            fig2 = disegna_forecast_plotly(t_arr_sim, c_sim, soglia_target, t_new_ore, dose_new, storia_completa)
            st.plotly_chart(fig2, use_container_width=True)

    if salva and id_paz and id_paz != "":
        try:
            record = {
                'data_val': dt_val, 'farmaco': farmaco,
                'eta_giorni': gg, 'peso': peso,
                'creatinina': crea, 'crrt': bool(crrt), 'Sc': Sc, 'Qef': Qef,
                'dialisi_tipo': dialisi,
                'dose_mg': dose_mg, 'dose_mg_kg': dose_mg_kg,
                'tau_h': tau_h_display,
                'n_dosi': len(storia_completa), 'mic': mic,
                'n_livelli': len(livelli_completi),
                'CL_ind': CL_ind, 'V_ind': Vss_d, 'AUC24': AUC24_v,
                'AUC_MIC': auc_mic, 'Cmin_ss': Cmin_ss, 'modello': nm,
                'modello_dalba': modello_dalba,
                'storia_dosi': storia_completa,
                'livelli': livelli_completi,
                'data_prima_dose': t0_assoluto.strftime('%Y-%m-%d %H:%M'),
            }
            st.session_state['storico']['monitoraggi'].append(record)
            st.session_state['storico']['codice_paziente'] = id_paz
            st.session_state['storico']['sesso'] = sesso

            st.success("💾 Analisi aggiunta allo storico di sessione.")

            json_bytes = json.dumps(st.session_state['storico'], indent=2,
                                    ensure_ascii=False).encode('utf-8')
            st.download_button(
                "📥 Scarica storico paziente aggiornato (JSON)",
                data=json_bytes,
                file_name=f"TDMBuddy_storico_{id_paz}.json",
                mime="application/json",
                use_container_width=True,
                type="primary",
            )
        except Exception as e:
            st.warning(f"Errore nel salvataggio in sessione: {e}")

# ============================================================
# INTERFACCIA STREAMLIT
# ============================================================

st.set_page_config(page_title="TDM Buddy v5.3", layout="wide", page_icon="💊")

st.title("💊 TDM Buddy v5.3")
st.caption(
    "Monitoraggio Terapeutico Bayesiano — Daptomicina, Ciprofloxacina, Dalbavancin  |  "
    "MAP cumulativo 1-cmt / 2-cmt / 2-cmt allometrico"
)

st.sidebar.header("🗂️ Storico Paziente (locale)")

file_up = st.sidebar.file_uploader(
    "📤 Carica storico paziente (JSON)", type=["json"],
    help="Carica il file scaricato in una sessione precedente per continuare la terapia."
)
if file_up is not None:
    try:
        dati_caricati = json.load(file_up)
        assert 'monitoraggi' in dati_caricati
        st.session_state['storico'] = dati_caricati
        st.sidebar.success(f"✅ Storico caricato: '{dati_caricati.get('codice_paziente','?')}'")
    except Exception as e:
        st.sidebar.error(f"File non valido: {e}")

id_paz = st.sidebar.text_input(
    "Codice paziente", value=st.session_state['storico'].get('codice_paziente', '')
)
sesso = st.sidebar.selectbox(
    "Sesso biologico", ["Maschio", "Femmina"],
    index=0 if st.session_state['storico'].get('sesso', 'Maschio') == 'Maschio' else 1
)
st.session_state['storico']['codice_paziente'] = id_paz
st.session_state['storico']['sesso'] = sesso

n_rec = len(st.session_state['storico'].get('monitoraggi', []))
if n_rec > 0:
    st.sidebar.info(f"📋 {n_rec} analisi in questa sessione.")
    json_bytes_sb = json.dumps(st.session_state['storico'], indent=2, ensure_ascii=False).encode('utf-8')
    st.sidebar.download_button("📥 Scarica storico corrente", data=json_bytes_sb, file_name=f"TDMBuddy_{id_paz}.json", mime="application/json")
    if st.sidebar.button("🗑️ Cancella storico"):
        st.session_state['storico'] = storico_vuoto()
        st.rerun()

st.sidebar.markdown("---")
st.sidebar.markdown("### 📚 Modelli Integrati")
st.sidebar.markdown("""
**Daptomicina**
- Soraluce 2018
- Wei 2020 
- Ye 2022
- Ye 2025
- Lim 2021
- Bhavnani 2010

**Ciprofloxacina**
- Schlender 2018
- Balbas-Martinez 2019
- Forrest 1993
- Ambrose 2007

**Dalbavancin (Mariani / Cojutti)**
- Modello Mariani (IRCCS Gaslini).
- Cojutti 2021
- Cojutti 2023
""")

tab1, tab2 = st.tabs(["➕ Nuova Analisi TDM","📚 Storico Paziente"])

with tab1:
    if not id_paz:
        st.warning("Inserisci un codice paziente nella barra laterale per procedere.")
        st.stop()

    monitoraggi_tutti = st.session_state['storico'].get('monitoraggi', [])
    if monitoraggi_tutti:
        _ultimo = monitoraggi_tutti[-1]
        _gg_def = int(_ultimo.get('eta_giorni', 10950))
        _peso_def = float(_ultimo.get('peso', 70.0))
    else:
        _gg_def = 10950 
        _peso_def = 70.0

    col_e1, col_e2, col_e3, col_peso = st.columns(4)
    with col_e1: eta_anni_in = st.number_input("Età — Anni", 0, 100, min(int(_gg_def // 365), 100))
    with col_e2: eta_mesi_in = st.number_input("Età — Mesi", 0, 11, min(int((_gg_def % 365) // 30), 11))
    with col_e3: eta_giorni_in = st.number_input("Età — Giorni", 0, 30, int((_gg_def % 365) % 30))
    with col_peso: peso = st.number_input("Peso (kg)", 0.3, 300.0, _peso_def)

    gg  = int(round(eta_anni_in*365.25 + eta_mesi_in*30.44 + eta_giorni_in))
    eta = gg / 365.25

    farmaco = st.radio("💊 Farmaco da monitorare:", ["Daptomicina","Ciprofloxacina","Dalbavancin"], horizontal=True)
    modello_dalba = None

    if farmaco=="Dalbavancin":
        modello_dalba = st.selectbox("🧬 Modello PK di riferimento:", list(DALBA_MODELLI.keys()))

    st.markdown("### 🔗 Integrazione Storico Dosi Paziente")
    storia_prev_db, livelli_prev_db, t0_db = carica_storia_precedente(farmaco)
    usa_storico = False

    if storia_prev_db or livelli_prev_db:
        col_info, col_check = st.columns([3,1])
        with col_info: st.success(f"📂 **Trovato storico in sessione:** {len(storia_prev_db)} dosi e {len(livelli_prev_db)} livelli.")
        with col_check: usa_storico = st.checkbox("Continua terapia precedente", value=True)

    st.markdown("### 1. Parametri Clinici")
    col1,col2,col3 = st.columns(3)
    with col1:
        crea = st.number_input("Creatinina Sierica (mg/dL)", 0.05, 25.0, 1.0, 0.05)
        mic = st.number_input(f"MIC {farmaco} (mg/L)", 0.001, 16.0, 0.0625 if farmaco=="Dalbavancin" else 1.0, format="%.4f")
    with col2:
        crrt = st.checkbox("Paziente in CRRT?")
        Sc, Qef = (st.number_input("Sc", 0.01, 0.60, 0.20), st.number_input("Qef (L/h)", 0.5, 6.0, 2.0)) if crrt else (0.2, 2.0)
    with col3:
        if farmaco=="Daptomicina" and gg<6570:
            dialisi = st.selectbox("Funzione Renale / Dialisi:", ["Nessuna","Insufficienza Renale Lieve (GFR 60-90)","Insufficienza Renale Moderata (GFR 30-60)","Insufficienza Renale Grave (GFR 15-30)","ESRD (GFR <15)"] + (["Emodialisi (HD)","Dialisi Peritoneale (CCPD)"] if gg>=4380 else []))
        else:
            dialisi = "Nessuna"

    st.markdown("### 2. Nuove Somministrazioni e Livelli Plasmatici")
    dose_def = {"Daptomicina": 500.0, "Dalbavancin": 1500.0, "Ciprofloxacina": 400.0}[farmaco]

    t_ultima = t0_db + timedelta(hours=storia_prev_db[-1]['ore']) if (usa_storico and storia_prev_db) else datetime.now()
    data_def_dose = t_ultima.date() if usa_storico else date.today()
    ora_def_dose  = t_ultima.time() if usa_storico else datetime.strptime("08:00", "%H:%M").time()

    colA, colB = st.columns(2)
    with colA:
        df_dosi_edit = st.data_editor(pd.DataFrame({"Data (GG/MM/AAAA)": [data_def_dose.strftime("%d/%m/%Y")], "Ora (HH:MM)": [ora_def_dose.strftime("%H:%M")], "Dose (mg)": [dose_def]}), num_rows="dynamic", use_container_width=True)
    with colB:
        df_liv_edit = st.data_editor(pd.DataFrame({"Data (GG/MM/AAAA)": [""], "Ora (HH:MM)": [""], "Concentrazione (mg/L)": [None], "Nota (opz.)": [""]}), num_rows="dynamic", use_container_width=True)

    def _parse_data_ora(d_str, o_str):
        if not d_str or not o_str: return None
        try: d = datetime.strptime(str(d_str).strip(), "%d/%m/%Y").date()
        except:
            try: d = datetime.strptime(str(d_str).strip(), "%Y-%m-%d").date()
            except: return None
        try: t = datetime.strptime(str(o_str).strip()[:5], "%H:%M").time()
        except: return None
        return datetime.combine(d, t)

    nuove_somm = [{'dt': _parse_data_ora(r["Data (GG/MM/AAAA)"], r["Ora (HH:MM)"]), 'dose': float(r["Dose (mg)"])} for _, r in df_dosi_edit.iterrows() if _parse_data_ora(r.get("Data (GG/MM/AAAA)"), r.get("Ora (HH:MM)")) and pd.notna(r.get("Dose (mg)"))]
    nuovi_liv = [{'dt': _parse_data_ora(r["Data (GG/MM/AAAA)"], r["Ora (HH:MM)"]), 'conc': float(r["Concentrazione (mg/L)"]), 'nota': str(r.get("Nota (opz.)", ""))} for _, r in df_liv_edit.iterrows() if _parse_data_ora(r.get("Data (GG/MM/AAAA)"), r.get("Ora (HH:MM)")) and pd.notna(r.get("Concentrazione (mg/L)"))]

    t0_assoluto = t0_db if (usa_storico and storia_prev_db) else (min(s['dt'] for s in nuove_somm) if nuove_somm else datetime.combine(date.today(), datetime.strptime("08:00","%H:%M").time()))

    storia_corrente = costruisci_storia(nuove_somm, t0_assoluto, storia_prev_db if usa_storico else [])
    livelli_correnti = costruisci_livelli(nuovi_liv, t0_assoluto, livelli_prev_db if usa_storico else [])
    dt_val = (max((l['dt'] for l in nuovi_liv), default=None) or max((s['dt'] for s in nuove_somm), default=datetime.now())).strftime('%Y-%m-%d %H:%M')

    if st.button("🚀 Avvia Analisi PK Bayesiana", type="primary", use_container_width=True):
        if not storia_corrente: st.error("Nessuna dose valida.")
        elif len(livelli_correnti) == 0: st.error("Inserisci almeno un livello plasmatico.")
        else:
            with st.spinner("Calcolo MAP Bayesiano in corso…"):
                esegui_analisi(farmaco, id_paz, sesso, gg, peso, dt_val, crea, crrt, Sc, Qef, dialisi, mic, livelli_correnti, storia_corrente, storia_prev_db if usa_storico else [], livelli_prev_db if usa_storico else [], t0_assoluto, salva=True, modello_dalba=modello_dalba)

with tab2:
    monitoraggi = st.session_state['storico'].get('monitoraggi', [])
    if not monitoraggi:
        st.info("Nessuna analisi TDM in questa sessione.")
    else:
        st.dataframe(pd.DataFrame(monitoraggi)[['data_val','farmaco','eta_giorni','peso','dose_mg','n_dosi','n_livelli','CL_ind','AUC24']], use_container_width=True)
        scelta = st.selectbox("Rivisualizza analisi:", ["— Scegli —"] + [f"#{i} — {r['data_val']} [{r.get('farmaco','?')}]" for i, r in enumerate(monitoraggi)])
        if scelta != "— Scegli —":
            idx_sel = int(scelta.split("#")[1].split(" ")[0])
            riga = monitoraggi[idx_sel]
            if st.button("🗑️ Rimuovi questo record"):
                st.session_state['storico']['monitoraggi'].pop(idx_sel)
                st.rerun()
            st.divider()
            esegui_analisi(str(riga.get('farmaco', 'Daptomicina')), id_paz, sesso, int(riga.get('eta_giorni', gg)), float(riga.get('peso', peso)), riga['data_val'], float(riga['creatinina']), bool(riga['crrt']), float(riga.get('Sc', 0.20)), float(riga.get('Qef', 2.0)), str(riga.get('dialisi_tipo', 'Nessuna')), float(riga['mic']), riga.get('livelli', []), riga['storia_dosi'], [], [], datetime.strptime(riga['data_prima_dose'], '%Y-%m-%d %H:%M') if riga.get('data_prima_dose') else datetime.now(), salva=False, modello_dalba=riga.get('modello_dalba'))
