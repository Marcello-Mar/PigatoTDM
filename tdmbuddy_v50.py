"""
TDM Buddy v6.3
==============

FUNZIONALITÀ PRINCIPALI
  - MAP Bayesiano 1-compartimento (Dapto) e 2-compartimento (Cipro/Dalba)
  - Integrazione storia dosi cumulativa e Bayesian Updating
  - Archivio clinico SQLite / Gestione JSON State
  - Grafici PK Interattivi via Plotly (zoom, pan, toggle tracce)
  - Motore blindato: gestione sicura dei tipi (isinstance, .get) per evitare crash
"""

import streamlit as st
import numpy as np
from scipy.optimize import minimize
import plotly.graph_objects as go
import pandas as pd
from datetime import datetime, date, timedelta
import json

# ============================================================
# STORAGE
# ============================================================

def storico_vuoto():
    return {'codice_paziente': '', 'sesso': 'Maschio', 'monitoraggi': []}

if 'storico' not in st.session_state:
    st.session_state['storico'] = storico_vuoto()
if 'farmaco_attivo' not in st.session_state:
    st.session_state['farmaco_attivo'] = None

# ============================================================
# PARAMETRI PK E TARGETS
# ============================================================

SORALUCE = {'CLNR': 0.160, 'CLR_theta': 0.366, 'Clcr_ref': 49.0, 'V': 12.30, 'omega_CL': 0.367, 'omega_V': 0.278, 'sigma': 0.123, 'Sc_default': 0.20, 'nome': 'Soraluce et al. 2018 [Adulto ICU]'}
YE2025 = {'CL':19.6, 'SD':5.5, 'V':0.17, 'nome': 'Ye L, et al. 2025 [<1 anno, PBPK]'}
WEI = {
    '3-6m'  : {'CL':19.7,'SD':5.5, 'V':0.17,'nome':'Wei et al. 2020 [3-6 mesi]'},
    '7-12m' : {'CL':19.6,'SD':5.8, 'V':0.17,'nome':'Wei et al. 2020 [7-12 mesi]'},
    '13-24m': {'CL':21.8,'SD':3.0, 'V':0.11,'nome':'Wei et al. 2020 [13-24 mesi]'},
    '2-6y'  : {'CL':20.0,'SD':6.1, 'V':0.11,'nome':'Wei et al. 2020 [2-6 anni]'},
    '7-11y' : {'CL':17.0,'SD':5.7, 'V':0.11,'nome':'Wei et al. 2020 [7-11 anni]'},
    '12-17y': {'CL':11.0,'SD':2.1, 'V':0.11,'nome':'Wei et al. 2020 [12-17 anni]'},
}
YE2022 = {'12-17y': {'mild':0.98,'moderate':1.08,'severe':1.24,'esrd':1.33}, '7-11y': {'mild':1.00,'moderate':1.19,'severe':1.49,'esrd':1.73}, '2-6y': {'mild':1.02,'moderate':1.28,'severe':1.74,'esrd':2.20}, '1-2y': {'mild':1.04,'moderate':1.29,'severe':1.72,'esrd':2.13}}
LIM = {'Vss':5.90,'CLren':0.023,'CLnr':0.158,'CL_HD':0.619, 'omega_V1':0.112,'omega_knr':0.088,'sigma':0.085,'peso_ref':40.0, 'nome': 'Lim SY, et al. 2021 [Adolescente HD/CCPD]'}
PD_DAPTO = {'AUC_MIC_staph':666,'AUC_MIC_entero':143,'Cmin_tox':24.3}

CIPRO = {'CL_kg':0.50,'V1_kg':0.55,'Q_kg':0.30,'V2_kg':1.90, 'omega_CL':0.37,'omega_V1':0.32,'omega_Q':0.25,'omega_V2':0.30, 'sigma':0.15,'fu':0.67, 'nome':'Schlender 2018 / Balbas-Martinez 2019'}
CIPRO_PED = {(91, 365): 1.40, (365, 730): 1.30, (730,2190): 1.20, (2190,4380):1.10, (4380,6570):1.00}
CIPRO_FRAZ_RENAL = 0.70
PD_CIPRO = {'AUC_MIC_target':125,'AUC_alta':400}
IRC_MAP = {"Lieve (GFR 60-90)": 'mild', "Moderata (GFR 30-60)": 'moderate', "Grave (GFR 15-30)": 'severe', "ESRD (GFR <15)": 'esrd'}

DALBA_PARAMS = {'CL_typ70': 0.0305, 'V1_typ70': 4.81, 'Q_typ70': 0.0081, 'V2_typ70': 4.81, 'omega_CL': 0.198, 'sigma': 0.375, 't_inf': 0.5, 'fu': 0.07, 'nome': "Modello Mariani (IRCCS Gaslini)"}
COJUTTI_PARAMS = {'CL_typ70': 0.106, 'V1_typ70': 17.4, 'Q_typ70': 0.103, 'V2_typ70': 15.1, 'omega_CL': 0.362, 'sigma': 0.18, 't_inf': 0.5, 'fu': 0.07, 'nome': "Cojutti PG, et al. AAC 2021 (Adulti)"}
DALBA_MODELLI = {'Mariani (pediatrico, allometrico — Gaslini)': DALBA_PARAMS, 'Cojutti 2021 (adulti, osteoarticolare — AAC)': COJUTTI_PARAMS}
DALBA_MIC_EUCAST = {'MIC90 osteoarticolare (0.0625)': 0.0625, 'EUCAST attuale S.aureus (0.125)': 0.125, 'USCAST S.aureus (0.25)': 0.25}
DALBA_TARGETS = {'ABSSSI semplice (fAUC24h/MIC >= 27.1)': 27.1, 'Infezione complessa (fAUC24h/MIC >= 111.1)': 111.1}

def dalba_soglia_totale(target, mic, f): return (target * mic) / (24.0 * f)
def dalba_soglia_trough(mic, f): return (4.0 * mic) / f

# ============================================================
# MOTORI PK
# ============================================================

def pk_1cmt(CL, V, storia, t_obs, t_inf=0.5):
    if CL <= 0 or V <= 0: return 0.0
    ke = CL / V; C = 0.0
    for s in storia:
        dt = t_obs - s.get('ore', 0)
        if dt <= 0: continue
        D = s.get('dose', 0)
        if dt <= t_inf: C += (D/(V*ke*t_inf))*(1.0-np.exp(-ke*dt))
        else: C += (D/(V*ke*t_inf))*(1.0-np.exp(-ke*t_inf))*np.exp(-ke*(dt-t_inf))
    return max(C, 0.0)

def pk_2cmt(CL, V1, Q, V2, storia, t_obs, t_inf=0.5):
    if CL<=0 or V1<=0 or Q<=0 or V2<=0: return 0.0
    k10=CL/V1; k12=Q/V1; k21=Q/V2
    s_=k10+k12+k21; d_=max(s_**2-4*k10*k21,0)
    alp=(s_+np.sqrt(d_))/2; bet=(s_-np.sqrt(d_))/2
    if abs(alp-bet)<1e-10: return pk_1cmt(CL,V1+V2,storia,t_obs,t_inf)
    A=(k21-alp)/((bet-alp)*V1); B=(k21-bet)/((alp-bet)*V1); C=0.0
    for s in storia:
        dt=t_obs-s.get('ore', 0)
        if dt <= 0: continue
        D=s.get('dose', 0); k0=D/t_inf
        if dt<=t_inf: C += k0 * ((A/alp)*(1-np.exp(-alp*dt)) + (B/bet)*(1-np.exp(-bet*dt)))
        else: C += k0 * ((A/alp)*(1-np.exp(-alp*t_inf))*np.exp(-alp*(dt-t_inf)) + (B/bet)*(1-np.exp(-bet*t_inf))*np.exp(-bet*(dt-t_inf)))
    return max(C,0.0)

def metriche_esposizione_24h(pk_func, params, storia, t_inf, n_punti=300):
    t_ultima_dose = storia[-1].get('ore', 0)
    t_grid = np.linspace(t_ultima_dose, t_ultima_dose + 24, n_punti)
    c_grid = np.array([pk_func(*params, storia, t, t_inf) for t in t_grid])
    auc = float(np.trapezoid(c_grid, x=t_grid)) if hasattr(np, 'trapezoid') else float(np.trapz(c_grid, x=t_grid))
    return auc, float(c_grid.max()), float(c_grid.min())

def tau_mediano(storia):
    if len(storia) < 2: return None
    ore = sorted(s.get('ore', 0) for s in storia)
    diffs = sorted([ore[i+1]-ore[i] for i in range(len(ore)-1)])
    return diffs[len(diffs)//2] if len(diffs) % 2 == 1 else (diffs[len(diffs)//2-1]+diffs[len(diffs)//2])/2.0

def calcola_clcr(eta_a, peso, crea, sesso):
    if crea<=0: return 60.0
    clcr=((140-eta_a)*peso)/(72*crea)
    return max(clcr*0.85 if sesso=="Femmina" else clcr, 1.0)

# ============================================================
# PRIOR DI POPOLAZIONE (STRUTTURATI A DIZIONARIO)
# ============================================================

def prior_dapto(id_paziente, sesso, gg, peso, crea, crrt, Sc, Qef, dialisi):
    eta = gg/365.25; note = ""
    if gg<365:
        p=YE2025; CL_pop=(p['CL']*peso)/1000; V_pop=p['V']*peso
        oCL=p['SD']/p['CL']; oV=0.30; sig=0.12; t_inf=1.0; nm=p['nome']
        if crrt: CL_pop+=Sc*Qef; note+=f"CRRT: CLEC={Sc*Qef:.3f} L/h. "
    elif 4380<=gg<6570 and dialisi in ("Emodialisi (HD)","Dialisi Peritoneale (CCPD)"):
        sc=peso/LIM['peso_ref']
        CL_pop=(LIM['CLren']+LIM['CLnr'])*sc; V_pop=LIM['Vss']*sc; nm=LIM['nome']
        oCL, oV, sig, t_inf = LIM['omega_knr'], LIM['omega_V1'], LIM['sigma'], 0.5
        if dialisi=="Emodialisi (HD)": note+="Somministrare DOPO HD."
    elif 91<=gg<6570:
        ag = '12-17y'
        if gg<213: ag='3-6m'
        elif gg<365: ag='7-12m'
        elif gg<730: ag='13-24m'
        elif gg<2190: ag='2-6y'
        elif gg<4380: ag='7-11y'
        p=WEI[ag]; CL_pop=(p['CL']*peso)/1000; V_pop=p['V']*peso
        oCL=p['SD']/p['CL']; oV=0.25; sig=0.12; t_inf=1.0 if gg<730 else 0.5; nm=p['nome']
        age_ye = ag if ag in YE2022 else ('2-6y' if ag in ('13-24m','2-6y') else '1-2y')
        if dialisi in IRC_MAP and age_ye in YE2022: CL_pop /= YE2022[age_ye][IRC_MAP[dialisi]]
        if crrt: CL_pop+=Sc*Qef
    else:
        s=SORALUCE; clcr=calcola_clcr(eta,peso,crea,sesso)
        CL_pop=s['CLNR']+s['CLR_theta']*(clcr/s['Clcr_ref']); V_pop=s['V']
        nm=f"{s['nome']} (ClCr={clcr:.0f} mL/min)"
        oCL, oV, sig, t_inf = s['omega_CL'], s['omega_V'], s['sigma'], 0.5
        if crrt: CL_pop+=Sc*Qef; note+=f"CRRT: CLEC={Sc*Qef:.3f} L/h."
    
    CL_std, V_std = CL_pop, V_pop

    monitoraggi = st.session_state['storico'].get('monitoraggi', [])
    candidati = [m for m in monitoraggi if isinstance(m, dict) and str(m.get('farmaco','')).startswith("Dapto")]
    if candidati:
        last_tdm = candidati[-1]
        cl_hist = float(last_tdm.get('CL_ind', CL_pop))
        cr_hist = float(last_tdm.get('creatinina', 1.0))
        if cr_hist > 0 and crea > 0: CL_pop = cl_hist * (cr_hist / crea)
        else: CL_pop = cl_hist
        V_pop = float(last_tdm.get('V_ind', V_pop))
        oCL *= 0.5; oV *= 0.5
        note += f"⚠️ [BAYESIAN UPDATING] Usa storico: CL_adj={CL_pop:.2f}, V={V_pop:.2f}. "

    return {'CL_p':CL_pop, 'V_p':V_pop, 'Q_p':0.0, 'V2_p':0.0, 'oCL':oCL, 'oV1':oV, 'oQ':0.0, 'oV2':0.0, 'sig':sig, 't_inf':t_inf, 'nm':nm, 'note':note, 'CL_std':CL_std, 'V_std':V_std, 'Q_std':0.0, 'V2_std':0.0}

def prior_cipro(id_paziente, sesso, gg, peso, crea, crrt, dialisi):
    eta=gg/365.25; note=""; c=CIPRO; fc=1.0
    for (lo,hi),f in CIPRO_PED.items():
        if lo<=gg<hi: fc=f; break
    CL_pop=c['CL_kg']*peso*fc; V1_pop=c['V1_kg']*peso; Q_pop=c['Q_kg']*peso; V2_pop=c['V2_kg']*peso
    oCL, oV1, oQ, oV2, sig = c['omega_CL'], c['omega_V1'], c['omega_Q'], c['omega_V2'], c['sigma']

    clcr=calcola_clcr(eta,peso,crea,sesso)
    if clcr>=30:
        fr=min(clcr/120.0,1.0)
        CL_pop=(c['CL_kg']*peso*fc*CIPRO_FRAZ_RENAL*fr)+(c['CL_kg']*peso*fc*(1-CIPRO_FRAZ_RENAL))
    else:
        CL_pop=(c['CL_kg']*peso*fc*CIPRO_FRAZ_RENAL*0.20)+(c['CL_kg']*peso*fc*(1-CIPRO_FRAZ_RENAL)*0.63)

    nome=f"{c['nome']} (ClCr={clcr:.0f} mL/min)"
    CL_std, V1_std, Q_std, V2_std = CL_pop, V1_pop, Q_pop, V2_pop

    monitoraggi = st.session_state['storico'].get('monitoraggi', [])
    candidati = [m for m in monitoraggi if isinstance(m, dict) and str(m.get('farmaco','')).startswith("Cipro")]
    if candidati:
        last_tdm = candidati[-1]
        cl_hist = float(last_tdm.get('CL_ind', CL_pop))
        cr_hist = float(last_tdm.get('creatinina', 1.0))
        if cr_hist > 0 and crea > 0: CL_pop = cl_hist * (cr_hist / crea)
        else: CL_pop = cl_hist
        V1_pop = float(last_tdm.get('V_ind', V1_pop))
        Q_pop = float(last_tdm.get('Q_ind', Q_pop)) if last_tdm.get('Q_ind') else Q_pop
        V2_pop = float(last_tdm.get('V2_ind', V2_pop)) if last_tdm.get('V2_ind') else V2_pop
        
        oCL *= 0.5; oV1 *= 0.5
        note += f"⚠️ [BAYESIAN UPDATING] Usa storico: CL_adj={CL_pop:.2f}, V1={V1_pop:.2f}. "

    return {'CL_p':CL_pop, 'V_p':V1_pop, 'Q_p':Q_pop, 'V2_p':V2_pop, 'oCL':oCL, 'oV1':oV1, 'oQ':oQ, 'oV2':oV2, 'sig':sig, 't_inf':1.0, 'nm':nome, 'note':note, 'CL_std':CL_std, 'V_std':V1_std, 'Q_std':Q_std, 'V2_std':V2_std}

def prior_dalba(id_paziente, peso, modello_scelto):
    p = DALBA_MODELLI.get(modello_scelto, list(DALBA_MODELLI.values())[0])
    note = ""

    if str(modello_scelto).startswith('Mariani'):
        fw = peso / 70.0
        CL_pop = p['CL_typ70'] * (fw ** 0.75)
        V1_pop = p['V1_typ70'] * fw
        Q_pop  = p['Q_typ70']  * (fw ** 0.75)
        V2_pop = p['V2_typ70'] * fw
    else:
        CL_pop, V1_pop, Q_pop, V2_pop = p['CL_typ70'], p['V1_typ70'], p['Q_typ70'], p['V2_typ70']

    CL_std, V1_std, Q_std, V2_std = CL_pop, V1_pop, Q_pop, V2_pop
    oCL = p['omega_CL']

    monitoraggi = st.session_state['storico'].get('monitoraggi', [])
    candidati = [m for m in monitoraggi if isinstance(m, dict) and str(m.get('farmaco','')).startswith("Dalba")]
    if candidati:
        last_tdm = candidati[-1]
        CL_pop = float(last_tdm.get('CL_ind', CL_pop))
        oCL *= 0.5
        note += f"⚠️ [BAYESIAN UPDATING] Usa storico per fittare: CL_adj={CL_pop:.4f} L/h. "

    return {'CL_p':CL_pop, 'V_p':V1_pop, 'Q_p':Q_pop, 'V2_p':V2_pop, 'oCL':oCL, 'oV1':0.0, 'oQ':0.0, 'oV2':0.0, 'sig':p['sigma'], 't_inf':p['t_inf'], 'nm':p['nome'], 'note':note, 'CL_std':CL_std, 'V_std':V1_std, 'Q_std':Q_std, 'V2_std':V2_std}

# ============================================================
# MAP BAYESIANO
# ============================================================

def map_1cmt(CL_p, V_p, oCL, oV, sig, storia, livelli, t_inf):
    if len(livelli) == 0: return CL_p, V_p, "Prior di popolazione"
    def f(p):
        CL_, V_ = np.exp(p[0]), np.exp(p[1])
        pr = (((p[0]-np.log(CL_p))**2)/oCL**2 + ((p[1]-np.log(V_p))**2)/oV**2)
        lk = sum(((np.log(max(l.get('conc',0),1e-6)) - np.log(max(pk_1cmt(CL_, V_, storia, l.get('ore',0), t_inf),1e-6)))**2) / sig**2 for l in livelli)
        return 0.5*(pr+lk)
    bnd = [(np.log(CL_p*0.05), np.log(CL_p*8)), (np.log(V_p*0.15), np.log(V_p*5))]
    r = minimize(f, [np.log(CL_p), np.log(V_p)], method='L-BFGS-B', bounds=bnd, options={'ftol':1e-12})
    return float(np.exp(r.x[0])), float(np.exp(r.x[1])), f"MAP Bayesiano 1-cmt su {len(livelli)} livelli"

def map_2cmt(CL_p, V1_p, Q_p, V2_p, oCL, oV1, oQ, oV2, sig, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0: return CL_p, V1_p, Q_p, V2_p, "Prior di popolazione"

    OQ_eff  = 0.10 + (oQ - 0.10)*min(1.0, n/6.0)
    OV2_eff = 0.10 + (oV2 - 0.10)*min(1.0, n/6.0)

    def f(p):
        CL_, V1_, Q_, V2_ = np.exp(p[0]), np.exp(p[1]), np.exp(p[2]), np.exp(p[3])
        pr = (((p[0]-np.log(CL_p))**2)/oCL**2 + ((p[1]-np.log(V1_p))**2)/oV1**2 + ((p[2]-np.log(Q_p))**2)/OQ_eff**2 + ((p[3]-np.log(V2_p))**2)/OV2_eff**2)
        lk = sum(((np.log(max(l.get('conc',0),1e-6)) - np.log(max(pk_2cmt(CL_, V1_, Q_, V2_, storia, l.get('ore',0), t_inf),1e-6)))**2) / sig**2 for l in livelli)
        return 0.5*(pr+lk)

    bnd = [(np.log(CL_p*0.05), np.log(CL_p*10)), (np.log(V1_p*0.1), np.log(V1_p*6)), (np.log(Q_p*0.3), np.log(Q_p*3)), (np.log(V2_p*0.3), np.log(V2_p*3))]
    r = minimize(f, [np.log(CL_p), np.log(V1_p), np.log(Q_p), np.log(V2_p)], method='L-BFGS-B', bounds=bnd, options={'ftol':1e-12})
    return float(np.exp(r.x[0])), float(np.exp(r.x[1])), float(np.exp(r.x[2])), float(np.exp(r.x[3])), f"MAP Bayesiano 2-cmt su {n} livelli"

def map_dalba(CL_p, V1_p, Q_p, V2_p, oCL, sig, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0: return CL_p, "Prior di popolazione"
    def f(p):
        CL_ = np.exp(p[0])
        prior = ((p[0] - np.log(CL_p)) ** 2) / oCL ** 2
        lik = sum(((np.log(max(l.get('conc',0),1e-6)) - np.log(max(pk_2cmt(CL_, V1_p, Q_p, V2_p, storia, l.get('ore',0), t_inf),1e-6))) ** 2) / sig ** 2 for l in livelli)
        return 0.5 * (prior + lik)
    r = minimize(f, [np.log(CL_p)], method='L-BFGS-B', bounds=[(np.log(CL_p*0.05), np.log(CL_p*10.0))], options={'ftol': 1e-12})
    return float(np.exp(r.x[0])), f"MAP Bayesiano su CL, {n} livelli (V1/Q/V2 fissi)"

def trova_tempo_soglia(CL, V1, Q, V2, storia, soglia, t_inf, t_partenza=0.0, t_max=2000.0, passo=1.0):
    if pk_2cmt(CL, V1, Q, V2, storia, t_partenza, t_inf) < soglia: return t_partenza
    t_prev = t_partenza; t = t_partenza + passo
    while t <= t_max:
        if pk_2cmt(CL, V1, Q, V2, storia, t, t_inf) < soglia:
            lo, hi = t_prev, t
            for _ in range(60):
                mid = (lo + hi) / 2.0
                if pk_2cmt(CL, V1, Q, V2, storia, mid, t_inf) >= soglia: lo = mid
                else: hi = mid
            return hi
        t_prev = t; t += passo
    return None

# ============================================================
# GRAFICI PLOTLY INTERATTIVI
# ============================================================

def disegna_grafico_plotly(farmaco, storia, livelli, CL_ind, V_ind, Q_ind, V2_ind, CL_std, V_std, Q_std, V2_std, t_inf, t_ultimo_storico=-1, t_fine_extra=12):
    t_ultimo_livello = max((l.get('ore',0) for l in livelli), default=0)
    t_fine = max(t_ultimo_livello, max((s.get('ore',0) for s in storia), default=0) + 12) + t_fine_extra
    t_arr = np.linspace(0, t_fine, 800)

    if str(farmaco).startswith("Dapto"):
        C_ind = [pk_1cmt(CL_ind, V_ind, storia, t, t_inf) for t in t_arr]
        C_pop = [pk_1cmt(CL_std, V_std, storia, t, t_inf) for t in t_arr]
    else:
        C_ind = [pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia, t, t_inf) for t in t_arr]
        C_pop = [pk_2cmt(CL_std, V_std, Q_std, V2_std, storia, t, t_inf) for t in t_arr]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=t_arr, y=C_pop, mode='lines', line=dict(dash='dash', color='#94A3B8', width=2), name='Popolazione Standard (Baseline)'))
    fig.add_trace(go.Scatter(x=t_arr, y=C_ind, mode='lines', line=dict(color='#0F172A', width=3), name='Profilo Individuale Paziente (MAP)'))

    max_c = max(max(C_ind, default=100), max(C_pop, default=100))

    if t_ultimo_storico >= 0:
        fig.add_vrect(x0=0, x1=t_ultimo_storico, fillcolor="#6EA8D8", opacity=0.1, layer="below", line_width=0, name="Terapia precedente")

    x_dosi_vecchie, y_dosi_vecchie, x_dosi_nuove, y_dosi_nuove = [], [], [], []
    for s in storia:
        ore_s = s.get('ore', 0)
        if ore_s <= t_ultimo_storico:
            x_dosi_vecchie.extend([ore_s, ore_s, None]); y_dosi_vecchie.extend([0, max_c, None])
        else:
            x_dosi_nuove.extend([ore_s, ore_s, None]); y_dosi_nuove.extend([0, max_c, None])

    if x_dosi_vecchie: fig.add_trace(go.Scatter(x=x_dosi_vecchie, y=y_dosi_vecchie, mode='lines', line=dict(color='#6EA8D8', dash='dot', width=1), name='Dosi Pregresse', hoverinfo='none'))
    if x_dosi_nuove: fig.add_trace(go.Scatter(x=x_dosi_nuove, y=y_dosi_nuove, mode='lines', line=dict(color='#475569', dash='dot', width=1), name='Nuove Dosi', hoverinfo='none'))

    if livelli:
        x_liv = [l.get('ore',0) for l in livelli]
        y_liv = [l.get('conc',0) for l in livelli]
        testi = [f"Misurazione<br>t={l.get('ore',0):.1f}h<br>C={l.get('conc',0):.2f} mg/L<br>{l.get('nota','')}" for l in livelli]
        fig.add_trace(go.Scatter(x=x_liv, y=y_liv, mode='markers', marker=dict(symbol='star', size=14, color='#16A34A', line=dict(width=1, color='white')), name='Livelli misurati', text=testi, hoverinfo='text'))

    if str(farmaco).startswith("Dapto"): fig.add_hline(y=24.3, line_dash="dashdot", line_color="#DC2626", annotation_text="Soglia tossicità CPK (24.3 mg/L)", annotation_position="top right")

    fig.update_layout(xaxis_title="Ore dall'inizio della terapia (t=0)", yaxis_title=f"Concentrazione plasmatica {farmaco} (mg/L)", plot_bgcolor='#FAFAFA', paper_bgcolor='#FAFAFA', hovermode="x unified", legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1), margin=dict(l=40, r=40, t=40, b=40))
    return fig

def disegna_forecast_plotly(t_arr_sim, c_sim, soglia_target, t_new_ore, dose_new, storia_completa):
    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(x=t_arr_sim, y=c_sim, mode='lines', line=dict(color='#0F172A', width=2.5), name='Profilo proiettato'))
    fig2.add_hline(y=soglia_target, line_dash="dashdot", line_color="#DC2626", annotation_text=f"Soglia target ({soglia_target:.2f} mg/L)", annotation_position="top right")
    fig2.add_vline(x=t_new_ore, line_dash="dash", line_color="#16A34A", annotation_text=f"Nuova dose ({dose_new:.0f} mg)")

    x_dosi, y_dosi = [], []
    max_c_sim = max(c_sim, default=100)
    for s in storia_completa:
        ore_s = s.get('ore', 0)
        x_dosi.extend([ore_s, ore_s, None])
        y_dosi.extend([0, max_c_sim, None])
    
    if x_dosi: fig2.add_trace(go.Scatter(x=x_dosi, y=y_dosi, mode='lines', line=dict(color='#94A3B8', dash='dot', width=1), name='Dosi passate', hoverinfo='none'))
    fig2.update_layout(xaxis_title="Ore dall'inizio terapia", yaxis_title="Dalbavancina (mg/L)", plot_bgcolor='#FAFAFA', paper_bgcolor='#FAFAFA', hovermode="x unified", legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1), margin=dict(l=40, r=40, t=40, b=40))
    return fig2

# ============================================================
# UI PRINCIPALE - LAYOUT E SIDEBAR
# ============================================================

st.set_page_config(page_title="TDM Buddy v6.3", layout="wide", page_icon="💊")

st.title("💊 TDM Buddy v6.3")
st.caption("Piattaforma Clinica di Farmacocinetica Bayesiana (Daptomicina, Ciprofloxacina, Dalbavancina)")

# SIDEBAR: Gestione Dati e Sessione
st.sidebar.header("🛠️ Gestione Sessione")

if st.sidebar.button("🔄 Nuovo Paziente / Svuota Dati", type="primary", use_container_width=True):
    st.session_state['storico'] = storico_vuoto()
    st.session_state['farmaco_attivo'] = None
    st.rerun()

st.sidebar.markdown("---")
file_up = st.sidebar.file_uploader("📥 Carica storico paziente (JSON)", type=["json"])
if file_up is not None:
    file_hash = hash(file_up.getvalue())
    if st.session_state.get('uploaded_file_hash') != file_hash:
        try:
            dati_caricati = json.loads(file_up.getvalue())
            st.session_state['storico'] = dati_caricati
            st.session_state['uploaded_file_hash'] = file_hash
            
            # Seleziona automaticamente il farmaco dal file
            if 'farmaco' in dati_caricati and dati_caricati['farmaco'] in FARMACI:
                st.session_state['farmaco_attivo'] = dati_caricati['farmaco']
            elif dati_caricati.get('monitoraggi'):
                last_f = dati_caricati['monitoraggi'][-1].get('farmaco')
                if last_f in FARMACI: st.session_state['farmaco_attivo'] = last_f

            st.sidebar.success(f"✅ Storico caricato: {len(dati_caricati.get('monitoraggi', []))} TDM trovati.")
            st.rerun()
        except Exception as e: 
            st.sidebar.error("File non valido.")

id_paz = st.sidebar.text_input("Codice Paziente (no dati sensibili)", value=st.session_state['storico'].get('codice_paziente', ''))
sesso = st.sidebar.selectbox("Sesso biologico", ["Maschio", "Femmina"], index=0 if st.session_state['storico'].get('sesso', 'Maschio') == 'Maschio' else 1)
st.session_state['storico']['codice_paziente'] = id_paz
st.session_state['storico']['sesso'] = sesso

n_rec = len(st.session_state['storico'].get('monitoraggi', []))
if n_rec > 0:
    st.sidebar.markdown("---")
    st.sidebar.info(f"📋 {n_rec} analisi salvate nella sessione corrente.")
    ultimo_farmaco = st.session_state['storico']['monitoraggi'][-1].get('farmaco', 'Multi')
    
    df_sb = pd.DataFrame(st.session_state['storico']['monitoraggi'])
    cols_sb = [c for c in ['data_val', 'farmaco', 'creatinina', 'mic'] if c in df_sb.columns]
    st.sidebar.dataframe(df_sb[cols_sb].tail(4), use_container_width=True, hide_index=True)
    
    st.sidebar.download_button("📤 Scarica Storico Aggiornato", data=json.dumps(st.session_state['storico'], indent=2, ensure_ascii=False).encode('utf-8'), file_name=f"TDMBuddy_{ultimo_farmaco}_{id_paz}.json", mime="application/json", type="primary", use_container_width=True)

# ============================================================
# SELEZIONE MOLECOLA - STEP 1 (Bloccante)
# ============================================================

if st.session_state['farmaco_attivo'] is None:
    st.markdown("## Passo 1 — Seleziona la molecola")
    st.caption("La scelta è vincolante per questa sessione. Per cambiare molecola dovrai ricaricare o svuotare i dati.")

    cols = st.columns(3)
    for col, farmaco_scelta in zip(cols, FARMACI):
        info = FARMACO_INFO[farmaco_scelta]
        with col:
            with st.container(border=True):
                st.markdown(f"<h3 style='text-align: center; margin-bottom: 0;'>{info['icona']} {farmaco_scelta}</h3>", unsafe_allow_html=True)
                st.markdown(f"<p style='text-align: center; color: gray; margin-top: 0;'>{info['descrizione']}</p>", unsafe_allow_html=True)
                if st.button("Seleziona", key=f"sel_{farmaco_scelta}", use_container_width=True, type="primary"):
                    st.session_state['farmaco_attivo'] = farmaco_scelta
                    st.session_state['storico']['farmaco'] = farmaco_scelta
                    st.rerun()
    st.stop()

farmaco = st.session_state['farmaco_attivo']
col_badge, col_cambia = st.columns([5, 1])
with col_badge: st.info(f"{FARMACO_INFO[farmaco]['icona']} **Molecola in corso: {farmaco}**")
with col_cambia:
    if st.button("🔄 Cambia molecola", use_container_width=True):
        st.session_state['farmaco_attivo'] = None
        st.session_state['storico'] = storico_vuoto()
        st.rerun()


# ============================================================
# MAIN PAGE - TABS
# ============================================================

tab1, tab2, tab3 = st.tabs(["➕ Analisi TDM", "📚 Storico Sessione", "ℹ️ Riferimenti Scientifici"])

with tab3:
    st.markdown("""
    ### Modelli Clinici Integrati
    **Daptomicina**
    - Soraluce A, et al. *Int J Antimicrob Agents.* 2018;52(2):158-165.
    - Wei XC, et al. *J Clin Pharmacol.* 2020;60(6):768-774.
    - Ye L, et al. *Front Pharmacol.* 2022;13:838599.
    - Ye L, Zhou H, et al. *Eur J Pharm Sci.* 2025;208:107058.
    - Lim SY, et al. *J Pediatr Pharmacol Ther.* 2021;26(2):123-132.
    - Bhavnani SM, et al. *Clin Infect Dis.* 2010;50(12):1568-1574.

    **Ciprofloxacina**
    - Schlender JF, et al. *Clin Pharmacokinet.* 2018;57(12):1613-1634.
    - Balbas-Martinez V, et al. *Eur J Pharm Sci.* 2019;128:171-179.
    - Forrest A, et al. *Antimicrob Agents Chemother.* 1993;37(5):1065-1072.
    - Ambrose PG, et al. *Clin Infect Dis.* 2007;44(S2):S79-S86.

    **Dalbavancina**
    - Modello Mariani: interno IRCCS Gaslini (scaling allometrico).
    - Cojutti PG, et al. *Antimicrob Agents Chemother.* 2021;65(5):e02260-20.
    """)

with tab1:
    if not id_paz:
        st.warning("👈 Inserisci il codice paziente nella barra laterale sinistra per procedere.")
        st.stop()

    monitoraggi_farmaco = [m for m in st.session_state['storico'].get('monitoraggi', []) if isinstance(m, dict) and str(m.get('farmaco','')).startswith(farmaco[:8])]
    _ultimo = monitoraggi_farmaco[-1] if monitoraggi_farmaco else {}

    _gg_def = int(_ultimo.get('eta_giorni', 10950))
    _peso_def = float(_ultimo.get('peso', 70.0))
    def_crea = float(_ultimo.get('creatinina', 1.0))
    def_mic = float(_ultimo.get('mic', 0.0625 if farmaco=="Dalbavancina" else 1.0))

    if monitoraggi_farmaco: st.success(f"📋 Dati storici recuperati dall'ultima analisi di {farmaco}.")

    # 1. Dati Paziente
    with st.container(border=True):
        st.markdown("#### 1. Dati Biometrici & Modello")
        col_e1, col_e2, col_e3, col_peso = st.columns(4)
        with col_e1: eta_anni_in = st.number_input("Età — Anni", 0, 100, min(int(_gg_def // 365), 100))
        with col_e2: eta_mesi_in = st.number_input("Età — Mesi", 0, 11, min(int((_gg_def % 365) // 30), 11))
        with col_e3: eta_giorni_in = st.number_input("Età — Giorni", 0, 30, int((_gg_def % 365) % 30))
        with col_peso: peso = st.number_input("Peso (kg)", 0.3, 300.0, _peso_def)

        gg = int(round(eta_anni_in*365.25 + eta_mesi_in*30.44 + eta_giorni_in))
        eta = gg / 365.25

        modello_dalba = None
        if farmaco=="Dalbavancina":
            opzioni_mod = list(DALBA_MODELLI.keys())
            mod_salvato = _ultimo.get('modello_dalba')
            idx_mod = opzioni_mod.index(mod_salvato) if mod_salvato in opzioni_mod else 0
            modello_dalba = st.selectbox("🧬 Modello PK di riferimento:", opzioni_mod, index=idx_mod)

    # 2. Parametri Clinici
    with st.container(border=True):
        st.markdown("#### 2. Dati Clinici Attuali")
        col1,col2,col3 = st.columns(3)
        with col1:
            crea = st.number_input("Creatinina Sierica (mg/dL)", 0.05, 25.0, def_crea, 0.05)
            mic = st.number_input(f"MIC {farmaco} (mg/L)", 0.001, 16.0, def_mic, format="%.4f")
        with col2:
            crrt = st.checkbox("Paziente in CRRT?")
            Sc, Qef = (st.number_input("Sc", 0.01, 0.60, 0.20), st.number_input("Qef (L/h)", 0.5, 6.0, 2.0)) if crrt else (0.2, 2.0)
        with col3:
            if farmaco=="Daptomicina" and gg<6570:
                dialisi = st.selectbox("Funzione Renale / Dialisi:", ["Nessuna","Lieve (GFR 60-90)","Moderata (GFR 30-60)","Grave (GFR 15-30)","ESRD (GFR <15)"] + (["Emodialisi (HD)","Dialisi Peritoneale (CCPD)"] if gg>=4380 else []))
            else: dialisi = "Nessuna"

    # 3. Timeline
    st.markdown("### 3. Timeline Dosi e Prelievi")
    st.info("💡 **Le tabelle mostrano l'INTERA storia clinica del paziente.** Modifica le righe esistenti o aggiungine di nuove per la valutazione odierna.")

    dose_def = {"Daptomicina": 500.0, "Dalbavancina": 1500.0, "Ciprofloxacina": 400.0}.get(farmaco, 500.0)

    dosi_list = []
    t0_db = datetime.now() 
    t_ultimo_storico = -1
    
    if monitoraggi_farmaco and 'storia_dosi' in _ultimo and isinstance(_ultimo['storia_dosi'], list) and len(_ultimo['storia_dosi']) > 0:
        t0_db = datetime.strptime(_ultimo.get('data_prima_dose', t0_db.strftime('%Y-%m-%d %H:%M')), '%Y-%m-%d %H:%M')
        for s in _ultimo['storia_dosi']:
            dt = t0_db + timedelta(hours=s.get('ore', 0))
            dosi_list.append({"Data (GG/MM/AAAA)": dt.strftime("%d/%m/%Y"), "Ora (HH:MM)": dt.strftime("%H:%M"), "Dose (mg)": s.get('dose', 0.0)})
        
        t_ultimo_storico = _ultimo['storia_dosi'][-1].get('ore', 0)
        
        if farmaco != "Dalbavancina":
            tau = int(_ultimo.get('tau_h', 24)) if _ultimo.get('tau_h') else 24
            next_dt = t0_db + timedelta(hours=t_ultimo_storico + tau)
            dosi_list.append({"Data (GG/MM/AAAA)": next_dt.strftime("%d/%m/%Y"), "Ora (HH:MM)": next_dt.strftime("%H:%M"), "Dose (mg)": _ultimo['storia_dosi'][-1].get('dose', dose_def)})
    else:
        dosi_list.append({"Data (GG/MM/AAAA)": date.today().strftime("%d/%m/%Y"), "Ora (HH:MM)": "08:00", "Dose (mg)": dose_def})

    liv_list = []
    if monitoraggi_farmaco and 'livelli' in _ultimo and isinstance(_ultimo['livelli'], list):
        for l in _ultimo['livelli']:
            dt = t0_db + timedelta(hours=l.get('ore', 0))
            liv_list.append({"Data (GG/MM/AAAA)": dt.strftime("%d/%m/%Y"), "Ora (HH:MM)": dt.strftime("%H:%M"), "Conc. (mg/L)": l.get('conc', 0.0), "Nota (opz.)": l.get('nota', '')})
    liv_list.append({"Data (GG/MM/AAAA)": "", "Ora (HH:MM)": "", "Conc. (mg/L)": None, "Nota (opz.)": ""})

    colA, colB = st.columns(2)
    with colA: 
        st.markdown("**💉 Somministrazioni**")
        df_dosi_edit = st.data_editor(pd.DataFrame(dosi_list), num_rows="dynamic", use_container_width=True)
    with colB: 
        st.markdown("**🔬 Livelli Plasmatici**")
        df_liv_edit = st.data_editor(pd.DataFrame(liv_list), num_rows="dynamic", use_container_width=True)

    def _parse_data_ora(d_str, o_str):
        if not d_str or not o_str: return None
        try: d = datetime.strptime(str(d_str).strip(), "%d/%m/%Y").date()
        except:
            try: d = datetime.strptime(str(d_str).strip(), "%Y-%m-%d").date()
            except: return None
        try: t = datetime.strptime(str(o_str).strip()[:5], "%H:%M").time()
        except: return None
        return datetime.combine(d, t)

    def safe_float(val):
        try: return float(val)
        except (ValueError, TypeError): return None

    storia_completa = []
    for _, r in df_dosi_edit.iterrows():
        dt_obj = _parse_data_ora(r.get("Data (GG/MM/AAAA)"), r.get("Ora (HH:MM)"))
        dose_val = safe_float(r.get("Dose (mg)"))
        if dt_obj and dose_val is not None: storia_completa.append({'dt': dt_obj, 'dose': dose_val})
    
    livelli_completi = []
    for _, r in df_liv_edit.iterrows():
        dt_obj = _parse_data_ora(r.get("Data (GG/MM/AAAA)"), r.get("Ora (HH:MM)"))
        conc_val = safe_float(r.get("Conc. (mg/L)"))
        if dt_obj and conc_val is not None: livelli_completi.append({'dt': dt_obj, 'conc': conc_val, 'nota': str(r.get("Nota (opz.)", ""))})

    st.markdown("<br>", unsafe_allow_html=True)
    if st.button("🚀 AVVIA FITTING MAP BAYESIANO", type="primary", use_container_width=True):
        if not storia_completa: st.error("Nessuna dose valida.")
        elif len(livelli_completi) == 0: st.error("Inserisci almeno un livello plasmatico.")
        else:
            t0_assoluto = min(s['dt'] for s in storia_completa)
            for s in storia_completa: s['ore'] = (s['dt'] - t0_assoluto).total_seconds() / 3600.0
            storia_completa.sort(key=lambda x: x['ore'])
            
            for l in livelli_completi: l['ore'] = (l['dt'] - t0_assoluto).total_seconds() / 3600.0
            livelli_completi.sort(key=lambda x: x['ore'])

            dt_val = (max((l['dt'] for l in livelli_completi), default=None) or max((s['dt'] for s in storia_completa), default=datetime.now())).strftime('%Y-%m-%d %H:%M')

            with st.spinner("Calcolo MAP Bayesiano in corso…"):
                dose_mg = storia_completa[-1]['dose']; dose_mg_kg = dose_mg / peso
                tau_h_display = tau_mediano(storia_completa)  

                if farmaco == "Daptomicina":
                    prior = prior_dapto(id_paz, sesso, gg, peso, crea, crrt, Sc, Qef, dialisi)
                    CL_ind, V_ind, tipo_fit = map_1cmt(prior['CL_p'], prior['V_p'], prior['oCL'], prior['oV1'], prior['sig'], storia_completa, livelli_completi, prior['t_inf'])
                    Q_ind = V2_ind = 0.0
                    Vss_d = V_ind
                    Vss_std = prior['V_std']
                    
                    AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_1cmt, (CL_ind,V_ind), storia_completa, prior['t_inf'])
                    auc_mic = AUC24_v / mic; emivita = np.log(2)*V_ind/CL_ind
                    ok_staph, ok_ent, tox = auc_mic>=PD_DAPTO['AUC_MIC_staph'], auc_mic>=PD_DAPTO['AUC_MIC_entero'], Cmin_ss>=PD_DAPTO['Cmin_tox']
                
                elif farmaco == "Dalbavancina":
                    prior = prior_dalba(id_paz, peso, modello_dalba)
                    CL_ind, tipo_fit = map_dalba(prior['CL_p'], prior['V_p'], prior['Q_p'], prior['V2_p'], prior['oCL'], prior['sig'], storia_completa, livelli_completi, prior['t_inf'])
                    V_ind, Q_ind, V2_ind = prior['V_p'], prior['Q_p'], prior['V2_p']
                    Vss_d = V_ind + V2_ind
                    Vss_std = prior['V_std'] + prior['V2_std']
                    
                    AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_2cmt, (CL_ind,V_ind,Q_ind,V2_ind), storia_completa, prior['t_inf'])
                    auc_mic = (AUC24_v * prior['fu']) / mic  # fAUC24/MIC
                    ok_target, tox = auc_mic >= DALBA_TARGETS['ABSSSI semplice (fAUC24h/MIC >= 27.1)'], False  
                
                else: 
                    prior = prior_cipro(id_paz, sesso, gg, peso, crea, crrt, dialisi)
                    CL_ind, V_ind, Q_ind, V2_ind, tipo_fit = map_2cmt(prior['CL_p'], prior['V_p'], prior['Q_p'], prior['V2_p'], prior['oCL'], prior['oV1'], prior['oQ'], prior['oV2'], prior['sig'], storia_completa, livelli_completi, prior['t_inf'])
                    Vss_d = V_ind + V2_ind
                    Vss_std = prior['V_std'] + prior['V2_std']
                    
                    AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_2cmt, (CL_ind,V_ind,Q_ind,V2_ind), storia_completa, prior['t_inf'])
                    auc_mic = AUC24_v/mic; emivita = np.log(2)*Vss_d/CL_ind
                    ok_target, tox = auc_mic>=PD_CIPRO['AUC_MIC_target'], AUC24_v>=PD_CIPRO['AUC_alta']

                st.markdown("---")
                st.subheader("🔬 Risultati PK/PD")
                st.caption(f"**Modello Popolazione Base:** {prior['nm']}")
                if prior['note']: st.info(prior['note'])

                c1,c2,c3 = st.columns(3)
                c1.metric("Clearance (CL) Paziente", f"{CL_ind:.3f} L/h", delta=f"Normale pop: {prior['CL_std']:.3f} L/h", delta_color="off")
                c2.metric("Vol. Distribuzione (Vss)", f"{Vss_d:.1f} L", delta=f"Normale pop: {Vss_std:.1f} L", delta_color="off")
                c3.metric("Emivita (T½)", f"{emivita:.1f} h")
                
                st.markdown("<br>", unsafe_allow_html=True)
                c4,c5,c6 = st.columns(3)
                c4.metric("AUC₂₄ stimata", f"{AUC24_v:.0f} mg·h/L")

                if farmaco=="Daptomicina":
                    c5.metric("AUC₂₄/MIC (≥ 666)",f"{auc_mic:.0f}", delta="✅ TARGET" if ok_staph else "⚠️ SOTTO TARGET", delta_color="normal" if ok_staph else "inverse")
                    c6.metric(f"Cmin SS (CPK ≥ 24.3)",f"{Cmin_ss:.1f} mg/L", delta="⚠️ RISCHIO CPK" if tox else "✅ Sicuro", delta_color="inverse" if tox else "normal")
                elif farmaco=="Dalbavancina":
                    c5.metric("fAUC₂₄h/MIC (≥ 27.1)",f"{auc_mic:.0f}", delta="✅ TARGET" if ok_target else "⚠️ SOTTO TARGET", delta_color="normal" if ok_target else "inverse")
                    c6.metric("Cmax / Cmin previste",f"{Cmax_ss:.1f} / {Cmin_ss:.1f} mg/L")
                else:
                    c5.metric("AUC₂₄/MIC (≥ 125)",f"{auc_mic:.0f}", delta="✅ TARGET" if ok_target else "⚠️ SOTTO TARGET", delta_color="normal" if ok_target else "inverse")
                    c6.metric(f"AUC₂₄ (Soglia ≥ 400)",f"{AUC24_v:.0f} mg·h/L", delta="⚠️ AUC ELEVATA" if tox else "✅ Esposizione sicura", delta_color="inverse" if tox else "normal")

                st.markdown("### 📈 Profilo Farmacocinetico Interattivo")
                t_cambio_val = -1
                if monitoraggi_farmaco and 'storia_dosi' in _ultimo:
                    for l in livelli_completi:
                        if l['ore'] == _ultimo['livelli'][-1].get('ore', -1) if _ultimo.get('livelli') else -1:
                            t_cambio_val = l['ore']
                
                fig = disegna_grafico_plotly(farmaco, storia_completa, livelli_completi, CL_ind, V_ind, Q_ind, V2_ind, prior['CL_std'], prior['V_std'], prior['Q_std'], prior['V2_std'], prior['t_inf'], t_cambio_val)
                st.plotly_chart(fig, use_container_width=True)

                def E(ok): return "TARGET RAGGIUNTO" if ok else "INSUFFICIENTE"
                
                if farmaco=="Dalbavancina":
                    pd_block = f"  Efficacia (fAUC24h/MIC >= 27.1): {auc_mic:.0f}  ->  {E(ok_target)}\n  Cmax prevista: {Cmax_ss:.1f} mg/L\n  Cmin prevista: {Cmin_ss:.1f} mg/L\n"
                elif farmaco=="Daptomicina":
                    pd_block = f"  Efficacia vs S. aureus/MRSA (AUC24/MIC >= 666): {auc_mic:.0f}  ->  {E(ok_staph)}\n  Efficacia vs E. faecium (AUC24/MIC >= 143): {auc_mic:.0f}  ->  {E(ok_ent)}\n  Sicurezza muscolare (Cmin_ss >= 24.3): {Cmin_ss:.1f} mg/L  ->  {'ATTENZIONE: MONITORARE CPK' if tox else 'SICURO'}\n"
                else:
                    pd_block = f"  Efficacia gram-negativi (AUC24/MIC >= 125): {auc_mic:.0f}  ->  {E(ok_target)}\n  Esposizione sistemica: AUC24 = {AUC24_v:.0f} mg*h/L\n  Rischio neurotossicita' (AUC24 >= 400): {'MONITORARE: AUC ELEVATA' if tox else 'SICURO'}\n"

                referto = f"""=== TDM BUDDY v6.3 — REFERTO {farmaco.upper()} ===
Data: {dt_val} | Paziente: {id_paz} | Sesso: {sesso} | Eta': {eta:.1f} anni | Peso: {peso} kg

[MODELLO APPLICATO]
{prior['nm']}

[DATI CLINICI]
Creatinina sierica: {crea} mg/dL
Dialisi / CRRT: {dialisi}
Dose più recente: {dose_mg:.0f} mg ({dose_mg_kg:.1f} mg/kg)
MIC patogeno: {mic} mg/L

[PARAMETRI PK INDIVIDUALI (MAP Bayesiano)]
Clearance (CL): {CL_ind:.3f} L/h
Volume di distribuzione (Vss): {Vss_d:.2f} L
Emivita terminale (T1/2): {emivita:.1f} h

[VALUTAZIONE PK/PD]
AUC24 stimata: {AUC24_v:.0f} mg*h/L
{pd_block}"""

                st.markdown("### 📋 Referto per Cartella Clinica")
                st.info("Copia il testo sottostante cliccando sull'icona in alto a destra nel riquadro.")
                st.code(referto, language='markdown')

                if farmaco == "Dalbavancina":
                    with st.container(border=True):
                        st.markdown("#### 🔮 Forecast — Prossima Somministrazione")

                        colF1, colF2, colF3 = st.columns(3)
                        with colF1:
                            target_scelto = st.selectbox("Target PK/PD di riferimento:", list(DALBA_TARGETS.keys()) + ["Trough (100% fT > 4×MIC)"])
                        with colF2:
                            mic_scelta_label = st.selectbox("MIC / breakpoint:", list(DALBA_MIC_EUCAST.keys()) + ["Personalizzato"])
                            if mic_scelta_label == "Personalizzato": mic_forecast = st.number_input("MIC (mg/L)", 0.001, 4.0, 0.125, step=0.001, format="%.4f")
                            else: mic_forecast = DALBA_MIC_EUCAST[mic_scelta_label]
                        with colF3:
                            f_unbound = st.number_input("Frazione libera (f)", 0.01, 0.10, DALBA_PARAMS['fu'], 0.01)

                        if target_scelto == "Trough (100% fT > 4×MIC)":
                            soglia_target = dalba_soglia_trough(mic_forecast, f_unbound)
                            desc_target = f"trough 4×MIC libera → soglia totale = {soglia_target:.2f} mg/L"
                        else:
                            soglia_target = dalba_soglia_totale(DALBA_TARGETS.get(target_scelto, 27.1), mic_forecast, f_unbound)
                            desc_target = f"{target_scelto} → soglia totale equivalente = {soglia_target:.2f} mg/L"

                        st.caption(f"📐 Soglia calcolata: **{desc_target}**")

                        t_valutazione = max(max((l['ore'] for l in livelli_completi), default=0), storia_completa[-1]['ore'])
                        c_attuale = pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia_completa, t_valutazione, prior['t_inf'])
                        t_crossing = trova_tempo_soglia(CL_ind, V_ind, Q_ind, V2_ind, storia_completa, soglia_target, prior['t_inf'], t_partenza=t_valutazione, t_max=t_valutazione + 2000, passo=2.0)

                        if t_crossing is None: st.warning("⚠️ Soglia non raggiunta entro 83 giorni.")
                        elif t_crossing <= t_valutazione + 0.5: st.error(f"🔴 Concentrazione già sotto soglia ({c_attuale:.2f} mg/L).")
                        else: st.success(f"🗓️ La concentrazione scenderà sotto soglia il **{(t0_assoluto + timedelta(hours=t_crossing)).strftime('%d/%m/%Y %H:%M')}** ({(t_crossing - t_valutazione) / 24.0:.1f} giorni).")

                        st.markdown("##### Simula nuova somministrazione ipotetica")
                        col_sim1, col_sim2, col_sim3 = st.columns(3)
                        with col_sim1: simula = st.checkbox("Aggiungi dose al grafico Forecast")
                        with col_sim2:
                            if simula:
                                data_new = st.date_input("Data nuova dose", value=(t0_assoluto + timedelta(hours=t_valutazione+24)).date(), format="DD/MM/YYYY")
                                ora_new = st.time_input("Ora dose", value=datetime.strptime("10:00","%H:%M").time())
                        with col_sim3:
                            if simula: dose_new = st.number_input("Dose ipotetica (mg)", 100.0, 1500.0, 1500.0, 50.0)

                        if simula:
                            t_new_ore = (datetime.combine(data_new, ora_new) - t0_assoluto).total_seconds() / 3600.0
                            storia_sim = sorted(storia_completa + [{'ore': t_new_ore, 'dose': dose_new}], key=lambda x: x['ore'])
                            t_crossing_sim = trova_tempo_soglia(CL_ind, V_ind, Q_ind, V2_ind, storia_sim, soglia_target, prior['t_inf'], t_partenza=t_new_ore, t_max=t_new_ore+2000, passo=2.0)
                            
                            if t_crossing_sim: st.success(f"✅ Con dose di {dose_new:.0f} mg, concentrazione sopra soglia fino al **{(t0_assoluto + timedelta(hours=t_crossing_sim)).strftime('%d/%m/%Y')}**.")
                            fig2 = disegna_forecast_plotly(np.linspace(0, t_crossing_sim + 48 if t_crossing_sim else t_new_ore + 500, 1000), [pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia_sim, t, prior['t_inf']) for t in np.linspace(0, t_crossing_sim + 48 if t_crossing_sim else t_new_ore + 500, 1000)], soglia_target, t_new_ore, dose_new, storia_completa)
                            st.plotly_chart(fig2, use_container_width=True)

                # Salvataggio in sessione
                record = {
                    'data_val': dt_val, 'farmaco': farmaco, 'eta_giorni': gg, 'peso': peso,
                    'creatinina': crea, 'crrt': bool(crrt), 'Sc': Sc, 'Qef': Qef, 'dialisi_tipo': dialisi,
                    'dose_mg': dose_mg, 'dose_mg_kg': dose_mg_kg, 'tau_h': tau_h_display,
                    'n_dosi': len(storia_completa), 'mic': mic, 'n_livelli': len(livelli_completi),
                    'CL_ind': CL_ind, 'V_ind': Vss_d, 'AUC24': AUC24_v, 'AUC_MIC': auc_mic, 'Cmin_ss': Cmin_ss, 
                    'modello': prior['nm'], 'modello_dalba': modello_dalba, 'storia_dosi': storia_completa,
                    'livelli': livelli_completi, 'data_prima_dose': t0_assoluto.strftime('%Y-%m-%d %H:%M'),
                }
                st.session_state['storico']['monitoraggi'].append(record)
                st.session_state['storico']['codice_paziente'] = id_paz
                st.session_state['storico']['sesso'] = sesso
                st.success("💾 Analisi validata e aggiunta allo storico di sessione. Ricordati di scaricare il JSON nella barra laterale sinistra per conservarlo.")

with tab2:
    monitoraggi = st.session_state['storico'].get('monitoraggi', [])
    if not monitoraggi: st.info("Nessuna analisi TDM in questa sessione.")
    else:
        df_sto = pd.DataFrame(monitoraggi)
        cols_vis = [c for c in ['data_val','farmaco','creatinina','mic','dose_mg','n_dosi','n_livelli','CL_ind','AUC24'] if c in df_sto.columns]
        st.dataframe(df_sto[cols_vis], use_container_width=True)
        
        colS1, colS2 = st.columns(2)
        with colS1:
            scelta = st.selectbox("Seleziona per rimuovere un record errato:", ["— Scegli —"] + [f"#{i} — {r['data_val']} [{r.get('farmaco','?')}]" for i, r in enumerate(monitoraggi)])
        with colS2:
            st.write("")
            st.write("")
            if scelta != "— Scegli —":
                idx_sel = int(scelta.split("#")[1].split(" ")[0])
                if st.button("🗑️ Rimuovi questo record"):
                    st.session_state['storico']['monitoraggi'].pop(idx_sel)
                    st.rerun()
