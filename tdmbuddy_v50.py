"""
TDM Buddy v6.1
==============

FUNZIONALITÀ PRINCIPALI
  - MAP Bayesiano 1-compartimento (Dapto) e 2-compartimento (Cipro/Dalba)
  - Interfaccia riprogettata: simmetria perfetta delle card iniziali.
  - Fix caricamento JSON: i widget si ripopolano automaticamente senza bloccare la UI.
  - Archivio clinico SQLite / Gestione JSON State (scoped a molecola e paziente).
  - Nessun salvataggio cloud: massima privacy, solo export JSON locale.
"""

import streamlit as st
import numpy as np
from scipy.optimize import minimize
import plotly.graph_objects as go
import pandas as pd
from datetime import datetime, date, timedelta
import json

# ============================================================
# COSTANTI E RIFERIMENTI
# ============================================================

FARMACI = ["Daptomicina", "Ciprofloxacina", "Dalbavancina"]

FARMACO_INFO = {
    "Daptomicina": {"icona": "🩸", "descrizione": "Lipopeptide"},
    "Ciprofloxacina": {"icona": "🧬", "descrizione": "Fluorochinolone"},
    "Dalbavancina": {"icona": "💉", "descrizione": "Lipoglicopeptide"},
}

RIFERIMENTI_PER_FASCIA = {
    "Daptomicina": [
        ("PK — < 1 anno", "Ye L, Zhou H, et al. Eur J Pharm Sci. 2025;208:107058"),
        ("PK — 1–17 anni", "Wei XC, et al. J Clin Pharmacol. 2020;60(6):768-774 (correzione IRC: Ye L, et al. Front Pharmacol. 2022)"),
        ("PK — Adolescenti in HD/CCPD", "Lim SY, et al. J Pediatr Pharmacol Ther. 2021;26(2):123-132"),
        ("PK — ≥ 18 anni", "Soraluce A, et al. Int J Antimicrob Agents. 2018;52(2):158-165"),
        ("Target S. aureus/MRSA (AUC/MIC ≥ 666)", "Louie A, et al. Antimicrob Agents Chemother. 2001; Canut A, et al. Eur J Clin Microbiol Infect Dis. 2012"),
        ("Target E. faecium (AUC/MIC ≥ 143)", "Safdar N, et al. Antimicrob Agents Chemother. 2004; Wei XC, et al. 2020"),
        ("Target tossicità (Cmin ≥ 24.3 mg/L)", "Bhavnani SM, et al. Clin Infect Dis. 2010;50(12):1568-1574"),
    ],
    "Ciprofloxacina": [
        ("PK — Lifespan", "Schlender JF, et al. Clin Pharmacokinet. 2018; Balbas-Martinez V, et al. Eur J Pharm Sci. 2019"),
        ("IIV", "Forrest A, et al. Antimicrob Agents Chemother. 1993;37(5):1065-1072"),
        ("Target gram-negativi (AUC/MIC ≥ 125)", "Forrest A, et al. Antimicrob Agents Chemother. 1993;37(5):1073-1081"),
    ],
    "Dalbavancina": [
        ("PK — Pediatrico (Mariani)", "Dati interni IRCCS Gaslini, non pubblicati"),
        ("PK — Adulti (Cojutti)", "Cojutti PG, et al. Antimicrob Agents Chemother. 2021;65(5):e02260-20"),
        ("Target efficacia", "Cojutti PG, et al. Antimicrob Agents Chemother. 2021;65(5):e02260-20"),
        ("Breakpoint MIC", "EUCAST Breakpoint Tables v.16.0 (2026) / USCAST"),
    ],
}

# ============================================================
# STORAGE
# ============================================================

def storico_vuoto(farmaco=None):
    return {'codice_paziente': '', 'sesso': 'Maschio', 'farmaco': farmaco, 'monitoraggi': []}

if 'storico' not in st.session_state:
    st.session_state['storico'] = storico_vuoto()
if 'farmaco_attivo' not in st.session_state:
    st.session_state['farmaco_attivo'] = None

def pulisci_widget_paziente():
    chiavi = ['input_codice_paziente', 'input_sesso', 'editor_dosi', 'editor_livelli', 'input_unita_eta', 'input_valore_eta', 'input_peso', 'input_crea', 'input_mic', 'input_crrt', 'input_sc', 'input_qef', 'input_dialisi']
    for k in chiavi: st.session_state.pop(k, None)

# ============================================================
# PARAMETRI PK
# ============================================================

SORALUCE = {'CLNR': 0.160, 'CLR_theta': 0.366, 'Clcr_ref': 49.0, 'V': 12.30, 'omega_CL': 0.367, 'omega_V': 0.278, 'sigma': 0.123, 'Sc_default': 0.20, 'nome': 'Soraluce 2018'}
YE2025 = {'CL': 19.6, 'SD': 5.5, 'V': 0.17, 'nome': 'Ye L 2025'}
WEI = {'3-6m': {'CL':19.7,'SD':5.5, 'V':0.17,'nome':'Wei 2020'}, '7-12m': {'CL':19.6,'SD':5.8, 'V':0.17,'nome':'Wei 2020'}, '13-24m': {'CL':21.8,'SD':3.0, 'V':0.11,'nome':'Wei 2020'}, '2-6y': {'CL':20.0,'SD':6.1, 'V':0.11,'nome':'Wei 2020'}, '7-11y': {'CL':17.0,'SD':5.7, 'V':0.11,'nome':'Wei 2020'}, '12-17y': {'CL':11.0,'SD':2.1, 'V':0.11,'nome':'Wei 2020'}}
YE2022 = {'12-17y': {'mild':0.98,'moderate':1.08,'severe':1.24,'esrd':1.33}, '7-11y': {'mild':1.00,'moderate':1.19,'severe':1.49,'esrd':1.73}, '2-6y': {'mild':1.02,'moderate':1.28,'severe':1.74,'esrd':2.20}, '1-2y': {'mild':1.04,'moderate':1.29,'severe':1.72,'esrd':2.13}}
LIM = {'Vss':5.90,'CLren':0.023,'CLnr':0.158,'CL_HD':0.619, 'omega_V1':0.112,'omega_knr':0.088,'sigma':0.085,'peso_ref':40.0, 'nome': 'Lim SY 2021'}
PD_DAPTO = {'AUC_MIC_staph': 666, 'AUC_MIC_entero': 143, 'Cmin_tox': 24.3}

CIPRO = {'CL_kg':0.50, 'V1_kg':0.55, 'Q_kg':0.30, 'V2_kg':1.90, 'omega_CL':0.37, 'omega_V1':0.32, 'omega_Q':0.25, 'omega_V2':0.30, 'sigma':0.15, 'fu':0.67, 'nome': 'Schlender 2018 / Balbas-Martinez 2019'}
CIPRO_PED = {(91,365):1.40, (365,730):1.30, (730,2190):1.20, (2190,4380):1.10, (4380,6570):1.00}
CIPRO_FRAZ_RENAL = 0.70
PD_CIPRO = {'AUC_MIC_target': 125}
IRC_MAP = {"Lieve (GFR 60-90)":'mild', "Moderata (GFR 30-60)":'moderate', "Grave (GFR 15-30)":'severe', "ESRD (GFR <15)":'esrd'}

DALBA_PARAMS = {'CL_typ70': 0.0305, 'V1_typ70': 4.81, 'Q_typ70': 0.0081, 'V2_typ70': 4.81, 'omega_CL': 0.198, 'sigma': 0.375, 't_inf': 0.5, 'fu': 0.07, 'nome': "Modello Mariani (Gaslini)"}
COJUTTI_PARAMS = {'CL_typ70': 0.106, 'V1_typ70': 17.4, 'Q_typ70': 0.103, 'V2_typ70': 15.1, 'omega_CL': 0.362, 'sigma': 0.18, 't_inf': 0.5, 'fu': 0.07, 'nome': "Modello Cojutti 2021 (Adulti)"}
DALBA_MODELLI = {'Mariani — pediatrico, allometrico (Gaslini)': DALBA_PARAMS, 'Cojutti 2021 — adulti, osteoarticolare': COJUTTI_PARAMS}
DALBA_MIC_OPZIONI = {'MIC90 osteoarticolare (0.0625 mg/L)': 0.0625, 'EUCAST S. aureus (≤0.125 mg/L)': 0.125, 'USCAST S. aureus (≤0.25 mg/L)': 0.25}
DALBA_TARGETS = {'Target minimo (≥ 27.1)': 27.1, 'Target intermedio (≥ 53.3)': 53.3, 'Target restrittivo (≥ 111.1)': 111.1}

def dalba_soglia_totale(target_auc_mic, mic, f): return (target_auc_mic * mic) / (24.0 * f)
def dalba_soglia_trough_c0mic(mic): return 200.0 * mic

# ============================================================
# MOTORI PK
# ============================================================

def pk_1cmt(CL, V, storia, t_obs, t_inf=0.5):
    if CL <= 0 or V <= 0: return 0.0
    ke = CL / V; C = 0.0
    for s in storia:
        dt = t_obs - s['ore']
        if dt <= 0: continue
        D = s['dose']
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
        dt=t_obs-s['ore']
        if dt <= 0: continue
        D=s['dose']; k0=D/t_inf
        if dt<=t_inf: C += k0 * ((A/alp)*(1-np.exp(-alp*dt)) + (B/bet)*(1-np.exp(-bet*dt)))
        else: C += k0 * ((A/alp)*(1-np.exp(-alp*t_inf))*np.exp(-alp*(dt-t_inf)) + (B/bet)*(1-np.exp(-bet*t_inf))*np.exp(-bet*(dt-t_inf)))
    return max(C,0.0)

def metriche_esposizione_24h(pk_func, params, storia, t_inf, n_punti=300):
    t_ultima_dose = storia[-1]['ore']
    t_grid = np.linspace(t_ultima_dose, t_ultima_dose + 24, n_punti)
    c_grid = np.array([pk_func(*params, storia, t, t_inf) for t in t_grid])
    auc = float(np.sum((c_grid[1:] + c_grid[:-1]) / 2.0 * np.diff(t_grid)))
    return auc, float(c_grid.max()), float(c_grid.min())

def tau_mediano(storia):
    if len(storia) < 2: return None
    ore = sorted(s['ore'] for s in storia)
    diffs = sorted(ore[i+1]-ore[i] for i in range(len(ore)-1))
    return diffs[len(diffs)//2] if len(diffs) % 2 == 1 else (diffs[len(diffs)//2-1]+diffs[len(diffs)//2])/2.0

def calcola_clcr(eta_a, peso, crea, sesso):
    if crea <= 0: return 60.0
    clcr = ((140-eta_a)*peso)/(72*crea)
    return max(clcr*0.85 if sesso == "Femmina" else clcr, 1.0)

# ============================================================
# PRIOR DI POPOLAZIONE
# ============================================================

def prior_dapto(sesso, gg, peso, crea, crrt, Sc, Qef, dialisi):
    eta = gg/365.25; note = ""
    if gg < 365:
        p = YE2025; CL_pop = (p['CL']*peso)/1000; V_pop = p['V']*peso
        oCL = p['SD']/p['CL']; oV = 0.30; sig = 0.12; t_inf = 1.0; nm = p['nome']
        if crrt: CL_pop += Sc*Qef; note += f"CRRT: CLEC={Sc*Qef:.3f} L/h. "
    elif 4380 <= gg < 6570 and dialisi in ("Emodialisi (HD)","Dialisi Peritoneale (CCPD)"):
        sc = peso/LIM['peso_ref']
        CL_pop = (LIM['CLren']+LIM['CLnr'])*sc; V_pop = LIM['Vss']*sc; nm = LIM['nome']
        oCL, oV, sig, t_inf = LIM['omega_knr'], LIM['omega_V1'], LIM['sigma'], 0.5
    elif 91 <= gg < 6570:
        ag = '12-17y'
        if gg < 213: ag='3-6m'
        elif gg < 365: ag='7-12m'
        elif gg < 730: ag='13-24m'
        elif gg < 2190: ag='2-6y'
        elif gg < 4380: ag='7-11y'
        p = WEI[ag]; CL_pop=(p['CL']*peso)/1000; V_pop=p['V']*peso
        oCL=p['SD']/p['CL']; oV=0.25; sig=0.12; t_inf = 1.0 if gg<730 else 0.5; nm = p['nome']
        age_ye = ag if ag in YE2022 else ('2-6y' if ag in ('13-24m','2-6y') else '1-2y')
        if dialisi in IRC_MAP and age_ye in YE2022:
            CL_pop /= YE2022[age_ye][IRC_MAP[dialisi]]
        if crrt: CL_pop += Sc*Qef
    else:
        s = SORALUCE; clcr = calcola_clcr(eta, peso, crea, sesso)
        CL_pop = s['CLNR'] + s['CLR_theta']*(clcr/s['Clcr_ref']); V_pop = s['V']
        nm = f"{s['nome']} (ClCr={clcr:.0f} mL/min)"
        oCL, oV, sig, t_inf = s['omega_CL'], s['omega_V'], s['sigma'], 0.5
        if crrt: CL_pop += Sc*Qef

    return {'CL_p':CL_pop, 'V_p':V_pop, 'Q_p':0.0, 'V2_p':0.0, 'oCL':oCL, 'oV1':oV, 'oQ':0.0, 'oV2':0.0, 'sig':sig, 't_inf':t_inf, 'nm':nm, 'note':note}

def prior_cipro(sesso, gg, peso, crea, crrt, dialisi):
    eta = gg/365.25; note=""; c=CIPRO; fc=1.0
    for (lo,hi),f in CIPRO_PED.items():
        if lo <= gg < hi: fc = f; break

    clcr = calcola_clcr(eta, peso, crea, sesso)
    if clcr >= 30: CL_pop = (c['CL_kg']*peso*fc*CIPRO_FRAZ_RENAL*min(clcr/120.0,1.0)) + (c['CL_kg']*peso*fc*(1-CIPRO_FRAZ_RENAL))
    else: CL_pop = (c['CL_kg']*peso*fc*CIPRO_FRAZ_RENAL*0.20) + (c['CL_kg']*peso*fc*(1-CIPRO_FRAZ_RENAL)*0.63)

    return {'CL_p':CL_pop, 'V_p':c['V1_kg']*peso, 'Q_p':c['Q_kg']*peso, 'V2_p':c['V2_kg']*peso,
            'oCL':c['omega_CL'], 'oV1':c['omega_V1'], 'oQ':c['omega_Q'], 'oV2':c['omega_V2'],
            'sig':c['sigma'], 't_inf':1.0, 'nm':f"{c['nome']} (ClCr={clcr:.0f})", 'note':note}

def prior_dalba(peso, modello_scelto):
    p = DALBA_MODELLI[modello_scelto]
    if modello_scelto.startswith('Mariani'):
        fw = peso/70.0
        return {'CL_p':p['CL_typ70']*(fw**0.75), 'V_p':p['V1_typ70']*fw, 'Q_p':p['Q_typ70']*(fw**0.75), 'V2_p':p['V2_typ70']*fw, 'oCL':p['omega_CL'], 'oV1':0.0, 'oQ':0.0, 'oV2':0.0, 'sig':p['sigma'], 't_inf':p['t_inf'], 'nm':p['nome'], 'note':"Scaling allometrico applicato."}
    return {'CL_p':p['CL_typ70'], 'V_p':p['V1_typ70'], 'Q_p':p['Q_typ70'], 'V2_p':p['V2_typ70'], 'oCL':p['omega_CL'], 'oV1':0.0, 'oQ':0.0, 'oV2':0.0, 'sig':p['sigma'], 't_inf':p['t_inf'], 'nm':p['nome'], 'note':"Nessuno scaling applicato."}

# ============================================================
# MAP BAYESIANO
# ============================================================

def map_1cmt(CL_p, V_p, oCL, oV, sig, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0: return CL_p, V_p, "Prior di popolazione"
    def f(p):
        CL_, V_ = np.exp(p[0]), np.exp(p[1])
        pr = (((p[0]-np.log(CL_p))**2)/oCL**2 + ((p[1]-np.log(V_p))**2)/oV**2)
        lk = sum(((np.log(max(l['conc'],1e-6)) - np.log(max(pk_1cmt(CL_,V_,storia,l['ore'],t_inf),1e-6)))**2)/sig**2 for l in livelli)
        return 0.5*(pr+lk)
    r = minimize(f, [np.log(CL_p),np.log(V_p)], method='L-BFGS-B', bounds=[(np.log(CL_p*0.05), np.log(CL_p*8)), (np.log(V_p*0.15), np.log(V_p*5))], options={'ftol':1e-12})
    return float(np.exp(r.x[0])), float(np.exp(r.x[1])), f"MAP Bayesiano 1-cmt su {n} livelli"

def map_2cmt(CL_p, V1_p, Q_p, V2_p, oCL, oV1, oQ, oV2, sig, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0: return CL_p, V1_p, Q_p, V2_p, "Prior di popolazione"
    OQ_eff  = 0.10 + (oQ  - 0.10)*min(1.0, n/6.0)
    OV2_eff = 0.10 + (oV2 - 0.10)*min(1.0, n/6.0)
    def f(p):
        CL_,V1_,Q_,V2_ = np.exp(p[0]),np.exp(p[1]),np.exp(p[2]),np.exp(p[3])
        pr = (((p[0]-np.log(CL_p))**2)/oCL**2 + ((p[1]-np.log(V1_p))**2)/oV1**2 + ((p[2]-np.log(Q_p))**2)/OQ_eff**2 + ((p[3]-np.log(V2_p))**2)/OV2_eff**2)
        lk = sum(((np.log(max(l['conc'],1e-6)) - np.log(max(pk_2cmt(CL_,V1_,Q_,V2_,storia,l['ore'],t_inf),1e-6)))**2)/sig**2 for l in livelli)
        return 0.5*(pr+lk)
    r = minimize(f, [np.log(CL_p),np.log(V1_p),np.log(Q_p),np.log(V2_p)], method='L-BFGS-B', bounds=[(np.log(CL_p*0.05),np.log(CL_p*10)), (np.log(V1_p*0.1),np.log(V1_p*6)), (np.log(Q_p*0.3),np.log(Q_p*3)), (np.log(V2_p*0.3),np.log(V2_p*3))], options={'ftol':1e-12})
    return float(np.exp(r.x[0])),float(np.exp(r.x[1])),float(np.exp(r.x[2])),float(np.exp(r.x[3])), f"MAP Bayesiano 2-cmt su {n} livelli"

def map_dalba(CL_p, V1_p, Q_p, V2_p, oCL, sig, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0: return CL_p, "Prior di popolazione"
    def f(p):
        CL_ = np.exp(p[0])
        prior = ((p[0]-np.log(CL_p))**2)/oCL**2
        lik = sum(((np.log(max(l['conc'],1e-6)) - np.log(max(pk_2cmt(CL_,V1_p,Q_p,V2_p,storia,l['ore'],t_inf),1e-6)))**2)/sig**2 for l in livelli)
        return 0.5*(prior+lik)
    r = minimize(f, [np.log(CL_p)], method='L-BFGS-B', bounds=[(np.log(CL_p*0.05),np.log(CL_p*10.0))], options={'ftol':1e-12})
    return float(np.exp(r.x[0])), f"MAP Bayesiano su CL, {n} livelli (V1/Q/V2 fissi)"

def trova_tempo_soglia(CL, V1, Q, V2, storia, soglia, t_inf, t_partenza=0.0, t_max=2000.0, passo=1.0):
    if pk_2cmt(CL,V1,Q,V2,storia,t_partenza,t_inf) < soglia: return t_partenza
    t_prev = t_partenza; t = t_partenza + passo
    while t <= t_max:
        if pk_2cmt(CL,V1,Q,V2,storia,t,t_inf) < soglia:
            lo, hi = t_prev, t
            for _ in range(60):
                mid = (lo+hi)/2.0
                if pk_2cmt(CL,V1,Q,V2,storia,mid,t_inf) >= soglia: lo = mid
                else: hi = mid
            return hi
        t_prev = t; t += passo
    return None

# ============================================================
# GRAFICI PLOTLY 
# ============================================================

def disegna_grafico(farmaco, storia, livelli, CL_ind, V_ind, Q_ind, V2_ind, CL_std, V_std, Q_std, V2_std, t_inf):
    t_ultimo_livello = max((l['ore'] for l in livelli), default=0)
    t_fine = max(t_ultimo_livello, max((s['ore'] for s in storia), default=0)+12) + 12
    t_arr = np.linspace(0, t_fine, 800)

    if farmaco == "Daptomicina":
        C_ind = [pk_1cmt(CL_ind,V_ind,storia,t,t_inf) for t in t_arr]
        C_std = [pk_1cmt(CL_std,V_std,storia,t,t_inf) for t in t_arr]
    else:
        C_ind = [pk_2cmt(CL_ind,V_ind,Q_ind,V2_ind,storia,t,t_inf) for t in t_arr]
        C_std = [pk_2cmt(CL_std,V_std,Q_std,V2_std,storia,t,t_inf) for t in t_arr]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=t_arr, y=C_std, mode='lines', line=dict(dash='dash', color='#94A3B8', width=2), name='Popolazione standard'))
    fig.add_trace(go.Scatter(x=t_arr, y=C_ind, mode='lines', line=dict(color='#1E3A5F', width=3), name='Paziente (MAP)'))

    max_c = max(max(C_ind, default=100), max(C_std, default=100))
    x_dosi, y_dosi = [], []
    for s in storia: x_dosi += [s['ore'], s['ore'], None]; y_dosi += [0, max_c, None]
    if x_dosi: fig.add_trace(go.Scatter(x=x_dosi, y=y_dosi, mode='lines', line=dict(color='#475569', dash='dot', width=1), name='Somministrazioni', hoverinfo='none'))

    if livelli:
        x_liv = [l['ore'] for l in livelli]; y_liv = [l['conc'] for l in livelli]
        testi = [f"t={l['ore']:.1f}h<br>C={l['conc']:.3f} mg/L" + (f"<br>{l['nota']}" if l.get('nota') else "") for l in livelli]
        fig.add_trace(go.Scatter(x=x_liv, y=y_liv, mode='markers', marker=dict(symbol='star', size=14, color='#16A34A', line=dict(width=1, color='white')), name='Livelli misurati', text=testi, hoverinfo='text'))

    if farmaco == "Daptomicina": fig.add_hline(y=24.3, line_dash="dashdot", line_color="#DC2626", annotation_text="Soglia CPK")

    fig.update_layout(xaxis_title="Ore dall'inizio della terapia (t=0)", yaxis_title=f"{farmaco} (mg/L)", plot_bgcolor='#FFFFFF', paper_bgcolor='#FFFFFF', hovermode="x unified", legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1), margin=dict(l=40, r=40, t=40, b=40))
    return fig

def disegna_forecast(t_arr, c_arr, soglia, t_new_ore, dose_new, storia):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=t_arr, y=c_arr, mode='lines', line=dict(color='#1E3A5F', width=2.5), name='Proiezione'))
    fig.add_hline(y=soglia, line_dash="dashdot", line_color="#DC2626", annotation_text=f"Soglia ({soglia:.2f} mg/L)")
    fig.add_vline(x=t_new_ore, line_dash="dash", line_color="#16A34A", annotation_text=f"Nuova dose ({dose_new:.0f} mg)")

    max_c = max(c_arr, default=100)
    x_dosi, y_dosi = [], []
    for s in storia: x_dosi += [s['ore'], s['ore'], None]; y_dosi += [0, max_c, None]
    if x_dosi: fig.add_trace(go.Scatter(x=x_dosi, y=y_dosi, mode='lines', line=dict(color='#94A3B8', dash='dot', width=1), name='Dosi pregresse', hoverinfo='none'))

    fig.update_layout(xaxis_title="Ore dall'inizio", yaxis_title="Dalbavancina (mg/L)", plot_bgcolor='#FFFFFF', paper_bgcolor='#FFFFFF', hovermode="x unified", legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1), margin=dict(l=40, r=40, t=40, b=40))
    return fig

# ============================================================
# INTERFACCIA
# ============================================================

st.set_page_config(page_title="TDM Buddy v6.1", layout="wide", page_icon="💊")

st.title("💊 TDM Buddy v6.1")
st.caption("Monitoraggio terapeutico bayesiano — Daptomicina, Ciprofloxacina, Dalbavancina")

# ------------------------------------------------------------
# STEP 1 — Selezione molecola (Layout Semplificato e Simmetrico)
# ------------------------------------------------------------

if st.session_state['farmaco_attivo'] is None:
    st.markdown("## Passo 1 — Seleziona la molecola")
    st.caption("La scelta è vincolante per questa sessione. Per cambiare molecola dovrai ricominciare.")

    cols = st.columns(3)
    for col, farmaco in zip(cols, FARMACI):
        info = FARMACO_INFO[farmaco]
        with col:
            with st.container(border=True):
                st.markdown(f"<h3 style='text-align: center; margin-bottom: 0;'>{info['icona']} {farmaco}</h3>", unsafe_allow_html=True)
                st.markdown(f"<p style='text-align: center; color: gray; margin-top: 0;'>{info['descrizione']}</p>", unsafe_allow_html=True)
                if st.button("Seleziona", key=f"sel_{farmaco}", use_container_width=True, type="primary"):
                    st.session_state['farmaco_attivo'] = farmaco
                    st.session_state['storico'] = storico_vuoto(farmaco)
                    st.rerun()
    st.stop()

farmaco = st.session_state['farmaco_attivo']

# Badge e Cambio
col_badge, col_cambia = st.columns([5, 1])
with col_badge: st.info(f"{FARMACO_INFO[farmaco]['icona']} **Molecola in corso: {farmaco}**")
with col_cambia:
    if st.button("🔄 Cambia molecola", use_container_width=True):
        pulisci_widget_paziente()
        st.session_state['farmaco_attivo'] = None
        st.session_state['storico'] = storico_vuoto()
        st.session_state.pop('risultati', None)
        st.rerun()

# ------------------------------------------------------------
# SIDEBAR — Anagrafica e JSON
# ------------------------------------------------------------

st.sidebar.header("🗂️ Storico paziente (locale)")

if st.sidebar.button("🔄 Svuota Dati / Nuovo Paziente", type="primary", use_container_width=True):
    pulisci_widget_paziente()
    st.session_state['storico'] = storico_vuoto(farmaco)
    st.session_state.pop('risultati', None)
    st.rerun()

st.sidebar.markdown("---")
file_up = st.sidebar.file_uploader("📤 Carica storico (JSON)", type=["json"])
if file_up is not None:
    file_hash = hash(file_up.getvalue())
    if st.session_state.get('uploaded_file_hash') != file_hash:
        try:
            dati_caricati = json.loads(file_up.getvalue())
            farmaco_file = dati_caricati.get('farmaco')
            if farmaco_file != farmaco:
                st.sidebar.error(f"❌ Errore: file di {farmaco_file or 'altra molecola'}.")
            else:
                pulisci_widget_paziente()
                st.session_state['storico'] = dati_caricati
                st.session_state['uploaded_file_hash'] = file_hash
                st.session_state['input_codice_paziente'] = dati_caricati.get('codice_paziente', '')
                st.session_state['input_sesso'] = dati_caricati.get('sesso', 'Maschio')
                st.rerun()
        except Exception:
            st.sidebar.error("File non valido o corrotto.")

id_paz = st.sidebar.text_input("Codice paziente (no nomi veri)", key="input_codice_paziente")
sesso = st.sidebar.selectbox("Sesso biologico", ["Maschio", "Femmina"], key="input_sesso")
st.session_state['storico']['codice_paziente'] = id_paz
st.session_state['storico']['sesso'] = sesso

n_rec = len(st.session_state['storico'].get('monitoraggi', []))
if n_rec > 0:
    st.sidebar.markdown("---")
    st.sidebar.info(f"📋 {n_rec} analisi in questa sessione.")
    st.sidebar.download_button("📥 Scarica storico aggiornato", data=json.dumps(st.session_state['storico'], indent=2, ensure_ascii=False).encode('utf-8'), file_name=f"TDMBuddy_{farmaco}_{id_paz or 'paz'}.json", mime="application/json", type="primary", use_container_width=True)

if not id_paz:
    st.warning("👈 Inserisci un codice paziente nella barra laterale per procedere.")
    st.stop()

# ------------------------------------------------------------
# STEP 2 — Dati biometrici e clinici
# ------------------------------------------------------------

monitoraggi_farmaco = st.session_state['storico'].get('monitoraggi', [])
_ultimo = monitoraggi_farmaco[-1] if monitoraggi_farmaco else None

_gg_def   = int(_ultimo.get('eta_giorni')) if _ultimo else None
_peso_def = float(_ultimo.get('peso', 70.0))
_crea_def = float(_ultimo.get('creatinina', 1.0))
_mic_def  = float(_ultimo.get('mic', 0.0625 if farmaco=="Dalbavancina" else 1.0))

st.markdown("## Passo 2 — Dati biometrici e clinici")

with st.expander(f"📚 Modelli di riferimento applicati — {farmaco}", expanded=False):
    for fascia, rif in RIFERIMENTI_PER_FASCIA[farmaco]: st.caption(f"**{fascia}:** {rif}")

with st.container(border=True):
    col_u, col_v, col_peso = st.columns(3)
    unita_eta = col_u.selectbox("Unità età", ["Anni", "Mesi", "Giorni"], key="input_unita_eta")
    valore_eta_def = round(_gg_def / 365.25, 2) if _gg_def is not None else None
    valore_eta = col_v.number_input(f"Età ({unita_eta.lower()})", min_value=0.0, value=valore_eta_def, key="input_valore_eta")
    peso = col_peso.number_input("Peso (kg)", 0.3, 300.0, value=_peso_def, key="input_peso")

    if valore_eta is None: gg = None
    elif unita_eta == "Anni": gg = int(round(valore_eta * 365.25))
    elif unita_eta == "Mesi": gg = int(round(valore_eta * 30.44))
    else: gg = int(round(valore_eta))
    eta = gg / 365.25 if gg is not None else None

    modello_dalba = None
    if farmaco == "Dalbavancina":
        idx_def = list(DALBA_MODELLI.keys()).index(_ultimo['modello_dalba']) if _ultimo and _ultimo.get('modello_dalba') in DALBA_MODELLI else 0
        modello_dalba = st.selectbox("🧬 Modello PK di riferimento:", list(DALBA_MODELLI.keys()), index=idx_def, key="input_modello_dalba")

with st.container(border=True):
    col1, col2, col3 = st.columns(3)
    crea = col1.number_input("Creatinina sierica (mg/dL)", 0.05, 25.0, value=_crea_def, step=0.05, key="input_crea")
    mic = col1.number_input("MIC di riferimento (mg/L)", 0.001, 16.0, value=_mic_def, step=0.001, format="%.4f", key="input_mic")
    
    crrt = col2.checkbox("Paziente in CRRT?", key="input_crrt")
    if crrt:
        Sc  = col2.number_input("Sc", 0.01, 0.60, 0.20, 0.01, key="input_sc")
        Qef = col2.number_input("Qef (L/h)", 0.5, 6.0, 2.0, 0.1, key="input_qef")
    else: Sc, Qef = 0.20, 2.0

    if farmaco == "Daptomicina" and gg is not None and gg < 6570:
        opt = ["Nessuna", "Lieve (GFR 60-90)", "Moderata (GFR 30-60)", "Grave (GFR 15-30)", "ESRD (GFR <15)"]
        if 4380 <= gg < 6570: opt += ["Emodialisi (HD)", "Dialisi Peritoneale (CCPD)"]
        dialisi = col3.selectbox("Funzione renale / dialisi:", opt, key="input_dialisi")
    else: dialisi = "Nessuna"

campi_mancanti = [k for k, v in [("età", gg), ("peso", peso), ("MIC", mic)] + ([("creatinina", crea)] if farmaco != "Dalbavancina" else []) if v is None]
if campi_mancanti: st.warning(f"⚠️ Compila: {', '.join(campi_mancanti)}.")

# ------------------------------------------------------------
# STEP 3 — Timeline dosi e livelli
# ------------------------------------------------------------

st.markdown("## Passo 3 — Timeline Somministrazioni e Prelievi")
st.info("💡 Aggiungi o modifica le righe. Inserisci tutte le dosi e i livelli del paziente.")

t0_db = datetime.now()
dosi_list = []
if _ultimo and 'storia_dosi' in _ultimo and _ultimo['storia_dosi']:
    t0_db = datetime.strptime(_ultimo['data_prima_dose'], '%Y-%m-%d %H:%M')
    for s in _ultimo['storia_dosi']: dosi_list.append({"Data": (t0_db + timedelta(hours=s['ore'])).strftime("%d/%m/%Y"), "Ora": (t0_db + timedelta(hours=s['ore'])).strftime("%H:%M"), "Dose (mg)": s['dose']})
else: dosi_list.append({"Data": "", "Ora": "", "Dose (mg)": None})

liv_list = []
if _ultimo and 'livelli' in _ultimo:
    for l in _ultimo['livelli']: liv_list.append({"Data": (t0_db + timedelta(hours=l['ore'])).strftime("%d/%m/%Y"), "Ora": (t0_db + timedelta(hours=l['ore'])).strftime("%H:%M"), "Conc. (mg/L)": l['conc'], "Nota": l.get('nota', '')})
liv_list.append({"Data": "", "Ora": "", "Conc. (mg/L)": None, "Nota": ""})

colA, colB = st.columns(2)
with colA:
    st.markdown("**💉 Somministrazioni**")
    df_dosi_edit = st.data_editor(pd.DataFrame(dosi_list), num_rows="dynamic", use_container_width=True, height=280, key="editor_dosi")
with colB:
    st.markdown("**🔬 Livelli plasmatici**")
    df_liv_edit = st.data_editor(pd.DataFrame(liv_list), num_rows="dynamic", use_container_width=True, height=280, key="editor_livelli")

def _parse_dt(d_str, o_str):
    d_str, o_str = str(d_str).strip(), str(o_str).strip()
    if not d_str or not o_str or d_str.lower() == 'nan' or o_str.lower() == 'nan': return None
    for fmt in ("%d/%m/%Y", "%Y-%m-%d"):
        try: d = datetime.strptime(d_str, fmt).date(); break
        except ValueError: continue
    else: return None
    try: t = datetime.strptime(o_str[:5], "%H:%M").time()
    except ValueError: return None
    return datetime.combine(d, t)

storia_completa = []
for _, r in df_dosi_edit.iterrows():
    dt_obj = _parse_dt(r.get("Data"), r.get("Ora"))
    try: dose_val = float(r.get("Dose (mg)"))
    except (TypeError, ValueError): dose_val = None
    if dt_obj is not None and dose_val is not None and not pd.isna(dose_val): storia_completa.append({'dt': dt_obj, 'dose': dose_val})

livelli_completi = []
for _, r in df_liv_edit.iterrows():
    dt_obj = _parse_dt(r.get("Data"), r.get("Ora"))
    try: conc_val = float(r.get("Conc. (mg/L)"))
    except (TypeError, ValueError): conc_val = None
    if dt_obj is not None and conc_val is not None and not pd.isna(conc_val): livelli_completi.append({'dt': dt_obj, 'conc': conc_val, 'nota': str(r.get("Nota", "") or "")})

# ------------------------------------------------------------
# AVVIO ANALISI E RISULTATI
# ------------------------------------------------------------

if st.button("🚀 Avvia il calcolo bayesiano", type="primary", use_container_width=True):
    if campi_mancanti: st.error("⚠️ Compila prima i campi obbligatori.")
    elif not storia_completa: st.error("Nessuna somministrazione valida in tabella.")
    elif not livelli_completi: st.error("Inserisci almeno un livello plasmatico valido.")
    else:
        t0_assoluto = min(s['dt'] for s in storia_completa)
        for s in storia_completa: s['ore'] = (s['dt'] - t0_assoluto).total_seconds() / 3600.0
        storia_completa.sort(key=lambda x: x['ore'])
        for l in livelli_completi: l['ore'] = (l['dt'] - t0_assoluto).total_seconds() / 3600.0
        livelli_completi.sort(key=lambda x: x['ore'])

        dt_val = (max((l['dt'] for l in livelli_completi), default=None) or max((s['dt'] for s in storia_completa), default=datetime.now())).strftime('%Y-%m-%d %H:%M')

        storia_completa = [{'ore': s['ore'], 'dose': float(s['dose'])} for s in storia_completa]
        livelli_completi = [{'ore': l['ore'], 'conc': float(l['conc']), 'nota': l.get('nota', '')} for l in livelli_completi]

        with st.spinner("Calcolo MAP Bayesiano in corso…"):
            dose_mg = storia_completa[-1]['dose']; dose_mg_kg = dose_mg / peso; tau_h_display = tau_mediano(storia_completa)
            ok_staph = ok_ent = ok_target = tox = None

            if farmaco == "Daptomicina":
                prior = prior_dapto(id_paz, sesso, gg, peso, crea, crrt, Sc, Qef, dialisi)
                CL_ind, V_ind, tipo_fit = map_1cmt(prior['CL_p'], prior['V_p'], prior['oCL'], prior['oV1'], prior['sig'], storia_completa, livelli_completi, prior['t_inf'])
                Q_ind = V2_ind = 0.0; Vss_d = V_ind
                AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_1cmt, (CL_ind, V_ind), storia_completa, prior['t_inf'])
                auc_mic = AUC24_v / mic; emivita = float(np.log(2)*V_ind/CL_ind)
                ok_staph = auc_mic >= PD_DAPTO['AUC_MIC_staph']; ok_ent = auc_mic >= PD_DAPTO['AUC_MIC_entero']; tox = Cmin_ss >= PD_DAPTO['Cmin_tox']

            elif farmaco == "Dalbavancina":
                prior = prior_dalba(id_paz, peso, modello_dalba)
                CL_ind, tipo_fit = map_dalba(prior['CL_p'], prior['V_p'], prior['Q_p'], prior['V2_p'], prior['oCL'], prior['sig'], storia_completa, livelli_completi, prior['t_inf'])
                V_ind, Q_ind, V2_ind = prior['V_p'], prior['Q_p'], prior['V2_p']; Vss_d = V_ind + V2_ind
                AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_2cmt, (CL_ind, V_ind, Q_ind, V2_ind), storia_completa, prior['t_inf'])
                auc_mic = AUC24_v / mic; emivita = float(np.log(2)*Vss_d/CL_ind)
                ok_target = auc_mic >= DALBA_TARGETS['ABSSSI semplice (fAUC24h/MIC >= 27.1)']; tox = False

            else:
                prior = prior_cipro(id_paz, sesso, gg, peso, crea, crrt, dialisi)
                CL_ind, V_ind, Q_ind, V2_ind, tipo_fit = map_2cmt(prior['CL_p'], prior['V_p'], prior['Q_p'], prior['V2_p'], prior['oCL'], prior['oV1'], prior['oQ'], prior['oV2'], prior['sig'], storia_completa, livelli_completi, prior['t_inf'])
                Vss_d = V_ind + V2_ind
                AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_2cmt, (CL_ind, V_ind, Q_ind, V2_ind), storia_completa, prior['t_inf'])
                auc_mic = AUC24_v / mic; emivita = float(np.log(2)*Vss_d/CL_ind)
                ok_target = auc_mic >= PD_CIPRO['AUC_MIC_target']; tox = None

            st.session_state['risultati'] = {
                'farmaco': farmaco, 'id_paz': id_paz, 'CL_ind': CL_ind, 'V_ind': V_ind, 'Q_ind': Q_ind, 'V2_ind': V2_ind,
                'Vss_d': Vss_d, 'AUC24_v': AUC24_v, 'Cmax_ss': Cmax_ss, 'Cmin_ss': Cmin_ss, 'auc_mic': auc_mic, 'emivita': emivita,
                'tipo_fit': tipo_fit, 'prior': prior, 'storia_completa': storia_completa, 'livelli_completi': livelli_completi,
                't0_assoluto': t0_assoluto, 'dt_val': dt_val, 'dose_mg': dose_mg, 'dose_mg_kg': dose_mg_kg, 'tau_h_display': tau_h_display,
                'mic': mic, 'ok_staph': ok_staph, 'ok_ent': ok_ent, 'ok_target': ok_target, 'tox': tox, 'sesso': sesso, 'eta': eta, 'gg': gg,
                'peso': peso, 'crea': crea, 'dialisi': dialisi, 'modello_dalba': modello_dalba, 'crrt': crrt, 'Sc': Sc, 'Qef': Qef, 'registrato': False,
            }
            st.info("🧮 Calcolo completato.")

# ------------------------------------------------------------
# VISUALIZZAZIONE RISULTATI
# ------------------------------------------------------------

if ('risultati' in st.session_state and st.session_state['risultati'].get('farmaco') == farmaco and st.session_state['risultati'].get('id_paz') == id_paz):
    R = st.session_state['risultati']
    st.markdown("---")
    st.subheader("🔬 Risultati PK/PD")
    
    c1, c2, c3 = st.columns(3)
    c1.metric("Clearance (CL) del paziente", f"{R['CL_ind']:.3f} L/h", delta=f"Normale: {R['prior']['CL_std']:.3f} L/h", delta_color="off")
    c2.metric("Volume di distribuzione (Vss)", f"{R['Vss_d']:.1f} L", delta=f"Normale: {(R['prior']['V_std'] + R['prior']['V2_std']):.1f} L", delta_color="off")
    c3.metric("Emivita (T½)", f"{R['emivita']:.1f} h")

    c4, c5, c6 = st.columns(3)
    c4.metric("AUC₂₄ stimata", f"{R['AUC24_v']:.0f} mg·h/L")

    if farmaco == "Daptomicina":
        c5.metric("AUC₂₄/MIC (≥ 666)", f"{R['auc_mic']:.0f}", delta="✅ Target raggiunto" if R['ok_staph'] else "⚠️ Sotto target", delta_color="normal" if R['ok_staph'] else "inverse")
        c6.metric("Cmin (Soglia CPK ≥ 24.3)", f"{R['Cmin_ss']:.1f} mg/L", delta="⚠️ Rischio CPK" if R['tox'] else "✅ Sicuro", delta_color="inverse" if R['tox'] else "normal")
    elif farmaco == "Dalbavancina":
        c5.metric("fAUC₂₄h/MIC (≥ 27.1)", f"{R['auc_mic']:.0f}", delta="✅ Target raggiunto" if R['ok_target'] else "⚠️ Sotto target", delta_color="normal" if R['ok_target'] else "inverse")
        c6.metric("Cmax / Cmin previste", f"{R['Cmax_ss']:.1f} / {R['Cmin_ss']:.1f} mg/L")
    else:
        c5.metric("AUC₂₄/MIC (≥ 125)", f"{R['auc_mic']:.0f}", delta="✅ Target raggiunto" if R['ok_target'] else "⚠️ Sotto target", delta_color="normal" if R['ok_target'] else "inverse")
        c6.metric("AUC₂₄ (Soglia neurotox)", f"{R['AUC24_v']:.0f} mg·h/L", help="Nessuna soglia convalidata.")

    st.markdown("### 📈 Profilo farmacocinetico")
    fig = disegna_grafico(farmaco, R['storia_completa'], R['livelli_completi'], R['CL_ind'], R['V_ind'], R['Q_ind'], R['V2_ind'], R['prior']['CL_std'], R['prior']['V_std'], R['prior']['Q_std'], R['prior']['V2_std'], R['prior']['t_inf'])
    st.plotly_chart(fig, use_container_width=True, theme=None)

    def _esito(ok): return "TARGET RAGGIUNTO" if ok else "INSUFFICIENTE"
    if farmaco == "Dalbavancina": pd_block = f"Efficacia fAUC24h/MIC: {R['auc_mic']:.0f} -> {_esito(R['ok_target'])}\nCmax prevista: {R['Cmax_ss']:.1f} mg/L\nCmin prevista: {R['Cmin_ss']:.1f} mg/L"
    elif farmaco == "Daptomicina": pd_block = f"Efficacia S. aureus: {R['auc_mic']:.0f} -> {_esito(R['ok_staph'])}\nEfficacia E. faecium: {R['auc_mic']:.0f} -> {_esito(R['ok_ent'])}\nSicurezza (Cmin >= 24.3): {R['Cmin_ss']:.1f} mg/L -> {'ATTENZIONE' if R['tox'] else 'SICURO'}"
    else: pd_block = f"Efficacia gram-negativi: {R['auc_mic']:.0f} -> {_esito(R['ok_target'])}\nEsposizione sistemica: {R['AUC24_v']:.0f} mg*h/L"

    referto = f"""=== TDM BUDDY v6.1 — REFERTO {farmaco.upper()} ===
Data: {R['dt_val']} | Paziente: {R['id_paz']} | Sesso: {R['sesso']} | Età: {R['eta']:.1f} anni | Peso: {R['peso']} kg

[MODELLO APPLICATO]
{R['prior']['nm']}

[DATI CLINICI]
Creatinina sierica: {R['crea']} mg/dL
Dialisi / CRRT: {R['dialisi']}
Dose più recente: {R['dose_mg']:.0f} mg ({R['dose_mg_kg']:.1f} mg/kg)
MIC patogeno: {R['mic']} mg/L

[PARAMETRI PK INDIVIDUALI — MAP Bayesiano]
Clearance (CL): {R['CL_ind']:.3f} L/h
Volume di distribuzione (Vss): {R['Vss_d']:.2f} L
Emivita terminale (T1/2): {R['emivita']:.1f} h

[VALUTAZIONE PK/PD]
AUC24 stimata: {R['AUC24_v']:.0f} mg*h/L
{pd_block}"""

    st.markdown("### 📋 Referto per cartella clinica")
    st.code(referto, language='text')

    if farmaco == "Dalbavancina":
        with st.container(border=True):
            st.markdown("#### 🔮 Forecast — pianificazione della prossima dose")
            colF1, colF2, colF3 = st.columns(3)
            with colF1: target_scelto = st.selectbox("Target PK/PD:", list(DALBA_TARGETS.keys()) + ["Trough C0/MIC ≥ 200"], key="fc_target")
            with colF2:
                mic_label = st.selectbox("MIC / breakpoint:", list(DALBA_MIC_OPZIONI.keys()) + ["Personalizzata"], key="fc_mic_label")
                if mic_label == "Personalizzata": mic_fc = st.number_input("MIC personalizzata", 0.001, 4.0, 0.125, format="%.4f", key="fc_mic_custom")
                else: mic_fc = DALBA_MIC_OPZIONI[mic_label]
            with colF3: f_unb = st.number_input("Frazione libera (f)", 0.01, 0.10, DALBA_PARAMS['fu'], 0.01, key="fc_funb")

            if target_scelto == "Trough C0/MIC ≥ 200":
                soglia = dalba_soglia_trough_c0mic(mic_fc)
                desc = f"C0/MIC ≥ 200 → soglia totale = {soglia:.2f} mg/L"
            else:
                soglia = dalba_soglia_totale(DALBA_TARGETS[target_scelto], mic_fc, f_unb)
                desc = f"{target_scelto} → soglia totale = {soglia:.2f} mg/L"

            st.info(f"📐 Soglia calcolata: **{desc}**")
            t_val = max(max((l['ore'] for l in R['livelli_completi']), default=0), R['storia_completa'][-1]['ore'])
            c_att = pk_2cmt(R['CL_ind'], R['V_ind'], R['Q_ind'], R['V2_ind'], R['storia_completa'], t_val, R['prior']['t_inf'])
            t_cross = trova_tempo_soglia(R['CL_ind'], R['V_ind'], R['Q_ind'], R['V2_ind'], R['storia_completa'], soglia, R['prior']['t_inf'], t_partenza=t_val, t_max=t_val+2000, passo=2.0)

            if t_cross is None: st.warning("⚠️ Soglia non raggiunta entro ~83 giorni.")
            elif t_cross <= t_val + 0.5: st.error(f"🔴 Concentrazione già sotto soglia ({c_att:.2f} < {soglia:.2f} mg/L).")
            else: st.success(f"🗓️ La concentrazione scenderà sotto soglia il **{(R['t0_assoluto'] + timedelta(hours=t_cross)).strftime('%d/%m/%Y alle %H:%M')}** (tra circa {(t_cross-t_val)/24.0:.1f} giorni).")

            st.markdown("##### Simula una nuova somministrazione ipotetica")
            col_s1, col_s2, col_s3 = st.columns(3)
            with col_s1: simula = st.checkbox("Aggiungi dose ipotetica al forecast", key="fc_simula")
            with col_s2:
                if simula:
                    data_new = st.date_input("Data", value=(R['t0_assoluto'] + timedelta(hours=t_val+24)).date(), format="DD/MM/YYYY", key="fc_data_new")
                    ora_new = st.time_input("Ora", value=datetime.strptime("10:00","%H:%M").time(), key="fc_ora_new")
            with col_s3:
                if simula: dose_new = st.number_input("Dose ipotetica (mg)", 100.0, 1500.0, 1500.0, 50.0, key="fc_dose_new")

            if simula:
                t_new = (datetime.combine(data_new, ora_new) - R['t0_assoluto']).total_seconds()/3600.0
                storia_sim = sorted(R['storia_completa'] + [{'ore': t_new, 'dose': dose_new}], key=lambda x: x['ore'])
                t_cross_sim = trova_tempo_soglia(R['CL_ind'], R['V_ind'], R['Q_ind'], R['V2_ind'], storia_sim, soglia, R['prior']['t_inf'], t_partenza=t_new, t_max=t_new+2000, passo=2.0)
                if t_cross_sim: st.success(f"✅ Con {dose_new:.0f} mg, concentrazione sopra soglia fino al **{(R['t0_assoluto'] + timedelta(hours=t_cross_sim)).strftime('%d/%m/%Y')}**.")
                t_arr_sim = np.linspace(0, (t_cross_sim + 48) if t_cross_sim else (t_new + 500), 1000)
                c_arr_sim = [pk_2cmt(R['CL_ind'], R['V_ind'], R['Q_ind'], R['V2_ind'], storia_sim, t, R['prior']['t_inf']) for t in t_arr_sim]
                fig2 = disegna_forecast(t_arr_sim, c_arr_sim, soglia, t_new, dose_new, R['storia_completa'])
                st.plotly_chart(fig2, use_container_width=True, theme=None)

    st.markdown("---")
    if R.get('registrato'):
        st.success("✅ Questa valutazione è già stata registrata nello storico di sessione.")
    else:
        if st.button("💾 Registra questa valutazione nello storico", type="primary", use_container_width=True):
            record = {
                'data_val': R['dt_val'], 'farmaco': farmaco, 'eta_giorni': R['gg'], 'peso': R['peso'],
                'creatinina': R['crea'], 'crrt': bool(R['crrt']), 'Sc': R['Sc'], 'Qef': R['Qef'], 'dialisi_tipo': R['dialisi'],
                'dose_mg': R['dose_mg'], 'dose_mg_kg': R['dose_mg_kg'], 'tau_h': R['tau_h_display'],
                'n_dosi': len(R['storia_completa']), 'mic': R['mic'], 'n_livelli': len(R['livelli_completi']),
                'CL_ind': R['CL_ind'], 'V_ind': R['Vss_d'], 'Q_ind': R['Q_ind'], 'V2_ind': R['V2_ind'], 'AUC24': R['AUC24_v'], 
                'AUC_MIC': R['auc_mic'], 'Cmin_ss': R['Cmin_ss'], 'modello': R['prior']['nm'], 'modello_dalba': R['modello_dalba'],
                'storia_dosi': R['storia_completa'], 'livelli': R['livelli_completi'], 'data_prima_dose': R['t0_assoluto'].strftime('%Y-%m-%d %H:%M'),
            }
            st.session_state['storico']['monitoraggi'].append(record)
            st.session_state['storico']['codice_paziente'] = id_paz
            st.session_state['storico']['sesso'] = sesso
            st.session_state['storico']['farmaco'] = farmaco
            st.session_state['risultati']['registrato'] = True
            st.rerun()
