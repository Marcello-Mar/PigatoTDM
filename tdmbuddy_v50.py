"""
TDM Buddy v5.2
==============

FUNZIONALITÀ PRINCIPALI
  - MAP Bayesiano 1-compartimento (Daptomicina) e 2-compartimento (Ciprofloxacina)
  - Integrazione storia dosi: la TDM successiva recupera automaticamente
    tutte le dosi precedenti dallo storico di sessione e le fonde nella timeline, permettendo
    di stimare l'esposizione cumulativa anche dopo modifiche di posologia.
  - Archivio clinico SQLite con storico per paziente
  - Export PDF referto strutturato

MODELLI INTEGRATI — DAPTOMICINA
  [1] Soraluce A, Asín-Prieto E, Rodríguez-Gascón A, Barrasa H, Maynar J,
      Carcelero E, Soy D, Isla A. Population pharmacokinetics of daptomycin
      in critically ill patients. Int J Antimicrob Agents. 2018 Aug;52(2):
      158-165. doi:10.1016/j.ijantimicag.2018.03.008. PMID:29572042
  [2] Wei XC, Zhao MF, Li X, Xiao X. Pharmacokinetic/Pharmacodynamic
      Analysis of Daptomycin Against Staphylococcus aureus and Enterococcus
      faecium in Pediatric Patients by Monte Carlo Simulation. J Clin
      Pharmacol. 2020 Jun;60(6):768-774. doi:10.1002/jcph.1576.
  [3] Ye L, You X, Zhou J, Wu C, Ke M, Wu W, Huang P, Lin C. Physiologically
      based pharmacokinetic modeling of daptomycin dose optimization in
      pediatric patients with renal impairment. Front Pharmacol. 2022 Aug
      16;13:838599. doi:10.3389/fphar.2022.838599.
  [4] Ye L, Zhou H, Guo G, Chen M, Zhang J. Physiologically-based
      pharmacokinetic modeling to predict the exposure and to assess
      pharmacodynamics of daptomycin in infants within 1 year old. Eur J
      Pharm Sci. 2025 Mar;208:107058. doi:10.1016/j.ejps.2025.107058.
  [5] Lim SY, Lewis T, Woo S, Turman M, Bourne DWA, Burton ME, Rianthavorn P.
      Daptomycin Pharmacokinetics in Adolescents Undergoing Hemodialysis and
      Peritoneal Dialysis: A Case Series With Pharmacokinetic Modeling. J
      Pediatr Pharmacol Ther. 2021 Feb 15;26(2):123-132.
      doi:10.5863/1551-6776-26.2.123.
  [6] Bhavnani SM, Rubino CM, Ambrose PG, Drusano GL. Daptomycin exposure
      and the probability of elevations in the creatine phosphokinase level:
      data from a randomized trial of patients with bacteremia and
      endocarditis. Clin Infect Dis. 2010 Jun 15;50(12):1568-1574.
      doi:10.1086/652767.

MODELLI INTEGRATI — CIPROFLOXACINA
  [7] Schlender JF, Teutonico D, Coboeken K, Schnizler K, Eissing T,
      Willmann S, Jaehde U, Stass H. A Physiologically-Based Pharmacokinetic
      Model to Describe Ciprofloxacin Pharmacokinetics Over the Entire Span
      of Life. Clin Pharmacokinet. 2018 Dec;57(12):1613-1634.
      doi:10.1007/s40262-018-0661-6.
  [8] Balbas-Martinez V, Michelet R, Edginton AN, Meesters K, Trocóniz IF,
      Vermeulen A. Physiologically-Based Pharmacokinetic model for
      Ciprofloxacin in children with complicated Urinary Tract Infection. Eur
      J Pharm Sci. 2019 Feb 1;128:171-179. doi:10.1016/j.ejps.2018.11.033.
  [9] Forrest A, Ballow CH, Nix DE, Birmingham MC, Schentag JJ. Development
      of a population pharmacokinetic model and optimal sampling strategies
      for intravenous ciprofloxacin. Antimicrob Agents Chemother. 1993
      May;37(5):1065-1072. doi:10.1128/aac.37.5.1065.
  [10] Ambrose PG, Bhavnani SM, Rubino CM, Goss TF, Bhavnani SM, Drusano GL.
       Pharmacokinetics-pharmacodynamics of antimicrobial therapy: it's not
       just for mice anymore. Clin Infect Dis. 2007;44(Suppl 2):S79-S86.
       doi:10.1086/510079.
"""

import streamlit as st
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import pandas as pd
from datetime import datetime, date, timedelta
import json
from fpdf import FPDF

# ============================================================
# STORAGE — nessun database. Tutto vive in st.session_state
# (isolato per sessione/browser) e viene esportato/importato come
# file JSON che l'operatore scarica e gestisce localmente.
# Nessun dato tocca il disco del server in nessun momento.
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
# Fonte: "Population Pharmacokinetics of Dalbavancin in Pediatric
# Patients: A Real-World Therapeutic Drug Monitoring Study" (dati
# IRCCS Istituto Giannina Gaslini, Genova; modello non ancora
# pubblicato, fornito dall'utente). Modello 2-compartimenti,
# SAEM/nlmixr2-rxode2, 35 pazienti / 83 eventi di dose / 137 concentrazioni.
#
# Parametri tipici (individuo di riferimento 70 kg):
#   lCL = -3.49  -> CL_typ = exp(-3.49) = 0.0305 L/h   (RSE 1.54%, BSV 19.8% CV, shrink 20.9%)
#   lV1 =  1.57  -> V1_typ = exp( 1.57) = 4.81  L      (RSE 3.87%, nessuna IIV)
#   lQ  = -4.81  -> Q_typ  = exp(-4.81) = 0.0081 L/h   (RSE 3.24%, nessuna IIV)
#   lV2 =  1.57  -> V2_typ = exp( 1.57) = 4.81  L      (RSE 15.9%, nessuna IIV)
#   Errore residuo proporzionale = 0.375 (37.5%)
#
# Scaling allometrico (covariata a priori, esponenti fissi):
#   CL = CL_typ * (WT/70)^0.75 * exp(eta_CL)
#   V1 = V1_typ * (WT/70)
#   Q  = Q_typ  * (WT/70)^0.75
#   V2 = V2_typ * (WT/70)
# Nessuna covariata aggiuntiva (eGFR, albumina, eta') ha superato i
# criteri di inclusione: la IIV e' modellata SOLO su CL; V1/Q/V2 sono
# quindi valori deterministici (fissi) una volta noto il peso corporeo.
DALBA_PARAMS = {
    'CL_typ70': 0.0305, 'V1_typ70': 4.81, 'Q_typ70': 0.0081, 'V2_typ70': 4.81,
    'omega_CL': 0.198, 'sigma': 0.375, 't_inf': 0.5, 'fu': 0.07,
    'nome': ("Modello Mariani (dati IRCCS Gaslini, non ancora pubblicato, "
             "modello interno fornito dall'utente)"),
}

# Modello alternativo — Cojutti PG, et al. Population pharmacokinetics of
# dalbavancin and dosing consideration for optimal treatment of adult
# patients with staphylococcal osteoarticular infections. Antimicrob
# Agents Chemother. 2021;65(9):e00879-21. Popolazione: 15 pazienti adulti
# con infezioni osteoarticolari da stafilococco.
#
# Parametri tipici di popolazione (valori assoluti, NON normalizzati per
# peso — nel dato fornito non è presente un'equazione di scaling
# allometrico esplicita per questo modello, a differenza di Mariani).
# Applicati quindi come valori fissi di popolazione, appropriati per un
# paziente adulto vicino al peso medio della coorte originale; usare con
# cautela per pazienti con peso molto diverso dalla popolazione di
# derivazione (15 adulti con infezione osteoarticolare).
#   CL = 0.106 L/h, V1 = 17.4 L, Q = 0.103 L/h, V2 = 15.1 L
#   IIV_CL = 36.2% CV, errore residuo proporzionale = 18%
#   Emivita terminale calcolata dal modello: 11.3 giorni
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

# Target PK/PD dalbavancina — riferimenti citati nel documento sorgente:
#   fAUC24h/MIC >= 27.1   -> target ABSSSI semplice
#   fAUC24h/MIC >= 111.1  -> target infezioni complesse/osteoarticolari
#                            (Cojutti P, et al. Antimicrob Agents Chemother 2023)
#   100% fT > 4xMIC       -> target trough, attivita' battericida tempo-dipendente
# Soglie di concentrazione TOTALE equivalenti (f=0.07, 93% legame proteico):
#   C_total = (target_AUC/MIC * MIC) / (24 * f)
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
    """C_total target = (target * MIC) / (24 * f) -- soglia media/AUC-equivalente."""
    return (target_auc_mic * mic) / (24.0 * f)

def dalba_soglia_trough(mic, f):
    """Soglia di trough: concentrazione totale equivalente a 4x MIC libera."""
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
    """
    Calcola AUC, Cmax e Cmin numericamente sulle 24h successive
    all'ultima dose somministrata, indipendentemente dal fatto che il
    dosaggio sia regolare o irregolare (nessun concetto di 'tau' fisso
    richiesto). pk_func e' pk_1cmt o pk_2cmt; params sono i parametri
    posizionali richiesti da pk_func DOPO 'storia, t, t_inf' (es. (CL,V)
    per 1cmt, (CL,V1,Q,V2) per 2cmt).
    Restituisce (AUC24, Cmax, Cmin).
    """
    t_ultima_dose = storia[-1]['ore']
    t_grid = np.linspace(t_ultima_dose, t_ultima_dose + 24, n_punti)
    c_grid = np.array([pk_func(*params, storia, t, t_inf) for t in t_grid])
    # Trapezoidale manuale: np.trapz è stata rimossa in NumPy 2.0
    # (np.trapezoid la sostituisce, ma questa forma funziona su ogni versione)
    auc = float(np.sum((c_grid[1:] + c_grid[:-1]) / 2.0 * np.diff(t_grid)))
    return auc, float(c_grid.max()), float(c_grid.min())


def tau_mediano(storia):
    """Intervallo mediano tra dosi consecutive -- SOLO a scopo informativo/
    di visualizzazione nel referto; nessun calcolo dipende da questo valore."""
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
        elif gg<2190: ag='2-6y'
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
    return CL,V1,Q,V2,c['omega_CL'],c['omega_V1'],c['omega_Q'],c['sigma'],1.0,nm,note


def prior_dalba(peso, modello_scelto):
    """
    Prior dalbavancina -- supporta due modelli alternativi:
      - 'Mariani (pediatrico, allometrico — Gaslini)': scaling allometrico
        a esponenti fissi (0.75 su CL/Q, 1 su V1/V2), normalizzato a 70 kg.
      - 'Cojutti 2021 (adulti, osteoarticolare — AAC)': valori di
        popolazione fissi, nessuno scaling per peso disponibile nella
        pubblicazione originale.
    """
    p = DALBA_MODELLI[modello_scelto]

    if modello_scelto.startswith('Mariani'):
        fw = peso / 70.0
        CL_pop = p['CL_typ70'] * (fw ** 0.75)
        V1_pop = p['V1_typ70'] * fw
        Q_pop  = p['Q_typ70']  * (fw ** 0.75)
        V2_pop = p['V2_typ70'] * fw
        note = ("Scaling allometrico applicato (CL,Q: peso/70^0.75; V1,V2: peso/70). "
                "Nessuna covariata su eGFR/albumina/eta' ritenuta significativa nel "
                "modello originale -- per insufficienza renale grave (ClCr<30 mL/min) "
                "fare riferimento alla scheda tecnica DALVANCE (FDA 2025) per "
                "l'aggiustamento posologico, non incorporato nel modello PK.")
    else:
        CL_pop, V1_pop, Q_pop, V2_pop = p['CL_typ70'], p['V1_typ70'], p['Q_typ70'], p['V2_typ70']
        note = ("⚠️ Modello Cojutti 2021: nessuno scaling per peso disponibile nella "
                "pubblicazione originale -- valori di popolazione applicati come fissi, "
                "indipendentemente dal peso inserito. Derivato da 15 pazienti adulti "
                "con infezione osteoarticolare: usare con cautela per pazienti "
                "molto diversi da questa popolazione.")

    return CL_pop, V1_pop, Q_pop, V2_pop, p['omega_CL'], p['sigma'], p['t_inf'], p['nome'], note

# ============================================================
# MAP BAYESIANO
# ============================================================

def map_1cmt(CL_p, V_p, oCL, oV, sig, storia, livelli, t_inf):
    """
    MAP Bayesiano 1-compartimento su un numero QUALSIASI di livelli
    plasmatici misurati (0, 1, 2, ... N). Ogni elemento di 'livelli' e'
    un dict {'ore': <ore da t0>, 'conc': <mg/L misurata>}.
    """
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


def map_2cmt(CL_p, V1_p, Q_p, V2_p, oCL, oV1, oQ, sig, storia, livelli, t_inf):
    """
    MAP Bayesiano 2-compartimento su un numero QUALSIASI di livelli.
    Con pochi campioni (<=2) Q e V2 restano ancorati stretti al prior
    (identificabilita' limitata); con piu' campioni il vincolo si
    allenta automaticamente in proporzione (omega scala con 1/sqrt(n),
    fino al valore "libero" oQ/oV2 originale).
    """
    n = len(livelli)
    if n == 0:
        return CL_p, V1_p, Q_p, V2_p, "Prior di popolazione (nessun livello inserito)"

    # con pochi campioni Q/V2 sono poco identificabili: parto stretto
    # (0.10) e allento progressivamente con n, senza mai superare oQ/oV2
    fattore_allentamento = min(1.0, n / 6.0)
    OQ_eff  = 0.10 + (oQ  - 0.10) * fattore_allentamento
    OV2_eff = 0.10 + (oQ  - 0.10) * fattore_allentamento  # stessa cinetica di allentamento

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


def map_dalba(CL_pop, V1_pop, Q_pop, V2_pop, omega_CL, sigma,
              storia, livelli, t_inf):
    """
    MAP dalbavancina: stima SOLO log(CL) (unico parametro con IIV nel
    modello originale) su un numero QUALSIASI di livelli plasmatici
    misurati nel tempo -- essenziale per il TDM cumulativo di un
    farmaco long-acting come la dalbavancina, dove il paziente puo'
    avere decine di prelievi su mesi di follow-up. V1, Q, V2 restano
    fissi ai valori di popolazione (nessuna IIV ritenuta significativa
    su questi tre parametri nei modelli sorgente).
    """
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
    r = minimize(f, x0, method='L-BFGS-B', bounds=bnd,
                options={'maxiter': 2000, 'ftol': 1e-12})
    lbl = f"MAP Bayesiano su CL, {n} livell{'o' if n==1 else 'i'} plasmatic{'o' if n==1 else 'i'}; V1/Q/V2 fissi (allometria/popolazione)"
    return np.exp(r.x[0]), lbl


def trova_tempo_soglia(CL, V1, Q, V2, storia, soglia, t_inf, t_partenza=0.0, t_max=2000.0, passo=1.0):
    """
    Cerca (scansione + bisezione) il primo istante t > t_partenza in cui
    la concentrazione prevista scende sotto 'soglia'. Assume curva
    monotona decrescente dopo l'ultima dose (nessuna nuova somministrazione
    tra t_partenza e t_max). Restituisce None se non raggiunta entro t_max
    o gia' sotto soglia al momento della valutazione (in tal caso ritorna
    t_partenza stesso).
    """
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
    """
    Carica dosi E livelli plasmatici più recenti dallo storico di
    SESSIONE (nessun database) per il farmaco indicato. Restituisce
    (storia_dosi, livelli_precedenti, data_prima_dose_dt) oppure
    ([], [], None) se non presente.
    """
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
    """
    Fonde la storia precedente (già in ore assolute da t0_assoluto)
    con le nuove somministrazioni (datetime objects).
    nuove_somm: lista di dict {'dt': datetime, 'dose': float}
    Restituisce lista ordinata di {'ore': float, 'dose': float}.
    """
    storia = list(storia_precedente)
    for s in nuove_somm:
        ore = (s['dt'] - t0_assoluto).total_seconds() / 3600.0
        storia.append({'ore': ore, 'dose': s['dose']})
    storia.sort(key=lambda x: x['ore'])
    return storia


def costruisci_livelli(nuovi_liv, t0_assoluto, livelli_precedenti):
    """
    Fonde i livelli plasmatici precedenti (già in ore assolute da
    t0_assoluto) con i nuovi livelli misurati (datetime objects).
    nuovi_liv: lista di dict {'dt': datetime, 'conc': float, 'nota': str,
                              'creatinina': float|None}
    Restituisce lista ordinata di {'ore': float, 'conc': float, 'nota': str,
                                    'creatinina': float|None}.
    Nessun limite al numero di livelli: possono essere quanti il
    paziente ne ha effettivamente avuti nel corso della terapia.
    La creatinina per livello è puramente informativa/di registrazione:
    NON è usata come covariata in nessun calcolo PK (il valore di
    creatinina che guida i modelli resta quello inserito in "Parametri
    Clinici" per la valutazione corrente).
    """
    livelli = list(livelli_precedenti)
    for l in nuovi_liv:
        ore = (l['dt'] - t0_assoluto).total_seconds() / 3600.0
        livelli.append({'ore': ore, 'conc': l['conc'], 'nota': l.get('nota', ''),
                        'creatinina': l.get('creatinina')})
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
    pdf.cell(0,9,"TDM Buddy v5.2 — Referto di Monitoraggio Terapeutico Farmacologico",ln=1,align='C')
    pdf.set_font("Arial",size=9); pdf.ln(2)
    for riga in safe.split('\n'):
        pdf.multi_cell(0,5,txt=riga)
    return pdf.output(dest='S').encode('latin-1')

# ============================================================
# GRAFICO
# ============================================================

def disegna_grafico(farmaco, storia, livelli,
                    CL_ind, V_ind, Q_ind, V2_ind, CL_pop, V_pop, Q_pop, V2_pop,
                    t_inf, storia_precedente, t_fine_extra=12):
    t_ultimo_livello = max((l['ore'] for l in livelli), default=0)
    t_fine = max(t_ultimo_livello,
                 max(s['ore'] for s in storia) + 12) + t_fine_extra
    t_arr = np.linspace(0, t_fine, 800)

    if farmaco == "Daptomicina":
        C_ind = [pk_1cmt(CL_ind, V_ind, storia, t, t_inf) for t in t_arr]
        C_pop = [pk_1cmt(CL_pop, V_pop, storia, t, t_inf) for t in t_arr]
    else:
        C_ind = [pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia, t, t_inf) for t in t_arr]
        C_pop = [pk_2cmt(CL_pop, V_pop, Q_pop, V2_pop, storia, t, t_inf) for t in t_arr]

    fig, ax = plt.subplots(figsize=(13, 4.5))
    fig.patch.set_facecolor('#FAFAFA')
    ax.set_facecolor('#FAFAFA')

    # Ombra periodo storico vs nuovo
    if storia_precedente:
        t_cambio = max(s['ore'] for s in storia_precedente)
        ax.axvspan(0, t_cambio, alpha=0.06, color='#6EA8D8', zorder=0, label='Terapia precedente (storico sessione)')
        ax.axvline(t_cambio, color='#6EA8D8', ls=':', lw=1.2, alpha=0.7)
        ax.text(t_cambio+0.3, ax.get_ylim()[1]*0.95 if ax.get_ylim()[1]>0 else 5,
                'Cambio posologia', fontsize=7, color='#3B6EA8', va='top')

    ax.plot(t_arr, C_pop, color='#94A3B8', lw=1.8, ls='--', alpha=0.8,
            label='Profilo Popolazione (Prior)')
    ax.plot(t_arr, C_ind, color='#0F172A', lw=2.5,
            label='Profilo Individuale (MAP Bayesiano)')

    if farmaco == "Daptomicina":
        ax.axhline(24.3, color='#DC2626', ls='-.', lw=1.5,
                   label='Soglia tossicità CPK — Cmin ≥ 24.3 mg/L\n'
                         '(Bhavnani SM, et al. Clin Infect Dis 2010;50(12):1568-74)')

    # Marcatori dosi — distingui vecchie da nuove
    t_cambio_val = max(s['ore'] for s in storia_precedente) if storia_precedente else -1
    for s in storia:
        col = '#6EA8D8' if s['ore'] <= t_cambio_val else '#475569'
        ax.axvline(s['ore'], color=col, ls=':', alpha=0.35, lw=1.0)

    # Livelli misurati — quanti ce ne sono, tutti a scatter, con etichetta
    # solo per i primi (evita legende infinite con decine di prelievi)
    MAX_LABEL = 6
    for i, l in enumerate(livelli):
        etichetta = (f"Livello {i+1} = {l['conc']:.2f} mg/L (t={l['ore']:.0f}h)"
                    if i < MAX_LABEL else None)
        ax.scatter(l['ore'], l['conc'], color='#16A34A', s=140, zorder=7,
                   marker='*', edgecolors='white', linewidths=0.5,
                   label=etichetta)
    if len(livelli) > MAX_LABEL:
        ax.scatter([], [], color='#16A34A', marker='*', s=140,
                   label=f"... + altri {len(livelli)-MAX_LABEL} livelli")

    ax.set_xlabel("Ore dall'inizio della terapia (t = 0 prima dose)", fontsize=10)
    ax.set_ylabel(f"Concentrazione plasmatica {farmaco} (mg/L)", fontsize=10)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color('#CBD5E1')
    ax.grid(True, alpha=0.20, color='#CBD5E1')
    ax.legend(fontsize=8, loc='upper right', framealpha=0.92,
              edgecolor='#CBD5E1', fancybox=False)
    plt.tight_layout()
    return fig

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
    tau_h_display = tau_mediano(storia_completa)  # SOLO informativo nel referto

    # --- Prior + MAP (su un numero QUALSIASI di livelli, 0..N) ---
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
        V_ind, Q_ind, V2_ind = V1_p, Qp, V2_p   # fissi (nessuna IIV nel modello originale)
        Q_p = Qp; Vss_d = V_ind + V2_ind
        AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_2cmt, (CL_ind,V_ind,Q_ind,V2_ind), storia_completa, t_inf)
        auc_mic  = AUC24_v / mic
        emivita  = np.log(2)*Vss_d/CL_ind
        ok_target= auc_mic >= DALBA_TARGETS['ABSSSI semplice (fAUC24h/MIC >= 27.1)']
        tox      = False  # nessuna soglia di tossicita' da concentrazione validata per dalba
    else:
        CL_p,V1_p,Qp,V2_p,oCL,oV1,oQ,sig,t_inf,nm,note = prior_cipro(sesso,gg,peso,crea,crrt,dialisi)
        CL_ind,V_ind,Q_ind,V2_ind,tipo_fit = map_2cmt(CL_p,V1_p,Qp,V2_p,oCL,oV1,oQ,sig,storia_completa,livelli_completi,t_inf)
        Q_p=Qp; Vss_d=V_ind+V2_ind
        AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(pk_2cmt, (CL_ind,V_ind,Q_ind,V2_ind), storia_completa, t_inf)
        auc_mic  = AUC24_v/mic
        emivita  = np.log(2)*Vss_d/CL_ind
        ok_target= auc_mic>=PD_CIPRO['AUC_MIC_target']
        tox      = AUC24_v>=PD_CIPRO['AUC_alta']

    # ---- DISPLAY ----
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
              help="1-cmt (dapto): V; 2-cmt (cipro): V1+V2")
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
                  delta="⚠️ RISCHIO CPK — MONITORARE" if tox else "✅ Sicuro",
                  delta_color="inverse" if tox else "normal")
    elif farmaco=="Dalbavancin":
        c5.metric("fAUC₂₄h/MIC  (target ABSSSI ≥ 27.1)",f"{auc_mic:.0f}",
                  delta="✅ TARGET" if ok_target else "⚠️ SOTTO TARGET",
                  delta_color="normal" if ok_target else "inverse",
                  help="fAUC calcolata numericamente sulle 24h successive all'ultima dose")
        c6.metric("Cmax / Cmin previste",f"{Cmax_ss:.1f} / {Cmin_ss:.1f} mg/L",
                  help="Vedi sezione Forecast qui sotto per target su infezioni complesse e pianificazione redose")
    else:
        c5.metric("AUC₂₄/MIC  (target gram-neg ≥ 125)",f"{auc_mic:.0f}",
                  delta="✅ TARGET" if ok_target else "⚠️ SOTTO TARGET",
                  delta_color="normal" if ok_target else "inverse")
        c6.metric(f"AUC₂₄  (soglia neurotox ≥ 400 mg·h/L)",f"{AUC24_v:.0f} mg·h/L",
                  delta="⚠️ AUC ELEVATA — monitorare CNS" if tox else "✅ Esposizione sicura",
                  delta_color="inverse" if tox else "normal")

    # Grafico
    st.subheader("📈 Profilo Farmacocinetico — Timeline Completa")
    fig = disegna_grafico(farmaco, storia_completa, livelli_completi,
                          CL_ind, V_ind, Q_ind, V2_ind, CL_p, V_p if farmaco=="Daptomicina" else V1_p,
                          Q_p, V2_p, t_inf, storia_precedente)
    st.pyplot(fig); plt.close()

    # Referto
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
            f"  Nota: per target su infezioni complesse/osteoarticolari\n"
            f"    (fAUC24h/MIC >= 111.1) e per la pianificazione della\n"
            f"    prossima somministrazione vedere la sezione Forecast.\n"
            f"  Nota: nessuna soglia di tossicita' da concentrazione e'\n"
            f"    validata per dalbavancina (a differenza della Cmin_tox\n"
            f"    daptomicina); il monitoraggio di sicurezza segue la\n"
            f"    scheda tecnica DALVANCE.\n"
        )
    elif farmaco=="Daptomicina":
        pd_block=(
            f"  Efficacia vs S. aureus / MRSA     (AUC24/MIC >= 666):\n"
            f"    AUC24/MIC = {auc_mic:.0f}  ->  {E(ok_staph)}\n"
            f"    [Rif: Soraluce A, et al. Int J Antimicrob Agents. 2018;52(2):158-165]\n"
            f"  Efficacia vs E. faecium           (AUC24/MIC >= 143):\n"
            f"    AUC24/MIC = {auc_mic:.0f}  ->  {E(ok_ent)}\n"
            f"    [Rif: Wei XC, et al. J Clin Pharmacol. 2020;60(6):768-774]\n"
            f"  Sicurezza — tossicita' muscolare  (Cmin_ss >= 24.3 mg/L):\n"
            f"    Cmin_ss = {Cmin_ss:.1f} mg/L  ->  "
            f"{'ATTENZIONE: MONITORARE CPK' if tox else 'SICURO'}\n"
            f"    [Rif: Bhavnani SM, et al. Clin Infect Dis. 2010;50(12):1568-1574]\n"
        )
    else:
        pd_block=(
            f"  Efficacia gram-negativi           (AUC24/MIC >= 125, EUCAST/Ambrose 2007):\n"
            f"    AUC24/MIC = {auc_mic:.0f}  ->  {E(ok_target)}\n"
            f"    [Rif: Ambrose PG, et al. Clin Infect Dis. 2007;44(Suppl 2):S79-S86]\n"
            f"  Esposizione sistemica: AUC24 = {AUC24_v:.0f} mg*h/L\n"
            f"  Rischio neurotossicita' (AUC24 >= 400 mg*h/L):\n"
            f"    {'MONITORARE: AUC ELEVATA' if tox else 'SICURO'}\n"
            f"  Nota: cipro non ha soglia Cmin validata come daptomicina.\n"
            f"  Nota: MAP 2-cmt con <=2 campioni -> Q/V2 ancorati al prior;\n"
            f"    stima CL e V1 affidabile; Q e V2 invariati rispetto alla popolazione.\n"
        )

    referto = f"""\
=== TDM BUDDY v5.2 — REFERTO {farmaco.upper()} ===
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
  Intervallo mediano dosi: {(f'{tau_h_display:.1f} h' if tau_h_display else "N/D (una sola dose)")}
  Durata infusione:        {t_inf*60:.0f} min
  MIC patogeno:            {mic} mg/L
  Livelli plasmatici usati per il fit: {len(livelli_completi)}
{chr(10).join(f"    t={l['ore']:.1f}h  ->  {l['conc']:.4f} mg/L" + (f"  [creat.={l['creatinina']:.2f} mg/dL]" if l.get('creatinina') is not None else "") + (f"  ({l['nota']})" if l.get('nota') else "") for l in livelli_completi) if livelli_completi else "    (nessun livello inserito — solo prior di popolazione)"}

[PARAMETRI PK INDIVIDUALI STIMATI]
  Clearance (CL):               {CL_ind:.3f} L/h  ({(CL_ind/peso)*1000:.1f} mL/h/kg)
  Volume di distribuzione (Vss):{Vss_d:.2f} L     ({Vss_d/peso:.2f} L/kg)
  Emivita terminale (T1/2):     {emivita:.1f} h

[VALUTAZIONE PK/PD]
  AUC24 stimata (finestra 24h dopo l'ultima dose): {AUC24_v:.0f} mg*h/L

{pd_block}
[RIFERIMENTI BIBLIOGRAFICI COMPLETI]
  DAPTOMICINA
  [1] Soraluce A, Asin-Prieto E, Rodriguez-Gascon A, Barrasa H, Maynar J,
      Carcelero E, Soy D, Isla A. Population pharmacokinetics of daptomycin
      in critically ill patients. Int J Antimicrob Agents. 2018 Aug;52(2):
      158-165. doi:10.1016/j.ijantimicag.2018.03.008. PMID: 29572042.
  [2] Wei XC, Zhao MF, Li X, Xiao X. Pharmacokinetic/Pharmacodynamic
      Analysis of Daptomycin Against Staphylococcus aureus and Enterococcus
      faecium in Pediatric Patients by Monte Carlo Simulation. J Clin
      Pharmacol. 2020 Jun;60(6):768-774. doi:10.1002/jcph.1576.
  [3] Ye L, You X, Zhou J, Wu C, Ke M, Wu W, Huang P, Lin C. Physiologically
      based pharmacokinetic modeling of daptomycin dose optimization in
      pediatric patients with renal impairment. Front Pharmacol. 2022 Aug
      16;13:838599. doi:10.3389/fphar.2022.838599.
  [4] Ye L, Zhou H, Guo G, Chen M, Zhang J. Physiologically-based
      pharmacokinetic modeling to predict the exposure and to assess
      pharmacodynamics of daptomycin in infants within 1 year old. Eur J
      Pharm Sci. 2025 Mar;208:107058. doi:10.1016/j.ejps.2025.107058.
  [5] Lim SY, Lewis T, Woo S, Turman M, Bourne DWA, Burton ME, Rianthavorn
      P. Daptomycin Pharmacokinetics in Adolescents Undergoing Hemodialysis
      and Peritoneal Dialysis: A Case Series With Pharmacokinetic Modeling.
      J Pediatr Pharmacol Ther. 2021 Feb 15;26(2):123-132.
      doi:10.5863/1551-6776-26.2.123.
  [6] Bhavnani SM, Rubino CM, Ambrose PG, Drusano GL. Daptomycin exposure
      and the probability of elevations in the creatine phosphokinase level:
      data from a randomized trial. Clin Infect Dis. 2010 Jun 15;50(12):
      1568-1574. doi:10.1086/652767.
  CIPROFLOXACINA
  [7] Schlender JF, Teutonico D, Coboeken K, Schnizler K, Eissing T,
      Willmann S, Jaehde U, Stass H. A Physiologically-Based Pharmacokinetic
      Model to Describe Ciprofloxacin Pharmacokinetics Over the Entire Span
      of Life. Clin Pharmacokinet. 2018 Dec;57(12):1613-1634.
      doi:10.1007/s40262-018-0661-6.
  [8] Balbas-Martinez V, Michelet R, Edginton AN, Meesters K, Trocóniz IF,
      Vermeulen A. Physiologically-Based Pharmacokinetic model for
      Ciprofloxacin in children with complicated Urinary Tract Infection. Eur
      J Pharm Sci. 2019 Feb 1;128:171-179. doi:10.1016/j.ejps.2018.11.033.
  [9] Forrest A, Ballow CH, Nix DE, Birmingham MC, Schentag JJ. Development
      of a population pharmacokinetic model and optimal sampling strategies
      for intravenous ciprofloxacin. Antimicrob Agents Chemother. 1993
      May;37(5):1065-1072. doi:10.1128/aac.37.5.1065.
  [10] Ambrose PG, Bhavnani SM, Rubino CM, et al. Pharmacokinetics-
       pharmacodynamics of antimicrobial therapy: it is not just for mice
       anymore. Clin Infect Dis. 2007;44(Suppl 2):S79-S86.
       doi:10.1086/510079.
  DALBAVANCINA
  [11] Modello Mariani: "Population Pharmacokinetics of Dalbavancin in
       Pediatric Patients: A Real-World Therapeutic Drug Monitoring
       Study" -- dati IRCCS Istituto Giannina Gaslini, Genova. Modello
       interno, non ancora pubblicato (fornito dall'utente).
  [12] Cojutti PG, et al. Population pharmacokinetics of dalbavancin and
       dosing consideration for optimal treatment of adult patients with
       staphylococcal osteoarticular infections. Antimicrob Agents
       Chemother. 2023.
  [13] Cojutti PG, et al. Real-time TDM-based Bayesian forecasting to
       optimize dosing of dalbavancin in patients with prolonged
       treatment. Int J Antimicrob Agents. 2021.
  [14] EUCAST breakpoint tables -- Staphylococcus aureus, dalbavancin.
       European Committee on Antimicrobial Susceptibility Testing.
  [15] DALVANCE (dalbavancin) prescribing information. FDA, 2025.
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

    # ============================================================
    # PANNELLO FORECAST — solo Dalbavancin
    # ============================================================
    if farmaco == "Dalbavancin":
        st.markdown("---")
        st.subheader("🔮 Forecast — Pianificazione Prossima Somministrazione")
        st.caption(
            "Proietta la curva individuale (CL stimato via MAP, V1/Q/V2 allometrici) "
            "in avanti nel tempo per stimare quando la concentrazione scende sotto "
            "la soglia terapeutica scelta, ed eventualmente valutare l'effetto di "
            "una nuova dose ipotetica."
        )

        colF1, colF2, colF3 = st.columns(3)
        with colF1:
            target_scelto = st.selectbox(
                "Target PK/PD di riferimento:",
                list(DALBA_TARGETS.keys()) + ["Trough (100% fT > 4×MIC)"],
                help="Cojutti PG, et al. Antimicrob Agents Chemother 2023 (complessa); "
                     "target ABSSSI semplice per infezioni non complicate"
            )
        with colF2:
            mic_scelta_label = st.selectbox("MIC / breakpoint:",
                                            list(DALBA_MIC_EUCAST.keys()) + ["Personalizzato"])
            if mic_scelta_label == "Personalizzato":
                mic_forecast = st.number_input("MIC personalizzata (mg/L)", 0.001, 4.0, 0.125,
                                               step=0.001, format="%.4f")
            else:
                mic_forecast = DALBA_MIC_EUCAST[mic_scelta_label]
        with colF3:
            f_unbound = st.number_input(
                "Frazione libera (f)", 0.01, 0.10, DALBA_PARAMS['fu'], 0.01,
                help="Default 0.07 (93% legame proteico). Range sensibilita' "
                     "riportato nel documento originale: 0.01-0.07 (93-99% legame)."
            )

        if target_scelto == "Trough (100% fT > 4×MIC)":
            soglia_target = dalba_soglia_trough(mic_forecast, f_unbound)
            desc_target = f"trough 4×MIC libera → soglia concentrazione totale = {soglia_target:.2f} mg/L"
        else:
            target_val = DALBA_TARGETS[target_scelto]
            soglia_target = dalba_soglia_concentrazione_totale(target_val, mic_forecast, f_unbound)
            desc_target = f"{target_scelto} → soglia concentrazione totale equivalente = {soglia_target:.2f} mg/L"

        st.info(f"📐 Soglia calcolata: **{desc_target}**")

        # Tempo dall'ultima osservazione (ultimo livello o ultima dose, il più recente)
        t_ultimo_liv = max((l['ore'] for l in livelli_completi), default=0)
        t_valutazione = max(t_ultimo_liv, storia_completa[-1]['ore'])
        c_attuale = pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia_completa, t_valutazione, t_inf)

        st.markdown(f"**Concentrazione attuale stimata** (t = {t_valutazione:.1f} h dall'inizio terapia): "
                   f"{c_attuale:.2f} mg/L")

        t_crossing = trova_tempo_soglia(CL_ind, V_ind, Q_ind, V2_ind, storia_completa,
                                        soglia_target, t_inf, t_partenza=t_valutazione,
                                        t_max=t_valutazione + 2000, passo=2.0)

        if t_crossing is None:
            st.warning("⚠️ La soglia non viene raggiunta entro ~83 giorni dalla valutazione attuale "
                      "con la sola dose già somministrata — verificare i parametri inseriti.")
        elif t_crossing <= t_valutazione + 0.5:
            st.error(f"🔴 La concentrazione è **già sotto la soglia** target al momento della valutazione "
                    f"({c_attuale:.2f} < {soglia_target:.2f} mg/L). Considerare una nuova somministrazione.")
        else:
            data_crossing = t0_assoluto + timedelta(hours=t_crossing)
            giorni_rimanenti = (t_crossing - t_valutazione) / 24.0
            st.success(
                f"🗓️ La concentrazione scenderà sotto la soglia target il "
                f"**{data_crossing.strftime('%d/%m/%Y alle %H:%M')}** "
                f"(tra circa **{giorni_rimanenti:.1f} giorni**)."
            )

        # --- Simulazione nuova dose ipotetica ---
        st.markdown("#### 💉 Simula una nuova somministrazione ipotetica")
        col_sim1, col_sim2, col_sim3 = st.columns(3)
        with col_sim1:
            simula = st.checkbox("Aggiungi dose ipotetica al forecast")
        with col_sim2:
            if simula:
                data_new = st.date_input("Data nuova dose ipotetica",
                                         value=(t0_assoluto + timedelta(hours=t_valutazione+24)).date(),
                                         format="DD/MM/YYYY", key="dalba_new_date")
                ora_new = st.time_input("Ora", value=datetime.strptime("10:00","%H:%M").time(),
                                        key="dalba_new_time")
        with col_sim3:
            if simula:
                dose_new = st.number_input("Dose ipotetica (mg)", 100.0, 1500.0, 1500.0, 50.0,
                                           key="dalba_new_dose")

        if simula:
            dt_new = datetime.combine(data_new, ora_new)
            t_new_ore = (dt_new - t0_assoluto).total_seconds() / 3600.0
            storia_sim = storia_completa + [{'ore': t_new_ore, 'dose': dose_new}]
            storia_sim.sort(key=lambda x: x['ore'])

            t_crossing_sim = trova_tempo_soglia(CL_ind, V_ind, Q_ind, V2_ind, storia_sim,
                                                soglia_target, t_inf,
                                                t_partenza=t_new_ore, t_max=t_new_ore+2000, passo=2.0)
            if t_crossing_sim:
                data_crossing_sim = t0_assoluto + timedelta(hours=t_crossing_sim)
                giorni_sim = (t_crossing_sim - t_new_ore) / 24.0
                st.success(
                    f"✅ Con la dose ipotetica di {dose_new:.0f} mg il "
                    f"{dt_new.strftime('%d/%m/%Y')}, la concentrazione resterà sopra soglia fino al "
                    f"**{data_crossing_sim.strftime('%d/%m/%Y')}** "
                    f"(**{giorni_sim:.1f} giorni** di copertura dalla nuova dose)."
                )

            # Grafico esteso con dose ipotetica
            t_fine_sim = t_crossing_sim + 48 if t_crossing_sim else t_new_ore + 500
            t_arr_sim = np.linspace(0, t_fine_sim, 1000)
            c_sim = [pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia_sim, t, t_inf) for t in t_arr_sim]

            fig2, ax2 = plt.subplots(figsize=(13, 4.2))
            fig2.patch.set_facecolor('#FAFAFA'); ax2.set_facecolor('#FAFAFA')
            ax2.plot(t_arr_sim, c_sim, color='#0F172A', lw=2.2, label='Profilo proiettato (con dose ipotetica)')
            ax2.axhline(soglia_target, color='#DC2626', ls='-.', lw=1.5,
                       label=f'Soglia target ({soglia_target:.2f} mg/L)')
            ax2.axvline(t_new_ore, color='#16A34A', ls='--', lw=1.5, alpha=0.7,
                       label=f'Nuova dose ({dose_new:.0f} mg)')
            for s in storia_completa:
                ax2.axvline(s['ore'], color='#94A3B8', ls=':', alpha=0.3)
            ax2.set_xlabel("Ore dall'inizio terapia (t=0 prima dose assoluta)", fontsize=10)
            ax2.set_ylabel("Dalbavancina (mg/L)", fontsize=10)
            ax2.set_xlim(left=0); ax2.set_ylim(bottom=0)
            ax2.legend(fontsize=8, loc='upper right')
            ax2.grid(alpha=0.15)
            plt.tight_layout()
            st.pyplot(fig2); plt.close()

    if salva and id_paz and id_paz != "":
        try:
            record = {
                'data_val': dt_val, 'farmaco': farmaco,
                'eta_giorni': gg, 'peso': peso,
                'creatinina': crea, 'crrt': bool(crrt), 'Sc': Sc, 'Qef': Qef,
                'dialisi_tipo': dialisi,
                'dose_mg': dose_mg, 'dose_mg_kg': dose_mg_kg,
                'tau_h': tau_h_display,   # solo informativo, nessun calcolo dipende da questo
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

            st.success("💾 Analisi aggiunta allo storico di **questa sessione**. "
                      "Nulla è stato salvato sul server: scarica il file aggiornato "
                      "qui sotto per conservarlo e ricaricarlo la prossima volta.")

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

st.set_page_config(page_title="TDM Buddy v5.2", layout="wide", page_icon="💊")

st.title("💊 TDM Buddy v5.2")
st.caption(
    "Monitoraggio Terapeutico Bayesiano — Daptomicina, Ciprofloxacina, Dalbavancin  |  "
    "MAP 1-cmt / 2-cmt / 2-cmt allometrico — nessun dato salvato sul server"
)

st.warning(
    "⚠️ **Strumento a scopo educational / di supporto decisionale.** Nessun dato "
    "viene conservato sul server in nessun momento: l'intera cronologia clinica "
    "esiste solo nel file JSON che scarichi e gestisci localmente. **Non inserire "
    "nomi, cognomi o date di nascita esatte** — usa un codice paziente locale "
    "non identificativo ed età espressa in anni/mesi/giorni. Questo strumento "
    "non sostituisce il giudizio clinico e non è un dispositivo medico certificato."
)

# ---- SIDEBAR: gestione storico paziente (JSON, nessun DB) ----
st.sidebar.header("🗂️ Storico Paziente (locale)")

file_up = st.sidebar.file_uploader(
    "📤 Carica storico paziente (JSON)", type=["json"],
    help="Carica il file scaricato in una sessione precedente per continuare "
         "la terapia di questo paziente con lo storico dosi già integrato."
)
if file_up is not None:
    try:
        dati_caricati = json.load(file_up)
        assert 'monitoraggi' in dati_caricati
        st.session_state['storico'] = dati_caricati
        st.sidebar.success(
            f"✅ Storico caricato: {len(dati_caricati.get('monitoraggi', []))} "
            f"analisi pregresse per il codice "
            f"'{dati_caricati.get('codice_paziente','?')}'."
        )
    except Exception as e:
        st.sidebar.error(f"File non valido: {e}")

id_paz = st.sidebar.text_input(
    "Codice paziente (non identificativo)",
    value=st.session_state['storico'].get('codice_paziente', ''),
    help="Es. 'PZ-001' o iniziali. Non inserire nome e cognome completi."
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
    json_bytes_sb = json.dumps(st.session_state['storico'], indent=2,
                               ensure_ascii=False).encode('utf-8')
    st.sidebar.download_button(
        "📥 Scarica storico corrente (JSON)",
        data=json_bytes_sb,
        file_name=f"TDMBuddy_storico_{id_paz or 'paziente'}.json",
        mime="application/json",
        use_container_width=True,
    )
    if st.sidebar.button("🗑️ Cancella storico di sessione"):
        st.session_state['storico'] = storico_vuoto()
        st.rerun()
else:
    st.sidebar.caption("Nessuna analisi ancora in questa sessione. "
                       "Carica un JSON precedente o inizia una nuova analisi.")

st.sidebar.markdown("---")
st.sidebar.markdown("### 📚 Modelli Integrati")
st.sidebar.markdown("""
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

**Dalbavancin (Mariani / Cojutti)**
- Modello Mariani: interno IRCCS Gaslini (non pubblicato, fornito dall'utente).
- Cojutti PG, et al. *Antimicrob Agents Chemother.* 2021;65(9):e00879-21 (modello osteoarticolare adulti).
- Cojutti PG, et al. *Antimicrob Agents Chemother.* 2023 (target complessa).
- Cojutti PG, et al. *Int J Antimicrob Agents.* 2021 (TDM Bayesiano).
- EUCAST breakpoints — dalbavancin / S. aureus.
- DALVANCE (dalbavancin) — FDA prescribing information, 2025.
""")

# ---- TABS ----
tab1, tab2 = st.tabs(["➕ Nuova Analisi TDM","📚 Storico Paziente"])

with tab1:
    if not id_paz:
        st.warning("Inserisci un codice paziente (non identificativo) nella barra laterale per procedere.")
        st.stop()

    st.markdown("### 👤 Dati del paziente per questa valutazione")
    st.caption("Età e peso possono cambiare nel tempo (specialmente in pediatria): "
              "verificali e aggiornali per OGNI valutazione TDM, anche se il "
              "paziente è lo stesso di una sessione precedente.")

    monitoraggi_tutti = st.session_state['storico'].get('monitoraggi', [])
    if monitoraggi_tutti:
        _ultimo = monitoraggi_tutti[-1]
        _gg_def = int(_ultimo.get('eta_giorni', 10950))
        _peso_def = float(_ultimo.get('peso', 70.0))
    else:
        _gg_def = 10950  # ~30 anni, default adulto
        _peso_def = 70.0

    _anni_def  = _gg_def // 365
    _resto_def = _gg_def % 365
    _mesi_def  = _resto_def // 30
    _giorni_def = _resto_def % 30

    col_e1, col_e2, col_e3, col_peso = st.columns(4)
    with col_e1:
        eta_anni_in = st.number_input("Età — Anni", 0, 100, min(int(_anni_def), 100))
    with col_e2:
        eta_mesi_in = st.number_input("Età — Mesi", 0, 11, min(int(_mesi_def), 11))
    with col_e3:
        eta_giorni_in = st.number_input("Età — Giorni", 0, 30, int(_giorni_def))
    with col_peso:
        peso = st.number_input("Peso (kg)", 0.3, 300.0, _peso_def)

    gg  = int(round(eta_anni_in*365.25 + eta_mesi_in*30.44 + eta_giorni_in))
    eta = gg / 365.25
    st.caption(f"Età calcolata: **{gg} giorni** ({eta:.2f} anni)")

    farmaco = st.radio("💊 Farmaco da monitorare:",
                        ["Daptomicina","Ciprofloxacina","Dalbavancin"],
                        horizontal=True)
    modello_dalba = None  # rilevante solo per farmaco == "Dalbavancin"

    if farmaco=="Daptomicina":
        if gg<365:     st.info(f"👶 Lattante ({gg} gg) → Modello: Ye et al., Eur J Pharm Sci 2025;208:107058")
        elif gg<6570:  st.info(f"🧒 Pediatrico ({eta:.1f}a) → Modello: Wei et al., J Clin Pharmacol 2020;60(6):768-774")
        else:          st.info(f"👨 Adulto ({eta:.0f}a) → Modello: Soraluce et al., Int J Antimicrob Agents 2018;52(2):158-165")
    elif farmaco=="Dalbavancin":
        modello_dalba = st.selectbox(
            "🧬 Modello PK di riferimento per il forecast:",
            list(DALBA_MODELLI.keys()),
            help="Mariani: pediatrico, scaling allometrico, dati Gaslini non pubblicati. "
                 "Cojutti 2021: adulti con infezione osteoarticolare, valori fissi (nessuno "
                 "scaling per peso disponibile nella pubblicazione originale)."
        )
        if modello_dalba.startswith('Mariani'):
            st.info(f"🧬 {modello_dalba} — scaling allometrico su peso {peso} kg. "
                   f"Solo CL ha variabilità interindividuale stimabile via MAP; "
                   f"V1/Q/V2 fissi (allometria).")
        else:
            st.info(f"🧬 {modello_dalba} — parametri di popolazione fissi "
                   f"(nessuno scaling per peso). Solo CL stimabile via MAP.")
    else:
        if gg<91:      st.warning("⚠️ Neonato <3 mesi: Schlender 2018 non validato in questa fascia d'età.")
        elif gg<6570:  st.info(f"🧒 Pediatrico ({eta:.1f}a) → Schlender JF, et al., Clin Pharmacokinet 2018;57(12):1613-1634 (scaling)")
        else:          st.info(f"👨 Adulto ({eta:.0f}a) → Schlender 2018 / Balbas-Martinez V, et al., Eur J Pharm Sci 2019;128:171-179")

    # ---- SEZIONE STORICO DOSI ----
    st.markdown("### 🔗 Integrazione Storico Dosi Paziente")

    storia_prev_db, livelli_prev_db, t0_db = carica_storia_precedente(farmaco)

    usa_storico = False
    if storia_prev_db or livelli_prev_db:
        col_info, col_check = st.columns([3,1])
        with col_info:
            st.success(
                f"📂 **Trovato storico in sessione:** {len(storia_prev_db)} dosi e "
                f"{len(livelli_prev_db)} livelli pregressi "
                f"(prima dose assoluta: {t0_db.strftime('%d/%m/%Y %H:%M')})."
            )
        with col_check:
            usa_storico = st.checkbox("Continua terapia precedente", value=True,
                                       help="Se attivato, le nuove dosi/livelli vengono aggiunti "
                                            "alla timeline esistente e il MAP usa tutta la storia "
                                            "cumulativa (quanti dosi e livelli il paziente abbia "
                                            "effettivamente avuto, senza limiti).")
        if usa_storico:
            cE1, cE2 = st.columns(2)
            with cE1:
                df_prev_d = pd.DataFrame([
                    {'Ore da t=0': f"{s['ore']:.1f}", 'Dose (mg)': s['dose']}
                    for s in storia_prev_db]) if storia_prev_db else pd.DataFrame()
                with st.expander(f"📋 {len(storia_prev_db)} dosi pregresse (sola lettura)", expanded=False):
                    st.dataframe(df_prev_d, use_container_width=True)
            with cE2:
                df_prev_l = pd.DataFrame([
                    {'Ore da t=0': f"{l['ore']:.1f}", 'Conc. (mg/L)': l['conc']}
                    for l in livelli_prev_db]) if livelli_prev_db else pd.DataFrame()
                with st.expander(f"🔬 {len(livelli_prev_db)} livelli pregressi (sola lettura)", expanded=False):
                    st.dataframe(df_prev_l, use_container_width=True)
    else:
        st.info("ℹ️ Nessuno storico disponibile in questa sessione per questo farmaco. "
                "Questa sarà la prima analisi TDM.")

    # ---- PARAMETRI CLINICI ----
    st.markdown("### 1. Parametri Clinici")
    col1,col2,col3 = st.columns(3)
    with col1:
        crea = st.number_input("Creatinina Sierica (mg/dL)", 0.05, 25.0, 1.0, 0.05)
        if farmaco == "Dalbavancin":
            mic = st.number_input("MIC di riferimento (mg/L)", 0.008, 4.0, 0.0625,
                                  step=0.001, format="%.4f",
                                  help="Il target PK/PD dettagliato con selezione MIC EUCAST "
                                       "è disponibile nella sezione Forecast dopo l'analisi. "
                                       "Accetta fino a 4 decimali (es. 0.0625, 0.125).")
        else:
            mic  = st.number_input(f"MIC {farmaco} (mg/L)", 0.001, 16.0, 1.0,
                                   step=0.001, format="%.4f",
                                   help="Accetta fino a 4 decimali (es. 0.125, 0.0625).")
    with col2:
        crrt = st.checkbox("Paziente in CRRT?")
        if crrt:
            Sc  = st.number_input("Sieving Coefficient (Sc)", 0.01, 0.60, 0.20, 0.01,
                                   help="Sc medio daptomicina = 0.20 ± 0.05 (Soraluce 2018)")
            Qef = st.number_input("Flusso Effluente Qef (L/h)", 0.5, 6.0, 2.0, 0.1)
        else:
            Sc=SORALUCE['Sc_default']; Qef=2.0
    with col3:
        if farmaco=="Daptomicina" and gg<6570:
            opt=["Nessuna","Insufficienza Renale Lieve (GFR 60-90)",
                 "Insufficienza Renale Moderata (GFR 30-60)",
                 "Insufficienza Renale Grave (GFR 15-30)","ESRD (GFR <15)"]
            if 4380<=gg<6570: opt+=["Emodialisi (HD)","Dialisi Peritoneale (CCPD)"]
            dialisi=st.selectbox("Funzione Renale / Dialisi:",opt)
        else:
            dialisi="Nessuna"
            st.caption("📌 Cockcroft-Gault calcolato internamente da creatinina, età, peso, sesso.")

    # ---- NUOVE SOMMINISTRAZIONI E LIVELLI (liste libere, nessun limite) ----
    st.markdown("### 2. Nuove Somministrazioni e Livelli Plasmatici")
    st.caption(
        "Inserisci ESATTAMENTE le somministrazioni e i prelievi effettuati dal "
        "paziente — quanti sono, con qualsiasi intervallo. Nessuno schema fisso "
        "(intervallo regolare, singolo Cmax/Cmin) viene imposto: usa quante righe "
        "servono, anche irregolari nel tempo."
    )

    dose_def = {"Daptomicina": 500.0, "Dalbavancin": 1500.0, "Ciprofloxacina": 400.0}[farmaco]
    if farmaco == "Dalbavancin":
        st.caption("📌 Dalbavancina: dosaggio tipico 22.5 mg/kg (<6a) o 18 mg/kg (≥6a), "
                  "max 1500 mg/dose. Aggiungi una riga per ogni somministrazione "
                  "realmente avvenuta, con la sua data/ora esatta.")

    if usa_storico and storia_prev_db:
        t_ultima = t0_db + timedelta(hours=storia_prev_db[-1]['ore'])
        data_def_dose = (t_ultima + timedelta(hours=24)).date()
        ora_def_dose  = (t_ultima + timedelta(hours=24)).time()
    else:
        data_def_dose = date.today()
        ora_def_dose  = datetime.strptime("08:00", "%H:%M").time()

    colA, colB = st.columns(2)

    with colA:
        st.markdown("**💉 Nuove Somministrazioni**")
        df_dosi_init = pd.DataFrame({
            "Data (GG/MM/AAAA)": [data_def_dose.strftime("%d/%m/%Y")],
            "Ora (HH:MM)": [ora_def_dose.strftime("%H:%M")],
            "Dose (mg)": [dose_def],
        })
        df_dosi_edit = st.data_editor(
            df_dosi_init, num_rows="dynamic", use_container_width=True,
            key="tbl_dosi_edit",
            column_config={
                "Dose (mg)": st.column_config.NumberColumn(format="%.1f", min_value=0.1),
            },
        )

    with colB:
        st.markdown("**🔬 Nuovi Livelli Plasmatici Misurati**")
        st.caption("Aggiungi (o lascia vuota) una riga per ogni prelievo. "
                  "Nessun limite al numero di livelli inseribili. La creatinina "
                  "per livello è solo di registrazione — non entra nei calcoli PK "
                  "(quella usata dai modelli resta il valore in \"Parametri Clinici\").")
        df_liv_init = pd.DataFrame({
            "Data (GG/MM/AAAA)": [""],
            "Ora (HH:MM)": [""],
            "Concentrazione (mg/L)": [None],
            "Creatinina (mg/dL, opz.)": [None],
            "Nota (opz.)": [""],
        })
        df_liv_edit = st.data_editor(
            df_liv_init, num_rows="dynamic", use_container_width=True,
            key="tbl_livelli_edit",
            column_config={
                "Concentrazione (mg/L)": st.column_config.NumberColumn(format="%.4f", min_value=0.0),
                "Creatinina (mg/dL, opz.)": st.column_config.NumberColumn(format="%.2f", min_value=0.0),
            },
        )

    # ---- PARSING TABELLE -> LISTE DI EVENTI ----
    def _parse_data_ora(data_str, ora_str):
        data_str = str(data_str).strip()
        ora_str  = str(ora_str).strip()
        if not data_str or not ora_str or data_str.lower() == 'nan' or ora_str.lower() == 'nan':
            return None
        for fmt_d in ("%d/%m/%Y", "%Y-%m-%d"):
            try:
                d = datetime.strptime(data_str, fmt_d).date()
                break
            except ValueError:
                d = None
        if d is None:
            return None
        try:
            t = datetime.strptime(ora_str[:5], "%H:%M").time()
        except ValueError:
            return None
        return datetime.combine(d, t)

    nuove_somm = []
    errori_dosi = []
    for _, r in df_dosi_edit.iterrows():
        dt_row = _parse_data_ora(r.get("Data (GG/MM/AAAA)", ""), r.get("Ora (HH:MM)", ""))
        dose_val = r.get("Dose (mg)", None)
        if dt_row is None or dose_val is None or pd.isna(dose_val):
            continue
        nuove_somm.append({'dt': dt_row, 'dose': float(dose_val)})

    nuovi_liv = []
    for _, r in df_liv_edit.iterrows():
        dt_row = _parse_data_ora(r.get("Data (GG/MM/AAAA)", ""), r.get("Ora (HH:MM)", ""))
        conc_val = r.get("Concentrazione (mg/L)", None)
        if dt_row is None or conc_val is None or pd.isna(conc_val):
            continue
        crea_val = r.get("Creatinina (mg/dL, opz.)", None)
        crea_val = float(crea_val) if (crea_val is not None and not pd.isna(crea_val)) else None
        nuovi_liv.append({'dt': dt_row, 'conc': float(conc_val),
                          'nota': str(r.get("Nota (opz.)", "") or ""),
                          'creatinina': crea_val})

    if not nuove_somm:
        st.warning("⚠️ Inserisci almeno una nuova somministrazione con data, ora e dose valide.")

    # ---- COSTRUZIONE TIMELINE ----
    if usa_storico and storia_prev_db and t0_db:
        t0_assoluto = t0_db
    elif nuove_somm:
        t0_assoluto = min(s['dt'] for s in nuove_somm)
    else:
        t0_assoluto = datetime.combine(date.today(), datetime.strptime("08:00","%H:%M").time())

    storia_corrente = costruisci_storia(
        nuove_somm, t0_assoluto,
        storia_prev_db if (usa_storico and storia_prev_db) else []
    )
    livelli_correnti = costruisci_livelli(
        nuovi_liv, t0_assoluto,
        livelli_prev_db if (usa_storico and livelli_prev_db) else []
    )

    dt_val = (max((l['dt'] for l in nuovi_liv), default=None) or
             max((s['dt'] for s in nuove_somm), default=datetime.now())
             ).strftime('%Y-%m-%d %H:%M')

    # Anteprima timeline
    with st.expander("👁️ Anteprima timeline completa prima di lanciare l'analisi", expanded=False):
        n_prev_dosi = len(storia_prev_db) if (usa_storico and storia_prev_db) else 0
        df_tl = pd.DataFrame([
            {'Ore da t=0': f"{s['ore']:.1f}h",
             'Data/Ora assoluta': (t0_assoluto + timedelta(hours=s['ore'])).strftime('%d/%m/%Y %H:%M'),
             'Dose (mg)': s['dose'],
             'Origine': '📂 Storico' if i < n_prev_dosi else '✏️ Nuova'}
            for i, s in enumerate(storia_corrente)])
        st.markdown(f"**Dosi ({len(storia_corrente)} totali)**")
        st.dataframe(df_tl, use_container_width=True)

        n_prev_liv = len(livelli_prev_db) if (usa_storico and livelli_prev_db) else 0
        df_tl_liv = pd.DataFrame([
            {'Ore da t=0': f"{l['ore']:.1f}h",
             'Data/Ora assoluta': (t0_assoluto + timedelta(hours=l['ore'])).strftime('%d/%m/%Y %H:%M'),
             'Conc. (mg/L)': f"{l['conc']:.4f}",
             'Creatinina (mg/dL)': (f"{l['creatinina']:.2f}" if l.get('creatinina') is not None else '—'),
             'Nota': l.get('nota',''),
             'Origine': '📂 Storico' if i < n_prev_liv else '✏️ Nuovo'}
            for i, l in enumerate(livelli_correnti)])
        st.markdown(f"**Livelli ({len(livelli_correnti)} totali)**")
        st.dataframe(df_tl_liv, use_container_width=True)
        st.caption(f"t=0 assoluto: **{t0_assoluto.strftime('%d/%m/%Y %H:%M')}**")

    if st.button("🚀 Avvia Analisi PK Bayesiana",type="primary",use_container_width=True):
        if not storia_corrente:
            st.error("Nessuna dose valida in timeline. Verifica la tabella delle somministrazioni.")
        elif len(livelli_correnti) == 0:
            st.error("Inserisci almeno un livello plasmatico (in una qualsiasi riga della tabella).")
        elif any(l['ore'] < 0 for l in livelli_correnti):
            st.error("⚠️ Un prelievo ha data/ora antecedente al t=0 della terapia. Verificare le date.")
        else:
            with st.spinner("Calcolo MAP Bayesiano in corso…"):
                esegui_analisi(
                    farmaco, id_paz, sesso, gg, peso, dt_val,
                    crea, crrt, Sc, Qef, dialisi,
                    mic, livelli_correnti,
                    storia_corrente,
                    storia_prev_db if (usa_storico and storia_prev_db) else [],
                    livelli_prev_db if (usa_storico and livelli_prev_db) else [],
                    t0_assoluto, salva=True,
                    modello_dalba=modello_dalba)

# ---- TAB STORICO ----
with tab2:
    monitoraggi = st.session_state['storico'].get('monitoraggi', [])
    if not monitoraggi:
        st.info("Nessuna analisi TDM in questa sessione. Carica un file JSON "
               "precedente dalla barra laterale, oppure esegui una prima "
               "analisi nella tab 'Nuova Analisi TDM'.")
    else:
        st.subheader(f"📋 Storico TDM di sessione — {id_paz}")
        st.caption("Questi dati esistono solo nella sessione corrente del browser. "
                  "Scarica il JSON aggiornato dalla barra laterale per non perderli.")

        df_sto = pd.DataFrame(monitoraggi)
        cols_vis = [c for c in ['data_val','farmaco','eta_giorni','peso','dose_mg',
                                 'n_dosi','n_livelli','mic','CL_ind','V_ind',
                                 'AUC24','AUC_MIC','Cmin_ss','modello']
                    if c in df_sto.columns]
        df_sto_show = df_sto[cols_vis].copy()
        df_sto_show.insert(0, 'idx', range(len(df_sto_show)))
        st.dataframe(df_sto_show, use_container_width=True)

        opz = [f"#{i} — {r['data_val']} [{r.get('farmaco','?')}]"
               for i, r in enumerate(monitoraggi)]
        scelta = st.selectbox("Rivisualizza analisi:", ["— Scegli —"] + opz)

        if scelta != "— Scegli —":
            idx_sel = int(scelta.split("#")[1].split(" ")[0])
            riga = monitoraggi[idx_sel]

            if st.button("🗑️ Rimuovi questo record dallo storico di sessione", type="primary"):
                st.session_state['storico']['monitoraggi'].pop(idx_sel)
                st.warning("Record rimosso dalla sessione. Ricorda di scaricare "
                          "il JSON aggiornato dalla barra laterale se vuoi "
                          "conservare questa modifica.")
                st.rerun()

            st.divider()
            storia_s  = riga['storia_dosi']
            livelli_s = riga.get('livelli', [])
            t0_s      = (datetime.strptime(riga['data_prima_dose'], '%Y-%m-%d %H:%M')
                        if riga.get('data_prima_dose') else datetime.now())
            gg_st   = int(riga.get('eta_giorni', gg))
            peso_st = float(riga.get('peso', peso))

            esegui_analisi(
                str(riga.get('farmaco', 'Daptomicina')),
                id_paz, sesso, gg_st, peso_st,
                riga['data_val'], float(riga['creatinina']),
                bool(riga['crrt']),
                float(riga.get('Sc', 0.20)), float(riga.get('Qef', 2.0)),
                str(riga.get('dialisi_tipo', 'Nessuna')), float(riga['mic']),
                livelli_s,
                storia_s, [], [], t0_s, salva=False,
                modello_dalba=riga.get('modello_dalba'))
