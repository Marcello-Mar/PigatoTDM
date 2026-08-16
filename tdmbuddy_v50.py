"""
TDM Buddy v6.0
==============

Piattaforma di monitoraggio terapeutico bayesiano per Daptomicina,
Ciprofloxacina e Dalbavancina, con storico paziente gestito localmente
(nessun dato salvato sul server — solo file JSON scaricabili/ricaricabili
dall'operatore). Ogni file JSON è scoped a UN paziente E UNA molecola,
per evitare di mescolare per errore lo storico di farmaci diversi.

RIFERIMENTI BIBLIOGRAFICI

DAPTOMICINA
  Soraluce A, Asín-Prieto E, Rodríguez-Gascón A, Barrasa H, Maynar J,
  Carcelero E, Soy D, Isla A. Population pharmacokinetics of daptomycin
  in critically ill patients. Int J Antimicrob Agents. 2018;52(2):
  158-165. doi:10.1016/j.ijantimicag.2018.03.008. PMID:29572042
  Wei XC, Zhao MF, Li X, Xiao X. Pharmacokinetic/Pharmacodynamic
  Analysis of Daptomycin Against S. aureus and E. faecium in Pediatric
  Patients by Monte Carlo Simulation. J Clin Pharmacol. 2020;60(6):
  768-774. doi:10.1002/jcph.1576
  Ye L, You X, Zhou J, Wu C, Ke M, Wu W, Huang P, Lin C. PBPK modeling
  of daptomycin dose optimization in pediatric patients with renal
  impairment. Front Pharmacol. 2022;13:838599. doi:10.3389/fphar.2022.838599
  Ye L, Zhou H, Guo G, Chen M, Zhang J. PBPK modeling to predict the
  exposure and assess pharmacodynamics of daptomycin in infants within
  1 year old. Eur J Pharm Sci. 2025;208:107058. doi:10.1016/j.ejps.2025.107058
  Lim SY, Lewis T, Woo S, Turman M, Bourne DWA, Burton ME, Rianthavorn P.
  Daptomycin Pharmacokinetics in Adolescents Undergoing Hemodialysis and
  Peritoneal Dialysis. J Pediatr Pharmacol Ther. 2021;26(2):123-132.
  doi:10.5863/1551-6776-26.2.123
  Bhavnani SM, Rubino CM, Ambrose PG, Drusano GL. Daptomycin exposure
  and the probability of elevations in the creatine phosphokinase
  level. Clin Infect Dis. 2010;50(12):1568-1574. doi:10.1086/652767
  [target tossicità: Cmin >= 24.3 mg/L]
  Louie A, Kaw P, Liu W, Jumbe N, Miller MH, Drusano GL. Pharmacodynamics
  of daptomycin in a murine-thigh model of Staphylococcus aureus
  infection. Antimicrob Agents Chemother. 2001;45:845-851.
  [target efficacia S. aureus/MRSA: AUC/MIC >= 666, modello murino]
  Canut A, Isla A, Betriu C, Gascón AR. Pharmacokinetic-pharmacodynamic
  evaluation of daptomycin, tigecycline, and linezolid versus vancomycin
  for the treatment of MRSA infections in four western European
  countries. Eur J Clin Microbiol Infect Dis. 2012;31(9):2227-2235.
  [adozione clinica del target 666 in adulti critici, citato in Soraluce 2018]
  Safdar N, Andes D, Craig WA. In vivo pharmacodynamic activity of
  daptomycin. Antimicrob Agents Chemother. 2004;48:63-68.
  [target efficacia E. faecium: AUC/MIC libero 12.9, convertito a
  AUC/MIC totale >= 143 in Wei 2020 assumendo legame proteico 91%]

CIPROFLOXACINA
  Schlender JF, Teutonico D, Coboeken K, Schnizler K, Eissing T,
  Willmann S, Jaehde U, Stass H. A PBPK Model to Describe Ciprofloxacin
  Pharmacokinetics Over the Entire Span of Life. Clin Pharmacokinet.
  2018;57(12):1613-1634. doi:10.1007/s40262-018-0661-6
  Balbas-Martinez V, Michelet R, Edginton AN, Meesters K, Trocóniz IF,
  Vermeulen A. PBPK model for Ciprofloxacin in children with
  complicated Urinary Tract Infection. Eur J Pharm Sci. 2019;128:
  171-179. doi:10.1016/j.ejps.2018.11.033
  Forrest A, Ballow CH, Nix DE, Birmingham MC, Schentag JJ. Development
  of a population pharmacokinetic model and optimal sampling strategies
  for intravenous ciprofloxacin. Antimicrob Agents Chemother. 1993;
  37(5):1065-1072. doi:10.1128/aac.37.5.1065 [modello PK]
  Forrest A, Nix DE, Ballow CH, Goss TF, Birmingham MC, Schentag JJ.
  Pharmacodynamics of intravenous ciprofloxacin in seriously ill
  patients. Antimicrob Agents Chemother. 1993;37(5):1073-1081.
  doi:10.1128/AAC.37.5.1073 [target efficacia AUC/MIC >= 125, 74
  adulti acutamente critici — companion paper del precedente]
  Nessuna soglia AUC quantitativa validata per la neurotossicità da
  fluorochinoloni: non riportata nell'app per questo motivo.

DALBAVANCINA
  Modello Mariani: dati IRCCS Istituto Giannina Gaslini, Genova.
  Modello interno non ancora pubblicato, fornito dall'utente.
  Cojutti PG, Rinaldi M, Zamparini E, Rossi N, Tedeschi S, Conti M,
  Pea F, Viale P. Population pharmacokinetics of dalbavancin and
  dosing consideration for optimal treatment of adult patients with
  staphylococcal osteoarticular infections. Antimicrob Agents
  Chemother. 2021;65(5):e02260-20. doi:10.1128/AAC.02260-20
  (target fAUC24h/MIC 27.1 / 53.3 / 111.1 e breakpoint EUCAST/USCAST
  entrambi citati in questo stesso lavoro)
  EUCAST Clinical Breakpoint Tables v.16.0, valide dal 2026-01-01.
  European Committee on Antimicrobial Susceptibility Testing.
  www.eucast.org (breakpoint dalbavancina/S. aureus: S ≤ 0.125 mg/L,
  valore stabile nelle tabelle EUCAST dal 2016; verificare l'ultima
  versione ufficiale alla data della prescrizione).

CORREZIONI v5.7.1 (Gemini) -> v6.0:
  - Rimosso il meccanismo di "Bayesian Updating" che rileggeva lo
    storico per spostare il prior: la stessa informazione era già
    conteggiata nella verosimiglianza (storico precompilato in
    tabella), quindi il prior veniva aggiornato due volte con lo
    stesso dato, producendo una falsa precisione.
  - "Dalbavancin" -> "Dalbavancina" (coerenza con italiano)
  - JSON ora scoped a paziente+molecola (non più multi-farmaco in un
    unico file)
  - Interfaccia riorganizzata a step guidati con selezione farmaco
    prominente
  - Grafici Plotly: risolto bug testo bianco su sfondo chiaro
    (Streamlit theme override non disattivato)
  - Target PK/PD dalbavancina rietichettati con fonte precisa
    (Cojutti 2021, citazione corretta: 65(5):e02260-20, non 65(9))
  - MIC EUCAST/USCAST rietichettate correttamente (0.125 = EUCAST
    attuale, 0.25 = USCAST, non "EUCAST precedente")
"""

import streamlit as st
import numpy as np
from scipy.optimize import minimize
import plotly.graph_objects as go
import pandas as pd
from datetime import datetime, date, timedelta
import json

# ============================================================
# COSTANTI GENERALI
# ============================================================

FARMACI = ["Daptomicina", "Ciprofloxacina", "Dalbavancina"]

FARMACO_INFO = {
    "Daptomicina": {
        "icona": "🩸",
        "descrizione": "Lipopeptide — Soraluce, Wei, Ye, Lim",
    },
    "Ciprofloxacina": {
        "icona": "🧬",
        "descrizione": "Fluorochinolone — Schlender, Balbas-Martinez",
    },
    "Dalbavancina": {
        "icona": "💉",
        "descrizione": "Lipoglicopeptide — Mariani, Cojutti",
    },
}

# Riferimenti bibliografici completi (modello PK, target di efficacia,
# target di tossicità/sicurezza, con il setting/popolazione in cui
# ciascuno è stato validato), mostrati in piccolo in cima al Passo 2.
RIFERIMENTI_PER_FASCIA = {
    "Daptomicina": [
        ("PK — < 1 anno", "Ye L, Zhou H, et al. *Eur J Pharm Sci.* 2025;208:107058"),
        ("PK — 1–17 anni", "Wei XC, et al. *J Clin Pharmacol.* 2020;60(6):768-774 "
                           "(correzione insufficienza renale: Ye L, et al. *Front Pharmacol.* "
                           "2022;13:838599)"),
        ("PK — Adolescenti in HD/CCPD", "Lim SY, et al. *J Pediatr Pharmacol Ther.* 2021;26(2):123-132"),
        ("PK — ≥ 18 anni", "Soraluce A, et al. *Int J Antimicrob Agents.* 2018;52(2):158-165 "
                           "(adulti in terapia intensiva)"),
        ("Target efficacia — S. aureus/MRSA, AUC/MIC ≥ 666",
         "Louie A, et al. *Antimicrob Agents Chemother.* 2001;45:845-851 (modello murino, "
         "thigh infection); adottato clinicamente in adulti critici da Canut A, et al. "
         "*Eur J Clin Microbiol Infect Dis.* 2012;31(9):2227-2235"),
        ("Target efficacia — E. faecium, AUC/MIC ≥ 143",
         "Safdar N, Andes D, Craig WA. *Antimicrob Agents Chemother.* 2004;48:63-68 "
         "(modello murino); conversione a concentrazione totale riportata in "
         "Wei XC, et al. *J Clin Pharmacol.* 2020;60(6):768-774"),
        ("Target tossicità — Cmin ≥ 24.3 mg/L",
         "Bhavnani SM, et al. *Clin Infect Dis.* 2010;50(12):1568-1574 (trial randomizzato, "
         "adulti con batteriemia ed endocardite)"),
    ],
    "Ciprofloxacina": [
        ("PK — tutte le età (scaling pediatrico incluso)",
         "Schlender JF, et al. *Clin Pharmacokinet.* 2018;57(12):1613-1634; "
         "Balbas-Martinez V, et al. *Eur J Pharm Sci.* 2019;128:171-179"),
        ("Variabilità individuale (IIV)", "Forrest A, et al. *Antimicrob Agents Chemother.* 1993;37(5):1065-1072"),
        ("Target efficacia — gram-negativi, AUC/MIC ≥ 125",
         "Forrest A, Nix DE, Ballow CH, et al. *Antimicrob Agents Chemother.* 1993;37(5):1073-1081 "
         "(74 adulti acutamente critici, ciprofloxacina EV 200mg q12h–400mg q8h)"),
        ("Target tossicità", "Nessuna soglia AUC quantitativa validata per la neurotossicità "
                             "da fluorochinoloni: il rischio è documentato soprattutto in "
                             "insufficienza renale non corretta (case report), non riportiamo "
                             "quindi un valore soglia numerico."),
    ],
    "Dalbavancina": [
        ("PK — pediatrico, scaling allometrico (modello Mariani)",
         "Dati interni IRCCS Gaslini, non ancora pubblicato"),
        ("PK — adulti, infezione osteoarticolare (modello Cojutti)",
         "Cojutti PG, et al. *Antimicrob Agents Chemother.* 2021;65(5):e02260-20"),
        ("Target efficacia — fAUC24h/MIC 27.1 / 53.3 / 111.1",
         "Cojutti PG, et al. *Antimicrob Agents Chemother.* 2021;65(5):e02260-20 "
         "(Monte Carlo simulation, adulti con infezione osteoarticolare, n=15)"),
        ("Breakpoint MIC", "EUCAST Breakpoint Tables v.16.0, valide dal 2026-01-01 (eucast.org)"),
        ("Target tossicità", "Nessuna soglia di concentrazione validata per la tossicità da "
                             "dalbavancina; il monitoraggio di sicurezza segue la scheda "
                             "tecnica DALVANCE."),
    ],
}

# ============================================================
# STORAGE — nessun database, solo session_state + JSON scaricabile.
# Ogni sessione/file è scoped a UNA molecola.
# ============================================================

def storico_vuoto(farmaco=None):
    return {'codice_paziente': '', 'sesso': 'Maschio', 'farmaco': farmaco, 'monitoraggi': []}

# Chiavi dei widget legati al singolo paziente/valutazione. Streamlit
# mantiene lo stato di un widget associato alla sua key anche dopo che
# session_state a monte è stato resettato — il parametro "value=" viene
# ignorato una volta che quella key esiste già. Per evitare che dati di
# un paziente (o di una molecola precedente) restino "a galla" nei
# widget dopo un reset, queste chiavi vanno cancellate esplicitamente.
CHIAVI_WIDGET_PAZIENTE = [
    'input_codice_paziente', 'input_sesso', 'editor_dosi', 'editor_livelli',
    'input_unita_eta', 'input_valore_eta', 'input_peso',
    'input_modello_dalba', 'input_crea', 'input_mic', 'input_crrt',
    'input_sc', 'input_qef', 'input_dialisi',
    'fc_target', 'fc_mic_label', 'fc_mic_custom', 'fc_funb', 'fc_simula',
    'fc_data_new', 'fc_ora_new', 'fc_dose_new',
]

def pulisci_widget_paziente():
    for k in CHIAVI_WIDGET_PAZIENTE:
        st.session_state.pop(k, None)

if 'storico' not in st.session_state:
    st.session_state['storico'] = storico_vuoto()
if 'farmaco_attivo' not in st.session_state:
    st.session_state['farmaco_attivo'] = None

# ============================================================
# PARAMETRI PK — DAPTOMICINA
# ============================================================

SORALUCE = {
    'CLNR': 0.160, 'CLR_theta': 0.366, 'Clcr_ref': 49.0, 'V': 12.30,
    'omega_CL': 0.367, 'omega_V': 0.278, 'sigma': 0.123, 'Sc_default': 0.20,
    'nome': 'Soraluce A, et al. Int J Antimicrob Agents 2018;52(2):158-165 [Adulto ICU]',
}
YE2025 = {'CL': 19.6, 'SD': 5.5, 'V': 0.17,
          'nome': 'Ye L, Zhou H, et al. Eur J Pharm Sci 2025;208:107058 [<1 anno]'}
WEI = {
    '3-6m'  : {'CL':19.7,'SD':5.5, 'V':0.17,'nome':'Wei XC, et al. J Clin Pharmacol 2020;60(6):768-774 [3-6 mesi]'},
    '7-12m' : {'CL':19.6,'SD':5.8, 'V':0.17,'nome':'Wei XC, et al. J Clin Pharmacol 2020;60(6):768-774 [7-12 mesi]'},
    '13-24m': {'CL':21.8,'SD':3.0, 'V':0.11,'nome':'Wei XC, et al. J Clin Pharmacol 2020;60(6):768-774 [13-24 mesi]'},
    '2-6y'  : {'CL':20.0,'SD':6.1, 'V':0.11,'nome':'Wei XC, et al. J Clin Pharmacol 2020;60(6):768-774 [2-6 anni]'},
    '7-11y' : {'CL':17.0,'SD':5.7, 'V':0.11,'nome':'Wei XC, et al. J Clin Pharmacol 2020;60(6):768-774 [7-11 anni]'},
    '12-17y': {'CL':11.0,'SD':2.1, 'V':0.11,'nome':'Wei XC, et al. J Clin Pharmacol 2020;60(6):768-774 [12-17 anni]'},
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
    'nome': 'Lim SY, et al. J Pediatr Pharmacol Ther 2021;26(2):123-132 [Adolescente HD/CCPD]',
}
PD_DAPTO = {'AUC_MIC_staph': 666, 'AUC_MIC_entero': 143, 'Cmin_tox': 24.3}

# ============================================================
# PARAMETRI PK — CIPROFLOXACINA
# ============================================================

CIPRO = {
    'CL_kg':0.50, 'V1_kg':0.55, 'Q_kg':0.30, 'V2_kg':1.90,
    'omega_CL':0.37, 'omega_V1':0.32, 'omega_Q':0.25, 'omega_V2':0.30,
    'sigma':0.15, 'fu':0.67,
    'nome': 'Schlender JF, et al. Clin Pharmacokinet 2018;57(12):1613-1634 / '
            'Balbas-Martinez V, et al. Eur J Pharm Sci 2019;128:171-179',
}
CIPRO_PED = {(91,365):1.40, (365,730):1.30, (730,2190):1.20, (2190,4380):1.10, (4380,6570):1.00}
CIPRO_FRAZ_RENAL = 0.70
PD_CIPRO = {'AUC_MIC_target': 125}
IRC_MAP = {"Lieve (GFR 60-90)":'mild', "Moderata (GFR 30-60)":'moderate',
           "Grave (GFR 15-30)":'severe', "ESRD (GFR <15)":'esrd'}

# ============================================================
# PARAMETRI PK — DALBAVANCINA
# ============================================================

DALBA_PARAMS = {
    'CL_typ70': 0.0305, 'V1_typ70': 4.81, 'Q_typ70': 0.0081, 'V2_typ70': 4.81,
    'omega_CL': 0.198, 'sigma': 0.375, 't_inf': 0.5, 'fu': 0.07,
    'nome': "Modello Mariani (IRCCS Gaslini, dati non pubblicati, scaling allometrico)",
}
COJUTTI_PARAMS = {
    'CL_typ70': 0.106, 'V1_typ70': 17.4, 'Q_typ70': 0.103, 'V2_typ70': 15.1,
    'omega_CL': 0.362, 'sigma': 0.18, 't_inf': 0.5, 'fu': 0.07,
    'nome': "Cojutti PG, et al. Antimicrob Agents Chemother 2021;65(5):e02260-20 "
            "(adulti, osteoarticolare, n=15, valori di popolazione fissi)",
}
DALBA_MODELLI = {
    'Mariani — pediatrico, allometrico (Gaslini)': DALBA_PARAMS,
    'Cojutti 2021 — adulti, osteoarticolare (AAC)': COJUTTI_PARAMS,
}

# MIC di riferimento: EUCAST (0.125) e USCAST (0.25) sono DUE ENTI
# DIVERSI con soglie diverse — 0.25 non è una versione "precedente"
# di EUCAST. Fonte per entrambe: Cojutti PG, et al. AAC 2021;65(5):
# e02260-20, che le cita esplicitamente insieme.
DALBA_MIC_OPZIONI = {
    'MIC90 osteoarticolare — Cojutti 2021 (0.0625 mg/L)': 0.0625,
    'EUCAST S. aureus, breakpoint attuale (≤0.125 mg/L)': 0.125,
    'USCAST S. aureus (≤0.25 mg/L)': 0.25,
}

# Target fAUC24h/MIC — TUTTI E TRE dallo stesso Monte Carlo simulation
# di Cojutti PG, et al. AAC 2021;65(5):e02260-20. Etichette brevi per
# la tendina (il dettaglio completo è mostrato sotto come didascalia).
DALBA_TARGETS = {
    'Target minimo (≥ 27.1)': 27.1,
    'Target intermedio (≥ 53.3)': 53.3,
    'Target restrittivo (≥ 111.1)': 111.1,
}
DALBA_TARGETS_DESC = {
    'Target minimo (≥ 27.1)': 'infezioni non complicate',
    'Target intermedio (≥ 53.3)': 'infezioni moderate',
    'Target restrittivo (≥ 111.1)': 'infezioni complesse / osteoarticolari',
}

def dalba_soglia_totale(target_auc_mic, mic, f):
    """Soglia di concentrazione totale equivalente al target fAUC24h/MIC scelto."""
    return (target_auc_mic * mic) / (24.0 * f)

def dalba_soglia_trough_c0mic(mic):
    """
    Target di trough C0/MIC >= 200 (concentrazione TOTALE, rapporto
    diretto, non richiede la frazione libera). Riportato in letteratura
    TDM sulla dalbavancina come equivalente a C0 ~ 25 mg/L al
    breakpoint EUCAST di 0.125 mg/L. Verificare la fonte primaria
    prima dell'uso clinico definitivo.
    """
    return 200.0 * mic

# ============================================================
# MOTORI PK
# ============================================================

def pk_1cmt(CL, V, storia, t_obs, t_inf=0.5):
    if CL <= 0 or V <= 0:
        return 0.0
    ke = CL / V
    C = 0.0
    for s in storia:
        dt = t_obs - s['ore']
        if dt <= 0:
            continue
        D = s['dose']
        if dt <= t_inf:
            C += (D/(V*ke*t_inf))*(1.0-np.exp(-ke*dt))
        else:
            C += (D/(V*ke*t_inf))*(1.0-np.exp(-ke*t_inf))*np.exp(-ke*(dt-t_inf))
    return max(C, 0.0)


def pk_2cmt(CL, V1, Q, V2, storia, t_obs, t_inf=0.5):
    if CL<=0 or V1<=0 or Q<=0 or V2<=0:
        return 0.0
    k10=CL/V1; k12=Q/V1; k21=Q/V2
    s_=k10+k12+k21; d_=max(s_**2-4*k10*k21,0)
    alp=(s_+np.sqrt(d_))/2; bet=(s_-np.sqrt(d_))/2
    if abs(alp-bet)<1e-10:
        return pk_1cmt(CL,V1+V2,storia,t_obs,t_inf)
    A=(k21-alp)/((bet-alp)*V1); B=(k21-bet)/((alp-bet)*V1); C=0.0
    for s in storia:
        dt=t_obs-s['ore']
        if dt <= 0:
            continue
        D=s['dose']; k0=D/t_inf
        if dt<=t_inf:
            C += k0 * ((A/alp)*(1-np.exp(-alp*dt)) + (B/bet)*(1-np.exp(-bet*dt)))
        else:
            C += k0 * ((A/alp)*(1-np.exp(-alp*t_inf))*np.exp(-alp*(dt-t_inf)) +
                      (B/bet)*(1-np.exp(-bet*t_inf))*np.exp(-bet*(dt-t_inf)))
    return max(C,0.0)


def metriche_esposizione_24h(pk_func, params, storia, t_inf, n_punti=300):
    """AUC, Cmax, Cmin numerici sulle 24h dopo l'ultima dose (nessuna
    formula chiusa a 'tau' fisso, valida sia per dosaggio regolare che
    irregolare)."""
    t_ultima_dose = storia[-1]['ore']
    t_grid = np.linspace(t_ultima_dose, t_ultima_dose + 24, n_punti)
    c_grid = np.array([pk_func(*params, storia, t, t_inf) for t in t_grid])
    auc = float(np.sum((c_grid[1:] + c_grid[:-1]) / 2.0 * np.diff(t_grid)))
    return auc, float(c_grid.max()), float(c_grid.min())


def tau_mediano(storia):
    if len(storia) < 2:
        return None
    ore = sorted(s['ore'] for s in storia)
    diffs = sorted(ore[i+1]-ore[i] for i in range(len(ore)-1))
    n = len(diffs)
    return diffs[n//2] if n % 2 == 1 else (diffs[n//2-1]+diffs[n//2])/2.0


def calcola_clcr(eta_a, peso, crea, sesso):
    if crea <= 0:
        return 60.0
    clcr = ((140-eta_a)*peso)/(72*crea)
    return max(clcr*0.85 if sesso == "Femmina" else clcr, 1.0)

# ============================================================
# PRIOR DI POPOLAZIONE
# ============================================================

def prior_dapto(sesso, gg, peso, crea, crrt, Sc, Qef, dialisi):
    eta = gg/365.25
    note = ""
    if gg < 365:
        p = YE2025
        CL_pop = (p['CL']*peso)/1000; V_pop = p['V']*peso
        oCL = p['SD']/p['CL']; oV = 0.30; sig = 0.12; t_inf = 1.0; nm = p['nome']
        if crrt:
            CLEC = Sc*Qef; CL_pop += CLEC
            note += f"CRRT: CLEC={CLEC:.3f} L/h. "
    elif 4380 <= gg < 6570 and dialisi in ("Emodialisi (HD)","Dialisi Peritoneale (CCPD)"):
        sc = peso/LIM['peso_ref']
        CL_pop = (LIM['CLren']+LIM['CLnr'])*sc; V_pop = LIM['Vss']*sc; nm = LIM['nome']
        oCL, oV, sig, t_inf = LIM['omega_knr'], LIM['omega_V1'], LIM['sigma'], 0.5
        if dialisi == "Emodialisi (HD)":
            note += "⚠️ Somministrare DOPO la sessione di emodialisi per minimizzare le perdite. "
    elif 91 <= gg < 6570:
        if   gg < 213:  ag='3-6m'
        elif gg < 365:  ag='7-12m'
        elif gg < 730:  ag='13-24m'
        elif gg < 2190: ag='2-6y'
        elif gg < 4380: ag='7-11y'
        else:           ag='12-17y'
        p = WEI[ag]
        CL_pop=(p['CL']*peso)/1000; V_pop=p['V']*peso
        oCL=p['SD']/p['CL']; oV=0.25; sig=0.12
        t_inf = 1.0 if gg<730 else 0.5; nm = p['nome']
        age_ye = ag if ag in YE2022 else ('2-6y' if ag in ('13-24m','2-6y') else '1-2y')
        if dialisi in IRC_MAP and age_ye in YE2022:
            fattore = YE2022[age_ye][IRC_MAP[dialisi]]
            CL_pop /= fattore
            note += f"Correzione IRC (Ye L, et al. Front Pharmacol 2022;13:838599): CL×1/{fattore:.2f}. "
        if crrt:
            CLEC = Sc*Qef; CL_pop += CLEC
            note += f"CRRT: CLEC={CLEC:.3f} L/h. "
    else:
        s = SORALUCE
        clcr = calcola_clcr(eta, peso, crea, sesso)
        CL_pop = s['CLNR'] + s['CLR_theta']*(clcr/s['Clcr_ref']); V_pop = s['V']
        nm = f"{s['nome']} (ClCr={clcr:.0f} mL/min)"
        oCL, oV, sig, t_inf = s['omega_CL'], s['omega_V'], s['sigma'], 0.5
        if crrt:
            CLEC = Sc*Qef; CL_pop += CLEC
            note += f"CRRT: CLEC={CLEC:.3f} L/h. "

    return {'CL_p':CL_pop, 'V_p':V_pop, 'Q_p':0.0, 'V2_p':0.0,
            'oCL':oCL, 'oV1':oV, 'oQ':0.0, 'oV2':0.0, 'sig':sig, 't_inf':t_inf,
            'nm':nm, 'note':note}


def prior_cipro(sesso, gg, peso, crea, crrt, dialisi):
    eta = gg/365.25; note=""; c=CIPRO; fc=1.0
    for (lo,hi),f in CIPRO_PED.items():
        if lo <= gg < hi:
            fc = f
            note += f"Scaling pediatrico (Schlender 2018): CL×{fc:.2f}. "
            break

    clcr = calcola_clcr(eta, peso, crea, sesso)
    if clcr >= 30:
        CL_pop = (c['CL_kg']*peso*fc*CIPRO_FRAZ_RENAL*min(clcr/120.0,1.0)) + \
                 (c['CL_kg']*peso*fc*(1-CIPRO_FRAZ_RENAL))
        if clcr < 90:
            note += f"IRC: ClCr={clcr:.0f} mL/min → CL ridotta. "
    else:
        CL_pop = (c['CL_kg']*peso*fc*CIPRO_FRAZ_RENAL*0.20) + \
                 (c['CL_kg']*peso*fc*(1-CIPRO_FRAZ_RENAL)*0.63)
        note += f"⚠️ IRC grave (ClCr={clcr:.0f} mL/min): CL significativamente ridotta. "

    V1_pop = c['V1_kg']*peso; Q_pop = c['Q_kg']*peso; V2_pop = c['V2_kg']*peso
    nome = f"{c['nome']} (ClCr={clcr:.0f} mL/min)"

    return {'CL_p':CL_pop, 'V_p':V1_pop, 'Q_p':Q_pop, 'V2_p':V2_pop,
            'oCL':c['omega_CL'], 'oV1':c['omega_V1'], 'oQ':c['omega_Q'], 'oV2':c['omega_V2'],
            'sig':c['sigma'], 't_inf':1.0, 'nm':nome, 'note':note}


def prior_dalba(peso, modello_scelto):
    p = DALBA_MODELLI[modello_scelto]
    if modello_scelto.startswith('Mariani'):
        fw = peso/70.0
        CL_pop = p['CL_typ70']*(fw**0.75); V1_pop = p['V1_typ70']*fw
        Q_pop  = p['Q_typ70']*(fw**0.75);  V2_pop = p['V2_typ70']*fw
        note = ("Scaling allometrico applicato (CL,Q: peso/70^0.75; V1,V2: peso/70). "
                "Nessuna covariata renale/età ritenuta significativa nel modello originale.")
    else:
        CL_pop, V1_pop, Q_pop, V2_pop = p['CL_typ70'], p['V1_typ70'], p['Q_typ70'], p['V2_typ70']
        note = ("⚠️ Nessuno scaling per peso disponibile nella pubblicazione originale: "
                "valori di popolazione fissi, derivati da 15 adulti con infezione "
                "osteoarticolare. Usare con cautela per pazienti molto diversi da questa coorte.")

    return {'CL_p':CL_pop, 'V_p':V1_pop, 'Q_p':Q_pop, 'V2_p':V2_pop,
            'oCL':p['omega_CL'], 'oV1':0.0, 'oQ':0.0, 'oV2':0.0,
            'sig':p['sigma'], 't_inf':p['t_inf'], 'nm':p['nome'], 'note':note}

# ============================================================
# MAP BAYESIANO (N livelli arbitrari)
# ============================================================

def map_1cmt(CL_p, V_p, oCL, oV, sig, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0:
        return CL_p, V_p, "Prior di popolazione (nessun livello inserito)"
    def f(p):
        CL_, V_ = np.exp(p[0]), np.exp(p[1])
        pr = (((p[0]-np.log(CL_p))**2)/oCL**2 + ((p[1]-np.log(V_p))**2)/oV**2)
        lk = sum(((np.log(max(l['conc'],1e-6)) -
                  np.log(max(pk_1cmt(CL_,V_,storia,l['ore'],t_inf),1e-6)))**2)/sig**2
                for l in livelli)
        return 0.5*(pr+lk)
    bnd = [(np.log(CL_p*0.05), np.log(CL_p*8)), (np.log(V_p*0.15), np.log(V_p*5))]
    r = minimize(f, [np.log(CL_p),np.log(V_p)], method='L-BFGS-B', bounds=bnd,
                options={'maxiter':2000,'ftol':1e-12})
    return float(np.exp(r.x[0])), float(np.exp(r.x[1])), f"MAP Bayesiano 1-cmt su {n} livell{'o' if n==1 else 'i'}"


def map_2cmt(CL_p, V1_p, Q_p, V2_p, oCL, oV1, oQ, oV2, sig, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0:
        return CL_p, V1_p, Q_p, V2_p, "Prior di popolazione (nessun livello inserito)"
    OQ_eff  = 0.10 + (oQ  - 0.10)*min(1.0, n/6.0)
    OV2_eff = 0.10 + (oV2 - 0.10)*min(1.0, n/6.0)
    def f(p):
        CL_,V1_,Q_,V2_ = np.exp(p[0]),np.exp(p[1]),np.exp(p[2]),np.exp(p[3])
        pr = (((p[0]-np.log(CL_p))**2)/oCL**2 + ((p[1]-np.log(V1_p))**2)/oV1**2 +
             ((p[2]-np.log(Q_p))**2)/OQ_eff**2 + ((p[3]-np.log(V2_p))**2)/OV2_eff**2)
        lk = sum(((np.log(max(l['conc'],1e-6)) -
                  np.log(max(pk_2cmt(CL_,V1_,Q_,V2_,storia,l['ore'],t_inf),1e-6)))**2)/sig**2
                for l in livelli)
        return 0.5*(pr+lk)
    bnd = [(np.log(CL_p*0.05),np.log(CL_p*10)), (np.log(V1_p*0.1),np.log(V1_p*6)),
           (np.log(Q_p*0.3),np.log(Q_p*3)), (np.log(V2_p*0.3),np.log(V2_p*3))]
    r = minimize(f, [np.log(CL_p),np.log(V1_p),np.log(Q_p),np.log(V2_p)],
                method='L-BFGS-B', bounds=bnd, options={'maxiter':2000,'ftol':1e-12})
    nota_id = ("Q e V2 ancorati stretti al prior (pochi campioni)" if n<=2
              else f"Q e V2 con vincolo allentato ({n} campioni)")
    return (float(np.exp(r.x[0])),float(np.exp(r.x[1])),float(np.exp(r.x[2])),float(np.exp(r.x[3])),
            f"MAP Bayesiano 2-cmt su {n} livelli — {nota_id}")


def map_dalba(CL_p, V1_p, Q_p, V2_p, oCL, sig, storia, livelli, t_inf):
    n = len(livelli)
    if n == 0:
        return CL_p, "Prior di popolazione (nessun livello inserito)"
    def f(p):
        CL_ = np.exp(p[0])
        prior = ((p[0]-np.log(CL_p))**2)/oCL**2
        lik = sum(((np.log(max(l['conc'],1e-6)) -
                   np.log(max(pk_2cmt(CL_,V1_p,Q_p,V2_p,storia,l['ore'],t_inf),1e-6)))**2)/sig**2
                 for l in livelli)
        return 0.5*(prior+lik)
    r = minimize(f, [np.log(CL_p)], method='L-BFGS-B',
                bounds=[(np.log(CL_p*0.05),np.log(CL_p*10.0))], options={'maxiter':2000,'ftol':1e-12})
    return float(np.exp(r.x[0])), f"MAP Bayesiano su CL, {n} livell{'o' if n==1 else 'i'} (V1/Q/V2 fissi)"


def trova_tempo_soglia(CL, V1, Q, V2, storia, soglia, t_inf, t_partenza=0.0, t_max=2000.0, passo=1.0):
    if pk_2cmt(CL,V1,Q,V2,storia,t_partenza,t_inf) < soglia:
        return t_partenza
    t_prev = t_partenza; t = t_partenza + passo
    while t <= t_max:
        if pk_2cmt(CL,V1,Q,V2,storia,t,t_inf) < soglia:
            lo, hi = t_prev, t
            for _ in range(60):
                mid = (lo+hi)/2.0
                if pk_2cmt(CL,V1,Q,V2,storia,mid,t_inf) >= soglia:
                    lo = mid
                else:
                    hi = mid
            return hi
        t_prev = t; t += passo
    return None

# ============================================================
# GRAFICI PLOTLY (sfondo chiaro, testo scuro — sempre leggibili
# indipendentemente dal tema del browser/Streamlit dell'operatore)
# ============================================================

_PLOTLY_LAYOUT_BASE = dict(
    plot_bgcolor='#FFFFFF', paper_bgcolor='#FFFFFF',
    font=dict(color='#1F2937', size=13),
    xaxis=dict(color='#1F2937', gridcolor='#E5E7EB', linecolor='#9CA3AF'),
    yaxis=dict(color='#1F2937', gridcolor='#E5E7EB', linecolor='#9CA3AF'),
    hovermode="x unified",
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1,
               font=dict(color='#1F2937')),
    margin=dict(l=40, r=40, t=40, b=40),
)


def disegna_grafico(farmaco, storia, livelli, CL_ind, V_ind, Q_ind, V2_ind,
                    CL_std, V_std, Q_std, V2_std, t_inf, t_ultimo_storico=-1):
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
    fig.add_trace(go.Scatter(x=t_arr, y=C_std, mode='lines',
                             line=dict(dash='dash', color='#94A3B8', width=2),
                             name='Popolazione standard'))
    fig.add_trace(go.Scatter(x=t_arr, y=C_ind, mode='lines',
                             line=dict(color='#1E3A5F', width=3),
                             name='Paziente (MAP Bayesiano)'))

    max_c = max(max(C_ind, default=100), max(C_std, default=100))

    if t_ultimo_storico >= 0:
        fig.add_vrect(x0=0, x1=t_ultimo_storico, fillcolor="#6EA8D8", opacity=0.10,
                     layer="below", line_width=0)

    x_vecchie, y_vecchie, x_nuove, y_nuove = [], [], [], []
    for s in storia:
        if s['ore'] <= t_ultimo_storico:
            x_vecchie += [s['ore'], s['ore'], None]; y_vecchie += [0, max_c, None]
        else:
            x_nuove += [s['ore'], s['ore'], None]; y_nuove += [0, max_c, None]
    if x_vecchie:
        fig.add_trace(go.Scatter(x=x_vecchie, y=y_vecchie, mode='lines',
                                 line=dict(color='#6EA8D8', dash='dot', width=1),
                                 name='Dosi pregresse', hoverinfo='none'))
    if x_nuove:
        fig.add_trace(go.Scatter(x=x_nuove, y=y_nuove, mode='lines',
                                 line=dict(color='#475569', dash='dot', width=1),
                                 name='Nuove dosi', hoverinfo='none'))

    if livelli:
        x_liv = [l['ore'] for l in livelli]; y_liv = [l['conc'] for l in livelli]
        testi = [f"t={l['ore']:.1f}h<br>C={l['conc']:.3f} mg/L" +
                (f"<br>{l['nota']}" if l.get('nota') else "") for l in livelli]
        fig.add_trace(go.Scatter(x=x_liv, y=y_liv, mode='markers',
                                 marker=dict(symbol='star', size=14, color='#16A34A',
                                            line=dict(width=1, color='white')),
                                 name='Livelli misurati', text=testi, hoverinfo='text'))

    if farmaco == "Daptomicina":
        fig.add_hline(y=24.3, line_dash="dashdot", line_color="#DC2626",
                     annotation_text="Soglia tossicità CPK (Bhavnani 2010)",
                     annotation_font_color='#DC2626')

    fig.update_layout(xaxis_title="Ore dall'inizio della terapia (t = 0 prima dose)",
                      yaxis_title=f"{farmaco} plasmatica (mg/L)", **_PLOTLY_LAYOUT_BASE)
    return fig


def disegna_forecast(t_arr, c_arr, soglia, t_new_ore, dose_new, storia):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=t_arr, y=c_arr, mode='lines',
                             line=dict(color='#1E3A5F', width=2.5), name='Proiezione'))
    fig.add_hline(y=soglia, line_dash="dashdot", line_color="#DC2626",
                 annotation_text=f"Soglia target ({soglia:.2f} mg/L)",
                 annotation_font_color='#DC2626')
    fig.add_vline(x=t_new_ore, line_dash="dash", line_color="#16A34A",
                 annotation_text=f"Nuova dose ({dose_new:.0f} mg)",
                 annotation_font_color='#16A34A')

    max_c = max(c_arr, default=100)
    x_dosi, y_dosi = [], []
    for s in storia:
        x_dosi += [s['ore'], s['ore'], None]; y_dosi += [0, max_c, None]
    if x_dosi:
        fig.add_trace(go.Scatter(x=x_dosi, y=y_dosi, mode='lines',
                                 line=dict(color='#94A3B8', dash='dot', width=1),
                                 name='Dosi già somministrate', hoverinfo='none'))

    fig.update_layout(xaxis_title="Ore dall'inizio della terapia",
                      yaxis_title="Dalbavancina plasmatica (mg/L)", **_PLOTLY_LAYOUT_BASE)
    return fig

# ============================================================
# INTERFACCIA
# ============================================================

st.set_page_config(page_title="TDM Buddy v6.0", layout="wide", page_icon="💊")

st.title("💊 TDM Buddy v6.0")
st.caption("Monitoraggio terapeutico bayesiano — Daptomicina, Ciprofloxacina, Dalbavancina")

st.warning(
    "⚠️ **Strumento esclusivamente a scopo educational.** Non fornisce supporto alle "
    "decisioni cliniche e non è idoneo a orientare la terapia di pazienti reali: un "
    "utilizzo con finalità di supporto decisionale ne configurerebbe la natura di "
    "dispositivo medico, per cui questo strumento non è né validato né certificato. "
    "Nessun dato viene conservato sul server: l'intera cronologia esiste solo nel file "
    "JSON che scarichi e gestisci localmente. Non inserire nomi, cognomi o date di "
    "nascita esatte — usa un codice paziente non identificativo. Ogni file JSON è "
    "specifico per un paziente E per una molecola: non caricare lo storico di un "
    "farmaco diverso da quello selezionato."
)

# ------------------------------------------------------------
# STEP 1 — Selezione molecola (gate obbligatorio, card grandi)
# ------------------------------------------------------------

if st.session_state['farmaco_attivo'] is None:
    st.markdown("## Passo 1 — Seleziona la molecola da monitorare")
    st.caption("La scelta è vincolante per questa sessione: per cambiare molecola dovrai "
              "ricominciare (i dati non ancora scaricati andranno persi).")

    cols = st.columns(3)
    for col, farmaco in zip(cols, FARMACI):
        info = FARMACO_INFO[farmaco]
        with col:
            with st.container(border=True, height=165):
                st.markdown(f"### {info['icona']} {farmaco}")
                st.caption(info['descrizione'])
                if st.button(f"Seleziona {farmaco}", key=f"sel_{farmaco}",
                            use_container_width=True, type="primary"):
                    st.session_state['farmaco_attivo'] = farmaco
                    st.session_state['storico'] = storico_vuoto(farmaco)
                    st.rerun()
    st.stop()

farmaco = st.session_state['farmaco_attivo']

# ------------------------------------------------------------
# Badge molecola attiva + cambio
# ------------------------------------------------------------

col_badge, col_cambia = st.columns([5, 1])
with col_badge:
    st.info(f"{FARMACO_INFO[farmaco]['icona']} **Molecola in corso: {farmaco}**")
with col_cambia:
    if st.button("🔄 Cambia molecola", use_container_width=True):
        pulisci_widget_paziente()
        st.session_state['farmaco_attivo'] = None
        st.session_state['storico'] = storico_vuoto()
        st.session_state.pop('risultati', None)
        st.rerun()

# ------------------------------------------------------------
# SIDEBAR — anagrafica minima + storico JSON (scoped alla molecola)
# ------------------------------------------------------------

st.sidebar.header("🗂️ Storico paziente (locale)")
st.sidebar.caption(f"File compatibili: solo **{farmaco}**")

file_up = st.sidebar.file_uploader("📤 Carica storico (JSON)", type=["json"])
if file_up is not None:
    try:
        dati_caricati = json.load(file_up)
        farmaco_file = dati_caricati.get('farmaco')
        if farmaco_file != farmaco:
            st.sidebar.error(
                f"❌ Questo file è per **{farmaco_file or 'farmaco non specificato'}**, "
                f"non per **{farmaco}**. Cambia molecola (bottone in alto) per caricarlo, "
                f"oppure carica un file corretto."
            )
        else:
            pulisci_widget_paziente()
            st.session_state['storico'] = dati_caricati
            st.sidebar.success(
                f"✅ Storico caricato: {len(dati_caricati.get('monitoraggi', []))} "
                f"analisi per il codice '{dati_caricati.get('codice_paziente','?')}'."
            )
            st.rerun()
    except Exception:
        st.sidebar.error("File non valido o corrotto.")

id_paz = st.sidebar.text_input(
    "Codice paziente (non identificativo)",
    value=st.session_state['storico'].get('codice_paziente', ''),
    key="input_codice_paziente",
)
sesso = st.sidebar.selectbox(
    "Sesso biologico", ["Maschio", "Femmina"],
    index=0 if st.session_state['storico'].get('sesso', 'Maschio') == 'Maschio' else 1,
    key="input_sesso",
)
st.session_state['storico']['codice_paziente'] = id_paz
st.session_state['storico']['sesso'] = sesso
st.session_state['storico']['farmaco'] = farmaco

if st.sidebar.button("🆕 Nuovo paziente (stessa molecola)", use_container_width=True,
                     help="Azzera codice paziente, tabelle dosi/livelli e risultati "
                          "correnti. Non tocca file già scaricati."):
    pulisci_widget_paziente()
    st.session_state['storico'] = storico_vuoto(farmaco)
    st.session_state.pop('risultati', None)
    st.rerun()

n_rec = len(st.session_state['storico'].get('monitoraggi', []))
if n_rec > 0:
    st.sidebar.markdown("---")
    st.sidebar.info(f"📋 {n_rec} analisi in questa sessione.")
    df_sb = pd.DataFrame(st.session_state['storico']['monitoraggi'])
    cols_sb = [c for c in ['data_val', 'creatinina', 'mic'] if c in df_sb.columns]
    st.sidebar.dataframe(df_sb[cols_sb].tail(4), use_container_width=True, hide_index=True)
    st.sidebar.download_button(
        "📥 Scarica storico aggiornato",
        data=json.dumps(st.session_state['storico'], indent=2, ensure_ascii=False).encode('utf-8'),
        file_name=f"TDMBuddy_{farmaco}_{id_paz or 'paziente'}.json",
        mime="application/json", type="primary", use_container_width=True,
    )

if not id_paz:
    st.warning("👈 Inserisci un codice paziente nella barra laterale per procedere.")
    st.stop()

# ------------------------------------------------------------
# Recupero ultimo record per pre-compilazione (stesso farmaco, per
# costruzione: lo storico è già scoped a un'unica molecola)
# ------------------------------------------------------------

monitoraggi_farmaco = st.session_state['storico'].get('monitoraggi', [])
_ultimo = monitoraggi_farmaco[-1] if monitoraggi_farmaco else None

_gg_def   = int(_ultimo.get('eta_giorni')) if _ultimo else None
_peso_def = float(_ultimo.get('peso')) if _ultimo else None
_crea_def = float(_ultimo.get('creatinina')) if _ultimo else None
_mic_def  = float(_ultimo.get('mic')) if _ultimo else None

# ------------------------------------------------------------
# STEP 2 — Dati biometrici e clinici
# ------------------------------------------------------------

st.markdown("## Passo 2 — Dati biometrici e clinici")

with st.expander(f"📚 Modelli di riferimento applicati — {farmaco}", expanded=False):
    for fascia, rif in RIFERIMENTI_PER_FASCIA[farmaco]:
        st.caption(f"**{fascia}:** {rif}")

with st.container(border=True):
    st.markdown("#### Età e peso al momento di questa valutazione")
    st.caption("Possono cambiare rispetto all'ultima analisi: verificali sempre. "
              "Campi obbligatori, lasciati vuoti finché non li compili.")

    col_u, col_v, col_peso = st.columns(3)
    unita_default = "Anni"
    valore_eta_default = None
    if _gg_def is not None:
        valore_eta_default = round(_gg_def / 365.25, 2)

    with col_u:
        unita_eta = st.selectbox("Unità età", ["Anni", "Mesi", "Giorni"],
                                 key="input_unita_eta")
    max_per_unita = {"Anni": 100.0, "Mesi": 1200.0, "Giorni": 36500.0}
    with col_v:
        valore_eta = st.number_input(
            f"Età ({unita_eta.lower()})", min_value=0.0, max_value=max_per_unita[unita_eta],
            value=valore_eta_default, key="input_valore_eta",
            help="Inserisci l'età nell'unità di misura scelta a sinistra.")
    with col_peso:
        peso = st.number_input("Peso (kg)", 0.3, 300.0, value=_peso_def, key="input_peso")

    if valore_eta is None:
        gg = None
    elif unita_eta == "Anni":
        gg = int(round(valore_eta * 365.25))
    elif unita_eta == "Mesi":
        gg = int(round(valore_eta * 30.44))
    else:
        gg = int(round(valore_eta))
    eta = gg / 365.25 if gg is not None else None
    st.caption(f"Età calcolata: **{gg} giorni** ({eta:.2f} anni)" if gg is not None
              else "Età non ancora inserita.")

    modello_dalba = None
    if farmaco == "Dalbavancina":
        opzioni_modello = list(DALBA_MODELLI.keys())
        idx_default = 0
        if _ultimo and _ultimo.get('modello_dalba') in opzioni_modello:
            idx_default = opzioni_modello.index(_ultimo['modello_dalba'])
        modello_dalba = st.selectbox("🧬 Modello PK di riferimento:", opzioni_modello,
                                     index=idx_default, key="input_modello_dalba")

with st.container(border=True):
    st.markdown("#### Funzione renale, MIC, dialisi")
    col1, col2, col3 = st.columns(3)
    with col1:
        crea_help = None if farmaco != "Dalbavancina" else (
            "Non utilizzata dai modelli PK della dalbavancina (nessuna covariata "
            "renale validata) — raccolta solo a scopo di registrazione clinica.")
        crea = st.number_input("Creatinina sierica (mg/dL)", 0.05, 25.0, value=_crea_def,
                               step=0.05, key="input_crea", help=crea_help)
        if farmaco == "Dalbavancina":
            mic = st.number_input("MIC di riferimento (mg/L)", 0.008, 4.0, value=_mic_def,
                                  step=0.001, format="%.4f", key="input_mic",
                                  help="Il target PK/PD dettagliato con selezione MIC "
                                       "EUCAST/USCAST è disponibile nel pannello Forecast "
                                       "dopo l'analisi.")
        else:
            mic = st.number_input(f"MIC {farmaco} (mg/L)", 0.001, 16.0, value=_mic_def,
                                  step=0.001, format="%.4f", key="input_mic")
    with col2:
        crrt = st.checkbox("Paziente in CRRT?", key="input_crrt")
        if crrt:
            Sc  = st.number_input("Sieving coefficient (Sc)", 0.01, 0.60, 0.20, 0.01,
                                  key="input_sc")
            Qef = st.number_input("Flusso effluente Qef (L/h)", 0.5, 6.0, 2.0, 0.1,
                                  key="input_qef")
        else:
            Sc, Qef = 0.20, 2.0
    with col3:
        if farmaco == "Daptomicina" and gg is not None and gg < 6570:
            opt = ["Nessuna", "Lieve (GFR 60-90)", "Moderata (GFR 30-60)",
                   "Grave (GFR 15-30)", "ESRD (GFR <15)"]
            if 4380 <= gg < 6570:
                opt += ["Emodialisi (HD)", "Dialisi Peritoneale (CCPD)"]
            dialisi = st.selectbox("Funzione renale / dialisi:", opt, key="input_dialisi")
        else:
            dialisi = "Nessuna"
            st.caption("📌 ClCr calcolata con Cockcroft-Gault.")

campi_mancanti = []
if gg is None:
    campi_mancanti.append("età")
if peso is None:
    campi_mancanti.append("peso")
if mic is None:
    campi_mancanti.append("MIC")
if crea is None and farmaco != "Dalbavancina":
    campi_mancanti.append("creatinina")
if campi_mancanti:
    st.warning(f"⚠️ Campi ancora da compilare: {', '.join(campi_mancanti)}.")

# ------------------------------------------------------------
# STEP 3 — Timeline dosi e livelli
# ------------------------------------------------------------

st.markdown("## Passo 3 — Somministrazioni e livelli plasmatici")
st.info("💡 Le tabelle mostrano l'intera storia clinica del paziente per questa molecola. "
       "Modifica le righe esistenti solo se ti accorgi di un errore di trascrizione; "
       "aggiungi righe nuove per la valutazione odierna. Nessun limite al numero di "
       "dosi o livelli inseribili.")

if farmaco == "Dalbavancina":
    st.caption("📌 Dalbavancina: dosaggio tipico 22.5 mg/kg (<6 anni) o 18 mg/kg (≥6 anni), "
              "max 1500 mg/dose. Aggiungi una riga per ogni somministrazione realmente "
              "avvenuta, con la sua data/ora esatta.")

t0_db = datetime.now()
t_ultimo_storico = -1
dosi_list = []
if _ultimo and 'storia_dosi' in _ultimo and _ultimo['storia_dosi']:
    t0_db = datetime.strptime(_ultimo['data_prima_dose'], '%Y-%m-%d %H:%M')
    for s in _ultimo['storia_dosi']:
        dt = t0_db + timedelta(hours=s['ore'])
        dosi_list.append({"Data": dt.strftime("%d/%m/%Y"), "Ora": dt.strftime("%H:%M"),
                          "Dose (mg)": s['dose']})
    t_ultimo_storico = _ultimo['storia_dosi'][-1]['ore']
else:
    dosi_list.append({"Data": "", "Ora": "", "Dose (mg)": None})

liv_list = []
if _ultimo and 'livelli' in _ultimo:
    for l in _ultimo['livelli']:
        dt = t0_db + timedelta(hours=l['ore'])
        liv_list.append({"Data": dt.strftime("%d/%m/%Y"), "Ora": dt.strftime("%H:%M"),
                         "Conc. (mg/L)": l['conc'], "Nota": l.get('nota', '')})
liv_list.append({"Data": "", "Ora": "", "Conc. (mg/L)": None, "Nota": ""})

colA, colB = st.columns(2)
with colA:
    st.markdown("**💉 Somministrazioni**")
    df_dosi_edit = st.data_editor(
        pd.DataFrame(dosi_list), num_rows="dynamic", use_container_width=True,
        height=280,
        column_config={"Dose (mg)": st.column_config.NumberColumn(format="%.1f", min_value=0.1)},
        key="editor_dosi",
    )
with colB:
    st.markdown("**🔬 Livelli plasmatici**")
    df_liv_edit = st.data_editor(
        pd.DataFrame(liv_list), num_rows="dynamic", use_container_width=True,
        height=280,
        column_config={"Conc. (mg/L)": st.column_config.NumberColumn(format="%.4f", min_value=0.0)},
        key="editor_livelli",
    )

def _parse_data_ora(d_str, o_str):
    d_str, o_str = str(d_str).strip(), str(o_str).strip()
    if not d_str or not o_str or d_str.lower() == 'nan' or o_str.lower() == 'nan':
        return None
    d = None
    for fmt in ("%d/%m/%Y", "%Y-%m-%d"):
        try:
            d = datetime.strptime(d_str, fmt).date(); break
        except ValueError:
            continue
    if d is None:
        return None
    try:
        t = datetime.strptime(o_str[:5], "%H:%M").time()
    except ValueError:
        return None
    return datetime.combine(d, t)

storia_completa = []
for _, r in df_dosi_edit.iterrows():
    dt_obj = _parse_data_ora(r.get("Data"), r.get("Ora"))
    try:
        dose_val = float(r.get("Dose (mg)"))
    except (TypeError, ValueError):
        dose_val = None
    if dt_obj is not None and dose_val is not None and not pd.isna(dose_val):
        storia_completa.append({'dt': dt_obj, 'dose': dose_val})

livelli_completi = []
for _, r in df_liv_edit.iterrows():
    dt_obj = _parse_data_ora(r.get("Data"), r.get("Ora"))
    try:
        conc_val = float(r.get("Conc. (mg/L)"))
    except (TypeError, ValueError):
        conc_val = None
    if dt_obj is not None and conc_val is not None and not pd.isna(conc_val):
        livelli_completi.append({'dt': dt_obj, 'conc': conc_val, 'nota': str(r.get("Nota", "") or "")})

# ------------------------------------------------------------
# AVVIO ANALISI
# ------------------------------------------------------------

if st.button("🚀 Avvia il calcolo bayesiano", type="primary", use_container_width=True):
    if campi_mancanti:
        st.error(f"⚠️ Compila prima i campi obbligatori: {', '.join(campi_mancanti)}.")
    elif not storia_completa:
        st.error("Nessuna somministrazione valida in tabella.")
    elif not livelli_completi:
        st.error("Inserisci almeno un livello plasmatico valido.")
    else:
        t0_assoluto = min(s['dt'] for s in storia_completa)
        for s in storia_completa:
            s['ore'] = (s['dt'] - t0_assoluto).total_seconds() / 3600.0
        storia_completa.sort(key=lambda x: x['ore'])
        for l in livelli_completi:
            l['ore'] = (l['dt'] - t0_assoluto).total_seconds() / 3600.0
        livelli_completi.sort(key=lambda x: x['ore'])

        dt_val = (max((l['dt'] for l in livelli_completi), default=None) or
                 max((s['dt'] for s in storia_completa), default=datetime.now())
                 ).strftime('%Y-%m-%d %H:%M')

        # Rimuove la chiave 'dt' (oggetto datetime, non serializzabile in
        # JSON): da qui in avanti servono solo 'ore' e 'dose'/'conc'.
        storia_completa = [{'ore': s['ore'], 'dose': float(s['dose'])} for s in storia_completa]
        livelli_completi = [{'ore': l['ore'], 'conc': float(l['conc']), 'nota': l.get('nota', '')}
                            for l in livelli_completi]

        with st.spinner("Calcolo MAP Bayesiano in corso…"):
            dose_mg = storia_completa[-1]['dose']
            dose_mg_kg = dose_mg / peso
            tau_h_display = tau_mediano(storia_completa)
            ok_staph = ok_ent = ok_target = tox = None

            if farmaco == "Daptomicina":
                prior = prior_dapto(sesso, gg, peso, crea, crrt, Sc, Qef, dialisi)
                CL_ind, V_ind, tipo_fit = map_1cmt(prior['CL_p'], prior['V_p'], prior['oCL'],
                                                   prior['oV1'], prior['sig'], storia_completa,
                                                   livelli_completi, prior['t_inf'])
                Q_ind = V2_ind = 0.0
                Vss_d = V_ind
                AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(
                    pk_1cmt, (CL_ind, V_ind), storia_completa, prior['t_inf'])
                auc_mic = AUC24_v / mic
                emivita = float(np.log(2)*V_ind/CL_ind)
                ok_staph = auc_mic >= PD_DAPTO['AUC_MIC_staph']
                ok_ent   = auc_mic >= PD_DAPTO['AUC_MIC_entero']
                tox      = Cmin_ss >= PD_DAPTO['Cmin_tox']

            elif farmaco == "Dalbavancina":
                prior = prior_dalba(peso, modello_dalba)
                CL_ind, tipo_fit = map_dalba(prior['CL_p'], prior['V_p'], prior['Q_p'],
                                             prior['V2_p'], prior['oCL'], prior['sig'],
                                             storia_completa, livelli_completi, prior['t_inf'])
                V_ind, Q_ind, V2_ind = prior['V_p'], prior['Q_p'], prior['V2_p']
                Vss_d = V_ind + V2_ind
                AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(
                    pk_2cmt, (CL_ind, V_ind, Q_ind, V2_ind), storia_completa, prior['t_inf'])
                auc_mic = AUC24_v / mic
                emivita = float(np.log(2)*Vss_d/CL_ind)
                target_riferimento = list(DALBA_TARGETS.values())[0]
                ok_target = auc_mic >= target_riferimento
                tox = False

            else:
                prior = prior_cipro(sesso, gg, peso, crea, crrt, dialisi)
                CL_ind, V_ind, Q_ind, V2_ind, tipo_fit = map_2cmt(
                    prior['CL_p'], prior['V_p'], prior['Q_p'], prior['V2_p'],
                    prior['oCL'], prior['oV1'], prior['oQ'], prior['oV2'], prior['sig'],
                    storia_completa, livelli_completi, prior['t_inf'])
                Vss_d = V_ind + V2_ind
                AUC24_v, Cmax_ss, Cmin_ss = metriche_esposizione_24h(
                    pk_2cmt, (CL_ind, V_ind, Q_ind, V2_ind), storia_completa, prior['t_inf'])
                auc_mic = AUC24_v / mic
                emivita = float(np.log(2)*Vss_d/CL_ind)
                ok_target = auc_mic >= PD_CIPRO['AUC_MIC_target']
                tox = None  # nessuna soglia di tossicità AUC validata per ciprofloxacina

            # Tutto ciò che serve alla visualizzazione viene salvato in
            # session_state: il pannello Forecast (sotto) legge da qui,
            # cosi' i suoi controlli si possono toccare liberamente senza
            # far sparire il calcolo (st.button torna False ad ogni rerun
            # successivo al click).
            st.session_state['risultati'] = {
                'farmaco': farmaco, 'id_paz': id_paz,
                'CL_ind': CL_ind, 'V_ind': V_ind, 'Q_ind': Q_ind, 'V2_ind': V2_ind,
                'Vss_d': Vss_d, 'AUC24_v': AUC24_v, 'Cmax_ss': Cmax_ss, 'Cmin_ss': Cmin_ss,
                'auc_mic': auc_mic, 'emivita': emivita, 'tipo_fit': tipo_fit, 'prior': prior,
                'storia_completa': storia_completa, 'livelli_completi': livelli_completi,
                't0_assoluto': t0_assoluto, 'dt_val': dt_val, 'dose_mg': dose_mg,
                'dose_mg_kg': dose_mg_kg, 'tau_h_display': tau_h_display, 'mic': mic,
                'ok_staph': ok_staph, 'ok_ent': ok_ent, 'ok_target': ok_target, 'tox': tox,
                't_ultimo_storico': t_ultimo_storico, 'sesso': sesso, 'eta': eta, 'gg': gg,
                'peso': peso, 'crea': crea, 'dialisi': dialisi, 'modello_dalba': modello_dalba,
                'crrt': crrt, 'Sc': Sc, 'Qef': Qef,
                'registrato': False,
            }
            st.info("🧮 Calcolo completato. Controlla i risultati qui sotto — se vuoi "
                   "ricalcolare con dati diversi puoi farlo liberamente prima di "
                   "registrare l'occasione nello storico.")

# ------------------------------------------------------------
# VISUALIZZAZIONE PERSISTENTE DEI RISULTATI
# (fuori dal blocco del bottone: i controlli del pannello Forecast
# possono essere modificati liberamente senza far sparire il calcolo)
# ------------------------------------------------------------

if ('risultati' in st.session_state
        and st.session_state['risultati'].get('farmaco') == farmaco
        and st.session_state['risultati'].get('id_paz') == id_paz):

    R = st.session_state['risultati']
    CL_ind, V_ind, Q_ind, V2_ind = R['CL_ind'], R['V_ind'], R['Q_ind'], R['V2_ind']
    Vss_d, AUC24_v, Cmax_ss, Cmin_ss = R['Vss_d'], R['AUC24_v'], R['Cmax_ss'], R['Cmin_ss']
    auc_mic, emivita, tipo_fit, prior = R['auc_mic'], R['emivita'], R['tipo_fit'], R['prior']
    storia_completa, livelli_completi = R['storia_completa'], R['livelli_completi']
    t0_assoluto, dt_val = R['t0_assoluto'], R['dt_val']
    dose_mg, dose_mg_kg = R['dose_mg'], R['dose_mg_kg']
    tau_h_display, mic = R['tau_h_display'], R['mic']
    ok_staph, ok_ent, ok_target, tox = R['ok_staph'], R['ok_ent'], R['ok_target'], R['tox']
    t_ultimo_storico = R['t_ultimo_storico']
    sesso, eta, gg = R['sesso'], R['eta'], R['gg']
    peso, crea, dialisi = R['peso'], R['crea'], R['dialisi']
    crrt, Sc, Qef = R.get('crrt', False), R.get('Sc', 0.20), R.get('Qef', 2.0)
    modello_dalba_r = R.get('modello_dalba')

    st.markdown("---")
    st.subheader("🔬 Risultati PK/PD")
    st.caption(f"**Modello di popolazione:** {prior['nm']}")
    if prior['note']:
        st.info(prior['note'])

    c1, c2, c3 = st.columns(3)
    c1.metric("Clearance (CL)", f"{CL_ind:.3f} L/h",
             delta=f"Popolazione: {prior['CL_p']:.3f} L/h", delta_color="off")
    c2.metric("Volume (Vss)", f"{Vss_d:.1f} L",
             delta=f"Popolazione: {(prior['V_p'] + prior.get('V2_p',0)):.1f} L", delta_color="off")
    c3.metric("Emivita (T½)", f"{emivita:.1f} h")

    c4, c5, c6 = st.columns(3)
    c4.metric("AUC₂₄ stimata", f"{AUC24_v:.0f} mg·h/L")

    if farmaco == "Daptomicina":
        c5.metric("AUC₂₄/MIC (target ≥ 666)", f"{auc_mic:.0f}",
                 delta="✅ Target raggiunto" if ok_staph else "⚠️ Sotto target",
                 delta_color="normal" if ok_staph else "inverse")
        c6.metric("Cmin prevista (soglia CPK ≥ 24.3)", f"{Cmin_ss:.1f} mg/L",
                 delta="⚠️ Rischio CPK" if tox else "✅ Sicuro",
                 delta_color="inverse" if tox else "normal")
    elif farmaco == "Dalbavancina":
        c5.metric("fAUC₂₄h/MIC (target minimo ≥ 27.1)", f"{auc_mic:.0f}",
                 delta="✅ Target raggiunto" if ok_target else "⚠️ Sotto target",
                 delta_color="normal" if ok_target else "inverse")
        c6.metric("Cmax / Cmin previste", f"{Cmax_ss:.1f} / {Cmin_ss:.1f} mg/L",
                 help="Target su infezioni complesse e pianificazione redose "
                      "nel pannello Forecast qui sotto.")
    else:
        c5.metric("AUC₂₄/MIC (target ≥ 125)", f"{auc_mic:.0f}",
                 delta="✅ Target raggiunto" if ok_target else "⚠️ Sotto target",
                 delta_color="normal" if ok_target else "inverse")
        c6.metric("AUC₂₄ (esposizione sistemica)", f"{AUC24_v:.0f} mg·h/L",
                 help="Nessuna soglia AUC quantitativa è validata in letteratura per la "
                      "neurotossicità da fluorochinoloni. Il rischio è documentato "
                      "soprattutto in insufficienza renale non corretta: valutare "
                      "clinicamente, non solo sul valore di AUC.")

    st.markdown("### 📈 Profilo farmacocinetico")
    fig = disegna_grafico(farmaco, storia_completa, livelli_completi, CL_ind, V_ind,
                          Q_ind, V2_ind, prior['CL_p'], prior['V_p'], prior['Q_p'],
                          prior.get('V2_p', 0.0), prior['t_inf'], t_ultimo_storico)
    st.plotly_chart(fig, use_container_width=True, theme=None)

    def _esito(ok):
        return "target raggiunto" if ok else "insufficiente — rivalutare la posologia"

    if farmaco == "Dalbavancina":
        pd_block = (f"  Efficacia (fAUC24h/MIC >= 27.1, target minimo): {auc_mic:.0f} "
                   f"-> {_esito(ok_target)}\n"
                   f"  Cmax prevista (24h dopo l'ultima dose): {Cmax_ss:.1f} mg/L\n"
                   f"  Cmin prevista (24h dopo l'ultima dose): {Cmin_ss:.1f} mg/L\n"
                   f"  Nota: per il target restrittivo (infezioni complesse) e la\n"
                   f"  pianificazione della prossima dose vedere il pannello Forecast.\n")
    elif farmaco == "Daptomicina":
        pd_block = (f"  Efficacia vs S. aureus/MRSA (AUC24/MIC >= 666): {auc_mic:.0f} "
                   f"-> {_esito(ok_staph)}\n"
                   f"  Efficacia vs E. faecium (AUC24/MIC >= 143): {auc_mic:.0f} "
                   f"-> {_esito(ok_ent)}\n"
                   f"  Sicurezza muscolare (Cmin_ss >= 24.3 mg/L): {Cmin_ss:.1f} mg/L "
                   f"-> {'attenzione, monitorare CPK' if tox else 'sicuro'}\n")
    else:
        pd_block = (f"  Efficacia gram-negativi (AUC24/MIC >= 125): {auc_mic:.0f} "
                   f"-> {_esito(ok_target)}\n"
                   f"  Esposizione sistemica: AUC24 = {AUC24_v:.0f} mg*h/L\n"
                   f"  Nota sicurezza: nessuna soglia AUC quantitativa validata per la\n"
                   f"  neurotossicita' da fluorochinoloni. Il rischio e' documentato\n"
                   f"  soprattutto in insufficienza renale non corretta: valutare\n"
                   f"  clinicamente.\n")

    referto = f"""\
=== TDM BUDDY v6.0 — REFERTO {farmaco.upper()} ===
Data: {dt_val}  |  Codice paziente: {id_paz}  |  Sesso: {sesso}  |  Età: {eta:.1f} anni  |  Peso: {peso} kg

[MODELLO APPLICATO]
{prior['nm']}

[DATI CLINICI]
Creatinina sierica: {crea} mg/dL
Dialisi / CRRT: {dialisi}
Dose più recente: {dose_mg:.0f} mg ({dose_mg_kg:.1f} mg/kg)
Intervallo mediano tra dosi: {(f'{tau_h_display:.1f} h' if tau_h_display else 'N/D (una sola dose)')}
MIC patogeno: {mic} mg/L
Livelli plasmatici usati per il fit: {len(livelli_completi)}

[PARAMETRI PK INDIVIDUALI — MAP Bayesiano]
Algoritmo: {tipo_fit}
Clearance (CL): {CL_ind:.3f} L/h
Volume di distribuzione (Vss): {Vss_d:.2f} L
Emivita terminale (T1/2): {emivita:.1f} h

[VALUTAZIONE PK/PD]
AUC24 stimata (finestra 24h dopo l'ultima dose): {AUC24_v:.0f} mg*h/L
{pd_block}"""

    st.markdown("### 📋 Referto per cartella clinica")
    st.code(referto, language='text')

    # ---- FORECAST DALBAVANCINA ----
    if farmaco == "Dalbavancina":
        with st.container(border=True):
            st.markdown("#### 🔮 Forecast — pianificazione della prossima dose")
            st.caption("Proietta il profilo del paziente in avanti nel tempo per "
                      "stimare quando la concentrazione scenderà sotto la soglia "
                      "target scelta.")

            colF1, colF2, colF3 = st.columns(3)
            with colF1:
                target_scelto = st.selectbox(
                    "Target PK/PD:",
                    list(DALBA_TARGETS.keys()) + ["Trough C0/MIC ≥ 200"],
                    help="Tutti e tre i target fAUC24h/MIC sono dallo stesso "
                         "Monte Carlo simulation: Cojutti PG, et al. Antimicrob "
                         "Agents Chemother 2021;65(5):e02260-20.",
                    key="fc_target",
                )
            with colF2:
                mic_label = st.selectbox("MIC / breakpoint:",
                                         list(DALBA_MIC_OPZIONI.keys()) + ["Personalizzata"],
                                         key="fc_mic_label")
                if mic_label == "Personalizzata":
                    mic_fc = st.number_input("MIC personalizzata (mg/L)", 0.001, 4.0,
                                             0.125, step=0.001, format="%.4f", key="fc_mic_custom")
                else:
                    mic_fc = DALBA_MIC_OPZIONI[mic_label]
            with colF3:
                f_unb = st.number_input("Frazione libera (f)", 0.01, 0.10,
                                       DALBA_PARAMS['fu'], 0.01, key="fc_funb",
                                       help="Default 0.07 (93% legame proteico).")

            st.caption("*EUCAST S. aureus* = breakpoint attuale, tabelle v.16.0 "
                      "(2026); *USCAST* = ente statunitense, soglia diversa da "
                      "EUCAST — non una versione precedente.")
            if target_scelto in DALBA_TARGETS_DESC:
                st.caption(f"**{target_scelto}** → {DALBA_TARGETS_DESC[target_scelto]} "
                          f"(Cojutti PG, et al. *Antimicrob Agents Chemother.* "
                          f"2021;65(5):e02260-20).")

            if target_scelto == "Trough C0/MIC ≥ 200":
                soglia = dalba_soglia_trough_c0mic(mic_fc)
                desc = f"C0/MIC ≥ 200 → soglia totale = {soglia:.2f} mg/L"
            else:
                soglia = dalba_soglia_totale(DALBA_TARGETS[target_scelto], mic_fc, f_unb)
                desc = f"{target_scelto} → soglia totale equivalente = {soglia:.2f} mg/L"

            st.info(f"📐 Soglia calcolata: **{desc}**")

            t_val = max(max((l['ore'] for l in livelli_completi), default=0),
                       storia_completa[-1]['ore'])
            c_att = pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia_completa, t_val, prior['t_inf'])
            t_cross = trova_tempo_soglia(CL_ind, V_ind, Q_ind, V2_ind, storia_completa,
                                         soglia, prior['t_inf'], t_partenza=t_val,
                                         t_max=t_val+2000, passo=2.0)

            if t_cross is None:
                st.warning("⚠️ La soglia non viene raggiunta entro ~83 giorni "
                          "dalla valutazione attuale.")
            elif t_cross <= t_val + 0.5:
                st.error(f"🔴 Concentrazione già sotto soglia "
                        f"({c_att:.2f} < {soglia:.2f} mg/L).")
            else:
                data_cross = t0_assoluto + timedelta(hours=t_cross)
                st.success(f"🗓️ La concentrazione scenderà sotto soglia il "
                          f"**{data_cross.strftime('%d/%m/%Y alle %H:%M')}** "
                          f"(tra circa {(t_cross-t_val)/24.0:.1f} giorni).")

            st.markdown("##### Simula una nuova somministrazione ipotetica")
            col_s1, col_s2, col_s3 = st.columns(3)
            with col_s1:
                simula = st.checkbox("Aggiungi dose ipotetica al forecast", key="fc_simula")
            with col_s2:
                if simula:
                    data_new = st.date_input(
                        "Data nuova dose",
                        value=(t0_assoluto + timedelta(hours=t_val+24)).date(),
                        format="DD/MM/YYYY", key="fc_data_new")
                    ora_new = st.time_input("Ora", value=datetime.strptime("10:00","%H:%M").time(),
                                            key="fc_ora_new")
            with col_s3:
                if simula:
                    dose_new = st.number_input("Dose ipotetica (mg)", 100.0, 1500.0, 1500.0, 50.0,
                                               key="fc_dose_new")

            if simula:
                t_new = (datetime.combine(data_new, ora_new) - t0_assoluto).total_seconds()/3600.0
                storia_sim = sorted(storia_completa + [{'ore': t_new, 'dose': dose_new}],
                                    key=lambda x: x['ore'])
                t_cross_sim = trova_tempo_soglia(CL_ind, V_ind, Q_ind, V2_ind, storia_sim,
                                                 soglia, prior['t_inf'], t_partenza=t_new,
                                                 t_max=t_new+2000, passo=2.0)
                if t_cross_sim:
                    data_cs = t0_assoluto + timedelta(hours=t_cross_sim)
                    st.success(f"✅ Con la dose ipotetica di {dose_new:.0f} mg, "
                              f"la concentrazione resterà sopra soglia fino al "
                              f"**{data_cs.strftime('%d/%m/%Y')}**.")
                t_fine_sim = (t_cross_sim + 48) if t_cross_sim else (t_new + 500)
                t_arr_sim = np.linspace(0, t_fine_sim, 1000)
                c_arr_sim = [pk_2cmt(CL_ind, V_ind, Q_ind, V2_ind, storia_sim, t, prior['t_inf'])
                            for t in t_arr_sim]
                fig2 = disegna_forecast(t_arr_sim, c_arr_sim, soglia, t_new, dose_new,
                                        storia_completa)
                st.plotly_chart(fig2, use_container_width=True, theme=None)

    # ---- REGISTRAZIONE ESPLICITA NELLO STORICO ----
    # Separata dal calcolo: premere "Avvia il calcolo" può servire solo a
    # ricalcolare (correggere un refuso, provare una MIC diversa) senza
    # per questo creare una nuova "occasione clinica" nello storico.
    st.markdown("---")
    if R.get('registrato'):
        st.success("✅ Questa valutazione è già stata registrata nello storico di sessione. "
                  "Ricalcolando otterrai una nuova valutazione da registrare separatamente.")
    else:
        st.info("La valutazione qui sopra **non è ancora stata registrata** nello storico "
               "del paziente. Registrala solo quando corrisponde a un'occasione clinica "
               "reale (dosi e prelievi effettivamente avvenuti).")
        if st.button("💾 Registra questa valutazione nello storico", type="primary",
                     use_container_width=True):
            record = {
                'data_val': dt_val, 'farmaco': farmaco, 'eta_giorni': gg, 'peso': peso,
                'creatinina': crea, 'crrt': bool(crrt), 'Sc': Sc, 'Qef': Qef,
                'dialisi_tipo': dialisi,
                'dose_mg': dose_mg, 'dose_mg_kg': dose_mg_kg, 'tau_h': tau_h_display,
                'n_dosi': len(storia_completa), 'mic': mic, 'n_livelli': len(livelli_completi),
                'CL_ind': CL_ind, 'V_ind': Vss_d, 'AUC24': AUC24_v, 'AUC_MIC': auc_mic,
                'Cmin_ss': Cmin_ss, 'modello': prior['nm'], 'modello_dalba': modello_dalba_r,
                'storia_dosi': storia_completa, 'livelli': livelli_completi,
                'data_prima_dose': t0_assoluto.strftime('%Y-%m-%d %H:%M'),
            }
            st.session_state['storico']['monitoraggi'].append(record)
            st.session_state['storico']['codice_paziente'] = id_paz
            st.session_state['storico']['sesso'] = sesso
            st.session_state['storico']['farmaco'] = farmaco
            st.session_state['risultati']['registrato'] = True
            st.success("💾 Valutazione registrata nello storico di questa sessione. "
                      "Scarica il JSON dalla barra laterale per conservarla.")
            st.rerun()

# ------------------------------------------------------------
# STORICO DI SESSIONE
# ------------------------------------------------------------

st.markdown("---")
with st.expander(f"📚 Storico completo di sessione — {farmaco}", expanded=False):
    monitoraggi = st.session_state['storico'].get('monitoraggi', [])
    if not monitoraggi:
        st.info("Nessuna analisi TDM in questa sessione.")
    else:
        df_sto = pd.DataFrame(monitoraggi)
        cols_vis = [c for c in ['data_val','creatinina','mic','dose_mg','n_dosi',
                                'n_livelli','CL_ind','AUC24'] if c in df_sto.columns]
        st.dataframe(df_sto[cols_vis], use_container_width=True)

        colS1, colS2 = st.columns([3, 1])
        with colS1:
            scelta = st.selectbox(
                "Seleziona per rimuovere un record errato:",
                ["— Scegli —"] + [f"#{i} — {r['data_val']}" for i, r in enumerate(monitoraggi)])
        with colS2:
            st.write(""); st.write("")
            if scelta != "— Scegli —":
                idx_sel = int(scelta.split("#")[1].split(" ")[0])
                if st.button("🗑️ Rimuovi record", use_container_width=True):
                    st.session_state['storico']['monitoraggi'].pop(idx_sel)
                    st.rerun()
