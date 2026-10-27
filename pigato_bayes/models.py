"""
Modelli di popolazione (prior) per la TDM bayesiana.

Ogni funzione prior_* restituisce un dict:
  par:    valori tipici {'CL','V1','Q','V2'} in L/h e L (Q=V2=0 → 1 cmt)
  omega:  SD su scala log per parametro (0 = fisso)
  sigma_prop, sigma_add: errore residuo
  t_inf:  durata infusione predefinita (h)
  nome:   riferimento bibliografico del modello applicato
  note:   lista di avvertenze da mostrare
  fuori_range: True se il paziente è fuori dalla popolazione del modello

Stato della verifica sulle fonti (26/10/2026): vedi VERIFICA in fondo al file.
"""

from .pk_engine import omega_da_cv
from . import renal

# ============================================================
# DAPTOMICINA
# ============================================================

# Soraluce 2018, Tabella 3 (modello finale): CL = CLNR + θ·(Clcr/49) [+ CLEC],
# CLNR 0,160 L/h, θ 0,367 L/h, V 12,30 L, IIV CL 36,7% e V 27,8%,
# errore additivo su scala log 0,123. Clcr con Cockcroft-Gault.
SORALUCE = {'CLNR': 0.160, 'CLR_theta': 0.367, 'Clcr_ref': 49.0, 'V': 12.30,
            'omega_CL': omega_da_cv(36.7), 'omega_V': omega_da_cv(27.8), 'sigma': 0.123}

# CL in mL/h/kg, SD in mL/h/kg, V in L/kg. Bande 3–24 mesi: Wei 2020, valori
# coerenti con Bradley JS, et al. Pediatr Infect Dis J 2014;33:936-939
# (CL = dose/AUC: 18,6 / 18,3 / 21,3 mL/h/kg a 3–6 / 7–12 / 13–24 mesi).
WEI = {
    '3-6m':   {'CL': 19.7, 'SD': 5.5, 'V': 0.17},
    '7-12m':  {'CL': 19.6, 'SD': 5.8, 'V': 0.17},
    '13-24m': {'CL': 21.8, 'SD': 3.0, 'V': 0.11},
    '2-6y':   {'CL': 20.0, 'SD': 6.1, 'V': 0.11},
    '7-11y':  {'CL': 17.0, 'SD': 5.7, 'V': 0.11},
    '12-17y': {'CL': 11.0, 'SD': 2.1, 'V': 0.11},
}
# Ye L, et al. Front Pharmacol 2022;13:838599, Tabella 4: rapporto di AUC
# (IR/sano) per fascia d'età; la CL si divide per lo stesso fattore.
YE2022 = {
    '12-17y': {'mild': 0.98, 'moderate': 1.08, 'severe': 1.24, 'esrd': 1.33},
    '7-11y':  {'mild': 1.00, 'moderate': 1.19, 'severe': 1.49, 'esrd': 1.73},
    '2-6y':   {'mild': 1.02, 'moderate': 1.28, 'severe': 1.74, 'esrd': 2.20},
    '1-2y':   {'mild': 1.04, 'moderate': 1.29, 'severe': 1.72, 'esrd': 2.13},
}
LIM = {'Vss': 5.90, 'CLren': 0.0226, 'CLnr': 0.158, 'CL_HD': 0.619,
       'omega_V1': 0.112, 'omega_knr': 0.0878, 'sigma': 0.085, 'peso_ref': 40.0}
PD_DAPTO = {'AUC_MIC_staph': 666, 'AUC_MIC_entero': 143, 'Cmin_tox': 24.3}

IRC_OPZIONI = ["Nessuna", "Lieve (GFR 60-90)", "Moderata (GFR 30-60)",
               "Grave (GFR 15-30)", "ESRD (GFR <15)"]
IRC_MAP = {"Lieve (GFR 60-90)": 'mild', "Moderata (GFR 30-60)": 'moderate',
           "Grave (GFR 15-30)": 'severe', "ESRD (GFR <15)": 'esrd'}


def _banda_wei(gg):
    if gg < 213:
        return '3-6m'
    if gg < 365:
        return '7-12m'
    if gg < 730:
        return '13-24m'
    if gg < 2190:
        return '2-6y'
    if gg < 4380:
        return '7-11y'
    return '12-17y'


def _banda_ye(banda_wei):
    # Ye 2022 non ha fasce sotto 1 anno: 3–24 mesi → 1-2y (estrapolazione sotto l'anno)
    return {'3-6m': '1-2y', '7-12m': '1-2y', '13-24m': '1-2y'}.get(banda_wei, banda_wei)


def prior_dapto(sesso, gg, peso, crea, crrt, Sc, Qef, dialisi, **_):
    eta = gg / 365.25
    note = []
    fuori = False
    if 4380 <= gg < 6570 and dialisi in ("Emodialisi (HD)", "Dialisi Peritoneale (CCPD)"):
        sc = peso / LIM['peso_ref']
        CL = (LIM['CLren'] + LIM['CLnr']) * sc
        V = LIM['Vss'] * sc
        om = {'CL': LIM['omega_knr'], 'V1': LIM['omega_V1']}
        nome = 'Lim SY, et al. J Pediatr Pharmacol Ther 2021;26(2):123-132 [adolescente HD/CCPD]'
        note.append("Valori assoluti di Lim scalati linearmente sul peso (riferimento 40 kg, "
                    "assunto: coorte 29–64 kg).")
        if dialisi == "Emodialisi (HD)":
            note.append("La clearance intradialitica (CL_HD 0,619 L/h) non è modellata: "
                        "somministrare dopo la seduta.")
        return dict(par={'CL': CL, 'V1': V, 'Q': 0, 'V2': 0}, omega=om, sigma_prop=LIM['sigma'],
                    sigma_add=0.0, t_inf=0.5, nome=nome, note=note, fuori_range=False)

    if gg < 6570:
        if gg < 91:
            fuori = True
            note.append("Età < 3 mesi: nessun dato PK pubblicato; applicati i valori 3–6 mesi "
                        "(Bradley 2014 / Wei 2020) come estrapolazione.")
        b = _banda_wei(max(gg, 91))
        p = WEI[b]
        CL = p['CL'] * peso / 1000.0
        V = p['V'] * peso
        om = {'CL': omega_da_cv(100 * p['SD'] / p['CL']), 'V1': 0.25}
        nome = f"Wei XC, et al. J Clin Pharmacol 2020;60(6):768-774 [{b}]"
        if gg < 730:
            nome += " (dati originali: Bradley JS, et al. Pediatr Infect Dis J 2014;33:936-939)"
        t_inf = 1.0 if gg < 2190 else 0.5
        if dialisi in IRC_MAP:
            by = _banda_ye(b)
            f = YE2022[by][IRC_MAP[dialisi]]
            CL /= f
            note.append(f"Correzione renale Ye 2022 (fascia {by}): CL ÷ {f:.2f}."
                        + (" Fascia 1–2 anni applicata sotto l'anno (estrapolazione)." if gg < 365 else ""))
    else:
        clcr = renal.cockcroft_gault(eta, peso, crea, sesso) or 60.0
        CL = SORALUCE['CLNR'] + SORALUCE['CLR_theta'] * (clcr / SORALUCE['Clcr_ref'])
        V = SORALUCE['V']
        om = {'CL': SORALUCE['omega_CL'], 'V1': SORALUCE['omega_V']}
        nome = (f"Soraluce A, et al. Int J Antimicrob Agents 2018;52(2):158-165 "
                f"[adulto in terapia intensiva; ClCr CG {clcr:.0f} mL/min]")
        t_inf = 0.5
    if crrt:
        CLEC = Sc * Qef
        CL += CLEC
        note.append(f"CRRT: CL extracorporea = Sc × Qef = {CLEC:.3f} L/h aggiunta.")
    return dict(par={'CL': CL, 'V1': V, 'Q': 0, 'V2': 0}, omega=om, sigma_prop=0.12 if gg < 6570 else SORALUCE['sigma'],
                sigma_add=0.0, t_inf=t_inf, nome=nome, note=note, fuori_range=fuori)


# ============================================================
# CIPROFLOXACINA
# ============================================================
# Adulto: Forrest A, et al. Antimicrob Agents Chemother 1993;37:1065-1072,
#   Tabella 2 (74 adulti acutamente malati, IT2S): V1 0,69 L/kg (CV 26%),
#   Vp 0,51 L/kg (33%), CLd 38 L/h/1,73 m² (24%), CLT 17 L/h/1,73 m² (44%);
#   CLT (L/h) = (0,00145·CLCR + 0,167)·peso, CLCR in mL/min/1,73 m².
# Pediatria: Rajagopalan P, Gastonguay MR. J Clin Pharmacol 2003;43:698-710
#   (150 pazienti 0,27–16,9 anni): CL = 30,3·(WT/70)^0,75·(1+0,045·(età−2,5)),
#   Vc = 56,7·(WT/70), Vp = 89,8·(WT/70), Q = 37,5·(WT/70)^0,75.
#   IIV per parametro non riportata nell'abstract (intervallo 22,5–49,8%):
#   valori assunti sotto, da sostituire con la tabella del full text.

# Hirt D, et al. Eur J Clin Pharmacol 2021;77:1687-1695 (189 bambini 0,01–17,8 anni,
#   EV e orale; Tabella 2, modello finale costruzione + validazione):
#   CL = 9,07·(BW/20)^0,75·(eGFR/199)^0,337·(età/7,3)^0,0828 L/h
#   V2 (centrale) = 26,5·BW/20; Q = 3,9·(BW/20)^0,75; V3 = 17,4·BW/20
#   IIV solo su CL (ω 0,347); errore proporzionale 0,423. eGFR con Schwartz.
HIRT = {'CL': 9.07, 'Vc': 26.5, 'Q': 3.9, 'Vp': 17.4, 'egfr_ref': 199.0, 'egfr_exp': 0.337,
        'age_ref': 7.3, 'age_exp': 0.0828, 'omega': {'CL': 0.347}, 'sigma_prop': 0.423}
CIPRO_MODELLI = ['Hirt 2021 — pediatrico 0–18 anni (eGFR, età)', 'Rajagopalan 2003 — pediatrico 3 mesi–17 anni',
                 'Forrest 1993 — adulto acutamente malato']

RAJAGOPALAN = {'CL': 30.3, 'Vc': 56.7, 'Vp': 89.8, 'Q': 37.5, 'age_eff': 0.045,
               'omega': {'CL': 0.40, 'V1': 0.40, 'Q': 0.30, 'V2': 0.30},
               'sigma_prop': 0.27, 'sigma_add': 0.04}
FORREST = {'V1': 0.69, 'Vp': 0.51, 'CLd_173': 38.0,
           'omega': {'CL': omega_da_cv(44), 'V1': omega_da_cv(26), 'Q': omega_da_cv(24), 'V2': omega_da_cv(33)},
           'sigma_prop': 0.15}
PD_CIPRO = {'AUC_MIC_target': 125}
CIPRO_FRAZ_RENALE_100 = 0.145 / (0.145 + 0.167)   # quota renale di CLT a CLCR 100 (Forrest)


def prior_cipro(sesso, gg, peso, crea, crrt, Sc, Qef, altezza=None, egfr=None, modello=None, **_):
    eta = gg / 365.25
    note = []
    fuori = False
    if modello is None:
        modello = CIPRO_MODELLI[2] if eta >= 18 else CIPRO_MODELLI[0]
    if modello.startswith('Hirt'):
        h = HIRT
        fw = peso / 20.0
        eg = egfr if egfr else h['egfr_ref']
        CL = h['CL'] * fw ** 0.75 * (eg / h['egfr_ref']) ** h['egfr_exp'] * (max(eta, 0.01) / h['age_ref']) ** h['age_exp']
        par = {'CL': CL, 'V1': h['Vc'] * fw, 'Q': h['Q'] * fw ** 0.75, 'V2': h['Vp'] * fw}
        om = dict(h['omega'])
        sp, sa = h['sigma_prop'], 0.0
        nome = (f"Hirt D, et al. Eur J Clin Pharmacol 2021;77:1687-1695 [0,01–17,8 anni; "
                f"eGFR {eg:.0f} mL/min/1,73 m²]")
        if not egfr:
            note.append("eGFR non disponibile: usato il valore di riferimento del modello (199).")
        note.append("Variabilità individuale stimata solo su CL, come nel modello pubblicato.")
        if eta > 18:
            fuori = True
            note.append("Età oltre la popolazione del modello (0–18 anni).")
    elif modello.startswith('Forrest'):
        clcr = egfr if egfr else (renal.ckd_epi_2021(eta, crea, sesso) or 90.0)
        bsa = renal.bsa_mosteller(peso, altezza) or 1.73
        CL = (0.00145 * clcr + 0.167) * peso
        par = {'CL': CL, 'V1': FORREST['V1'] * peso, 'Q': FORREST['CLd_173'] * bsa / 1.73,
               'V2': FORREST['Vp'] * peso}
        om = dict(FORREST['omega'])
        sp, sa = FORREST['sigma_prop'], 0.0
        nome = (f"Forrest A, et al. Antimicrob Agents Chemother 1993;37:1065-1072 "
                f"[adulto acutamente malato; CLCR {clcr:.0f} mL/min/1,73 m²]")
        if eta < 18:
            fuori = True
            note.append("Modello dell'adulto applicato sotto i 18 anni.")
    else:
        if gg < 100:
            fuori = True
            note.append("Età < 3 mesi: fuori dalla popolazione di Rajagopalan (0,27–16,9 anni).")
        r = RAJAGOPALAN
        fw = peso / 70.0
        CL = r['CL'] * fw ** 0.75 * (1 + r['age_eff'] * (eta - 2.5))
        par = {'CL': CL, 'V1': r['Vc'] * fw, 'Q': r['Q'] * fw ** 0.75, 'V2': r['Vp'] * fw}
        if egfr is not None and egfr < 100:
            f = (1 - CIPRO_FRAZ_RENALE_100) + CIPRO_FRAZ_RENALE_100 * egfr / 100.0
            par['CL'] *= f
            note.append(f"eGFR {egfr:.0f} mL/min/1,73 m²: CL × {f:.2f} sulla quota renale "
                        f"(46% a CLCR 100, da Forrest 1993). Assunzione: il modello pediatrico "
                        f"non ha covariata renale.")
        om = dict(r['omega'])
        sp, sa = r['sigma_prop'], r['sigma_add']
        nome = "Rajagopalan P, Gastonguay MR. J Clin Pharmacol 2003;43:698-710 [0,27–16,9 anni]"
        note.append("Variabilità individuale per parametro assunta (abstract: 22,5–49,8%).")
    if crrt:
        CLEC = Sc * Qef
        par['CL'] += CLEC
        note.append(f"CRRT: CL extracorporea = Sc × Qef = {CLEC:.3f} L/h aggiunta.")
    return dict(par=par, omega=om, sigma_prop=sp, sigma_add=sa, t_inf=1.0, nome=nome,
                note=note, fuori_range=fuori)


# ============================================================
# DALBAVANCINA
# ============================================================

DALBA_MODELLI = {
    'Mariani — pediatrico, allometrico (Gaslini)': {
        'CL': 0.0305, 'V1': 4.81, 'Q': 0.0081, 'V2': 4.81, 'allometria': True,
        'omega': {'CL': 0.198}, 'sigma': 0.375,
        'nome': "Modello Mariani (IRCCS Gaslini, non pubblicato; 2 cmt, allometria, IIV su CL)",
    },
    'Cojutti 2021 — adulti, osteoarticolare': {
        'CL': 0.106, 'V1': 17.4, 'Q': 0.103, 'V2': 15.1, 'allometria': False,
        # Tabella 3 del paper: IIV CV% 36,21 / 44,27 / 130,46 / 62,34. Q resta fisso
        # (CV 130%, RSE 37%: non identificabile con campionamento sparso).
        'omega': {'CL': omega_da_cv(36.21), 'V1': omega_da_cv(44.27), 'V2': omega_da_cv(62.34)},
        'sigma': 0.18,
        'nome': ("Cojutti PG, et al. Antimicrob Agents Chemother 2021;65(5):e02260-20 "
                 "(adulti, osteoarticolare, n=15; nessuna covariata nel modello finale)"),
    },
    'Baiardi 2025 — adulti, trattamento prolungato, allometrico': {
        'CL': 0.0273, 'V1': 3.6, 'Q': 0.0225, 'V2': 6.4, 'allometria': True,
        'omega': {'CL': omega_da_cv(22.0), 'V1': omega_da_cv(17.3), 'Q': omega_da_cv(55.9),
                  'V2': omega_da_cv(30.1)},
        'sigma': 0.144,
        'nome': ("Baiardi G, et al. Antibiotics 2025;14:190 (adulti, n=30, ClCr>30 mL/min; "
                 "allometria, riferimento 70 kg assunto, mediana coorte 72 kg)"),
    },
}
DALBA_FU = 0.07
DALBA_MIC_OPZIONI = {
    'EUCAST v16.1 · 0,25 mg/L': 0.25,
    'Cojutti 2021 · 0,125 mg/L': 0.125,
}


def prior_dalba(peso, modello, **_):
    p = DALBA_MODELLI[modello]
    note = []
    if p['allometria']:
        fw = peso / 70.0
        par = {'CL': p['CL'] * fw ** 0.75, 'V1': p['V1'] * fw, 'Q': p['Q'] * fw ** 0.75, 'V2': p['V2'] * fw}
        note.append("Allometria: CL e Q × (peso/70)^0,75; V1 e V2 × peso/70.")
    else:
        par = {'CL': p['CL'], 'V1': p['V1'], 'Q': p['Q'], 'V2': p['V2']}
        note.append("Modello senza covariata di peso: valori di popolazione fissi.")
    return dict(par=par, omega=dict(p['omega']), sigma_prop=p['sigma'], sigma_add=0.0,
                t_inf=0.5, nome=p['nome'], note=note, fuori_range=False)


def dalba_soglia_auc(target_auc_mic, mic, fu):
    """Concentrazione totale media equivalente a fAUC24/MIC = target."""
    return (target_auc_mic * mic) / (24.0 * fu)


def dalba_soglia_4xmic(mic, fu):
    """100% fT > 4×MIC (Baiardi 2025): soglia totale = 4·MIC/fu."""
    return 4.0 * mic / fu


# ============================================================
# AMINOGLICOSIDI
# ============================================================
# Gentamicina — Llanos-Paez CC, et al. Antimicrob Agents Chemother
#   2017;61:e00205-17 (423 bambini oncologici, 0,2–18 anni; 2 cmt):
#   CL = 5,77·(GFRmat/100)·(Scr_attesa/Scr)^0,55
#   GFRmat = (FFM/70)^0,75 · PMA^3,33/(55,4^3,33 + PMA^3,33) · 112
#   V1 = 21,6·FFM/70; Q = 0,62·(FFM/70)^0,75; V2 = 13,8·FFM/70
#   BSV CL 16% (BOV 20,7%), V1 21,5%, V2 62,4%; residuo 27,5% + 0,04 mg/L
# Amikacina — Alhadab AA, et al. Antimicrob Agents Chemother 2018;62:e01781-17
#   (34 bambini oncologici, 1–18 anni; 2 cmt, allometria a 70 kg):
#   CL 11,1 L/h, V1 30,2 L, Q 4,26 L/h, V2 14,9 L; BSV CL 33%, V1 63,1%.
#   Residuo riportato "2,28% CV" (verosimilmente refuso): qui assunto 25% + 0,5 mg/L.
# Tobramicina FC — Hennig S, et al. Br J Clin Pharmacol 2008;65:502-510
#   (35 bambini con fibrosi cistica, 0,5–17,9 anni): CL 6,37 L/h/70 kg^0,75,
#   Vc 18,7 L/70 kg, Q 0,39 L/h, Vper 1,32 L; BSV CL 11,7% (BOV 6,5%),
#   Vc 11,7%, Vper 42%; residuo 19%.

AMINO_MODELLI = {
    'Gentamicina': ['Llanos-Paez 2017 — oncologia pediatrica'],
    'Amikacina': ['Alhadab 2018 — oncologia pediatrica',
                  'Llanos-Paez 2017 (gentamicina) — per analogia di classe'],
    'Tobramicina': ['Llanos-Paez 2017 (gentamicina) — per analogia di classe',
                    'Hennig 2008 — fibrosi cistica'],
}

# Soglie di default (PIGATO TDM, §4.4; tdm_data.json)
AMINO_SOGLIE = {
    # valle: finestra libera per il ridosaggio (PIGATO TDM §4.4); picco: intervallo pratico
    # pediatrico da tdm_data.json (genta/tobra 15–20 mg/L; amikacina 20–50 mg/L)
    'Gentamicina': {'valle': 1.0, 'cmax_mic': (8, 10), 'picco': (15, 20), 'tox_valle': 2.0},
    'Tobramicina': {'valle': 1.0, 'cmax_mic': (8, 10), 'picco': (15, 20), 'tox_valle': 2.0},
    'Amikacina':   {'valle': 5.0, 'cmax_mic': (8, 10), 'picco': (20, 50), 'tox_valle': 5.0},
}


def _llanos_paez(gg, peso, altezza, crea, sesso, eg_sett):
    eta = gg / 365.25
    ffm = renal.ffm(peso, altezza, eta, sesso)
    pma = renal.pma_settimane(gg, eg_sett)
    gfrmat = (ffm / 70.0) ** 0.75 * pma ** 3.33 / (55.4 ** 3.33 + pma ** 3.33) * 112.0
    scr_att = renal.creatinina_attesa_umol(eta, sesso)
    scr = crea * renal.MGDL_TO_UMOL
    CL = 5.77 * (gfrmat / 100.0) * (scr_att / scr) ** 0.55
    par = {'CL': CL, 'V1': 21.6 * ffm / 70.0, 'Q': 0.62 * (ffm / 70.0) ** 0.75, 'V2': 13.8 * ffm / 70.0}
    info = (f"FFM {ffm:.1f} kg (Al-Sallami); PMA {pma:.0f} sett.; GFR maturativo {gfrmat:.0f} mL/min; "
            f"creatinina attesa {scr_att:.0f} µmol/L vs misurata {scr:.0f} µmol/L")
    om = {'CL': omega_da_cv((16.0 ** 2 + 20.7 ** 2) ** 0.5), 'V1': omega_da_cv(21.5), 'V2': omega_da_cv(62.4)}
    return par, om, info


def prior_amino(farmaco, modello, sesso, gg, peso, altezza, crea, crrt, Sc, Qef, eg_sett=None, **_):
    eta = gg / 365.25
    note = []
    fuori = False
    if modello.startswith('Llanos-Paez'):
        if not altezza:
            raise ValueError("Il modello Llanos-Paez richiede l'altezza (massa magra).")
        par, om, info = _llanos_paez(gg, peso, altezza, crea, sesso, eg_sett)
        note.append(info)
        sp, sa = 0.275, 0.04
        nome = ("Llanos-Paez CC, et al. Antimicrob Agents Chemother 2017;61:e00205-17 "
                "[bambini oncologici 0,2–18 anni]")
        if farmaco != 'Gentamicina':
            nome += f" — applicato a {farmaco.lower()} per analogia di classe"
        if eta < 0.2 or eta > 18.5:
            fuori = True
            note.append("Età fuori dalla popolazione del modello (0,2–18 anni).")
    elif modello.startswith('Alhadab'):
        fw = peso / 70.0
        par = {'CL': 11.1 * fw ** 0.75, 'V1': 30.2 * fw, 'Q': 4.26 * fw ** 0.75, 'V2': 14.9 * fw}
        om = {'CL': omega_da_cv(33.0), 'V1': omega_da_cv(63.1)}
        sp, sa = 0.25, 0.5
        nome = ("Alhadab AA, et al. Antimicrob Agents Chemother 2018;62:e01781-17 "
                "[34 bambini oncologici 1–18 anni]")
        note.append("Nessuna covariata renale nel modello: con funzione renale ridotta il prior "
                    "sovrastima la CL e il livello sposta la stima. Residuo assunto (25% + 0,5 mg/L).")
        if eta < 1 or eta > 18.5:
            fuori = True
            note.append("Età fuori dalla popolazione del modello (1–18 anni).")
    else:  # Hennig 2008
        fw = peso / 70.0
        par = {'CL': 6.37 * fw ** 0.75, 'V1': 18.7 * fw, 'Q': 0.39, 'V2': 1.32}
        om = {'CL': omega_da_cv((11.7 ** 2 + 6.47 ** 2) ** 0.5), 'V1': omega_da_cv(11.66)}
        sp, sa = 0.19, 0.0
        nome = ("Hennig S, et al. Br J Clin Pharmacol 2008;65:502-510 "
                "[fibrosi cistica 0,5–17,9 anni]")
        note.append("Variabilità del modello molto bassa (pazienti FC stabili): il prior pesa molto; "
                    "fuori dalla FC usare il modello per analogia di classe.")
        if eta < 0.5 or eta > 18.5:
            fuori = True
            note.append("Età fuori dalla popolazione del modello (0,5–17,9 anni).")
    if crrt:
        CLEC = Sc * Qef
        par['CL'] += CLEC
        note.append(f"CRRT: CL extracorporea = Sc × Qef = {CLEC:.3f} L/h aggiunta.")
    return dict(par=par, omega=om, sigma_prop=sp, sigma_add=sa, t_inf=0.5, nome=nome,
                note=note, fuori_range=fuori)


VERIFICA = """
Verifica dei parametri sulle fonti primarie (26/10/2026)
- Lim 2021: CLren 22,6, CLnr 158, CL_HD 619 mL/h, Vss 5,9 L, BSV V1 11,2%, knr 8,78%,
  errore 8,5%: coincidono con la tabella del paper (PMC7887886).
- Ye 2022: fattori renali = rapporti di AUC della Tabella 4 (PMC9424659).
- Wei 2020, bande <2 anni: coerenti con Bradley 2014 (CL = dose/AUC).
- Soraluce 2018: Tabella 3 verificata sul PDF (θ renale 0,367; bootstrap 0,366).
- Wei 2020: Tabella 1 verificata sul PDF (CL per fascia). Il paper non riporta volumi:
  V 0,17 L/kg (<1 anno) e 0,11 L/kg (≥1 anno) restano valori assunti dalle versioni precedenti.
- Hirt 2021: Tabella 2 verificata sul PDF; modello pediatrico di default per la ciprofloxacina.
- Cojutti 2021: CL 0,106, V1 17,4, Q 0,103, V2 15,1; IIV 36/44/130/62%; residuo 18%.
- Baiardi 2025: CL 0,0273, V1 3,6, Q 0,0225, V2 6,4; IIV 22/17,3/55,9/30,1%; b 0,144.
- Forrest 1993: Tabella 2 e equazione di CLT (PMC187899).
- Rajagopalan 2003: modello strutturale dall'abstract; IIV per parametro da completare.
- Llanos-Paez 2017, Alhadab 2018, Hennig 2008: tabelle dei paper (PMC5527572,
  PMC5913936, PMC2291370).
- Dalbavancina, breakpoint EUCAST v16.1 Staphylococcus spp.: S ≤0,25 mg/L.
"""
