"""
Modulo TDM bayesiana: daptomicina, ciprofloxacina, dalbavancina, aminoglicosidi.

Flusso: scelta della molecola (vincolante per la sessione) → dati biometrici →
dosi e livelli → MAP bayesiano → risultati, proposta, referto → registrazione
esplicita nello storico (JSON per paziente × molecola).
"""

from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import streamlit as st

from . import charts, models, renal, storage, ui
from .pk_engine import conc, mappa_bayesiana, metriche_finestra, tempo_sotto_soglia

MODULO = 'bayes'
P = 'b'   # prefisso delle chiavi dei widget

FARMACI = {
    'Daptomicina': "Lipopeptide · Soraluce, Wei/Bradley, Ye, Lim",
    'Ciprofloxacina': "Fluorochinolone · Hirt, Rajagopalan, Forrest",
    'Dalbavancina': "Lipoglicopeptide · Mariani, Cojutti, Baiardi",
    'Gentamicina': "Aminoglicoside · Llanos-Paez",
    'Amikacina': "Aminoglicoside · Alhadab, Llanos-Paez",
    'Tobramicina': "Aminoglicoside · Llanos-Paez, Hennig",
}
AMINO = ('Gentamicina', 'Amikacina', 'Tobramicina')

WIDGETS = ['eta_unita', 'eta_val', 'peso', 'altezza', 'crea', 'eg', 'mic', 'crrt', 'sc', 'qef',
           'dialisi', 'modello', 'editor_dosi', 'editor_livelli', 'codice', 'sesso']


def _k(nome):
    return f"{P}_{nome}"


def _pulisci_widget():
    for w in WIDGETS:
        st.session_state.pop(_k(w), None)
    for k in list(st.session_state.keys()):
        if k.startswith(f"{P}_fc_") or k.startswith(f"{P}_am_"):
            st.session_state.pop(k, None)


def _widget_da_json(dati):
    _pulisci_widget()
    st.session_state[_k('codice')] = dati.get('codice_paziente', '') or ''
    st.session_state[_k('sesso')] = dati.get('sesso', 'Maschio') or 'Maschio'
    ult = dati['monitoraggi'][-1] if dati.get('monitoraggi') else None
    if not ult:
        return
    if ult.get('eta_giorni') is not None:
        st.session_state[_k('eta_unita')] = 'Anni'
        st.session_state[_k('eta_val')] = round(ult['eta_giorni'] / 365.25, 2)
    for campo, w in (('peso', 'peso'), ('altezza', 'altezza'), ('creatinina', 'crea'), ('mic', 'mic'),
                     ('eg_sett', 'eg')):
        if ult.get(campo) is not None:
            st.session_state[_k(w)] = float(ult[campo])
    validi = set(models.DALBA_MODELLI) | set(models.CIPRO_MODELLI) | {m for v in models.AMINO_MODELLI.values() for m in v}
    mod = ult.get('modello_scelto') or ult.get('modello_dalba')   # 'modello_dalba': file v6.0
    if mod in validi:
        st.session_state[_k('modello')] = mod


def _parse_dt(d, o):
    d, o = str(d).strip(), str(o).strip()
    if not d or not o or d.lower() in ('nan', 'none') or o.lower() in ('nan', 'none'):
        return None
    for fmt in ("%d/%m/%Y", "%Y-%m-%d"):
        try:
            dd = datetime.strptime(d, fmt).date()
            break
        except ValueError:
            dd = None
    if dd is None:
        return None
    try:
        tt = datetime.strptime(o[:5], "%H:%M").time()
    except ValueError:
        return None
    return datetime.combine(dd, tt)


def _tabella(righe, numeriche):
    """DataFrame per st.data_editor: colonne numeriche come float (celle vuote = NaN,
    non la scritta 'None'), testo come stringa."""
    df = pd.DataFrame(righe)
    for c in df.columns:
        if c in numeriche:
            df[c] = pd.to_numeric(df[c], errors='coerce').astype(float)
        else:
            df[c] = df[c].fillna("").astype(str)
    return df


def _num(v):
    try:
        f = float(v)
        return None if np.isnan(f) else f
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------

def _scelta_molecola():
    ui.passo("Passo 1")
    st.markdown("### Seleziona la molecola")
    st.caption("La scelta vale per la sessione: cambiando molecola i dati non scaricati si perdono.")
    nomi = list(FARMACI)
    for riga in (nomi[:3], nomi[3:]):
        cols = st.columns(3, vertical_alignment="bottom")
        for col, f in zip(cols, riga):
            with col, st.container(border=True):
                st.markdown(f"**{f}**")
                st.caption(FARMACI[f])
                if st.button(f"Seleziona {f.lower()}", key=f"{P}_sel_{f}", width="stretch",
                             type="primary"):
                    st.session_state[_k('farmaco')] = f
                    storage.reset(MODULO, f)
                    _pulisci_widget()
                    st.session_state.pop(_k('risultati'), None)
                    st.rerun()


def render():
    ui.intestazione("TDM bayesiana")
    st.caption(ui.AVVISO)
    farmaco = st.session_state.get(_k('farmaco'))
    if farmaco is None:
        _scelta_molecola()
        return

    st.info(f"Molecola in corso: **{farmaco}**")
    if st.sidebar.button("Cambia molecola", width="stretch", key=f"{P}_cambia"):
        _pulisci_widget()
        st.session_state.pop(_k('farmaco'), None)
        st.session_state.pop(_k('risultati'), None)
        storage.reset(MODULO)
        st.rerun()

    sto, codice, sesso = storage.sidebar_storico(MODULO, farmaco, P, _widget_da_json)
    if st.sidebar.button("Nuovo paziente (stessa molecola)", width="stretch", key=f"{P}_nuovo"):
        _pulisci_widget()
        storage.reset(MODULO, farmaco)
        st.session_state.pop(_k('risultati'), None)
        st.session_state.pop(f'{P}_ultimo_file', None)
        st.session_state.pop(f'{P}_msg', None)
        st.rerun()
    if not codice:
        st.warning("Inserisci un codice paziente nella barra laterale per procedere.")
        return

    ult = sto['monitoraggi'][-1] if sto.get('monitoraggi') else None

    # ---------------- Passo 2 ----------------
    ui.passo("Passo 2")
    st.markdown("### Dati biometrici e clinici")
    with st.container(border=True):
        a, b, c = st.columns(3, vertical_alignment="bottom")
        unita = a.selectbox("Unità età", ["Anni", "Mesi", "Giorni"], key=_k('eta_unita'))
        maxu = {"Anni": 100.0, "Mesi": 1200.0, "Giorni": 36500.0}[unita]
        eta_val = b.number_input(f"Età ({unita.lower()})", 0.0, maxu, value=None, key=_k('eta_val'))
        peso = c.number_input("Peso (kg)", 0.3, 300.0, value=None, format="%.1f", key=_k('peso'))
        gg = None
        if eta_val is not None:
            gg = int(round(eta_val * {"Anni": 365.25, "Mesi": 30.44, "Giorni": 1}[unita]))
        eta = gg / 365.25 if gg is not None else None
        a, b, c = st.columns(3, vertical_alignment="bottom")
        altezza = a.number_input("Altezza (cm)", 25.0, 250.0, value=None, format="%.0f", key=_k('altezza'),
                                 help="Lunghezza nel lattante.")
        crea = b.number_input("Creatinina (mg/dL)", 0.05, 25.0, value=None, step=0.01, format="%.2f",
                              key=_k('crea'))
        if farmaco == 'Dalbavancina':
            mic = c.number_input("MIC di riferimento (mg/L)", 0.001, 4.0, value=None, step=0.001,
                                 format="%.4f", key=_k('mic'))
        else:
            mic = c.number_input("MIC (mg/L)", 0.001, 64.0, value=None, step=0.001,
                                 format="%.4f", key=_k('mic'))

        a, b, c = st.columns(3, vertical_alignment="bottom")
        eg_sett = None
        if eta is not None and eta < 2:
            eg_sett = a.number_input("EG alla nascita (sett.)", 22.0, 42.0, value=None, key=_k('eg'),
                                     help="Per età postmestruale (aminoglicosidi) e formula del pretermine. "
                                          "Vuoto = 40 settimane.")
        dialisi = "Nessuna"
        if farmaco == 'Daptomicina' and gg is not None and gg < 6570:
            opz = list(models.IRC_OPZIONI)
            if 4380 <= gg < 6570:
                opz += ["Emodialisi (HD)", "Dialisi Peritoneale (CCPD)"]
            dialisi = b.selectbox("Funzione renale / dialisi", opz, key=_k('dialisi'),
                                  help="Correzione Ye 2022 per fascia di IRC.")
        crrt = c.checkbox("CRRT", key=_k('crrt'), help="Terapia sostitutiva renale continua in corso")
        Sc, Qef = 0.2, 2.0
        if crrt:
            a, b = st.columns(2, vertical_alignment="bottom")
            Sc = a.number_input("Sieving coefficient", 0.01, 1.0, 0.2 if farmaco not in AMINO else 0.85, 0.01,
                                key=_k('sc'))
            Qef = b.number_input("Flusso effluente (L/h)", 0.1, 6.0, 2.0, 0.1, key=_k('qef'))
        modello = None
        if farmaco == 'Dalbavancina':
            modello = st.selectbox("Modello PK", list(models.DALBA_MODELLI), key=_k('modello'))
        elif farmaco in AMINO:
            modello = st.selectbox("Modello PK", models.AMINO_MODELLI[farmaco], key=_k('modello'))
        elif farmaco == 'Ciprofloxacina':
            idx = 2 if (eta is not None and eta >= 18) else 0
            modello = st.selectbox("Modello PK", models.CIPRO_MODELLI, index=idx, key=_k('modello'),
                                   help="Hirt 2021 ha eGFR ed età come covariate; Rajagopalan 2003 stima CL "
                                        "più alte (~0,7 vs ~0,45 L/h/kg a 8 anni): con i livelli il MAP corregge.")

    eg = None
    if eta is not None and crea:
        eg = renal.egfr_scelta(eta, sesso, peso, altezza, crea, eg_sett)
        with st.container(border=True):
            st.markdown("**Funzione renale stimata**")
            st.caption(eg['nota'])
            ui.schede(ui.schede_egfr(eg))

    mancano = [n for n, v in (("età", gg), ("peso", peso), ("MIC", mic)) if v is None]
    if farmaco != 'Dalbavancina' and crea is None:
        mancano.append("creatinina")
    if farmaco in AMINO and modello and modello.startswith('Llanos') and not altezza:
        mancano.append("altezza (massa magra)")
    if farmaco == 'Ciprofloxacina' and modello and modello.startswith('Forrest') and not altezza:
        mancano.append("altezza (BSA)")
    if mancano:
        st.warning("Campi da compilare: " + ", ".join(mancano) + ".")

    # ---------------- Passo 3 ----------------
    ui.passo("Passo 3")
    st.markdown("### Somministrazioni e livelli")
    t_inf_def = 0.5 if farmaco not in ('Ciprofloxacina',) else 1.0
    if farmaco in AMINO:
        st.caption("Aminoglicosidi: inserisci il livello con data e ora reali del prelievo (es. 10 h "
                   "dall'inizio dell'infusione). Un solo livello basta: con un campione la stima "
                   "corregge soprattutto la clearance, il volume resta vicino al prior.")
    t0 = datetime.now()
    t_ult_storico = -1
    dosi_rows, liv_rows = [], []
    if ult and ult.get('storia_dosi'):
        t0 = datetime.strptime(ult['data_prima_dose'], '%Y-%m-%d %H:%M')
        for s in ult['storia_dosi']:
            dt = t0 + timedelta(hours=s['ore'])
            dosi_rows.append({"Data": dt.strftime("%d/%m/%Y"), "Ora": dt.strftime("%H:%M"),
                              "Dose (mg)": s['dose'], "Inf. (h)": s.get('tinf', t_inf_def)})
        t_ult_storico = ult['storia_dosi'][-1]['ore']
        for l in ult.get('livelli', []):
            dt = t0 + timedelta(hours=l['ore'])
            liv_rows.append({"Data": dt.strftime("%d/%m/%Y"), "Ora": dt.strftime("%H:%M"),
                             "mg/L": l['conc'], "Nota": l.get('nota', '')})
    else:
        dosi_rows.append({"Data": "", "Ora": "", "Dose (mg)": None, "Inf. (h)": t_inf_def})
    liv_rows.append({"Data": "", "Ora": "", "mg/L": None, "Nota": ""})

    ca, cb = st.columns(2, vertical_alignment="bottom")
    with ca:
        st.markdown("**Somministrazioni**")
        df_d = st.data_editor(_tabella(dosi_rows, ["Dose (mg)", "Inf. (h)"]), num_rows="dynamic", width="stretch",
                              height=280, key=_k('editor_dosi'),
                              column_config={"Dose (mg)": st.column_config.NumberColumn(format="%.1f", min_value=0.1),
                                             "Inf. (h)": st.column_config.NumberColumn(format="%.2f", min_value=0.05,
                                                                                            max_value=24.0)})
    with cb:
        st.markdown("**Livelli plasmatici**")
        df_l = st.data_editor(_tabella(liv_rows, ["mg/L"]), num_rows="dynamic", width="stretch",
                              height=280, key=_k('editor_livelli'),
                              column_config={"mg/L": st.column_config.NumberColumn(format="%.3f", min_value=0.0)})
    st.caption("Date nel formato GG/MM/AAAA, ore HH:MM. Le righe già registrate restano modificabili: "
               "correggile solo per errori di trascrizione.")

    storia, livelli = [], []
    for _, r in df_d.iterrows():
        dt = _parse_dt(r.get("Data"), r.get("Ora"))
        dose = _num(r.get("Dose (mg)"))
        tinf = _num(r.get("Inf. (h)")) or t_inf_def
        if dt and dose:
            storia.append({'dt': dt, 'dose': dose, 'tinf': tinf})
    for _, r in df_l.iterrows():
        dt = _parse_dt(r.get("Data"), r.get("Ora"))
        cc = _num(r.get("mg/L"))
        if dt and cc is not None:
            livelli.append({'dt': dt, 'conc': cc, 'nota': str(r.get("Nota") or "")})

    if st.button("Avvia il calcolo bayesiano", type="primary", width="stretch", key=f"{P}_calcola"):
        if mancano:
            st.error("Compila prima: " + ", ".join(mancano) + ".")
        elif not storia:
            st.error("Nessuna somministrazione valida in tabella.")
        elif not livelli:
            st.error("Inserisci almeno un livello plasmatico.")
        else:
            _calcola(farmaco, codice, sesso, gg, eta, peso, altezza, crea, eg_sett, eg, mic, crrt, Sc, Qef,
                     dialisi, modello, storia, livelli, t_ult_storico)

    R = st.session_state.get(_k('risultati'))
    if R and R['farmaco'] == farmaco and R['codice'] == codice:
        _mostra_risultati(R, sto)

    with st.expander(f"Storico di sessione — {farmaco}"):
        mon = sto.get('monitoraggi', [])
        if not mon:
            st.info("Nessuna valutazione registrata.")
        else:
            df = pd.DataFrame(mon)
            cols = [c for c in ['data_val', 'creatinina', 'egfr_valore', 'mic', 'dose_mg', 'n_dosi',
                                'n_livelli', 'CL_ind', 'AUC24'] if c in df.columns]
            st.dataframe(df[cols], width="stretch")
            sel = st.selectbox("Rimuovi un record errato", ["—"] + [f"#{i} — {r.get('data_val')}" for i, r in enumerate(mon)],
                               key=f"{P}_rimuovi_sel")
            if sel != "—" and st.button("Rimuovi record", key=f"{P}_rimuovi"):
                mon.pop(int(sel.split('#')[1].split(' ')[0]))
                st.rerun()

    with st.expander("Modelli e verifica delle fonti"):
        st.code(models.VERIFICA.strip(), language='text', wrap_lines=True)


def _prior(farmaco, sesso, gg, peso, altezza, crea, eg_sett, eg, crrt, Sc, Qef, dialisi, modello):
    if farmaco == 'Daptomicina':
        return models.prior_dapto(sesso, gg, peso, crea, crrt, Sc, Qef, dialisi)
    if farmaco == 'Ciprofloxacina':
        egv = eg['valore'] if eg else None
        return models.prior_cipro(sesso, gg, peso, crea, crrt, Sc, Qef, altezza=altezza, egfr=egv, modello=modello)
    if farmaco == 'Dalbavancina':
        return models.prior_dalba(peso, modello)
    return models.prior_amino(farmaco, modello, sesso, gg, peso, altezza, crea, crrt, Sc, Qef, eg_sett=eg_sett)


def _calcola(farmaco, codice, sesso, gg, eta, peso, altezza, crea, eg_sett, eg, mic, crrt, Sc, Qef,
             dialisi, modello, storia, livelli, t_ult_storico):
    t0 = min(s['dt'] for s in storia)
    dt_val = max(l['dt'] for l in livelli).strftime('%Y-%m-%d %H:%M')
    for s in storia:
        s['ore'] = (s['dt'] - t0).total_seconds() / 3600
    for l in livelli:
        l['ore'] = (l['dt'] - t0).total_seconds() / 3600
    storia = sorted([{'ore': s['ore'], 'dose': float(s['dose']), 'tinf': float(s['tinf'])} for s in storia],
                    key=lambda x: x['ore'])
    livelli = sorted([{'ore': l['ore'], 'conc': float(l['conc']), 'nota': l['nota']} for l in livelli],
                     key=lambda x: x['ore'])
    if any(l['ore'] < 0 for l in livelli):
        st.error("Un livello precede la prima dose: controlla date e ore.")
        return
    try:
        pr = _prior(farmaco, sesso, gg, peso, altezza, crea, eg_sett, eg, crrt, Sc, Qef, dialisi, modello)
    except ValueError as exc:
        st.error(str(exc))
        return
    om = dict(pr['omega'])
    n = len(livelli)
    # Q e V2 non identificabili con pochi campioni: ancorati al prior (come v5–v6)
    for k in ('Q', 'V2'):
        if om.get(k, 0) > 0:
            om[k] = 0.10 + (om[k] - 0.10) * min(1.0, n / 6.0)
    with st.spinner("Stima MAP in corso…"):
        par, descr = mappa_bayesiana(pr['par'], om, storia, livelli, pr['sigma_prop'], pr['sigma_add'], pr['t_inf'])
    t_last = storia[-1]['ore']
    taus = np.diff([s['ore'] for s in storia])
    tau_med = float(np.median(taus)) if len(taus) else None
    # AUC24 nelle 24 h dopo l'ultima dose, proseguendo il regime in atto
    # (stessa dose, intervallo mediano): con intervalli < 24 h la finestra
    # deve contenere anche le dosi successive.
    proiez = list(storia)
    if tau_med and tau_med < 24:
        k = 1
        while t_last + k * tau_med < t_last + 24 - 1e-6:
            proiez.append({'ore': t_last + k * tau_med, 'dose': storia[-1]['dose'], 'tinf': storia[-1]['tinf']})
            k += 1
    auc24, cmax24, cmin24 = metriche_finestra(par, proiez, t_last, 24.0, pr['t_inf'])
    vss = par['V1'] + par.get('V2', 0)
    # emivita terminale
    if par.get('Q') and par.get('V2'):
        k10, k12, k21 = par['CL'] / par['V1'], par['Q'] / par['V1'], par['Q'] / par['V2']
        s_ = k10 + k12 + k21
        beta = (s_ - np.sqrt(max(s_ ** 2 - 4 * k10 * k21, 0))) / 2
        thalf = float(np.log(2) / beta)
    else:
        thalf = float(np.log(2) * par['V1'] / par['CL'])
    st.session_state[_k('risultati')] = dict(
        farmaco=farmaco, codice=codice, sesso=sesso, gg=gg, eta=eta, peso=peso, altezza=altezza, crea=crea,
        eg_sett=eg_sett, egfr_valore=eg['valore'] if eg else None, egfr_formula=eg.get('etichetta') if eg else None,
        egfr_tutte=eg['tutte'] if eg else {}, mic=mic, crrt=crrt, Sc=Sc, Qef=Qef, dialisi=dialisi, modello=modello,
        storia=storia, livelli=livelli, t0=t0.strftime('%Y-%m-%d %H:%M'), dt_val=dt_val,
        t_ult_storico=t_ult_storico, prior=pr, par=par, descr=descr, auc24=auc24, cmax24=cmax24, cmin24=cmin24,
        vss=vss, thalf=thalf, tau_med=tau_med, registrato=False)


# ------------------------------------------------------------------

def _mostra_risultati(R, sto):
    farmaco, pr, par = R['farmaco'], R['prior'], R['par']
    storia, livelli, mic = R['storia'], R['livelli'], R['mic']
    t0 = datetime.strptime(R['t0'], '%Y-%m-%d %H:%M')
    st.markdown("---")
    ui.passo("Risultati")
    st.markdown("### Parametri individuali")
    st.caption(f"Modello: {pr['nome']}")
    st.caption(R['descr'])
    for nota in pr['note']:
        st.info(nota)
    if pr['fuori_range']:
        st.warning("Paziente fuori dalla popolazione del modello: la stima va letta con cautela.")

    ui.schede([
        {'k': "Clearance", 'v': ui.n(par['CL'], 3), 'u': "L/h", 's': f"prior {ui.n(pr['par']['CL'], 3)} L/h"},
        {'k': "Volume (Vss)", 'v': ui.n(R['vss'], 2), 'u': "L",
         's': f"prior {ui.n(pr['par']['V1'] + pr['par'].get('V2', 0), 2)} L"},
        {'k': "Emivita terminale", 'v': ui.n(R['thalf'], 1), 'u': "h"},
        {'k': "AUC 24 h dopo l'ultima dose", 'v': ui.n(R['auc24'], 0), 'u': "mg·h/L"},
    ])

    soglie = []
    if farmaco == 'Daptomicina':
        soglie = [(models.PD_DAPTO['Cmin_tox'], "Cmin 24,3 mg/L (CPK, Bhavnani 2010)", charts.DANGER)]
    fig = charts.profilo(farmaco, storia, livelli, par, pr['par'], pr['t_inf'], R['t_ult_storico'], soglie)
    st.plotly_chart(fig, width="stretch", theme=None)

    blocco = ""
    if farmaco == 'Daptomicina':
        blocco = _pannello_dapto(R)
    elif farmaco == 'Ciprofloxacina':
        blocco = _pannello_cipro(R)
    elif farmaco == 'Dalbavancina':
        blocco = _pannello_dalba(R, t0)
    else:
        blocco = _pannello_amino(R, t0)

    eg_txt = (f"{ui.n(R['egfr_valore'], 0)} mL/min/1,73 m² ({R['egfr_formula']})"
              if R.get('egfr_valore') else "N/D")
    tau_txt = f"{ui.n(R['tau_med'], 1)} h" if R['tau_med'] else "una sola dose"
    referto = (
        f"TDM BAYESIANA — {farmaco.upper()}\n"
        f"Data: {R['dt_val']} | Codice: {R['codice']} | Sesso: {R['sesso']} | Età: {ui.n(R['eta'], 1)} anni | "
        f"Peso: {ui.n(R['peso'], 1)} kg | Altezza: {ui.n(R['altezza'], 0)} cm\n\n"
        f"Modello: {pr['nome']}\n"
        f"Creatinina: {ui.n(R['crea'], 2)} mg/dL | eGFR: {eg_txt}"
        + (f" | CRRT Sc {ui.n(R['Sc'], 2)}, Qef {ui.n(R['Qef'], 1)} L/h" if R['crrt'] else "") + "\n"
        f"Dosi: {len(storia)} (ultima {ui.n(storia[-1]['dose'], 0)} mg = {ui.n(storia[-1]['dose'] / R['peso'], 1)} mg/kg; "
        f"intervallo mediano {tau_txt}) | Livelli: {len(livelli)} | MIC: {ui.n(mic, 4)} mg/L\n\n"
        f"Parametri individuali (MAP): CL {ui.n(par['CL'], 3)} L/h; Vss {ui.n(R['vss'], 2)} L; "
        f"t½ {ui.n(R['thalf'], 1)} h\n"
        f"{blocco}"
    )
    st.markdown("### Referto")
    st.code(referto, language='text', wrap_lines=True)

    st.markdown("---")
    if R.get('registrato'):
        st.success("Valutazione registrata nello storico di sessione. Scarica il JSON dalla barra laterale.")
    else:
        st.info("La valutazione non è ancora nello storico. Registrala solo se corrisponde a un'occasione "
                "clinica reale.")
        if st.button("Registra questa valutazione nello storico", type="primary", width="stretch",
                     key=f"{P}_registra"):
            rec = {
                'data_val': R['dt_val'], 'farmaco': farmaco, 'eta_giorni': R['gg'], 'peso': R['peso'],
                'altezza': R['altezza'], 'creatinina': R['crea'], 'eg_sett': R['eg_sett'],
                'egfr_valore': R['egfr_valore'], 'egfr_formula': R['egfr_formula'], 'egfr_tutte': R['egfr_tutte'],
                'crrt': bool(R['crrt']), 'Sc': R['Sc'], 'Qef': R['Qef'], 'dialisi_tipo': R['dialisi'],
                'mic': mic, 'modello_scelto': R['modello'], 'modello': pr['nome'],
                'dose_mg': storia[-1]['dose'], 'tau_h': R['tau_med'], 'n_dosi': len(storia), 'n_livelli': len(livelli),
                'CL_ind': par['CL'], 'V1_ind': par['V1'], 'Q_ind': par.get('Q'), 'V2_ind': par.get('V2'),
                'Vss': R['vss'], 'AUC24': R['auc24'], 'Cmin24': R['cmin24'], 'referto': referto,
                'storia_dosi': storia, 'livelli': livelli, 'data_prima_dose': R['t0'],
            }
            sto['monitoraggi'].append(rec)
            R['registrato'] = True
            st.rerun()


def _pannello_dapto(R):
    mic, auc = R['mic'], R['auc24']
    ratio = auc / mic
    pd_ = models.PD_DAPTO
    tox = R['cmin24'] >= pd_['Cmin_tox']
    ok_s = ratio >= pd_['AUC_MIC_staph']
    ok_e = ratio >= pd_['AUC_MIC_entero']
    ui.schede([
        {'k': "AUC24/MIC · S. aureus (≥ 666)", 'v': ui.n(ratio, 0),
         's': "target raggiunto" if ok_s else "sotto target", 'stato': 'ok' if ok_s else 'ko'},
        {'k': "AUC24/MIC · E. faecium (≥ 143)", 'v': ui.n(ratio, 0),
         's': "target raggiunto" if ok_e else "sotto target", 'stato': 'ok' if ok_e else 'ko'},
        {'k': "Cmin prevista (soglia CPK 24,3)", 'v': ui.n(R['cmin24'], 1), 'u': "mg/L",
         's': "rischio CPK" if tox else "sotto soglia", 'stato': 'ko' if tox else 'ok'},
    ])
    dose = R['storia'][-1]['dose']
    tau = R['tau_med'] or 24
    # Finestra tdm_data.json: AUC24/MIC ≥ 666 (S. aureus); AUC24 > 939 mg·h/L e
    # Cmin ≥ 24,3 mg/L come limiti di sicurezza. PK lineare: dose ∝ AUC e Cmin.
    f_min = pd_['AUC_MIC_staph'] * mic / auc
    f_max = min(939.0 / auc, 0.9 * pd_['Cmin_tox'] / max(R['cmin24'], 1e-6))
    if f_min > 1:
        azione, f = "AUMENTA", f_min
    elif f_max < 1:
        azione, f = "RIDUCI", max(f_max, f_min)
    else:
        azione, f = "MANTIENI", 1.0
    if azione != "MANTIENI" and abs(f - 1) < 0.10 and f_min <= 1.0 / 0.9:
        azione, f = "MANTIENI", 1.0   # correzioni < 10% non cliniche
    dose_new = dose * f
    if f_min > f_max:
        st.warning("Con questa MIC il target AUC24/MIC 666 non è raggiungibile senza superare AUC 939 mg·h/L "
                   "o Cmin 24,3 mg/L: valutare alternativa terapeutica.")
    testo = (f"{azione}: {ui.n(dose_new, 0)} mg ({ui.n(dose_new / R['peso'], 1)} mg/kg) ogni {ui.n(tau, 0)} h "
             f"(AUC24 prevista {ui.n(auc * f, 0)} mg·h/L, Cmin {ui.n(R['cmin24'] * f, 1)} mg/L)")
    (st.success if azione == "MANTIENI" else st.info)("Proposta — " + testo)
    return (f"AUC24/MIC {ui.n(ratio, 0)} (target ≥666 S. aureus; ≥143 E. faecium); AUC24 {ui.n(auc, 0)} mg·h/L "
            f"(limite 939); Cmin prevista {ui.n(R['cmin24'], 1)} mg/L (soglia CPK 24,3)\n"
            f"Proposta: {testo}\n")


def _pannello_cipro(R):
    mic, auc = R['mic'], R['auc24']
    ratio = auc / mic
    ok = ratio >= models.PD_CIPRO['AUC_MIC_target']
    ui.schede([
        {'k': "AUC24/MIC (≥ 125)", 'v': ui.n(ratio, 0), 's': "target raggiunto" if ok else "sotto target",
         'stato': 'ok' if ok else 'ko'},
        {'k': "Cmax/MIC (≥ 10)", 'v': ui.n(R['cmax24'] / mic, 1), 's': f"Cmax {ui.n(R['cmax24'], 2)} mg/L"},
        {'k': "Cmin (tenere < 2)", 'v': ui.n(R['cmin24'], 2), 'u': "mg/L",
         'stato': 'ko' if R['cmin24'] >= 2 else 'ok', 's': "oltre 2 mg/L" if R['cmin24'] >= 2 else "sotto 2 mg/L"},
    ])
    dose = R['storia'][-1]['dose']
    tau = R['tau_med'] or 24
    if not ok:
        f, azione = 125 * mic / auc, "AUMENTA"
    elif R['cmin24'] > 2:
        f, azione = 0.75, "RIDUCI"
    else:
        f, azione = 1.0, "MANTIENI"
    testo = (f"{azione}: {ui.n(dose * f, 0)} mg ({ui.n(dose * f / R['peso'], 1)} mg/kg) ogni {ui.n(tau, 0)} h "
             f"(AUC24/MIC prevista {ui.n(ratio * f, 0)})")
    (st.success if azione == "MANTIENI" else st.info)("Proposta — " + testo)
    st.caption("Nessuna soglia AUC di neurotossicità validata per i fluorochinoloni.")
    return (f"AUC24/MIC {ui.n(ratio, 0)} (target ≥125); Cmax {ui.n(R['cmax24'], 2)} mg/L "
            f"(Cmax/MIC {ui.n(R['cmax24'] / mic, 1)}); Cmin {ui.n(R['cmin24'], 2)} mg/L\n"
            f"Proposta: {testo}\n")


def _pannello_dalba(R, t0):
    pr, par, storia, livelli = R['prior'], R['par'], R['storia'], R['livelli']
    with st.container(border=True):
        st.markdown("#### Forecast — prossima dose")
        a, b, c = st.columns(3, vertical_alignment="bottom")
        target = a.selectbox("Target", ["fAUC24/MIC ≥ 111,1 (Cojutti 2021)", "100% fT > 4×MIC (Baiardi 2025)"],
                             key=f"{P}_fc_target")
        mic_lab = b.selectbox("MIC", list(models.DALBA_MIC_OPZIONI) + ["MIC del paziente"], key=f"{P}_fc_mic")
        mic_fc = R['mic'] if mic_lab == "MIC del paziente" else models.DALBA_MIC_OPZIONI[mic_lab]
        fu = c.number_input("Frazione libera", 0.01, 0.2, models.DALBA_FU, 0.01, key=f"{P}_fc_fu")
        if target.startswith("100%"):
            soglia = models.dalba_soglia_4xmic(mic_fc, fu)
        else:
            soglia = models.dalba_soglia_auc(111.1, mic_fc, fu)
        st.info(f"Soglia di concentrazione totale: {ui.n(soglia, 2)} mg/L (MIC {ui.n(mic_fc, 4)}, fu {ui.n(fu, 2)}).")
        t_val = max(livelli[-1]['ore'], storia[-1]['ore'])
        c_att = float(conc(par, storia, t_val, pr['t_inf'])[0])
        t_x = tempo_sotto_soglia(par, storia, soglia, t_val, t_val + 2400, passo=2.0, t_inf=pr['t_inf'])
        if t_x is None:
            esito = "soglia non raggiunta entro 100 giorni"
            st.warning(esito)
        elif t_x <= t_val + 0.5:
            esito = f"concentrazione già sotto soglia ({ui.n(c_att, 2)} mg/L)"
            st.error(esito)
        else:
            data_x = t0 + timedelta(hours=t_x)
            esito = f"sotto soglia il {data_x.strftime('%d/%m/%Y %H:%M')} (fra {ui.n((t_x - t_val) / 24, 1)} giorni)"
            st.success(f"La concentrazione scende {esito}.")
        sim = st.checkbox("Simula una dose aggiuntiva", key=f"{P}_fc_sim")
        if sim:
            a, b, c = st.columns(3, vertical_alignment="bottom")
            dnew = a.date_input("Data", value=(t0 + timedelta(hours=t_val + 24)).date(), format="DD/MM/YYYY",
                                key=f"{P}_fc_data")
            onew = b.time_input("Ora", value=datetime.strptime("10:00", "%H:%M").time(), key=f"{P}_fc_ora")
            mg = c.number_input("Dose (mg)", 50.0, 1500.0, 1500.0, 50.0, key=f"{P}_fc_dose")
            t_new = (datetime.combine(dnew, onew) - t0).total_seconds() / 3600
            s2 = sorted(storia + [{'ore': t_new, 'dose': mg, 'tinf': pr['t_inf']}], key=lambda x: x['ore'])
            tx2 = tempo_sotto_soglia(par, s2, soglia, t_new, t_new + 2400, passo=2.0, t_inf=pr['t_inf'])
            if tx2:
                st.success(f"Con {ui.n(mg, 0)} mg la concentrazione resta sopra soglia fino al "
                           f"{(t0 + timedelta(hours=tx2)).strftime('%d/%m/%Y')}.")
            tg = np.linspace(0, (tx2 or t_new + 500) + 48, 1000)
            st.plotly_chart(charts.forecast("Dalbavancina", tg, conc(par, s2, tg, pr['t_inf']), soglia, s2, [t_new]),
                            width="stretch", theme=None)
    return f"Forecast ({target}, MIC {ui.n(mic_fc, 4)}): soglia {ui.n(soglia, 2)} mg/L, {esito}\n"


def _simula_regime(par, storia, dose, tau, tinf, t_start, n_dosi=5):
    """Prosegue con n dosi uguali da t_start; metriche sull'ultimo intervallo."""
    s2 = list(storia) + [{'ore': t_start + i * tau, 'dose': dose, 'tinf': tinf} for i in range(n_dosi)]
    tl = t_start + (n_dosi - 1) * tau
    picco = float(conc(par, s2, tl + tinf + 0.5, tinf)[0])
    valle = float(conc(par, s2, tl + tau - 0.01, tinf)[0])
    auc, _, _ = metriche_finestra(par, s2, tl, tau, tinf)
    return picco, valle, auc * 24.0 / tau


def _arrotonda_amino(mg, farmaco):
    passo = 10 if farmaco == 'Amikacina' else (5 if mg >= 50 else 1)
    return max(passo, round(mg / passo) * passo)


def _pannello_amino(R, t0):
    farmaco, par, storia = R['farmaco'], R['par'], R['storia']
    soglie = models.AMINO_SOGLIE[farmaco]
    mic = R['mic']
    tinf = storia[-1].get('tinf', 0.5)
    with st.container(border=True):
        st.markdown("#### Previsioni e proposta di regime")
        a, b, c = st.columns(3, vertical_alignment="bottom")
        obiettivo = a.number_input("Obiettivo picco/MIC", 4.0, 20.0, float(soglie['cmax_mic'][1]), 1.0,
                                   key=f"{P}_am_obj")
        p_lo = b.number_input("Picco minimo (mg/L)", 1.0, 120.0, float(soglie['picco'][0]), 1.0, key=f"{P}_am_plo")
        p_hi = c.number_input("Picco massimo (mg/L)", 1.0, 150.0, float(soglie['picco'][1]), 1.0, key=f"{P}_am_phi",
                              help="Intervallo pratico pediatrico (tdm_data.json); FC tobramicina 30–40 mg/L.")
        d, e = st.columns(2, vertical_alignment="bottom")
        valle_max = d.number_input("Valle massima (mg/L)", 0.1, 10.0, float(soglie['valle']), 0.1,
                                   key=f"{P}_am_valle",
                                   help="Default PIGATO: finestra libera 1 mg/L (genta/tobra), 5 mg/L amikacina. "
                                        "In pediatria i target di valle sono più bassi (≤0,5 / ≤2,5).")
        tau_att = e.number_input("Intervallo attuale (h)", 6.0, 96.0, float(round(R['tau_med'] or 24)), 6.0,
                                 key=f"{P}_am_tau")
        lo_eff = max(obiettivo * mic, p_lo)
        hi_eff = max(p_hi, lo_eff)
        st.caption(f"Picco obiettivo: {ui.n(lo_eff, 1)}–{ui.n(hi_eff, 1)} mg/L (il più alto tra picco minimo e "
                   f"{ui.n(obiettivo, 0)}×MIC).")
        t_last = storia[-1]['ore']
        t_lib = tempo_sotto_soglia(par, storia, valle_max, t_last + tinf, t_last + 240, 0.1, tinf)
        dose_att = storia[-1]['dose']
        pk_att = _simula_regime(par, storia, dose_att, tau_att, tinf, t_last + tau_att)
        cmax_ult = float(conc(par, storia, t_last + tinf + 0.5, tinf)[0])
        ui.schede([
            {'k': "Picco dopo l'ultima dose", 'v': ui.n(cmax_ult, 1), 'u': "mg/L",
             's': "30 min dopo la fine dell'infusione"},
            {'k': "Picco/MIC", 'v': ui.n(cmax_ult / mic, 1), 'stato': 'ok' if cmax_ult / mic >= obiettivo else 'ko',
             's': f"obiettivo ≥ {ui.n(obiettivo, 0)}"},
            {'k': f"Sotto {ui.n(valle_max, 1)} mg/L dopo", 'v': ui.n(t_lib - t_last, 1) if t_lib else "> 240", 'u': "h",
             's': (t0 + timedelta(hours=t_lib)).strftime('%d/%m/%Y %H:%M') if t_lib else None},
            {'k': "AUC24 allo stato stazionario", 'v': ui.n(pk_att[2], 0), 'u': "mg·h/L", 's': "regime attuale"},
        ])

        if lo_eff <= pk_att[0] <= hi_eff:
            dose_obj = dose_att
        else:
            mid = (lo_eff + hi_eff) / 2 if hi_eff > lo_eff else lo_eff * 1.1
            dose_obj = _arrotonda_amino(dose_att * mid / max(pk_att[0], 1e-6), farmaco)
        righe, scelta = [], None
        for tau in (24, 36, 48, 72):
            for dose in sorted({dose_att, dose_obj}, key=lambda x: (x != dose_obj, x)):
                pk = _simula_regime(par, storia, dose, tau, tinf, t_last + max(tau_att, (t_lib - t_last) if t_lib else tau_att))
                ok = pk[1] < valle_max and lo_eff * 0.95 <= pk[0] <= hi_eff * 1.05
                righe.append({"Dose (mg)": dose, "mg/kg": round(dose / R['peso'], 1), "Ogni (h)": tau,
                              "Picco": round(pk[0], 1), "Picco/MIC": round(pk[0] / mic, 1),
                              "Valle": round(pk[1], 2), "AUC24": round(pk[2]), "OK": "sì" if ok else "no"})
                if ok and scelta is None:
                    scelta = righe[-1]
        st.dataframe(pd.DataFrame(righe), width="stretch", hide_index=True)
        st.caption("Picco e valle in mg/L, AUC24 in mg·h/L, allo stato stazionario del regime indicato.")
        if scelta:
            st.success(f"Proposta: {ui.n(scelta['Dose (mg)'], 0)} mg ({ui.n(scelta['mg/kg'], 1)} mg/kg) ogni "
                       f"{scelta['Ogni (h)']} h, dalla prossima dose quando il livello è sotto "
                       f"{ui.n(valle_max, 1)} mg/L.")
            prop = (f"Proposta: {ui.n(scelta['Dose (mg)'], 0)} mg ({ui.n(scelta['mg/kg'], 1)} mg/kg) ogni "
                    f"{scelta['Ogni (h)']} h (picco previsto {ui.n(scelta['Picco'], 1)} mg/L, "
                    f"picco/MIC {ui.n(scelta['Picco/MIC'], 1)}, valle {ui.n(scelta['Valle'], 2)} mg/L)")
        else:
            st.warning("Nessuna combinazione soddisfa insieme picco/MIC e valle: valutare intervallo più lungo, "
                       "obiettivo diverso o alternativa terapeutica.")
            prop = "Nessun regime della tabella soddisfa insieme i due criteri."
    return (f"Picco dopo l'ultima dose {ui.n(cmax_ult, 1)} mg/L (picco/MIC {ui.n(cmax_ult / mic, 1)}; obiettivo ≥ "
            f"{ui.n(obiettivo, 0)}); sotto {ui.n(valle_max, 1)} mg/L dopo "
            f"{ui.n(t_lib - t_last, 1) + ' h' if t_lib else '> 240 h'}\n{prop}\n")
