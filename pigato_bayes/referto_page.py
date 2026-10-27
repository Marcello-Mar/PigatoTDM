"""
Modulo referto di consulenza TDM (32 farmaci di tdm_data.json).

Dati biometrici → funzione renale (formula per età) → terapia in atto →
livelli → sede (penetrazione PIGATO) → patogeno e MIC (misurata o breakpoint
EUCAST in empirica) → proposta (riduci / mantieni / aumenta / sospendi) che
l'operatore conferma o modifica → testo del referto → storico JSON.
"""

import hashlib
import json
import math
from datetime import date, datetime
from pathlib import Path

import pandas as pd
import streamlit as st

from . import renal, rules, storage, ui

MODULO = 'referto'
P = 'r'
DATA = Path(__file__).parent / 'data'


@st.cache_data
def _dati():
    tdm = json.loads((DATA / 'tdm_targets.json').read_text(encoding='utf-8'))
    pen = json.loads((DATA / 'penetration.json').read_text(encoding='utf-8'))
    bp = json.loads((DATA / 'breakpoints.json').read_text(encoding='utf-8'))
    return tdm, pen, bp


def _k(x):
    return f"{P}_{x}"


SEDE_PROSA = {'cns': 'SNC/liquor', 'lung': 'polmone', 'bone': 'osso', 'skin': 'cute e tessuti molli',
              'abdomen': 'addome/peritoneo', 'eye': 'occhio', 'urine': 'urine'}

WIDGETS_PAZ =['eta_unita', 'eta_val', 'peso', 'altezza', 'crea', 'eg', 'rrt', 'codice', 'sesso']


def _pulisci(tutto=True):
    for k in list(st.session_state.keys()):
        if k.startswith(f"{P}_") and k not in (f"{P}_farmaco", f"{P}_ultimo_file", f"{P}_msg"):
            if tutto or k not in (_k('codice'), _k('sesso')):
                st.session_state.pop(k, None)


def _widget_da_json(dati):
    _pulisci()
    st.session_state[_k('codice')] = dati.get('codice_paziente', '') or ''
    st.session_state[_k('sesso')] = dati.get('sesso', 'Maschio') or 'Maschio'
    ult = dati['monitoraggi'][-1] if dati.get('monitoraggi') else None
    if not ult:
        return
    if ult.get('eta_giorni') is not None:
        st.session_state[_k('eta_unita')] = 'Anni'
        st.session_state[_k('eta_val')] = round(ult['eta_giorni'] / 365.25, 2)
    for campo, w in (('peso', 'peso'), ('altezza', 'altezza'), ('creatinina', 'crea'), ('eg_sett', 'eg')):
        if ult.get(campo) is not None:
            st.session_state[_k(w)] = float(ult[campo])
    for campo in ('sede', 'gruppo', 'modalita'):
        if ult.get(campo) is not None:
            st.session_state[_k(campo)] = ult[campo]


def _agenti():
    tdm, _, _ = _dati()
    return sorted(tdm['agents'], key=lambda a: a['it'].lower())


def _scelta_farmaco():
    ui.passo("Passo 1")
    st.markdown("### Farmaco della consulenza")
    ag = _agenti()
    gruppi = sorted({a['g_it'] for a in ag})
    g = st.selectbox("Classe", ["Tutte"] + gruppi, key=f"{P}_sel_classe")
    lista = [a for a in ag if g == "Tutte" or a['g_it'] == g]
    nome = st.selectbox("Farmaco", [a['it'] for a in lista], key=f"{P}_sel_nome")
    a = next(x for x in lista if x['it'] == nome)
    st.caption(f"Indicazione: {a.get('ind_it') or a['ind']}")
    if st.button("Inizia la consulenza", type="primary", key=f"{P}_inizia"):
        st.session_state[_k('farmaco')] = a['id']
        storage.reset(MODULO, a['it'])
        _pulisci()
        st.rerun()


def _fmt_ratio(r):
    if r is None:
        return "N/D"
    return ui.n(r, 0) if r >= 10 else (ui.n(r, 1) if r >= 1 else ui.n(r, 2))


def render():
    ui.intestazione("Referto TDM")
    st.caption(ui.AVVISO)
    fid = st.session_state.get(_k('farmaco'))
    if not fid:
        _scelta_farmaco()
        return
    tdm, PEN, BP = _dati()
    ag = next(a for a in tdm['agents'] if a['id'] == fid)
    regola = rules.REGOLE.get(fid, {'tipo': 'info', 'analita': ag['it'].lower()})
    nome = ag['it']

    st.info(f"Consulenza in corso: **{nome}** · {ag['g_it']}")
    if st.sidebar.button("Cambia farmaco", width="stretch", key=f"{P}_cambia"):
        _pulisci()
        st.session_state.pop(_k('farmaco'), None)
        storage.reset(MODULO)
        st.rerun()

    sto, codice, sesso = storage.sidebar_storico(MODULO, nome, P, _widget_da_json)
    if st.sidebar.button("Nuovo paziente (stesso farmaco)", width="stretch", key=f"{P}_nuovo"):
        _pulisci()
        storage.reset(MODULO, nome)
        st.session_state.pop(f'{P}_ultimo_file', None)
        st.session_state.pop(f'{P}_msg', None)
        st.rerun()
    if not codice:
        st.warning("Inserisci un codice paziente nella barra laterale per procedere.")
        return

    with st.expander("Target e regole di aggiustamento (tdm_data.json)"):
        st.markdown("**Target:** " + "; ".join(ag.get('tgt_it') or ag['tgt']))
        st.markdown("**Tossicità:** " + (ag.get('tox_it') or ag['tox']))
        if ag.get('dec_it'):
            st.markdown("**Riduzione:** " + "; ".join(ag['dec_it']))
        if ag.get('inc_it'):
            st.markdown("**Aumento:** " + "; ".join(ag['inc_it']))
        st.markdown("**Tempi:** " + (ag.get('tim_it') or ag['tim']))
        st.markdown("**Determinanti:** " + (ag.get('drv_it') or ag['drv']))
        st.caption("Fonte: " + (ag.get('src') or tdm['src']))

    # ---------------- Paziente ----------------
    ui.passo("Passo 2")
    st.markdown("### Paziente e funzione renale")
    with st.container(border=True):
        a, b, c = st.columns(3, vertical_alignment="bottom")
        unita = a.selectbox("Unità età", ["Anni", "Mesi", "Giorni"], key=_k('eta_unita'))
        eta_val = b.number_input(f"Età ({unita.lower()})", 0.0, {"Anni": 100.0, "Mesi": 1200.0, "Giorni": 36500.0}[unita],
                                 value=None, key=_k('eta_val'))
        peso = c.number_input("Peso (kg)", 0.3, 300.0, value=None, format="%.1f", key=_k('peso'))
        a, b, c = st.columns(3, vertical_alignment="bottom")
        altezza = a.number_input("Altezza (cm)", 25.0, 250.0, value=None, format="%.0f", key=_k('altezza'))
        crea = b.number_input("Creatinina (mg/dL)", 0.05, 25.0, value=None, step=0.01, format="%.2f", key=_k('crea'))
        gg = int(round(eta_val * {"Anni": 365.25, "Mesi": 30.44, "Giorni": 1}[unita])) if eta_val is not None else None
        eta = gg / 365.25 if gg is not None else None
        rrt = c.selectbox("Depurazione extrarenale", ["Nessuna", "CRRT", "Emodialisi intermittente",
                                                      "Dialisi peritoneale"], key=_k('rrt'))
        eg_sett = None
        if eta is not None and eta < 1:
            eg_sett = st.number_input("EG alla nascita (sett.)", 22.0, 42.0, value=None, key=_k('eg'))
    eg = renal.egfr_scelta(eta, sesso, peso, altezza, crea, eg_sett) if (eta is not None and crea) else None
    if eg:
        ui.schede(ui.schede_egfr(eg))
        st.caption(eg['nota'])

    # ---------------- Terapia ----------------
    ui.passo("Passo 3")
    st.markdown("### Terapia in atto")
    tipo = regola['tipo']
    unita_dose = "MUI" if fid == 'colistin' else "mg"
    with st.container(border=True):
        a, b, c, d = st.columns(4, vertical_alignment="bottom")
        modalita = a.selectbox("Somministrazione", ["Infusione continua", "Infusione estesa", "Intermittente"],
                               index=0 if tipo in ('css', 'fauc') else 2, key=_k('modalita'))
        if modalita == "Infusione continua":
            dose_die = b.number_input(f"Dose giornaliera ({unita_dose}/die)", 0.0, 100000.0, value=None,
                                      step=0.5 if unita_dose == "MUI" else 50.0,
                                      key=_k('dose_die'))
            dose_som, tau, tinf = None, 24.0, 24.0
        else:
            dose_som = b.number_input(f"Dose per somministrazione ({unita_dose})", 0.0, 20000.0, value=None,
                                      step=0.5 if unita_dose == "MUI" else 10.0,
                                      key=_k('dose_som'))
            tau = c.number_input("Intervallo (h)", 1.0, 168.0, 8.0 if tipo == 'css' else 24.0, 1.0, key=_k('tau'))
            tinf = d.number_input("Durata infusione (h)", 0.05, 24.0, 3.0 if modalita == "Infusione estesa" else 0.5,
                                  0.05, key=_k('tinf'))
            dose_die = dose_som * 24.0 / tau if dose_som else None
        carico = st.checkbox("Previa dose di carico", value=(modalita == "Infusione continua"), key=_k('carico'))
        if fid == 'colistin':
            st.caption("Colistimetato: inserisci la dose in MUI al posto dei mg (le etichette del referto useranno MUI).")
        if dose_die and peso:
            st.caption(f"Dose giornaliera {ui.n(dose_die, 0)} {unita_dose}/die = {ui.n(dose_die / peso, 1)} {unita_dose}/kg/die")

    # ---------------- Livelli ----------------
    ui.passo("Passo 4")
    st.markdown("### Livelli")
    lv = _input_livelli(tipo, regola, modalita, tau, tinf, nome)

    # ---------------- Sede, patogeno, MIC ----------------
    ui.passo("Passo 5")
    st.markdown("### Sede d'infezione, patogeno, MIC")
    pen_mol = PEN['molecules'].get(ag.get('pen') or '')
    with st.container(border=True):
        a, b = st.columns(2, vertical_alignment="bottom")
        sedi = ["sistemico"] + list(PEN['site_labels'])
        etich = {"sistemico": "Sangue / sistemica", **PEN['site_labels']}
        sede = a.selectbox("Sede", sedi, format_func=lambda s: etich[s], key=_k('sede'))
        kp, kp_info = None, None
        if sede != "sistemico":
            kp, kp_info = _penetrazione(pen_mol, sede, eta, PEN)
            kp = b.number_input("Rapporto sito/plasma (Kp)", 0.0, 200.0, value=kp, step=0.05, key=_k(f'kp_{sede}'),
                                help="Kp armonizzato PIGATO = C(sito, totale)/C(plasma, totale). Modificabile.")
            if kp_info:
                st.caption(kp_info)
        a, b = st.columns(2, vertical_alignment="bottom")
        gruppi_bp = BP['farmaci'].get(fid, {})
        gruppo = a.selectbox("Gruppo di patogeni", ["—"] + list(gruppi_bp),
                             format_func=lambda g: BP['gruppi'].get(g, g), key=_k('gruppo'))
        patogeno = b.text_input("Patogeno (testo del referto)",
                                value=BP['gruppi'].get(gruppo, "") if gruppo != "—" else "", key=_k(f'patogeno_{gruppo}'))
        c, d = st.columns(2, vertical_alignment="center")
        empirica = c.radio("Terapia", ["Mirata (MIC misurata)", "Empirica (breakpoint EUCAST)"], key=_k('empirica'),
                           horizontal=True).startswith("Empirica")
        mic_def = None
        if empirica and gruppo != "—":
            bpg = gruppi_bp.get(gruppo, {})
            r_bp = bpg.get('R_men') if (sede == 'cns' and bpg.get('R_men') is not None) else bpg.get('R')
            mic_def = float(r_bp) if r_bp is not None else None
            if sede == 'cns' and bpg.get('R_men') is not None:
                st.caption(f"Sede SNC: applicato il breakpoint EUCAST per meningite (R > {ui.n(bpg['R_men'], 3)}).")
            st.caption(f"EUCAST v16.1, {BP['gruppi'][gruppo]}: S ≤ {ui.n(bpg.get('S'), 3) if bpg.get('S') else '— (solo I)'}"
                       f", R > {ui.n(bpg.get('R'), 3) if bpg.get('R') else 'N/D'}"
                       + (f". {bpg['nota']}" if bpg.get('nota') else "")
                       + ". In empirica si usa la massima MIC ancora S/I.")
        mic = d.number_input("MIC (mg/L)", 0.0005, 512.0, value=mic_def, step=0.001, format="%.4f",
                             key=_k(f"mic_{gruppo}_{empirica}_{sede == 'cns'}"))
        mic_le = d.checkbox("Riportata come ≤", value=empirica, key=_k(f"micle_{empirica}"))
        if tipo in ('range', 'smx') and regola.get('target'):
            st.caption("Per questo farmaco il target è una concentrazione assoluta: la MIC entra solo nel commento.")

    # ---------------- Valutazione ----------------
    ui.passo("Passo 6")
    st.markdown("### Valutazione e proposta")
    val = _valuta(tipo, regola, lv, mic, kp, sede, dose_die, BP, gruppo, fid)
    if val.get('errore'):
        st.info(val['errore'])
        return
    _mostra_valutazione(val, regola)

    rene_txt = rules.descrivi_rene(eg['valore'] if eg else None)
    nuova_dose = rules.arrotonda_dose(dose_die * (1 + val['pct'] / 100)) if (dose_die and val['pct'] not in (None, -100)) else None
    firma = hashlib.md5(json.dumps([val['azione'], val['pct'], nuova_dose, mic, kp, sede, lv.get('conc'),
                                    dose_die, modalita], default=str).encode()).hexdigest()[:8]
    with st.container(border=True):
        a, b = st.columns(2, vertical_alignment="bottom")
        azioni = ["RIDUCI", "MANTIENI", "AUMENTA", "SOSPENDI", "ALLUNGA INTERVALLO", "DOSE AGGIUNTIVA"]
        azione = a.selectbox("Azione (confermare o modificare)", azioni,
                             index=azioni.index(val['azione']) if val['azione'] in azioni else 1, key=_k(f'azione_{firma}'))
        nd = b.number_input(f"Nuova dose giornaliera ({unita_dose}/die)", 0.0, 100000.0,
                            value=float(nuova_dose) if nuova_dose else (float(dose_die) if dose_die else None),
                            step=50.0, key=_k(f'nuova_{firma}'))
        tim_def = _ricontrollo_default(ag, eta)
        ricontrollo = st.text_input("Ricontrollo", value=tim_def, key=_k('ricontrollo'))
        regime_txt = st.text_input("Nuovo regime (testo libero, facoltativo)",
                                   value=_regime_default(nd, modalita, tau, unita_dose, azione),
                                   key=_k(f'regime_{firma}_{azione}_{nd}'))
        commento = st.text_area("Commento (modificabile)", height=120,
                                value=_commento(val, azione, nome, regime_txt, rene_txt, regola, sede, etich),
                                key=_k(f'commento_{firma}_{azione}_{nd}_{regime_txt}'))

    referto = _componi_referto(nome, etich, sede, patogeno, azione, regola, lv, modalita, dose_die, unita_dose,
                               carico, peso, eg, rrt, kp, val, mic, mic_le, empirica, commento, ricontrollo, pen_mol, eta)
    st.markdown("### Referto")
    st.code(referto, language='text', wrap_lines=True)
    st.download_button("Scarica il referto (.txt)", data=referto.encode('utf-8'),
                       file_name=f"Referto_TDM_{nome.replace(' ', '_').replace('/', '-')}_{codice}_{date.today():%Y%m%d}.txt",
                       mime="text/plain", key=_k('dl_txt'))

    reg_key = _k(f'registrato_{firma}_{azione}')
    if st.session_state.get(reg_key):
        st.success("Consulenza registrata nello storico. Scarica il JSON dalla barra laterale.")
    elif st.button("Registra la consulenza nello storico", type="primary", width="stretch", key=_k('registra')):
        sto['monitoraggi'].append({
            'data_val': datetime.now().strftime('%Y-%m-%d %H:%M'), 'farmaco': nome, 'id_farmaco': fid,
            'eta_giorni': gg, 'peso': peso, 'altezza': altezza, 'creatinina': crea, 'eg_sett': eg_sett,
            'egfr_valore': eg['valore'] if eg else None, 'egfr_formula': eg.get('etichetta') if eg else None,
            'rrt': rrt, 'modalita': modalita, 'dose_die': dose_die, 'dose_somministrazione': dose_som,
            'intervallo_h': tau, 'infusione_h': tinf, 'carico': carico,
            'livelli': lv.get('registro', []), 'sede': sede, 'kp': kp, 'gruppo': gruppo, 'patogeno': patogeno,
            'empirica': empirica, 'mic': mic, 'mic_le': mic_le,
            'rapporto_plasma': val.get('r_plasma'), 'rapporto_sito': val.get('r_sito'),
            'azione_proposta': val['azione'], 'azione': azione, 'nuova_dose_die': nd,
            'nuovo_regime': regime_txt, 'ricontrollo': ricontrollo, 'commento': commento, 'referto': referto,
        })
        st.session_state[reg_key] = True
        st.rerun()

    mon = sto.get('monitoraggi', [])
    with st.expander(f"Consulenze registrate — {nome} ({len(mon)})"):
        for i, r in enumerate(mon):
            st.markdown(f"**#{i} · {r.get('data_val')} · {r.get('azione')}**")
            st.code(r.get('referto', ''), language='text', wrap_lines=True)
        if mon:
            sel = st.selectbox("Rimuovi una consulenza", ["—"] + [f"#{i}" for i in range(len(mon))], key=_k('rm_sel'))
            if sel != "—" and st.button("Rimuovi", key=_k('rm')):
                mon.pop(int(sel[1:]))
                st.rerun()


# ------------------------------------------------------------------

def _tabella_livelli(chiave, etichetta, n_righe=1):
    df0 = pd.DataFrame([{"Data": date.today().strftime('%d/%m/%Y'), "Conc. (mg/L)": float('nan')}] * n_righe)
    df = st.data_editor(df0,
                        num_rows="dynamic", width="stretch", key=_k(chiave),
                        column_config={"Conc. (mg/L)": st.column_config.NumberColumn(format="%.2f", min_value=0.0)})
    out = []
    for _, r in df.iterrows():
        try:
            c = float(r.get("Conc. (mg/L)"))
        except (TypeError, ValueError):
            continue
        if math.isnan(c):
            continue
        out.append({'data': str(r.get("Data") or ""), 'conc': c})
    if etichetta:
        st.caption(etichetta)
    return out


def _input_livelli(tipo, regola, modalita, tau, tinf, nome):
    lv = {'registro': []}
    with st.container(border=True):
        if tipo in ('css', 'info'):
            lab = "Css" if modalita == "Infusione continua" else "Cmin (valle)"
            st.markdown(f"**{regola['analita'].capitalize()} — {lab}**")
            righe = _tabella_livelli('liv_primo', "Più righe = più prelievi: la valutazione usa l'ultimo; tutti "
                                                  "compaiono nel referto.")
            lv.update(tipo_livello=lab, serie=righe, conc=righe[-1]['conc'] if righe else None)
            lv['registro'] += [{'analita': regola['analita'], 'tipo': lab, **x} for x in righe]
            if regola.get('comp'):
                cmp = regola['comp']
                cc = st.number_input(f"{cmp['nome'].capitalize()} (mg/L, facoltativo)", 0.0, 1000.0, value=None,
                                     step=0.1, key=_k('comp'))
                lv['comp'] = cc
                if cc is not None:
                    lv['registro'].append({'analita': cmp['nome'], 'tipo': lab, 'conc': cc})
        elif tipo == 'range':
            lab = regola.get('livello', 'Cmin')
            st.markdown(f"**{regola['analita'].capitalize()} — {lab}**")
            righe = _tabella_livelli('liv_primo', "La valutazione usa l'ultimo prelievo.")
            lv.update(tipo_livello=lab, serie=righe, conc=righe[-1]['conc'] if righe else None)
            lv['registro'] += [{'analita': regola['analita'], 'tipo': lab, **x} for x in righe]
        elif tipo in ('auc', 'fq', 'fauc'):
            opz = ["Infusione continua (Css)", "Due livelli (picco e valle)"] if tipo != 'fq' else \
                ["Due livelli (picco e valle)", "Solo picco", "Solo valle"]
            if tipo == 'auc':
                opz.append("Solo valle")
            mod = st.radio("Livelli disponibili", opz, horizontal=True,
                           index=0 if (modalita == "Infusione continua" or tipo == 'fq') else 1, key=_k('mod_liv'))
            lv['mod'] = mod
            if mod.startswith("Infusione continua"):
                css = st.number_input("Css (mg/L)", 0.0, 1000.0, value=None, step=0.1, key=_k('css'))
                lv.update(conc=css, tipo_livello='Css', auc=css * 24 if css is not None else None)
                if css is not None:
                    lv['registro'].append({'analita': regola['analita'], 'tipo': 'Css', 'conc': css})
            else:
                a, b, c, d = st.columns(4, vertical_alignment="bottom")
                picco = t1 = valle = t2 = None
                if mod != "Solo valle":
                    picco = a.number_input("Picco (mg/L)", 0.0, 1000.0, value=None, step=0.1, key=_k('picco'))
                    t1 = b.number_input("Ora del picco (h dall'inizio dose)", 0.0, 168.0, value=float(tinf + 0.5), step=0.25,
                                        key=_k('t1'))
                if mod != "Solo picco":
                    valle = c.number_input("Valle (mg/L)", 0.0, 1000.0, value=None, step=0.1, key=_k('valle'))
                    t2 = d.number_input("Ora della valle (h dall'inizio dose)", 0.0, 168.0, value=float(max(tau - 0.5, 0.5)),
                                        step=0.25, key=_k('t2'))
                lv.update(picco=picco, t1=t1, valle=valle, t2=t2, tau=tau, tinf=tinf)
                if picco and valle and t2 > t1:
                    r = rules.auc_da_due_livelli(picco, t1, valle, t2, tau, tinf)
                    if r:
                        lv.update(auc=r[0], cmax_est=r[1], cmin_est=r[2], ke=r[3])
                        st.caption(f"AUC24 {ui.n(r[0], 0)} mg·h/L; Cmax a fine infusione {ui.n(r[1], 1)}, "
                                   f"Cmin a fine intervallo {ui.n(r[2], 2)} mg/L; t½ {ui.n(math.log(2) / r[3], 1)} h "
                                   "(Sawchuk-Zaske, 1 compartimento).")
                for x, tl in ((picco, 'picco'), (valle, 'valle')):
                    if x is not None:
                        lv['registro'].append({'analita': regola['analita'], 'tipo': tl, 'conc': x})
        elif tipo == 'amino':
            st.caption(f"Con un solo livello random (es. a 10 h) usa il modulo TDM bayesiana per {nome.lower()}: "
                       "stima la clearance e propone dose e intervallo.")
            a, b = st.columns(2, vertical_alignment="bottom")
            picco = a.number_input("Picco 30 min dopo fine infusione (mg/L)", 0.0, 300.0, value=None, step=0.1,
                                   key=_k('picco'))
            valle = b.number_input("Valle prima della dose (mg/L)", 0.0, 100.0, value=None, step=0.1, key=_k('valle'))
            lv.update(picco=picco, valle=valle)
            for x, tl in ((picco, 'picco'), (valle, 'valle')):
                if x is not None:
                    lv['registro'].append({'analita': regola['analita'], 'tipo': tl, 'conc': x})
        elif tipo == 'smx':
            a, b = st.columns(2, vertical_alignment="bottom")
            tmp = a.number_input("Trimetoprim, picco (mg/L)", 0.0, 100.0, value=None, step=0.1, key=_k('tmp'))
            smx = b.number_input("Sulfametossazolo, picco (mg/L)", 0.0, 1000.0, value=None, step=1.0, key=_k('smx'))
            lv.update(tmp=tmp, smx=smx)
            for x, tl in ((tmp, 'trimetoprim'), (smx, 'sulfametossazolo')):
                if x is not None:
                    lv['registro'].append({'analita': tl, 'tipo': 'picco', 'conc': x})
    return lv


def _penetrazione(pen_mol, sede, eta, PEN):
    if not pen_mol:
        return None, "Nessun dato di penetrazione nel dataset PIGATO per questa molecola."
    s = pen_mol['sites'].get(sede)
    if not s or (not s.get('kp') and not s.get('ped_kp')):
        return None, "Nessun dato di penetrazione per questa sede (No data)."
    ped = eta is not None and eta < 18 and s.get('ped_kp')
    kpd = s['ped_kp'] if ped else s['kp']
    q = PEN['quality_it'].get(s.get('q'), s.get('q'))
    fonte = "pediatrica" if ped else "adulto"
    info = (f"Kp tipico {ui.n(kpd['typ'], 2)} (intervallo {ui.n(kpd['min'], 2)}–{ui.n(kpd['max'], 2)}; base "
            f"{kpd.get('basis')}; dato {fonte}; evidenza: {q}). "
            + ((s.get('ped_text_it') if ped else s.get('text_it')) or ''))
    return float(kpd['typ']), info


def _valuta(tipo, regola, lv, mic, kp, sede, dose_die, BP, gruppo, fid):
    v = {'azione': 'MANTIENI', 'pct': 0, 'righe': [], 'tox': None}
    urinario = sede == 'urine'
    if tipo in ('css', 'info'):
        c = lv.get('conc')
        if c is None or not mic:
            return {'errore': "Inserisci il livello e la MIC per la valutazione."}
        v['r_plasma'] = c / mic
        if kp is not None:
            v['c_sito'] = c * kp
            v['r_sito'] = v['c_sito'] / mic
        rif = min(v['r_plasma'], v['r_sito']) if v.get('r_sito') is not None else v['r_plasma']
        v['r_decisione'] = rif
        v['base'] = 'sito' if (v.get('r_sito') is not None and v['r_sito'] < v['r_plasma']) else 'plasma'
        if tipo == 'css':
            v['azione'], v['pct'] = rules.classifica_fasce(rif, rules.FASCE_CSS)
        else:
            v['azione'], v['pct'] = 'MANTIENI', 0
            v['righe'].append(regola.get('nota', ''))
        if regola.get('tox') and c > regola['tox']:
            v['tox'] = f"{regola['analita']} {ui.n(c, 1)} mg/L oltre la soglia di neurotossicità di {ui.n(regola['tox'], 1)} mg/L"
            if v['azione'] not in ('RIDUCI',):
                v['azione'], v['pct'] = 'RIDUCI', -25
        if regola.get('comp') and lv.get('comp') is not None:
            ct = regola['comp']['CT']
            v['comp'] = (regola['comp']['nome'], lv['comp'], ct, lv['comp'] / ct)
        if urinario and v.get('r_sito') is not None:
            v['righe'].append("Il valore urinario è un indice direzionale, non un rapporto tissutale: "
                              "la decisione usa il rapporto più basso, qui quello plasmatico.")
    elif tipo == 'range':
        c = lv.get('conc')
        if c is None:
            return {'errore': "Inserisci il livello."}
        v['azione'], v['pct'] = rules.classifica_fasce(c, regola['fasce'])
        v['c'] = c
        if mic:
            v['r_plasma'] = c / mic
            if kp is not None:
                v['c_sito'] = c * kp
                v['r_sito'] = v['c_sito'] / mic
        if fid == 'dalbavancin' and c < regola['target'][0]:
            v['azione'], v['pct'] = 'DOSE AGGIUNTIVA', 0
        if regola.get('tox') and c > regola['tox']:
            v['tox'] = f"{ui.n(c, 2)} mg/L oltre la soglia di {ui.n(regola['tox'], 2)} mg/L ({regola.get('tox_txt', '')})"
    elif tipo == 'auc':
        auc = lv.get('auc')
        if auc is None:
            valle = lv.get('valle') or lv.get('conc')
            if valle is None:
                return {'errore': "Inserisci Css oppure picco e valle."}
            v['c'] = valle
            if regola.get('tox_valle') and valle >= regola['tox_valle']:
                v['tox'] = f"valle {ui.n(valle, 1)} mg/L ≥ {ui.n(regola['tox_valle'], 1)} mg/L"
                v['azione'], v['pct'] = 'RIDUCI', -25
            v['righe'].append("Senza AUC la valutazione di efficacia non è possibile: aggiungere un picco o usare "
                              "il modulo bayesiano.")
        else:
            if not mic:
                return {'errore': "Inserisci la MIC."}
            v['auc'] = auc
            v['r_plasma'] = auc / mic
            lo, hi = regola['auc_mic']
            r = v['r_plasma']
            if r < lo:
                v['azione'], v['pct'] = 'AUMENTA', round((lo + hi) / 2 / r * 100 - 100)
            elif r > hi:
                v['azione'], v['pct'] = 'RIDUCI', round((lo + hi) / 2 / r * 100 - 100)
            if kp is not None and lv.get('conc') is not None:
                v['c_sito'] = lv['conc'] * kp
                v['r_sito'] = v['c_sito'] / mic
            cmin = lv.get('cmin_est') if lv.get('cmin_est') is not None else (lv.get('valle') or lv.get('conc'))
            if regola.get('tox_valle') and cmin is not None and cmin >= regola['tox_valle']:
                v['tox'] = f"Cmin {ui.n(cmin, 1)} mg/L ≥ {ui.n(regola['tox_valle'], 1)} mg/L"
                if v['azione'] != 'RIDUCI':
                    v['azione'], v['pct'] = 'RIDUCI', -25
            if regola.get('tox_auc') and auc > regola['tox_auc'] and fid == 'vancomycin':
                v['righe'].append(f"AUC24 {ui.n(auc, 0)} mg·h/L > 600: rischio di danno renale acuto.")
    elif tipo == 'amino':
        pc, vl = lv.get('picco'), lv.get('valle')
        if pc is None and vl is None:
            return {'errore': "Inserisci picco e/o valle, oppure usa il modulo bayesiano."}
        if pc is not None and mic:
            v['r_plasma'] = pc / mic
            if v['r_plasma'] < regola['picco_mic']:
                v['azione'], v['pct'] = 'AUMENTA', 25
            elif v['r_plasma'] > 12:
                v['azione'], v['pct'] = 'RIDUCI', -25
        if vl is not None and vl > regola['valle_max']:
            v['azione'], v['pct'] = 'ALLUNGA INTERVALLO', 0
            v['tox'] = (f"valle {ui.n(vl, 2)} mg/L > {ui.n(regola['valle_max'], 1)} mg/L"
                        + (" (soglia di nefrotossicità)" if vl > regola['valle_tox'] else ""))
    elif tipo == 'fq':
        pc = lv.get('picco')
        vl = lv.get('valle')
        if pc is None and vl is None:
            return {'errore': "Inserisci almeno un livello."}
        if pc is not None and mic:
            v['r_plasma'] = pc / mic
            if v['r_plasma'] < regola['cmax_mic']:
                v['azione'], v['pct'] = 'AUMENTA', 25
            if kp is not None:
                v['c_sito'] = pc * kp
                v['r_sito'] = v['c_sito'] / mic
        if lv.get('auc') and mic:
            v['auc'] = lv['auc']
            v['righe'].append(f"AUC24/MIC {ui.n(lv['auc'] / mic, 0)} (target > {regola['auc_mic']}).")
        if vl is not None and vl > regola['valle_max']:
            v['tox'] = f"valle {ui.n(vl, 2)} mg/L > {ui.n(regola['valle_max'], 0)} mg/L"
            v['azione'], v['pct'] = ('ALLUNGA INTERVALLO', 0) if fid == 'levofloxacin' else ('RIDUCI', -25)
    elif tipo == 'fauc':
        auc = lv.get('auc')
        if auc is None or not mic:
            return {'errore': "Inserisci Css (o picco e valle) e la MIC."}
        fauc_mic = auc * regola['fu'] / mic
        v['r_plasma'] = fauc_mic
        tgt = regola['target_per_gruppo'].get(gruppo)
        if tgt:
            if fauc_mic < tgt[0]:
                v['azione'], v['pct'] = 'AUMENTA', min(50, round(tgt[0] / fauc_mic * 100 - 100))
            v['righe'].append(f"fAUC24/MIC {ui.n(fauc_mic, 0)} (target {ui.n(tgt[0], 1)}"
                              + (f"–{ui.n(tgt[1], 0)}" if tgt[1] else "") + ").")
        else:
            v['righe'].append("Target fAUC/MIC definito solo per S. aureus, Enterobacterales e P. aeruginosa.")
    elif tipo == 'smx':
        tmp, smx = lv.get('tmp'), lv.get('smx')
        if tmp is None and smx is None:
            return {'errore': "Inserisci i picchi di trimetoprim e/o sulfametossazolo."}
        tox = []
        low = []
        if tmp is not None:
            if tmp > regola['tmp'][1]:
                tox.append(f"trimetoprim {ui.n(tmp, 1)} > {regola['tmp'][1]} mg/L")
            elif tmp < regola['tmp'][0]:
                low.append(f"trimetoprim {ui.n(tmp, 1)} < {regola['tmp'][0]} mg/L")
        if smx is not None:
            if smx > regola['smx'][1]:
                tox.append(f"sulfametossazolo {ui.n(smx, 0)} > {regola['smx'][1]} mg/L")
            elif smx < regola['smx'][0]:
                low.append(f"sulfametossazolo {ui.n(smx, 0)} < {regola['smx'][0]} mg/L")
        if tox:
            v['tox'] = "; ".join(tox)
            v['azione'], v['pct'] = 'RIDUCI', -25
        elif low:
            v['azione'], v['pct'] = 'AUMENTA', 25
            v['righe'].append("Target PCP: " + "; ".join(low) + ".")
    if regola.get('nota'):
        v['righe'].append(regola['nota'])
    return v


def _mostra_valutazione(v, regola):
    voci = []
    if v.get('r_plasma') is not None:
        lab = {"auc": "AUC24/MIC", "fauc": "fAUC24/MIC", "amino": "Picco/MIC", "fq": "Cmax/MIC"}.get(regola['tipo'], "C/MIC plasma")
        voci.append({'k': lab, 'v': _fmt_ratio(v['r_plasma'])})
    if v.get('auc') is not None:
        voci.append({'k': "AUC24", 'v': ui.n(v['auc'], 0), 'u': "mg·h/L"})
    if v.get('c_sito') is not None:
        voci.append({'k': "Concentrazione al sito", 'v': ui.n(v['c_sito'], 1), 'u': "mg/L"})
        voci.append({'k': "C sito/MIC", 'v': _fmt_ratio(v['r_sito'])})
    pct = f" {'+' if v['pct'] > 0 else ''}{v['pct']}%" if v['pct'] not in (0, -100) else ""
    voci.append({'k': "Proposta", 'v': v['azione'] + pct, 'sel': True,
                 'stato': 'ok' if v['azione'] == 'MANTIENI' else None})
    ui.schede(voci)
    if v.get('r_decisione') is not None and v.get('r_sito') is not None:
        st.caption(f"Decisione sul rapporto più basso ({v['base']}): ×{_fmt_ratio(v['r_decisione'])} la MIC.")
    if v.get('tox'):
        st.error("Tossicità: " + v['tox'])
    if v.get('comp'):
        nome, c, ct, r = v['comp']
        (st.success if r >= 1 else st.warning)(f"{nome.capitalize()} {ui.n(c, 1)} mg/L: Css/CT {ui.n(r, 1)} (CT {ui.n(ct, 0)} mg/L).")
    for r in v['righe']:
        if r:
            st.caption(r)


def _ricontrollo_default(ag, eta):
    tim = ag.get('tim', '')
    if 'β-lactam' in ag.get('g', '') or 'ECPA' in tim:
        return "fra 48 ore con eventuale nuova consulenza TDM" if (eta is not None and eta < 18) else \
            "fra 48–72 ore con eventuale nuova consulenza TDM"
    if 'twice weekly' in tim:
        return "fra 3–4 giorni (due volte a settimana)"
    return "alla prossima occasione utile con nuova consulenza TDM"


def _regime_default(nd, modalita, tau, unita, azione):
    if azione in ('SOSPENDI',) or not nd:
        return ""
    if modalita == "Infusione continua":
        return f"{ui.n_semplice(nd, 0)} {unita}/die in infusione continua"
    if azione == 'ALLUNGA INTERVALLO' and tau:
        per = nd * tau / 24.0
        return f"{ui.n_semplice(rules.arrotonda_dose(per) if per >= 50 else per, 0)} {unita} ogni {ui.n(tau + 12, 0)} h"
    n_som = 24.0 / tau if tau else 1
    per = nd / n_som if n_som else nd
    return f"{ui.n_semplice(rules.arrotonda_dose(per) if per >= 50 else per, 0 if per >= 10 else 1)} {unita} ogni {ui.n(tau, 0)} h ({ui.n_semplice(nd, 0)} {unita}/die)"


def _commento(v, azione, nome, regime, rene, regola, sede, etich):
    r = v.get('r_decisione') or v.get('r_plasma')
    dove = "al sito di interesse" if v.get('base') == 'sito' else "plasmatica"
    rene_frase = f", considerata anche la {rene}," if rene else ""
    farm = nome.lower()
    if regola['tipo'] == 'auc' and v.get('auc') is not None:
        lo, hi = regola['auc_mic']
        base = f"AUC24/MIC pari a {_fmt_ratio(v['r_plasma'])} (target {lo}–{hi})"
    elif regola['tipo'] == 'fauc' and v.get('r_plasma') is not None:
        base = f"fAUC24/MIC pari a {_fmt_ratio(v['r_plasma'])}"
    elif regola['tipo'] == 'amino' and v.get('r_plasma') is not None:
        base = f"Picco di {regola['analita']} pari a {_fmt_ratio(v['r_plasma'])} volte la MIC"
    elif regola['tipo'] in ('range',):
        base = f"Livelli di {regola['analita']} pari a {ui.n(v.get('c'), 2)} mg/L (target {ui.n(regola['target'][0], 2)}" + \
               (f"–{ui.n(regola['target'][1], 2)} mg/L)" if regola['target'][1] else " mg/L)")
    elif r is not None:
        base = f"Livelli nel range {'elevato' if azione == 'RIDUCI' else 'basso' if azione == 'AUMENTA' else 'target'} di concentrazione {dove} (×{_fmt_ratio(r)} la MIC)"
    else:
        base = f"Livelli di {regola['analita']} valutati"
    if azione == 'RIDUCI':
        t = f"{base} per cui{rene_frase} è possibile ridurre {farm}" + (f" a {regime}." if regime else ".")
    elif azione == 'AUMENTA':
        t = f"{base} per cui{rene_frase} si consiglia di aumentare {farm}" + (f" a {regime}." if regime else ".")
    elif azione == 'SOSPENDI':
        t = f"{base}: si consiglia di sospendere {farm} e ricontrollare il livello prima di riprendere a dose ridotta."
    elif azione == 'ALLUNGA INTERVALLO':
        t = f"{base}: si consiglia di allungare l'intervallo tra le dosi" + (f" ({regime})." if regime else ".")
    elif azione == 'DOSE AGGIUNTIVA':
        t = f"{base}: si consiglia una dose aggiuntiva di {farm}" + (f" ({regime})." if regime else ".")
    else:
        t = f"{base}: si consiglia di mantenere la posologia in atto."
    if v.get('tox'):
        t += f" Attenzione: {v['tox']}."
    return t


def _componi_referto(nome, etich, sede, patogeno, azione, regola, lv, modalita, dose_die, unita, carico, peso,
                     eg, rrt, kp, v, mic, mic_le, empirica, commento, ricontrollo, pen_mol, eta):
    sede_txt = etich[sede] if sede != 'sistemico' else "Sistemica"
    titolo = " — ".join(x for x in (nome.upper(), sede_txt, patogeno or None, azione) if x)
    # livelli
    reg = lv.get('registro', [])
    if reg:
        per_analita = {}
        for x in reg:
            per_analita.setdefault(x['analita'], []).append(
                f"{ui.n(x['conc'], 1 if x['conc'] >= 1 else 2)} mg/L" + (f" {x['tipo']}" if x.get('tipo') and x['tipo'] not in ('Css',) else "")
                + (f" del {x['data']}" if x.get('data') and sum(1 for y in reg if y['analita'] == x['analita']) > 1 else ""))
        visti = " e ".join(f"{a} ({', '.join(vv)})" for a, vv in per_analita.items())
        r1 = f"Visionati livelli di {visti}."
    else:
        r1 = "Livelli non inseriti."
    mod_txt = {"Infusione continua": "in infusione continua", "Infusione estesa": "in infusione estesa",
               "Intermittente": "in somministrazione intermittente"}[modalita]
    terapia = (f"- Terapia in atto: {nome.lower()} {ui.n_semplice(dose_die, 0)} {unita}/die"
               + (f" (pari a circa {ui.n(dose_die / peso, 0)} {unita}/kg/die)" if dose_die and peso else "")
               + f" {mod_txt}" + (" previa dose carico" if carico else "")) if dose_die else "- Terapia in atto: N/D"
    righe = [titolo, "", r1, terapia]
    if eg and eg.get('valore'):
        righe.append(f"- eGFR stimata: {ui.n(eg['valore'], 0)} mL/min/1,73 m² ({eg['etichetta']})"
                     + (f"; {rrt}" if rrt != "Nessuna" else ""))
    elif rrt != "Nessuna":
        righe.append(f"- Depurazione extrarenale: {rrt}")
    if v.get('c_sito') is not None and sede != 'sistemico':
        if sede == 'urine':
            righe.append(f"- Concentrazione nel sito di interesse (urine, indice urinario ×{ui.n(kp, 0)}): "
                         f"oltre {ui.n_semplice(v['c_sito'], 0)} mg/L")
        else:
            ped = " (dato pediatrico)" if (eta is not None and eta < 18 and pen_mol and
                                           pen_mol['sites'].get(sede, {}).get('ped_kp')) else ""
            righe.append(f"- Concentrazione nel sito di interesse ({SEDE_PROSA.get(sede, sede_txt.lower())}, penetrazione mediana "
                         f"{ui.n(kp * 100, 0)}%{ped}): stimata intorno a {ui.n(v['c_sito'], 1 if v['c_sito'] < 10 else 0)} mg/L")
    if v.get('auc') is not None:
        righe.append(f"- AUC24 stimata: {ui.n(v['auc'], 0)} mg·h/L")
    if v.get('comp'):
        nm, c, ct, r = v['comp']
        righe.append(f"- Livelli di {nm}: {'al di sopra' if r >= 1 else 'al di sotto'} del target di {ui.n(ct, 0)} mg/L")
    if mic:
        righe.append(f"- MIC del patogeno di interesse: {'≤' if mic_le else ''}{ui.n(mic, 3)} mg/L"
                     + (f" ({patogeno})" if patogeno else "")
                     + (" — terapia empirica, breakpoint EUCAST v16.1" if empirica else ""))
    righe += ["", "Commento:", "", commento.strip(), ""]
    if ricontrollo:
        righe.append(f"Ricontrolla dosaggio {ricontrollo}.")
    return "\n".join(righe)
