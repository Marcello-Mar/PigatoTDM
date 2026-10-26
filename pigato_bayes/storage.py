"""
Storico paziente: solo st.session_state + file JSON scaricabile/ricaricabile.
Un file = un paziente × una molecola × un modulo ('bayes' o 'referto').
Nessun dato viene scritto sul server.
"""

import json
from datetime import datetime, timedelta

import pandas as pd
import streamlit as st

VERSIONE_FORMATO = 2


def storico_vuoto(modulo, farmaco=None):
    return {'formato': VERSIONE_FORMATO, 'modulo': modulo, 'codice_paziente': '',
            'sesso': 'Maschio', 'farmaco': farmaco, 'monitoraggi': []}


def chiave(modulo):
    return f'storico_{modulo}'


def get(modulo):
    k = chiave(modulo)
    if k not in st.session_state:
        st.session_state[k] = storico_vuoto(modulo)
    return st.session_state[k]


def reset(modulo, farmaco=None):
    st.session_state[chiave(modulo)] = storico_vuoto(modulo, farmaco)


def carica(modulo, farmaco, file_up):
    """
    Valida e carica un JSON. Restituisce (ok, messaggio, dati).
    Accetta anche i file v6.0 di TDM Buddy (senza 'modulo': trattati come 'bayes').
    """
    try:
        dati = json.load(file_up)
    except Exception:
        return False, "File non valido o corrotto.", None
    mod_file = dati.get('modulo', 'bayes')
    if mod_file != modulo:
        return False, (f"Questo file appartiene al modulo '{mod_file}', non a '{modulo}'."), None
    if dati.get('farmaco') != farmaco:
        return False, (f"Questo file è per {dati.get('farmaco') or 'un farmaco non specificato'}, "
                       f"non per {farmaco}. Cambia molecola o carica il file corretto."), None
    dati.setdefault('monitoraggi', [])
    dati['modulo'] = modulo
    dati['formato'] = VERSIONE_FORMATO
    st.session_state[chiave(modulo)] = dati
    return True, (f"Storico caricato: {len(dati['monitoraggi'])} valutazioni per il codice "
                  f"'{dati.get('codice_paziente', '?')}'."), dati


def json_bytes(storico):
    return json.dumps(storico, indent=2, ensure_ascii=False, default=str).encode('utf-8')


def nome_file(storico, estensione='json'):
    f = (storico.get('farmaco') or 'farmaco').replace(' ', '_').replace('/', '-')
    c = storico.get('codice_paziente') or 'paziente'
    return f"PIGATO_{storico.get('modulo')}_{f}_{c}.{estensione}"


def csv_long(storico):
    """Una riga per evento (biometria, dose, livello, referto)."""
    righe = []
    for i, rec in enumerate(storico.get('monitoraggi', [])):
        eta_g = rec.get('eta_giorni')
        ctx = {'codice_paziente': storico.get('codice_paziente'), 'farmaco': storico.get('farmaco'),
               'modulo': storico.get('modulo'), 'sesso': storico.get('sesso'),
               'id_valutazione': i, 'data_valutazione': rec.get('data_val'),
               'eta_giorni': eta_g, 'eta_anni': round(eta_g / 365.25, 2) if eta_g is not None else None,
               'peso_kg': rec.get('peso'), 'altezza_cm': rec.get('altezza'),
               'creatinina_mgdl': rec.get('creatinina'), 'egfr': rec.get('egfr_valore'),
               'egfr_formula': rec.get('egfr_formula')}
        righe.append({**ctx, 'tipo_evento': 'biometria', 'ore_da_t0': None, 'data_ora_evento': rec.get('data_val'),
                      'valore': None, 'unita': '', 'nota': ''})
        t0 = None
        if rec.get('data_prima_dose'):
            try:
                t0 = datetime.strptime(rec['data_prima_dose'], '%Y-%m-%d %H:%M')
            except ValueError:
                t0 = None
        for d in rec.get('storia_dosi', []) or []:
            righe.append({**ctx, 'tipo_evento': 'dose', 'ore_da_t0': d.get('ore'),
                          'data_ora_evento': (t0 + timedelta(hours=d['ore'])).strftime('%Y-%m-%d %H:%M') if t0 else '',
                          'valore': d.get('dose'), 'unita': 'mg', 'nota': f"tinf {d.get('tinf', '')} h"})
        for l in rec.get('livelli', []) or []:
            righe.append({**ctx, 'tipo_evento': 'livello', 'ore_da_t0': l.get('ore'),
                          'data_ora_evento': (t0 + timedelta(hours=l['ore'])).strftime('%Y-%m-%d %H:%M') if t0 and l.get('ore') is not None else rec.get('data_val'),
                          'valore': l.get('conc'), 'unita': 'mg/L',
                          'nota': ' '.join(str(x) for x in (l.get('analita', ''), l.get('tipo', ''), l.get('nota', '')) if x)})
        if rec.get('azione'):
            righe.append({**ctx, 'tipo_evento': 'referto', 'ore_da_t0': None, 'data_ora_evento': rec.get('data_val'),
                          'valore': rec.get('nuova_dose_die'), 'unita': 'mg/die',
                          'nota': f"{rec.get('azione')} | sede {rec.get('sede', '')} | MIC {rec.get('mic', '')}"})
    cols = ['codice_paziente', 'farmaco', 'modulo', 'id_valutazione', 'data_valutazione', 'tipo_evento',
            'data_ora_evento', 'ore_da_t0', 'valore', 'unita', 'nota', 'sesso', 'eta_giorni', 'eta_anni',
            'peso_kg', 'altezza_cm', 'creatinina_mgdl', 'egfr', 'egfr_formula']
    return pd.DataFrame(righe, columns=cols)


def sidebar_storico(modulo, farmaco, key_prefix, widget_reset_cb=None):
    """Blocco comune nella sidebar: caricamento, codice, sesso, download."""
    sto = get(modulo)
    st.sidebar.markdown("### Storico paziente (locale)")
    st.sidebar.caption(f"File compatibili: {modulo} · {farmaco}")
    up = st.sidebar.file_uploader("Carica storico (JSON)", type=["json"], key=f"{key_prefix}_upload")
    if up is not None:
        fid = f"{modulo}_{up.name}_{up.size}"
        if st.session_state.get(f'{key_prefix}_ultimo_file') != fid:
            st.session_state[f'{key_prefix}_ultimo_file'] = fid
            ok, msg, dati = carica(modulo, farmaco, up)
            if ok:
                if widget_reset_cb:
                    widget_reset_cb(dati)
                st.session_state[f'{key_prefix}_msg'] = ('ok', msg)
                st.rerun()
            else:
                st.session_state[f'{key_prefix}_msg'] = ('err', msg)
    m = st.session_state.get(f'{key_prefix}_msg')
    if m:
        (st.sidebar.success if m[0] == 'ok' else st.sidebar.error)(m[1])
    sto = get(modulo)
    codice = st.sidebar.text_input("Codice paziente (non identificativo)", key=f"{key_prefix}_codice")
    sesso = st.sidebar.selectbox("Sesso biologico", ["Maschio", "Femmina"], key=f"{key_prefix}_sesso")
    sto['codice_paziente'] = codice
    sto['sesso'] = sesso
    sto['farmaco'] = farmaco
    sto['modulo'] = modulo
    n_rec = len(sto.get('monitoraggi', []))
    if n_rec:
        st.sidebar.info(f"{n_rec} valutazioni registrate in questa sessione.")
        st.sidebar.download_button("Scarica storico aggiornato (JSON)", data=json_bytes(sto),
                                   file_name=nome_file(sto), mime="application/json", type="primary",
                                   width="stretch", key=f"{key_prefix}_dl_json")
        st.sidebar.download_button("Scarica dataset long (CSV)", data=csv_long(sto).to_csv(index=False).encode('utf-8'),
                                   file_name=nome_file(sto, 'csv'), mime="text/csv",
                                   width="stretch", key=f"{key_prefix}_dl_csv")
    return sto, codice, sesso
