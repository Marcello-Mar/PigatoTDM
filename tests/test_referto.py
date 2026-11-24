"""
Referto di consulenza: logica di decisione e percorso completo (AppTest) sui 32 farmaci.
Eseguire con:  python -m pytest tests -q
"""

import json
import sys
from pathlib import Path

import pytest

RADICE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RADICE))

from pigato_bayes import referto_page as rp, rules  # noqa: E402

R = rules.REGOLE


def val(fid, lv, mic=None, kp=None, sede='sistemico', gruppo='—'):
    return rp._valuta(R[fid]['tipo'], R[fid], lv, mic, kp, sede, 1000, None, gruppo, fid)


# ------------------------------------------------------------------ decisioni
def test_css_decide_sul_rapporto_piu_basso():
    # piperacillina Css 40, MIC 8 → plasma ×5 (MANTIENI); polmone Kp 0,4 → sito ×2 → AUMENTA +25
    v = val('piperacillin-tazobactam', {'conc': 40.0}, mic=8, kp=0.4, sede='lung')
    assert v['r_plasma'] == pytest.approx(5) and v['r_sito'] == pytest.approx(2)
    assert (v['azione'], v['pct'], v['base']) == ('AUMENTA', 25, 'sito')


def test_css_tossicita_forza_riduzione():
    v = val('meropenem', {'conc': 70.0}, mic=8)            # ×8,75 → RIDUCI −25 e tox > 64,2
    assert v['azione'] == 'RIDUCI' and v['tox']
    v = val('meropenem', {'conc': 70.0}, mic=16)           # ×4,4 → MANTIENI, ma tox → RIDUCI
    assert (v['azione'], v['pct']) == ('RIDUCI', -25)


def test_range_linezolid_e_teicoplanina():
    assert val('linezolid', {'conc': 9.5})['azione'] == 'RIDUCI'
    assert val('linezolid', {'conc': 1.5})['pct'] == 25
    assert val('teicoplanin', {'conc': 35})['azione'] == 'SOSPENDI'


def test_dalbavancina_sotto_target():
    assert val('dalbavancin', {'conc': 6.0})['azione'] == 'DOSE AGGIUNTIVA'
    assert val('dalbavancin', {'conc': 12.0})['azione'] == 'MANTIENI'


def test_vancomicina_auc():
    v = val('vancomycin', {'auc': 300.0, 'conc': None}, mic=1)
    assert v['azione'] == 'AUMENTA' and v['pct'] == round(500 / 300 * 100 - 100)
    v = val('vancomycin', {'auc': 800.0, 'conc': None}, mic=1)
    assert v['azione'] == 'RIDUCI'


def test_vancomicina_continua_css_non_e_una_valle():
    # Css 22 → AUC 528, MIC 1: in target; la Css non va confrontata con la soglia di valle 20
    lv = {'mod': 'Infusione continua (Css)', 'conc': 22.0, 'auc': 22.0 * 24}
    v = val('vancomycin', lv, mic=1)
    assert v['azione'] == 'MANTIENI' and not v['tox']


def test_vancomicina_intermittente_valle_alta():
    lv = {'mod': 'Due livelli (picco e valle)', 'auc': 550.0, 'cmin_est': 21.0, 'valle': 21.0}
    v = val('vancomycin', lv, mic=1)
    assert v['azione'] == 'RIDUCI' and v['tox']


def test_aminoglicosidi():
    assert val('gentamicin', {'picco': 6.0, 'valle': 0.3}, mic=1)['azione'] == 'AUMENTA'
    assert val('gentamicin', {'picco': 10.0, 'valle': 0.3}, mic=1)['azione'] == 'MANTIENI'
    assert val('gentamicin', {'picco': 14.0, 'valle': 0.3}, mic=1)['azione'] == 'RIDUCI'
    v = val('amikacin', {'picco': 40.0, 'valle': 6.0}, mic=4)
    assert v['azione'] == 'ALLUNGA INTERVALLO' and 'nefrotossicità' in v['tox']


def test_fluorochinoloni():
    assert val('ciprofloxacin', {'picco': 3.0}, mic=0.5)['azione'] == 'AUMENTA'
    assert val('ciprofloxacin', {'picco': 6.0, 'valle': 2.5}, mic=0.5)['azione'] == 'RIDUCI'
    assert val('levofloxacin', {'picco': 12.0, 'valle': 2.4}, mic=1)['azione'] == 'ALLUNGA INTERVALLO'


def test_fosfomicina_e_cotrimossazolo():
    v = val('fosfomycin', {'auc': 50 * 24}, mic=32, gruppo='enterobacterales')   # 37,5 < 83,3
    assert v['azione'] == 'AUMENTA' and v['pct'] == 50
    assert val('sulfamethoxazole-trimethoprim', {'tmp': 20.0, 'smx': 150.0})['azione'] == 'RIDUCI'
    assert val('sulfamethoxazole-trimethoprim', {'tmp': 3.0, 'smx': 80.0})['azione'] == 'AUMENTA'


def test_errori_senza_dati():
    assert 'errore' in val('meropenem', {'conc': None}, mic=2)
    assert 'errore' in val('vancomycin', {'conc': None})
    assert 'errore' in val('gentamicin', {})


def test_regime_default():
    assert rp._regime_default(3000, "Intermittente", 8, "mg", "AUMENTA") == "1000 mg ogni 8 h (3000 mg/die)"
    assert "ogni 36 h" in rp._regime_default(400, "Intermittente", 24, "mg", "ALLUNGA INTERVALLO")
    assert rp._regime_default(3000, "Intermittente", 8, "mg", "SOSPENDI") == ""


# ------------------------------------------------------------------ percorso completo
def _harness(fid, livelli, radice):
    import sys
    sys.path.insert(0, radice)
    import streamlit as st
    from pigato_bayes import referto_page as rp
    rp._tabella_livelli = lambda chiave, etichetta, n_righe=1: [{'data': d, 'conc': c} for d, c in livelli]
    st.session_state.setdefault('r_farmaco', fid)
    rp.render()


AGENTI = [a['id'] for a in json.loads((RADICE / 'pigato_bayes/data/tdm_targets.json').read_text())['agents']]


def _imposta(at, etichetta, valore):
    for w in at.number_input:
        if w.label == etichetta:
            w.set_value(valore)
            at.run()
            return True
    return False


def _input_per_tipo(tipo):
    return {
        'css': {},
        'range': {},
        'info': {},
        'auc': {"Picco (mg/L)": 30.0, "Valle (mg/L)": 10.0},
        'fauc': {"Css (mg/L)": 60.0},
        'amino': {"Picco 30 min dopo fine infusione (mg/L)": 9.0, "Valle prima della dose (mg/L)": 0.4},
        'fq': {"Picco (mg/L)": 4.0, "Valle (mg/L)": 0.6},
        'smx': {"Trimetoprim, picco (mg/L)": 7.0, "Sulfametossazolo, picco (mg/L)": 150.0},
    }[tipo]


@pytest.mark.parametrize("fid", AGENTI)
@pytest.mark.parametrize("paziente", [
    {"Età (anni)": 8.0, "Peso (kg)": 26.0, "Altezza (cm)": 128.0, "Creatinina (mg/dL)": 0.4},
    {"Età (anni)": 62.0, "Peso (kg)": 78.0, "Altezza (cm)": 174.0, "Creatinina (mg/dL)": 2.4},
], ids=["bambino", "adulto_IR"])
def test_referto_tutti_i_farmaci(fid, paziente):
    from streamlit.testing.v1 import AppTest
    tipo = R.get(fid, {'tipo': 'info'})['tipo']
    livelli = ([('20/11/2026', 40.0), ('23/11/2026', 30.0)] if tipo == 'css' else
               [('20/11/2026', 12.0), ('23/11/2026', 10.0)] if tipo in ('range', 'info') else [])
    at = AppTest.from_function(_harness, args=(fid, livelli, str(RADICE)), default_timeout=60)
    at.run()
    at.sidebar.text_input[0].input("TEST01")
    at.run()
    [x for x in at.sidebar.selectbox if x.label == "Sesso biologico"][0].set_value("Maschio")
    at.run()
    for k, v in paziente.items():
        assert _imposta(at, k, v), k
    if tipo not in ('css', 'fauc'):
        dose = "Dose per somministrazione (MUI)" if fid == 'colistin' else "Dose per somministrazione (mg)"
        assert _imposta(at, dose, 3.0 if fid == 'colistin' else 400.0)
    else:
        assert _imposta(at, "Dose giornaliera (mg/die)", 6000.0)
    for k, v in _input_per_tipo(tipo).items():
        assert _imposta(at, k, v), (k, [w.label for w in at.number_input])
    gruppi = [s for s in at.selectbox if s.label == "Gruppo di patogeni"][0]
    if len(gruppi.options) > 1:
        gruppi.set_value(gruppi.options[1])
        at.run()
        [r for r in at.radio if r.label == "Terapia"][0].set_value("Empirica (breakpoint EUCAST)")
        at.run()
    if not [w for w in at.number_input if w.label == "MIC (mg/L)"][0].value:
        _imposta(at, "MIC (mg/L)", 1.0)
    [s for s in at.selectbox if s.label == "Sede"][0].set_value("lung")
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    referti = [c.value for c in at.code]
    assert referti, f"nessun referto per {fid}: {[i.value for i in at.info]}"
    testo = referti[0]
    assert "Commento:" in testo and "Ricontrolla dosaggio" in testo
    if tipo == "css":
        assert "Bibliografia:" in testo and "doi:10.1186/s13054-022-04050-9" in testo
    import re
    assert not re.search(r"\bNone\b|\bnan\b", testo)
    # registrazione nello storico e JSON
    [b for b in at.button if b.label.startswith("Registra")][0].click()
    at.run()
    assert not at.exception
    sto = at.session_state['storico_referto']
    assert len(sto['monitoraggi']) == 1 and sto['monitoraggi'][0]['id_farmaco'] == fid
    json.dumps(sto, default=str)


def test_cambio_azione_aggiorna_la_dose():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_function(_harness, args=('linezolid', [('23/11/2026', 1.5)], str(RADICE)), default_timeout=60)
    at.run()
    at.sidebar.text_input[0].input("TEST02")
    at.run()
    [x for x in at.sidebar.selectbox if x.label == "Sesso biologico"][0].set_value("Maschio")
    at.run()
    for k, v in {"Età (anni)": 40.0, "Peso (kg)": 70.0, "Altezza (cm)": 175.0, "Creatinina (mg/dL)": 0.9,
                 "Dose per somministrazione (mg)": 600.0, "Intervallo (h)": 12.0}.items():
        _imposta(at, k, v)
    nd = [w for w in at.number_input if w.label.startswith("Nuova dose")][0]
    assert nd.value == 1500.0                                  # 1200 mg/die +25%
    az = [s for s in at.selectbox if s.label.startswith("Azione")][0]
    az.set_value("MANTIENI")
    at.run()
    nd = [w for w in at.number_input if w.label.startswith("Nuova dose")][0]
    assert nd.value == 1200.0


# ------------------------------------------------------------------ storico JSON
def _harness_ricarica(radice, payload):
    import io
    import sys
    sys.path.insert(0, radice)
    import streamlit as st
    from pigato_bayes import referto_page as rp, storage
    if 'caricato' not in st.session_state:
        st.session_state['caricato'] = True
        ok, msg, dati = storage.carica('referto', 'Linezolid', io.BytesIO(payload))
        st.session_state['esito'] = (ok, msg)
        if ok:
            rp._widget_da_json(dati)
        st.session_state['r_farmaco'] = 'linezolid'
    rp.render()


def test_json_andata_e_ritorno():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_function(_harness, args=('linezolid', [('23/11/2026', 9.5)], str(RADICE)), default_timeout=60)
    at.run()
    at.sidebar.text_input[0].input("PZ-JSON")
    at.run()
    [x for x in at.sidebar.selectbox if x.label == "Sesso biologico"][0].set_value("Maschio")
    at.run()
    for k, v in {"Età (anni)": 6.0, "Peso (kg)": 21.0, "Altezza (cm)": 116.0, "Creatinina (mg/dL)": 0.35,
                 "Dose per somministrazione (mg)": 210.0}.items():
        _imposta(at, k, v)
    [s for s in at.selectbox if s.label == "Depurazione extrarenale"][0].set_value("CRRT")
    at.run()
    [b for b in at.button if b.label.startswith("Registra")][0].click()
    at.run()
    from pigato_bayes import storage
    payload = storage.json_bytes(at.session_state['storico_referto'])

    at2 = AppTest.from_function(_harness_ricarica, args=(str(RADICE), payload), default_timeout=60)
    at2.run()
    assert at2.session_state['esito'][0], at2.session_state['esito']
    assert not at2.exception
    valori = {w.label: w.value for w in at2.number_input}
    assert valori["Peso (kg)"] == 21.0 and valori["Creatinina (mg/dL)"] == 0.35
    assert valori["Età (anni)"] == pytest.approx(6.0, abs=0.01)
    assert [s for s in at2.selectbox if s.label == "Depurazione extrarenale"][0].value == "CRRT"
    assert len(at2.session_state['storico_referto']['monitoraggi']) == 1


def test_json_farmaco_sbagliato_rifiutato():
    import io
    from pigato_bayes import storage
    import streamlit as st
    dati = storage.storico_vuoto('referto', 'Meropenem')
    ok, msg, _ = storage.carica('referto', 'Linezolid', io.BytesIO(json.dumps(dati).encode()))
    assert not ok and 'Meropenem' in msg
    dati = storage.storico_vuoto('bayes', 'Linezolid')
    ok, msg, _ = storage.carica('referto', 'Linezolid', io.BytesIO(json.dumps(dati).encode()))
    assert not ok
    del st


def _harness_una_volta(radice):
    import sys
    sys.path.insert(0, radice)
    import streamlit as st
    from pigato_bayes import referto_page as rp
    if 'init' not in st.session_state:
        st.session_state['init'] = True
        st.session_state['r_farmaco'] = 'linezolid'
    rp.render()


def test_cambio_farmaco_chiede_conferma():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_function(_harness_una_volta, args=(str(RADICE),), default_timeout=60)
    at.run()
    [b for b in at.sidebar.button if b.label == "Cambia farmaco"][0].click()
    at.run()
    assert not at.exception
    assert at.session_state['r_farmaco'] == 'linezolid'        # niente è cancellato prima della conferma
    [b for b in at.button if b.label == "Annulla"][0].click()
    at.run()
    assert at.session_state['r_farmaco'] == 'linezolid'
    [b for b in at.sidebar.button if b.label == "Cambia farmaco"][0].click()
    at.run()
    [b for b in at.button if b.label == "Conferma"][0].click()
    at.run()
    assert not at.exception
    assert 'r_farmaco' not in at.session_state


def test_bibliografia_per_farmaco():
    agenti = {a['id']: a for a in json.loads((RADICE / 'pigato_bayes/data/tdm_targets.json').read_text())['agents']}
    b = rules.bibliografia('meropenem', 'css', agenti['meropenem'])
    assert len(b) == 2 and 'Crit Care. 2022;26:178' in b[0] and 'Expert Rev Clin Pharmacol' in b[1]
    assert any('Rybak' in x for x in rules.bibliografia('vancomycin', 'auc', agenti['vancomycin']))
    assert any('Bhavnani' in x for x in rules.bibliografia('daptomycin', 'auc', agenti['daptomycin']))
    assert any('ECPA' in agenti['linezolid']['tim'] for _ in [0]) and rules.bibliografia('linezolid', 'range', agenti['linezolid'])
    assert any('EUCAST' in x for x in rules.bibliografia('gentamicin', 'amino', agenti['gentamicin'], empirica=True))
    assert rules.bibliografia('gentamicin', 'amino', agenti['gentamicin']) == []


def test_regole_identiche_a_gatti_2022_tabella1():
    """Gatti M et al. Crit Care 2022;26:178, tabella 1."""
    c = rules.classifica_fasce
    for r, atteso in [(1.9, ('AUMENTA', 50)), (3, ('AUMENTA', 25)), (6, ('MANTIENI', 0)),
                      (9, ('RIDUCI', -25)), (11, ('RIDUCI', -50))]:
        assert c(r, rules.FASCE_CSS) == atteso
    lzd = R['linezolid']['fasce']
    assert [c(x, lzd) for x in (0.5, 1.5, 5, 10, 16)] == [('AUMENTA', 50), ('AUMENTA', 25), ('MANTIENI', 0),
                                                         ('RIDUCI', -25), ('RIDUCI', -50)]
    flu = R['fluconazole']['fasce']
    assert [c(x, flu)[0] for x in (8, 25, 40, 60)] == ['AUMENTA', 'MANTIENI', 'RIDUCI', 'RIDUCI']
    assert c(60, flu)[1] == -50
    assert c(0.5, R['voriconazole']['fasce'])[0] == 'ACCORCIA INTERVALLO' and R['voriconazole']['tau_accorcia'] == 8
    assert c(9, R['voriconazole']['fasce'])[0] == 'SOSPENDI'
    assert c(0.5, R['posaconazole']['fasce'])[0] == 'ACCORCIA INTERVALLO' and R['posaconazole']['tau_accorcia'] == 12
    assert c(0.3, R['aciclovir']['fasce'])[0] == 'ACCORCIA INTERVALLO' and R['aciclovir']['tau_accorcia'] == 6
    assert c(0.3, R['ganciclovir-ganciclovir']['fasce'])[0] == 'ACCORCIA INTERVALLO'
    assert c(5.05, R['isavuconazole']['fasce'])[0] == 'RIDUCI'
    for f, tox in (('piperacillin-tazobactam', 361), ('meropenem', 64.2), ('cefepime', 36)):
        assert R[f]['tox'] == tox


@pytest.mark.parametrize("tau,azione,regime", [(12, "ACCORCIA INTERVALLO", "200 mg ogni 8 h (600 mg/die)"),
                                               (8, "AUMENTA", "250 mg ogni 8 h (750 mg/die)")])
def test_voriconazolo_cmin_bassa(tau, azione, regime):
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_function(_harness, args=('voriconazole', [('23/11/2026', 0.6)], str(RADICE)), default_timeout=60)
    at.run()
    at.sidebar.text_input[0].input("TV")
    [x for x in at.sidebar.selectbox if x.label == "Sesso biologico"][0].set_value("Maschio")
    at.run()
    for k, v in {"Età (anni)": 12.0, "Peso (kg)": 40.0, "Altezza (cm)": 150.0, "Creatinina (mg/dL)": 0.5,
                 "Dose per somministrazione (mg)": 200.0, "Intervallo (h)": float(tau)}.items():
        _imposta(at, k, v)
    assert not at.exception
    az = [s for s in at.selectbox if s.label.startswith("Azione")][0].value
    assert az == azione
    testo = at.code[0].value
    assert regime in testo
