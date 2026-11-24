"""
Modulo TDM bayesiana: percorso completo per i 6 farmaci e casi limite.
Eseguire con:  python -m pytest tests -q
"""

import sys
from pathlib import Path

import pytest

RADICE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RADICE))

from pigato_bayes import models, pk_engine as pk  # noqa: E402


def _harness(radice, farmaco, modello, gg, peso, altezza, crea, mic, dosi, livelli, dialisi, crrt, eg_sett):
    import sys
    sys.path.insert(0, radice)
    from datetime import datetime, timedelta
    import streamlit as st
    from pigato_bayes import bayes_page as bp, renal, storage
    st.session_state.setdefault('b_farmaco', farmaco)
    sto = storage.get('bayes')
    sto['codice_paziente'] = 'TEST'
    if 'b_risultati' not in st.session_state:
        t0 = datetime(2026, 11, 20, 8, 0)
        storia = [{'dt': t0 + timedelta(hours=h), 'dose': d, 'tinf': ti} for h, d, ti in dosi]
        liv = [{'dt': t0 + timedelta(hours=h), 'conc': c, 'nota': ''} for h, c in livelli]
        eta = gg / 365.25
        eg = renal.egfr_scelta(eta, 'Maschio', peso, altezza, crea, eg_sett)
        bp._calcola(farmaco, 'TEST', 'Maschio', gg, eta, peso, altezza, crea, eg_sett, eg, mic, crrt, 0.85, 2.0,
                    dialisi, modello, storia, liv, -1)
    R = st.session_state.get('b_risultati')
    if R:
        bp._mostra_risultati(R, sto)


CASI = {
    'genta_10h_bambino': ('Gentamicina', models.AMINO_MODELLI['Gentamicina'][0], 5 * 365, 19.0, 110.0, 0.35, 1.0,
                          [(0, 133, 0.5)], [(10, 1.2)], "Nessuna", False, None),
    'genta_neonato_pretermine': ('Gentamicina', models.AMINO_MODELLI['Gentamicina'][0], 20, 1.6, 42.0, 0.7, 1.0,
                                 [(0, 8, 0.5)], [(10, 3.0)], "Nessuna", False, 30.0),
    'amika_10h': ('Amikacina', models.AMINO_MODELLI['Amikacina'][0], 8 * 365, 25.0, 128.0, 0.4, 4.0,
                  [(0, 500, 0.5)], [(10, 3.5)], "Nessuna", False, None),
    'amika_llanos': ('Amikacina', models.AMINO_MODELLI['Amikacina'][1], 8 * 365, 25.0, 128.0, 0.4, 4.0,
                     [(0, 500, 0.5)], [(10, 3.5)], "Nessuna", False, None),
    'tobra_fc': ('Tobramicina', models.AMINO_MODELLI['Tobramicina'][1], 12 * 365, 38.0, 148.0, 0.5, 1.0,
                 [(0, 380, 0.5)], [(1, 28.0), (10, 2.0)], "Nessuna", False, None),
    'tobra_crrt': ('Tobramicina', models.AMINO_MODELLI['Tobramicina'][0], 12 * 365, 38.0, 148.0, 2.5, 1.0,
                   [(0, 280, 0.5)], [(10, 9.0)], "Nessuna", True, None),
    'dapto_bambino': ('Daptomicina', None, 10 * 365, 30.0, 138.0, 0.45, 0.5,
                      [(0, 300, 0.5), (24, 300, 0.5), (48, 300, 0.5), (72, 300, 0.5)], [(73, 95.0), (95.5, 12.0)],
                      "Nessuna", False, None),
    'dapto_IR_moderata': ('Daptomicina', None, 4 * 365, 16.0, 102.0, 1.1, 1.0,
                          [(0, 160, 1.0), (24, 160, 1.0)], [(47.5, 20.0)], "Moderata (GFR 30-60)", False, None),
    'dapto_HD_adolescente': ('Daptomicina', None, 15 * 365, 45.0, 160.0, 6.0, 1.0,
                             [(0, 300, 0.5), (48, 300, 0.5)], [(47.5, 25.0)], "Emodialisi (HD)", False, None),
    'dapto_adulto': ('Daptomicina', None, 60 * 365, 80.0, 175.0, 1.0, 1.0,
                     [(0, 800, 0.5), (24, 800, 0.5)], [(24.9, 110.0), (47.5, 18.0)], "Nessuna", False, None),
    'cipro_hirt': ('Ciprofloxacina', models.CIPRO_MODELLI[0], 8 * 365, 26.0, 128.0, 0.5, 0.25,
                   [(0, 260, 1.0), (8, 260, 1.0), (16, 260, 1.0)], [(17.5, 3.0), (23.5, 0.6)],
                   "Nessuna", False, None),
    'cipro_forrest_adulto': ('Ciprofloxacina', models.CIPRO_MODELLI[2], 55 * 365, 75.0, 172.0, 1.6, 0.5,
                             [(0, 400, 1.0), (12, 400, 1.0)], [(13.5, 4.5), (23.5, 1.4)], "Nessuna", False, None),
    'cipro_rajagopalan': ('Ciprofloxacina', models.CIPRO_MODELLI[1], 4 * 365, 16.0, 102.0, 0.3, 0.5,
                          [(0, 160, 1.0), (8, 160, 1.0)], [(9.5, 2.5)], "Nessuna", False, None),
}
for i, m in enumerate(models.DALBA_MODELLI):
    CASI[f'dalba_{i}'] = ('Dalbavancina', m, 10 * 365, 30.0, 138.0, 0.45, 0.25,
                          [(0, 675, 0.5), (168, 675, 0.5)], [(168, 30.0), (336, 25.0)], "Nessuna", False, None)


@pytest.mark.parametrize("nome", list(CASI))
def test_bayes_percorso_completo(nome):
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_function(_harness, args=(str(RADICE), *CASI[nome]), default_timeout=120)
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    assert not at.error, [e.value for e in at.error]
    R = at.session_state['b_risultati']
    par = R['par']
    assert all(v > 0 for v in par.values() if v)
    # il MAP non deve restare incollato ai limiti (0,05–10 × prior)
    for k in ('CL', 'V1'):
        r = par[k] / R['prior']['par'][k]
        assert 0.06 < r < 9.5, (k, r)
    # le previsioni individuali devono avvicinarsi ai livelli più del prior
    storia, liv = R['storia'], R['livelli']
    t = [l['ore'] for l in liv]
    y = [l['conc'] for l in liv]
    import numpy as np
    err_ind = np.abs(np.log(pk.conc(par, storia, t, R['prior']['t_inf'])) - np.log(y)).sum()
    err_pop = np.abs(np.log(pk.conc(R['prior']['par'], storia, t, R['prior']['t_inf'])) - np.log(y)).sum()
    assert err_ind <= err_pop + 1e-6
    referto = [c.value for c in at.code if 'Proposta' in c.value or 'Forecast' in c.value]
    assert referto, "testo del referto mancante"
    assert "None" not in referto[-1] and "nan" not in referto[-1]
    btn = [b for b in at.button if b.label.startswith("Registra")]
    assert btn
    btn[0].click()
    at.run()
    assert not at.exception
    assert len(at.session_state['storico_bayes']['monitoraggi']) == 1


def test_regime_aminoglicoside_proposto_rispetta_i_criteri():
    """La riga scelta nella tabella deve davvero avere picco nel range e valle sotto soglia."""
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_function(_harness, args=(str(RADICE), *CASI['amika_10h']), default_timeout=120)
    at.run()
    tab = [d.value for d in at.dataframe if 'Picco/MIC' in d.value.columns][0]
    scelte = tab[tab['OK'] == 'sì']
    assert len(scelte), tab
    r = scelte.iloc[0]
    assert r['Valle'] < models.AMINO_SOGLIE['Amikacina']['valle']
    assert r['Picco'] >= 0.95 * max(10 * 4.0, 20)


def _harness_render(radice):
    import sys
    sys.path.insert(0, radice)
    import streamlit as st
    from pigato_bayes import bayes_page as bp
    st.session_state.setdefault('b_farmaco', 'Gentamicina')
    bp.render()


def test_pagina_bayes_si_apre_e_chiede_i_dati():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_function(_harness_render, args=(str(RADICE),), default_timeout=60)
    at.run()
    assert not at.exception
    assert any("sesso" in w.value for w in at.warning)          # gate: codice e sesso
    at.sidebar.text_input[0].input("T1")
    [x for x in at.sidebar.selectbox if x.label == "Sesso biologico"][0].set_value("Femmina")
    at.run()
    assert not at.exception
    assert any("Campi da compilare" in w.value for w in at.warning)
    [b for b in at.button if b.label.startswith("Avvia")][0].click()
    at.run()
    assert not at.exception and any("Compila prima" in e.value for e in at.error)
