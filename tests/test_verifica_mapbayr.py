"""
Regressione contro mapbayr (R): 72 casi di daptomicina, ciprofloxacina, dalbavancina e aminoglicosidi.

I valori attesi in tests/dati/riferimento_mapbayr.json sono stati calcolati il 24/11/2026 con
mrgsolve 2.0.1 e mapbayr 0.10.2 (CRAN), con i modelli scritti come equazioni differenziali e l'AUC per integrazione.
I prior attesi coincidono con quelli ricalcolati in R dalle formule delle fonti.
Gli omega sono quelli pubblicati, passati allo stimatore senza modifiche come fa l'app.

Tolleranze: parametri 0,01% (scarto osservato ≤ 0,0005%); AUC24 0,5% (errore dei trapezi su
400 punti, osservato ≤ 0,19%); Cmin 0,2%; tempo sotto soglia della dalbavancina 0,01 h.
"""
import json
import math
from pathlib import Path

import pytest

from pigato_bayes import models
from pigato_bayes import pk_engine as pk

RIF = json.loads((Path(__file__).parent / 'dati' / 'riferimento_mapbayr.json').read_text(encoding='utf-8'))
CASI = RIF['casi']


def _prior(c):
    return getattr(models, c['chiamata']['funzione'])(**c['chiamata']['args'])


def _map(c, pr):
    lim = (1e-6, 1e6) if c['limiti_attivi_pigato'] else (0.05, 10.0)
    p, _ = pk.mappa_bayesiana(pr['par'], c['omega'], c['storia'], c['livelli'], pr['sigma_prop'],
                              pr['sigma_add'], pr['t_inf'], limiti=lim)
    return p


@pytest.mark.parametrize('c', CASI, ids=[c['id'] for c in CASI])
def test_prior(c):
    pr = _prior(c)
    for k, v in c['prior_atteso'].items():
        assert pr['par'][k] == pytest.approx(v, rel=1e-9, abs=1e-12), k


@pytest.mark.parametrize('c', CASI, ids=[c['id'] for c in CASI])
def test_map_contro_mapbayr(c):
    pr = _prior(c)
    p = _map(c, pr)
    for k, atteso in c['atteso_mapbayr']['par'].items():
        assert p[k] == pytest.approx(atteso, rel=1e-4), k


@pytest.mark.parametrize('c', CASI, ids=[c['id'] for c in CASI])
def test_auc_cmin_contro_mapbayr(c):
    pr = _prior(c)
    p = _map(c, pr)
    f = c['finestra']
    auc, _, cmin = pk.metriche_finestra(p, f['dosi'], f['start'], f['len'], pr['t_inf'])
    assert auc * f['fattore'] == pytest.approx(c['atteso_mapbayr']['AUC24'], rel=5e-3)
    assert cmin == pytest.approx(c['atteso_mapbayr']['Cmin'], rel=2e-3, abs=1e-9)


DALBA = [c for c in CASI if c.get('soglia')]


@pytest.mark.parametrize('c', DALBA, ids=[c['id'] for c in DALBA])
def test_tempo_sotto_soglia_dalbavancina(c):
    pr = _prior(c)
    p = _map(c, pr)
    t = pk.tempo_sotto_soglia(p, c['storia'], c['soglia'], c['t_val'], c['t_val'] + 2400, passo=2.0,
                              t_inf=pr['t_inf'])
    assert t == pytest.approx(c['atteso_mapbayr']['t_sotto_soglia'], abs=0.01)


# ---- casi limite segnalati dalla verifica

def test_livello_zero_bloccato():
    pr = models.prior_dapto('Maschio', 3300, 30.0, 0.5, False, 0.2, 2.0, 'Nessuna')
    st = [{'ore': 24 * i, 'dose': 270, 'tinf': 0.5} for i in range(3)]
    with pytest.raises(ValueError, match='maggiori di zero'):
        pk.mappa_bayesiana(pr['par'], pr['omega'], st, [{'ore': 71.8, 'conc': 0.0}], pr['sigma_prop'])


def test_avviso_parametro_al_limite():
    c = next(c for c in CASI if c['id'] == 'D23')       # picco 5 mg/L contro ~60 attesi: V a 10×
    pr = _prior(c)
    p, _, avvisi = pk.mappa_bayesiana(pr['par'], c['omega'], c['storia'], c['livelli'], pr['sigma_prop'],
                                      pr['sigma_add'], pr['t_inf'], con_avvisi=True)
    assert p['V1'] == pytest.approx(10 * pr['par']['V1'], rel=1e-6)
    assert any('limite superiore' in a and 'volume centrale' in a for a in avvisi)


def test_nessun_avviso_caso_regolare():
    c = next(c for c in CASI if c['id'] == 'D15')
    pr = _prior(c)
    _, _, avvisi = pk.mappa_bayesiana(pr['par'], c['omega'], c['storia'], c['livelli'], pr['sigma_prop'],
                                      pr['sigma_add'], pr['t_inf'], con_avvisi=True)
    assert avvisi == []


def test_beta_stabile_parametri_estremi():
    """Con k10·k21 << s² la vecchia formula (s − √d)/2 dava β = 0 e divisione per zero."""
    st = [{'ore': 0, 'dose': 1500, 'tinf': 0.5}]
    c = pk.conc_2cmt(1.06e-7, 1.74e7, 0.103, 1.51e-5, st, [1.0, 100.0])
    assert all(math.isfinite(x) for x in c)
    assert c[0] == pytest.approx(1500 / 1.74e7, rel=1e-3)
