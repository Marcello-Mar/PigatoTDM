"""
Verifica numerica contro valori calcolati a mano o noti dalle fonti.
Eseguire con:  python -m pytest tests -q
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pigato_bayes import models, pk_engine as pk, renal, rules  # noqa: E402

approx = pytest.approx


# ------------------------------------------------------------------ rene
def test_schwartz_e_ckid():
    assert renal.schwartz_bedside(140, 0.5, 10) == approx(0.413 * 140 / 0.5)
    k = 39.0 * 1.008 ** (10 - 12)                     # maschio < 12 anni
    assert renal.ckid_u25(140, 0.5, 10, "Maschio") == approx(k * 1.40 / 0.5)
    assert renal.ckid_u25(165, 0.8, 15, "Femmina") == approx(36.1 * 1.023 ** 3 * 1.65 / 0.8)
    assert renal.ckid_u25(180, 1.0, 20, "Maschio") == approx(50.8 * 1.80)
    assert renal.schwartz_neonatale(50, 0.4, 0.1) == approx(0.31 * 50 / 0.4)


def test_ckd_epi_2021_valori_noti():
    # Inker 2021: uomo 50 anni, Scr 1,0 → ~92; donna 50 anni, Scr 0,7 → ~105
    assert renal.ckd_epi_2021(50, 1.0, "Maschio") == approx(91.7, abs=0.5)
    assert renal.ckd_epi_2021(50, 0.7, "Femmina") == approx(105.3, abs=0.5)
    assert renal.ckd_epi_2021(17, 1.0, "Maschio") is None


def test_cockcroft_gault():
    assert renal.cockcroft_gault(60, 70, 1.0, "Maschio") == approx(80 * 70 / 72)
    assert renal.cockcroft_gault(60, 70, 1.0, "Femmina") == approx(80 * 70 / 72 * 0.85)


def test_formula_per_eta():
    assert renal.formula_raccomandata(0.05)[0] == "neonatale"
    assert renal.formula_raccomandata(0.5, eg_sett=32)[0] == "pretermine"
    assert renal.formula_raccomandata(10)[0] == "ckid"
    assert renal.formula_raccomandata(40)[0] == "ckdepi"
    eg = renal.egfr_scelta(10, "Maschio", 32, 140, 0.5)
    assert eg['chiave'] == "ckid"
    assert eg['assoluta'] == approx(eg['valore'] * math.sqrt(32 * 140 / 3600) / 1.73)


def test_ffm_e_creatinina_attesa():
    bmi = 70 / 1.75 ** 2
    assert renal.ffm(70, 175, 40, "Maschio") == approx(9270 * 70 / (6680 + 216 * bmi))
    assert renal.ffm(70, 175, 40, "Femmina") == approx(9270 * 70 / (8780 + 244 * bmi))
    # a 18 anni il fattore pediatrico di Al-Sallami tende a 1 (maschio)
    assert renal.ffm(70, 175, 17.99, "Maschio") == approx(renal.ffm(70, 175, 18, "Maschio"), rel=0.02)
    assert renal.creatinina_attesa_umol(1, "Maschio") == approx(-2.37330 + 23.93581, abs=1e-6)
    assert renal.creatinina_attesa_umol(30, "Femmina") == 69.5


# ------------------------------------------------------------------ motore PK
def test_1cmt_infusione_analitica():
    CL, V, D, ti = 5.0, 20.0, 100.0, 1.0
    ke = CL / V
    c_fine = D / ti / CL * (1 - math.exp(-ke * ti))
    st = [{'ore': 0, 'dose': D, 'tinf': ti}]
    assert pk.conc_1cmt(CL, V, st, 1.0)[0] == approx(c_fine)
    assert pk.conc_1cmt(CL, V, st, 5.0)[0] == approx(c_fine * math.exp(-ke * 4))
    assert pk.conc_1cmt(CL, V, st, 0.0)[0] == 0.0


def test_auc_infinito_uguale_dose_su_cl():
    st = [{'ore': 0, 'dose': 500, 'tinf': 0.5}]
    for par in ({'CL': 2.0, 'V1': 15.0}, {'CL': 2.0, 'V1': 10.0, 'Q': 3.0, 'V2': 25.0}):
        t = np.linspace(0, 400, 200001)
        c = pk.conc(par, st, t, 0.5)
        assert np.trapezoid(c, t) == approx(500 / 2.0, rel=1e-3)


def test_2cmt_continuo_con_1cmt_quando_q_tende_a_zero():
    st = [{'ore': 0, 'dose': 100, 'tinf': 1}]
    a = pk.conc_2cmt(3.0, 12.0, 1e-7, 30.0, st, [0.5, 4, 10])
    b = pk.conc_1cmt(3.0, 12.0, st, [0.5, 4, 10])
    assert a == approx(b, rel=1e-4)


def test_sovrapposizione_stato_stazionario():
    CL, V, D, tau, ti = 4.0, 20.0, 200.0, 12.0, 1.0
    ke = CL / V
    st = [{'ore': i * tau, 'dose': D, 'tinf': ti} for i in range(40)]
    t = 39 * tau + ti
    css_max = D / ti / CL * (1 - math.exp(-ke * ti)) / (1 - math.exp(-ke * tau))
    assert pk.conc_1cmt(CL, V, st, t)[0] == approx(css_max, rel=1e-6)
    auc, cmax, cmin = pk.metriche_finestra({'CL': CL, 'V1': V}, st, 39 * tau, tau, ti, n=4000)
    assert auc == approx(D / CL, rel=1e-3)


def test_omega_da_cv():
    assert pk.omega_da_cv(30) == approx(math.sqrt(math.log(1.09)))


def test_map_senza_livelli_restituisce_prior():
    prior = {'CL': 2.0, 'V1': 10.0}
    p, _ = pk.mappa_bayesiana(prior, {'CL': 0.3, 'V1': 0.3}, [{'ore': 0, 'dose': 100}], [], 0.1)
    assert p == prior


def test_map_recupera_parametri_veri():
    vero = {'CL': 4.0, 'V1': 12.0}
    st = [{'ore': i * 24, 'dose': 300, 'tinf': 0.5} for i in range(3)]
    t_obs = [49.0, 52.0, 58.0, 66.0]
    liv = [{'ore': t, 'conc': float(pk.conc(vero, st, t)[0])} for t in t_obs]
    p, _ = pk.mappa_bayesiana({'CL': 2.0, 'V1': 15.0}, {'CL': 0.5, 'V1': 0.5}, st, liv, 0.01)
    assert p['CL'] == approx(4.0, rel=0.02)
    assert p['V1'] == approx(12.0, rel=0.05)


def test_map_un_livello_si_sposta_nella_direzione_giusta():
    prior = {'CL': 3.0, 'V1': 15.0}
    st = [{'ore': 0, 'dose': 300, 'tinf': 0.5}]
    c_prior = float(pk.conc(prior, st, 10)[0])
    p_alto, _ = pk.mappa_bayesiana(prior, {'CL': 0.3, 'V1': 0.2}, st, [{'ore': 10, 'conc': 2 * c_prior}], 0.2)
    p_basso, _ = pk.mappa_bayesiana(prior, {'CL': 0.3, 'V1': 0.2}, st, [{'ore': 10, 'conc': 0.5 * c_prior}], 0.2)
    assert p_alto['CL'] < prior['CL'] < p_basso['CL']


def test_tempo_sotto_soglia_analitico():
    CL, V = 1.0, 10.0
    st = [{'ore': 0, 'dose': 100, 'tinf': 0.01}]
    c0 = float(pk.conc_1cmt(CL, V, st, 0.01)[0])
    atteso = 0.01 + math.log(c0 / 2.0) / (CL / V)
    t = pk.tempo_sotto_soglia({'CL': CL, 'V1': V}, st, 2.0, 0.5, 200, passo=0.5, t_inf=0.01)
    assert t == approx(atteso, abs=1e-3)


# ------------------------------------------------------------------ modelli di popolazione
def test_prior_ai_valori_di_riferimento():
    # Soraluce: ClCr CG = 49 → CL = CLNR + θ
    eta = 140 - 49 * 72 / 70                      # età che dà CG 49 con 70 kg, Scr 1
    pr = models.prior_dapto("Maschio", eta * 365.25, 70, 1.0, False, 0, 0, "Nessuna")
    assert pr['par']['CL'] == approx(0.160 + 0.367, rel=1e-3)
    assert pr['par']['V1'] == 12.30
    # Wei 2–6 anni: 20 mL/h/kg
    pr = models.prior_dapto("Maschio", 4 * 365.25, 16, 0.3, False, 0, 0, "Nessuna")
    assert pr['par']['CL'] == approx(0.020 * 16)
    # Hirt: 20 kg, eGFR 199, 7,3 anni → 9,07 L/h
    pr = models.prior_cipro("Maschio", 7.3 * 365.25, 20, 0.4, False, 0, 0, egfr=199,
                            modello=models.CIPRO_MODELLI[0])
    assert pr['par']['CL'] == approx(9.07, rel=1e-6)
    assert pr['par']['V1'] == approx(26.5)
    # Rajagopalan: 70 kg, 2,5 anni → 30,3 L/h
    pr = models.prior_cipro("Maschio", 2.5 * 365.25, 70, 0.4, False, 0, 0, egfr=None,
                            modello=models.CIPRO_MODELLI[1])
    assert pr['par']['CL'] == approx(30.3)
    # Forrest: CLT = (0,00145·CLCR + 0,167)·peso
    pr = models.prior_cipro("Maschio", 50 * 365.25, 70, 1.0, False, 0, 0, altezza=175, egfr=100,
                            modello=models.CIPRO_MODELLI[2])
    assert pr['par']['CL'] == approx((0.145 + 0.167) * 70)
    # dalbavancina, 70 kg = valori tipici
    for nome, m in models.DALBA_MODELLI.items():
        pr = models.prior_dalba(70, nome)
        assert pr['par']['CL'] == approx(m['CL'])
        assert pr['par']['V2'] == approx(m['V2'])


def test_llanos_paez_adulto_tipico():
    # FFM 70 kg, PMA molto alta, creatinina = attesa → CL = 5,77 × 1,12
    gg = 40 * 365.25
    eta = 40
    # cerca peso/altezza con FFM ≈ 70 per il maschio adulto
    peso, alt = 90.0, 190.0
    f = renal.ffm(peso, alt, eta, "Maschio")
    scr = renal.creatinina_attesa_umol(eta, "Maschio") / renal.MGDL_TO_UMOL
    par, _, _ = models._llanos_paez(gg, peso, alt, scr, "Maschio", None)
    pma = renal.pma_settimane(gg)
    mat = pma ** 3.33 / (55.4 ** 3.33 + pma ** 3.33)
    assert par['CL'] == approx(5.77 * 1.12 * (f / 70) ** 0.75 * mat, rel=1e-9)
    assert par['V1'] == approx(21.6 * f / 70)


def test_alhadab_hennig_70kg():
    pr = models.prior_amino('Amikacina', models.AMINO_MODELLI['Amikacina'][0], "Maschio",
                            10 * 365.25, 70, 170, 0.6, False, 0, 0)
    assert pr['par'] == approx({'CL': 11.1, 'V1': 30.2, 'Q': 4.26, 'V2': 14.9})
    pr = models.prior_amino('Tobramicina', models.AMINO_MODELLI['Tobramicina'][1], "Maschio",
                            10 * 365.25, 70, 170, 0.6, False, 0, 0)
    assert pr['par'] == approx({'CL': 6.37, 'V1': 18.7, 'Q': 0.39, 'V2': 1.32})


def test_crrt_aggiunge_clearance_extracorporea():
    a = models.prior_dapto("Maschio", 50 * 365.25, 70, 1.0, False, 0, 0, "Nessuna")
    b = models.prior_dapto("Maschio", 50 * 365.25, 70, 1.0, True, 0.5, 2.0, "Nessuna")
    assert b['par']['CL'] - a['par']['CL'] == approx(1.0)


def test_soglie_dalbavancina():
    assert models.dalba_soglia_auc(111.1, 0.25, 0.07) == approx(16.53, abs=0.01)
    assert models.dalba_soglia_4xmic(0.25, 0.07) == approx(14.29, abs=0.01)


# ------------------------------------------------------------------ regole del referto
@pytest.mark.parametrize("r,az,pct", [(1.9, 'AUMENTA', 50), (3.9, 'AUMENTA', 25), (4.0, 'MANTIENI', 0),
                                      (7.99, 'MANTIENI', 0), (9.0, 'RIDUCI', -25), (12, 'RIDUCI', -50)])
def test_fasce_ecpa(r, az, pct):
    assert rules.classifica_fasce(r, rules.FASCE_CSS) == (az, pct)


def test_sawchuk_zaske_contro_simulazione():
    CL, V, D, tau, ti = 6.0, 20.0, 1000.0, 12.0, 1.0
    par = {'CL': CL, 'V1': V}
    st = [{'ore': i * tau, 'dose': D, 'tinf': ti} for i in range(30)]
    t0 = 29 * tau
    t1, t2 = 2.0, 11.5
    c1, c2 = (float(pk.conc(par, st, t0 + t)[0]) for t in (t1, t2))
    auc24, cmax, cmin, ke = rules.auc_da_due_livelli(c1, t1, c2, t2, tau, ti)
    assert ke == approx(CL / V)
    assert cmax == approx(float(pk.conc(par, st, t0 + ti)[0]), rel=1e-6)
    assert auc24 == approx(D / CL * 24 / tau, rel=0.01)


def test_arrotonda_dose():
    assert rules.arrotonda_dose(16100) == 16000
    assert rules.arrotonda_dose(1130) == 1250
    assert rules.arrotonda_dose(612) == 600
    assert rules.arrotonda_dose(73) == 75


def test_regole_coerenti_con_tdm_data():
    """Le fasce di azione devono seguire le regole 'dec'/'inc' di tdm_data.json."""
    tei = rules.REGOLE['teicoplanin']['fasce']
    assert rules.classifica_fasce(32, tei)[0] == 'SOSPENDI'          # "If Cmin >30: hold"
    assert rules.classifica_fasce(8, tei)[0] == 'AUMENTA'            # "If Cmin <10: increase"
    lzd = rules.REGOLE['linezolid']['fasce']
    assert rules.classifica_fasce(0.5, lzd) == ('AUMENTA', 50)
    assert rules.classifica_fasce(10, lzd) == ('RIDUCI', -25)
    assert rules.classifica_fasce(16, lzd) == ('RIDUCI', -50)
    assert rules.REGOLE['levofloxacin']['valle_max'] == 2          # "Space to q36–48h if Cmin > 2"


# ------------------------------------------------------------------ cistatina C
def test_ckid_u25_cistatina_pierce_tabella2():
    c = renal.ckid_u25_cys
    assert c(1.0, 15, "Maschio") == approx(87.2)
    assert c(1.0, 1, "Maschio") == approx(87.2 * 1.011 ** -14)          # ≈ 74,8 nel testo
    assert c(1.0, 1, "Maschio") == approx(74.8, abs=0.1)
    assert c(1.0, 18, "Maschio") == approx(77.1, abs=0.1)              # 87,2 × 0,960³
    assert c(1.0, 20, "Maschio") == approx(77.1)
    assert c(1.0, 12, "Femmina") == approx(79.9)
    assert c(1.0, 1, "Femmina") == approx(76.5, abs=0.1)
    assert c(1.0, 18, "Femmina") == approx(68.3, abs=0.1)
    assert c(2.0, 15, "Maschio") == approx(87.2 / 2)
    assert c(1.0, 0.5, "Maschio") is None


def test_ckd_epi_cistatina_2012():
    assert renal.ckd_epi_cys_2012(0.8, 40, "Maschio") == approx(133 * 0.996 ** 40)
    assert renal.ckd_epi_cys_2012(1.6, 40, "Femmina") == approx(133 * 2 ** -1.328 * 0.996 ** 40 * 0.932)
    assert renal.ckd_epi_cys_2012(0.6, 60, "Maschio") == approx(133 * 0.75 ** -0.499 * 0.996 ** 60)


def test_massa_muscolare_ridotta_usa_la_cistatina():
    # 12 anni, creatinina 0,20 (muscolo scarso): la creatinina dà ~270, la cistatina ~94
    eg = renal.egfr_scelta(12, "Maschio", 25, 140, 0.20, cys_mgl=0.9, muscolo_ridotto=True)
    assert eg['chiave'] == 'ckid_cys' and eg['valore'] == approx(renal.ckid_u25_cys(0.9, 12, "Maschio"))
    assert eg['tutte']['ckid'] > 2.5 * eg['valore']
    # la creatinina equivalente ridà, con CKiD U25 a creatinina, lo stesso eGFR
    assert renal.ckid_u25(140, eg['crea_modelli'], 12, "Maschio") == approx(eg['valore'], rel=1e-6)
    # senza la casella si resta sulla creatinina, anche se la cistatina c'è
    eg2 = renal.egfr_scelta(12, "Maschio", 25, 140, 0.20, cys_mgl=0.9)
    assert eg2['chiave'] == 'ckid' and eg2['crea_modelli'] == 0.20 and 'ckid_cys' in eg2['tutte']
    # massa ridotta senza cistatina: avviso, nessuna correzione
    eg3 = renal.egfr_scelta(12, "Maschio", 25, 140, 0.20, muscolo_ridotto=True)
    assert eg3['chiave'] == 'ckid' and any('cistatina' in a for a in eg3['avvisi'])
    # glucocorticoidi: avviso di sottostima solo se la stima è da cistatina
    eg4 = renal.egfr_scelta(12, "Maschio", 25, 140, 0.20, None, 0.9, True, True)
    assert any('Glucocorticoidi' in a for a in eg4['avvisi'])
    # adulto > 25 anni: CKD-EPI 2012 a cistatina, creatinina equivalente su CKD-EPI 2021
    eg5 = renal.egfr_scelta(40, "Femmina", 45, 160, 0.3, cys_mgl=1.1, muscolo_ridotto=True)
    assert eg5['chiave'] == 'ckdepi_cys'
    assert renal.ckd_epi_2021(40, eg5['crea_modelli'], "Femmina") == approx(eg5['valore'], rel=1e-6)
    # sotto 1 anno: nessuna formula a cistatina
    eg6 = renal.egfr_scelta(0.5, "Maschio", 7, 65, 0.15, cys_mgl=1.2, muscolo_ridotto=True)
    assert eg6['chiave'] == 'neonatale' and any('1 anno' in a for a in eg6['avvisi'])
