"""
Funzione renale e descrittori di taglia corporea.

Le equazioni eGFR e la scelta della formula raccomandata per età sono
le stesse di PIGATO TDM (sezione eGFR), così che le due app diano lo
stesso numero sullo stesso paziente.

Riferimenti
  Schwartz GJ, et al. J Am Soc Nephrol 2009;20:629-637 (bedside, k=0.413)
  Pierce CB, et al. Kidney Int 2021;99:948-956 (CKiD U25)
  Schwartz neonatale k=0.31 (<1 anno), come in PIGATO TDM
  Wilhelm-Bals, stima peso/creatinina nel pretermine, come in PIGATO TDM
  Cockcroft DW, Gault MH. Nephron 1976;16:31-41
  Inker LA, et al. N Engl J Med 2021;385:1737-1749 (CKD-EPI 2021)
  Pierce CB, et al. Kidney Int 2021;99:948-956, tab. 2 (CKiD U25 cistatina C)
  Inker LA, et al. N Engl J Med 2012;367:20-29, tab. 2 (CKD-EPI cistatina C 2012)
  Mosteller RD. N Engl J Med 1987;317:1098 (BSA)
  Al-Sallami HS, et al. Clin Pharmacokinet 2015;54:1169-1178 (FFM pediatrica)
  Janmahasatian S, et al. Clin Pharmacokinet 2005;44:1051-1065 (FFM adulto)
  Ceriotti F, et al. Clin Chem 2008;54:559-566 (creatinina attesa per età;
    forma analitica come implementata in clinPK::calc_creat, con
    interpolazione 15-18 anni da Johansen A, et al. Ther Drug Monit 2011)
"""

import math

MGDL_TO_UMOL = 88.4


def bsa_mosteller(peso_kg, altezza_cm):
    if not peso_kg or not altezza_cm:
        return None
    return math.sqrt(peso_kg * altezza_cm / 3600.0)


def schwartz_bedside(altezza_cm, scr_mgdl, eta_anni):
    if not (altezza_cm and scr_mgdl and eta_anni is not None) or not (1 <= eta_anni <= 18):
        return None
    return 0.413 * altezza_cm / scr_mgdl


def ckid_u25(altezza_cm, scr_mgdl, eta_anni, sesso):
    if not (altezza_cm and scr_mgdl and eta_anni is not None) or not (1 <= eta_anni <= 25):
        return None
    maschio = (sesso == "Maschio")
    if eta_anni < 12:
        k = (39.0 if maschio else 36.1) * (1.008 ** (eta_anni - 12))
    elif eta_anni < 18:
        k = 39.0 * (1.045 ** (eta_anni - 12)) if maschio else 36.1 * (1.023 ** (eta_anni - 12))
    else:
        k = 50.8 if maschio else 41.4
    return k * (altezza_cm / 100.0) / scr_mgdl


def schwartz_neonatale(altezza_cm, scr_mgdl, eta_anni):
    if not (altezza_cm and scr_mgdl and eta_anni is not None) or eta_anni >= 1:
        return None
    return 0.31 * altezza_cm / scr_mgdl


def pretermine_peso_creatinina(peso_kg, scr_mgdl, eg_sett):
    if not (peso_kg and scr_mgdl and eg_sett) or eg_sett >= 37:
        return None
    return 2.32 * (peso_kg * 1000.0) / (scr_mgdl * MGDL_TO_UMOL)


def cockcroft_gault(eta_anni, peso_kg, scr_mgdl, sesso):
    """mL/min, non indicizzata. Solo adulto."""
    if not (peso_kg and scr_mgdl and eta_anni is not None) or eta_anni < 18:
        return None
    v = ((140 - eta_anni) * peso_kg) / (72 * scr_mgdl)
    return v * 0.85 if sesso == "Femmina" else v


def ckd_epi_2021(eta_anni, scr_mgdl, sesso):
    if not (scr_mgdl and eta_anni is not None) or eta_anni < 18:
        return None
    f = (sesso == "Femmina")
    kap = 0.7 if f else 0.9
    al = -0.241 if f else -0.302
    return (142 * min(scr_mgdl / kap, 1) ** al * max(scr_mgdl / kap, 1) ** -1.2
            * 0.9938 ** eta_anni * (1.012 if f else 1.0))


def ckid_u25_cys(cys_mgl, eta_anni, sesso):
    """CKiD U25 a cistatina C (IFCC), 1–25 anni: eGFR = K / cisC."""
    if not (cys_mgl and eta_anni is not None) or not (1 <= eta_anni <= 25):
        return None
    if sesso == "Maschio":
        if eta_anni < 15:
            k = 87.2 * 1.011 ** (eta_anni - 15)
        elif eta_anni < 18:
            k = 87.2 * 0.960 ** (eta_anni - 15)
        else:
            k = 77.1
    else:
        if eta_anni < 12:
            k = 79.9 * 1.004 ** (eta_anni - 12)
        elif eta_anni < 18:
            k = 79.9 * 0.974 ** (eta_anni - 12)
        else:
            k = 68.3
    return k / cys_mgl


def ckd_epi_cys_2012(cys_mgl, eta_anni, sesso):
    """CKD-EPI 2012 a cistatina C, adulto."""
    if not (cys_mgl and eta_anni is not None) or eta_anni < 18:
        return None
    x = cys_mgl / 0.8
    return (133 * min(x, 1) ** -0.499 * max(x, 1) ** -1.328 * 0.996 ** eta_anni
            * (0.932 if sesso == "Femmina" else 1.0))


def formula_raccomandata(eta_anni, eg_sett=None, eta_giorni=None):
    """Stessa logica di PIGATO TDM. Restituisce (chiave, etichetta)."""
    if eg_sett and eg_sett < 37 and eta_anni is not None and eta_anni < 1:
        return "pretermine", "Pretermine (<37 sett.): stima peso/creatinina; formule standard non validate"
    if eta_anni is None:
        return None, "Età mancante"
    if eta_anni < 1 / 12:
        return "neonatale", "Neonato (<1 mese): Schwartz neonatale k=0,31, con cautela; seguire l'andamento"
    if eta_anni < 1:
        return "neonatale", "Lattante (<1 anno): Schwartz neonatale k=0,31; interpretare con cautela"
    if eta_anni < 25:
        return "ckid", "1–25 anni: CKiD U25 (KDIGO 2024); in alternativa Schwartz bedside"
    return "ckdepi", "Adulto: CKD-EPI 2021 senza coefficiente etnico"


ETICHETTE = {
    "schwartz": "Schwartz bedside 2009",
    "ckid": "CKiD U25",
    "neonatale": "Schwartz neonatale (k=0,31)",
    "pretermine": "Pretermine peso/creatinina",
    "cg": "Cockcroft-Gault",
    "ckdepi": "CKD-EPI 2021",
    "ckid_cys": "CKiD U25 cistatina C",
    "ckdepi_cys": "CKD-EPI 2012 cistatina C",
}
CISTATINA = ("ckid_cys", "ckdepi_cys")


def egfr_tutte(eta_anni, sesso, peso_kg, altezza_cm, scr_mgdl, eg_sett=None, cys_mgl=None):
    """Tutte le stime applicabili. Valori indicizzati in mL/min/1,73 m², tranne CG (mL/min)."""
    r = {
        "schwartz": schwartz_bedside(altezza_cm, scr_mgdl, eta_anni),
        "ckid": ckid_u25(altezza_cm, scr_mgdl, eta_anni, sesso),
        "neonatale": schwartz_neonatale(altezza_cm, scr_mgdl, eta_anni),
        "pretermine": pretermine_peso_creatinina(peso_kg, scr_mgdl, eg_sett),
        "cg": cockcroft_gault(eta_anni, peso_kg, scr_mgdl, sesso),
        "ckdepi": ckd_epi_2021(eta_anni, scr_mgdl, sesso),
        "ckid_cys": ckid_u25_cys(cys_mgl, eta_anni, sesso),
        "ckdepi_cys": ckd_epi_cys_2012(cys_mgl, eta_anni, sesso) if (eta_anni or 0) >= 25 else None,
    }
    return {k: v for k, v in r.items() if v is not None}


def _egfr_creatinina(chiave, eta_anni, sesso, peso_kg, altezza_cm, scr, eg_sett):
    return {"schwartz": lambda: schwartz_bedside(altezza_cm, scr, eta_anni),
            "ckid": lambda: ckid_u25(altezza_cm, scr, eta_anni, sesso),
            "neonatale": lambda: schwartz_neonatale(altezza_cm, scr, eta_anni),
            "pretermine": lambda: pretermine_peso_creatinina(peso_kg, scr, eg_sett),
            "ckdepi": lambda: ckd_epi_2021(eta_anni, scr, sesso)}[chiave]()


def creatinina_equivalente(chiave, eta_anni, sesso, peso_kg, altezza_cm, eg_sett, egfr_obiettivo):
    """
    Creatinina che, con la formula a creatinina `chiave`, darebbe `egfr_obiettivo`.
    Serve a passare ai modelli PK (che usano la creatinina) la funzione renale
    stimata con la cistatina quando la creatinina non è attendibile.
    """
    lo, hi = 0.02, 30.0
    if _egfr_creatinina(chiave, eta_anni, sesso, peso_kg, altezza_cm, lo, eg_sett) is None:
        return None
    for _ in range(80):                      # eGFR decresce con la creatinina
        mid = math.sqrt(lo * hi)
        if _egfr_creatinina(chiave, eta_anni, sesso, peso_kg, altezza_cm, mid, eg_sett) > egfr_obiettivo:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def egfr_scelta(eta_anni, sesso, peso_kg, altezza_cm, scr_mgdl, eg_sett=None, cys_mgl=None,
                muscolo_ridotto=False, glucocorticoidi=False):
    """
    Restituisce dict con la formula raccomandata, il valore indicizzato,
    la clearance assoluta (eGFR × BSA/1,73), gli avvisi e la creatinina da
    passare ai modelli PK ('crea_modelli'), oppure None se non calcolabile.

    Con massa muscolare ridotta la creatinina sovrastima il filtrato: se c'è la
    cistatina C la formula raccomandata diventa quella a cistatina (CKiD U25 fino
    a 25 anni, CKD-EPI 2012 oltre) e ai modelli si passa la creatinina equivalente.
    """
    tutte = egfr_tutte(eta_anni, sesso, peso_kg, altezza_cm, scr_mgdl, eg_sett, cys_mgl)
    chiave, nota = formula_raccomandata(eta_anni, eg_sett)
    avvisi = []
    chiave_cr = chiave
    if chiave not in tutte:
        # ripiego: prima stima disponibile nell'ordine di preferenza
        for k in ("ckid", "schwartz", "ckdepi", "neonatale", "pretermine"):
            if k in tutte:
                chiave = chiave_cr = k
                break
        else:
            chiave = chiave_cr = None
    crea_modelli = scr_mgdl
    if muscolo_ridotto:
        k_cys = "ckid_cys" if (eta_anni or 0) < 25 else "ckdepi_cys"
        if eta_anni is not None and eta_anni < 1:
            avvisi.append("Massa muscolare ridotta sotto 1 anno: nessuna formula a cistatina C validata "
                          "(CKiD U25 da 1 anno); la stima da creatinina va letta come limite alto.")
        elif k_cys in tutte:
            nota = ("Massa muscolare ridotta: la creatinina sovrastima il filtrato, si usa la cistatina C "
                    + ("(CKiD U25)" if k_cys == "ckid_cys" else "(CKD-EPI 2012)")
                    + " senza formule combinate.")
            chiave = k_cys
            if chiave_cr:
                eq = creatinina_equivalente(chiave_cr, eta_anni, sesso, peso_kg, altezza_cm, eg_sett, tutte[k_cys])
                if eq:
                    crea_modelli = eq
                    avvisi.append(f"Ai modelli PK si passa la creatinina equivalente alla stima da cistatina "
                                  f"({eq:.2f} mg/dL invece di {scr_mgdl:.2f})".replace(".", ","))
        else:
            avvisi.append("Massa muscolare ridotta: la stima da creatinina sovrastima il filtrato (anche di "
                          "più volte). Dosare la cistatina C; con farmaci a basso indice terapeutico dare più "
                          "peso ai livelli che all'eGFR.")
    if glucocorticoidi and chiave in CISTATINA:
        avvisi.append("Glucocorticoidi sistemici: aumentano la cistatina C indipendentemente dal filtrato, "
                      "quindi questa stima può sottostimare la funzione renale.")
    if chiave is None:
        return {"tutte": tutte, "chiave": None, "valore": None, "assoluta": None, "nota": nota,
                "avvisi": avvisi, "crea_modelli": crea_modelli}
    bsa = bsa_mosteller(peso_kg, altezza_cm)
    val = tutte[chiave]
    ass = val * bsa / 1.73 if bsa else None
    return {"tutte": tutte, "chiave": chiave, "etichetta": ETICHETTE[chiave],
            "valore": val, "assoluta": ass, "bsa": bsa, "nota": nota, "avvisi": avvisi,
            "crea_modelli": crea_modelli, "creatinina": scr_mgdl, "cistatina": cys_mgl, "muscolo_ridotto": bool(muscolo_ridotto),
            "glucocorticoidi": bool(glucocorticoidi)}


# ------------------------------------------------------------------
# Descrittori per i modelli degli aminoglicosidi
# ------------------------------------------------------------------

def bmi(peso_kg, altezza_cm):
    if not peso_kg or not altezza_cm:
        return None
    return peso_kg / (altezza_cm / 100.0) ** 2


def ffm(peso_kg, altezza_cm, eta_anni, sesso):
    """Fat-free mass (kg): Al-Sallami 2015 sotto i 18 anni, Janmahasatian 2005 nell'adulto."""
    b = bmi(peso_kg, altezza_cm)
    if b is None:
        return None
    if sesso == "Femmina":
        base = 9270.0 * peso_kg / (8780.0 + 244.0 * b)
        if eta_anni is not None and eta_anni < 18:
            return (1.11 + (1 - 1.11) / (1 + (max(eta_anni, 0.01) / 7.1) ** -1.1)) * base
        return base
    base = 9270.0 * peso_kg / (6680.0 + 216.0 * b)
    if eta_anni is not None and eta_anni < 18:
        return (0.88 + (1 - 0.88) / (1 + (max(eta_anni, 0.01) / 13.4) ** -12.7)) * base
    return base


def creatinina_attesa_umol(eta_anni, sesso):
    """Creatinina sierica media attesa per età (µmol/L), Ceriotti 2008 (forma clinPK)."""
    a = max(eta_anni, 0.05)
    if a <= 15:
        return -2.37330 - 12.91367 * math.log(a) + 23.93581 * math.sqrt(a)
    if a < 18:
        return 9.5471 * a - 87.847 if sesso == "Maschio" else 4.7137 * a - 15.347
    return 84.0 if sesso == "Maschio" else 69.5


def pma_settimane(eta_giorni, eg_sett=None):
    """Età postmestruale in settimane (EG alla nascita, default 40)."""
    return (eg_sett or 40.0) + eta_giorni / 7.0
