"""
Regole di interpretazione per il modulo referti (32 farmaci di tdm_data.json).

Le soglie numeriche vengono da tdm_data.json (Gatti et al. Crit Care
2022;26:178; Gatti & Pea, Expert Rev Clin Pharmacol 2023;16:1035, programma
ECPA, più le sintesi per agente citate in 'src'). Qui sono solo tradotte in
fasce eseguibili. La proposta è sempre un suggerimento: l'operatore la
conferma o la modifica prima di generare il referto.

Tipi di regola
  css   rapporto Css (o Cmin) / MIC con le fasce ECPA dei beta-lattamici
  range concentrazione assoluta (Cmin o Css) in fasce con azione
  auc   AUC24 da Css continua o da due livelli; target AUC/MIC
  amino picco/MIC e valle (monosomministrazione)
  fq    fluorochinoloni: Cmax/MIC, AUC/MIC, valle
  fauc  fosfomicina: fAUC24/MIC per patogeno
  smx   cotrimossazolo: picchi di TMP e SMX
  info  nessuna regola di concentrazione validata
"""

import math

# Fasce ECPA dei beta-lattamici (tdm_data.json, 'dec'/'inc')
FASCE_CSS = [
    (0, 2, 'AUMENTA', +50),
    (2, 4, 'AUMENTA', +25),
    (4, 8, 'MANTIENI', 0),
    (8, 10, 'RIDUCI', -25),
    (10, math.inf, 'RIDUCI', -50),
]

REGOLE = {
    # ---------------- beta-lattamici tradizionali ----------------
    'piperacillin-tazobactam': {'tipo': 'css', 'analita': 'piperacillina', 'tox': 361,
                                'comp': {'nome': 'tazobactam', 'CT': 4.0}},
    'meropenem': {'tipo': 'css', 'analita': 'meropenem', 'tox': 64.2,
                  'nota': 'VAP: restringere a Css/MIC 6–8.'},
    'ceftazidime': {'tipo': 'css', 'analita': 'ceftazidima', 'tox': None},
    'cefepime': {'tipo': 'css', 'analita': 'cefepime', 'tox': 36},
    'ampicillin': {'tipo': 'css', 'analita': 'ampicillina', 'tox': None},
    'ampicillin-sulbactam': {'tipo': 'css', 'analita': 'ampicillina', 'tox': None},
    # ---------------- nuovi beta-lattamici ----------------
    'avibactam-ceftazidime': {'tipo': 'css', 'analita': 'ceftazidima', 'tox': None,
                              'comp': {'nome': 'avibactam', 'CT': 4.0}},
    'meropenem-vaborbactam': {'tipo': 'css', 'analita': 'meropenem', 'tox': 64.2,
                              'nota': 'Vaborbactam: target AUC/MIC > 24 (non valutato in automatico).'},
    'cefiderocol': {'tipo': 'css', 'analita': 'cefiderocol', 'tox': None},
    'ceftobiprole': {'tipo': 'css', 'analita': 'ceftobiprolo', 'tox': None},
    'ceftaroline': {'tipo': 'css', 'analita': 'ceftarolina', 'tox': None,
                    'nota': 'Indice fT>MIC: fasce ECPA applicate a Cmin/Css per analogia; '
                            'nessuna soglia di tossicità validata.'},
    'ceftolozane-tazobactam': {'tipo': 'css', 'analita': 'ceftolozano', 'tox': None,
                               'comp': {'nome': 'tazobactam', 'CT': 4.0}},
    'avibactam-aztreonam': {'tipo': 'info', 'analita': 'aztreonam',
                            'nota': 'Nessuna regola di concentrazione validata: aztreonam 60% fT>MIC '
                                    '(MIC nominale 8) + avibactam 50% fT>CT 2,5 mg/L; dosaggio guidato '
                                    'dalla funzione renale.'},
    # ---------------- concentrazione assoluta ----------------
    'linezolid': {'tipo': 'range', 'analita': 'linezolid', 'livello': 'Cmin', 'target': (2, 8),
                  'fasce': [(0, 1, 'AUMENTA', +50), (1, 2, 'AUMENTA', +25), (2, 8, 'MANTIENI', 0),
                            (8, 15, 'RIDUCI', -25), (15, math.inf, 'RIDUCI', -50)],
                  'tox': 8, 'tox_txt': 'trombocitopenia'},
    'fluconazole': {'tipo': 'range', 'analita': 'fluconazolo', 'livello': 'Cmin', 'target': (10, 20),
                    'fasce': [(0, 10, 'AUMENTA', +25), (10, 30, 'MANTIENI', 0),
                              (30, 50, 'RIDUCI', -25), (50, math.inf, 'RIDUCI', -50)]},
    'voriconazole': {'tipo': 'range', 'analita': 'voriconazolo', 'livello': 'Cmin', 'target': (1, 3),
                     'fasce': [(0, 1, 'AUMENTA', +25), (1, 3.5, 'MANTIENI', 0),
                               (3.5, 8, 'RIDUCI', -25), (8, math.inf, 'SOSPENDI', -100)],
                     'tox': 3.5, 'tox_txt': 'epatotossicità e neurotossicità',
                     'nota': 'Cmin < 1: in alternativa accorciare a ogni 6–8 h. Cmin 3,5–8: −25/−50%.'},
    'posaconazole': {'tipo': 'range', 'analita': 'posaconazolo', 'livello': 'Cmin', 'target': (1, 3),
                     'fasce': [(0, 1, 'AUMENTA', +25), (1, 3, 'MANTIENI', 0), (3, math.inf, 'RIDUCI', -25)],
                     'tox': 3, 'tox_txt': 'pseudoiperaldosteronismo',
                     'nota': 'Cmin < 1: valutare frazionamento ogni 12 h; Cmin > 3: −25/−50%.'},
    'isavuconazole': {'tipo': 'range', 'analita': 'isavuconazolo', 'livello': 'Cmin', 'target': (1, 5.1),
                      'fasce': [(0, 1, 'AUMENTA', +25), (1, 5.1, 'MANTIENI', 0), (5.1, math.inf, 'RIDUCI', -25)],
                      'tox': 5.1, 'tox_txt': 'disturbi gastrointestinali'},
    'ganciclovir-ganciclovir': {'tipo': 'range', 'analita': 'ganciclovir', 'livello': 'Cmin', 'target': (0.7, 2),
                                'fasce': [(0, 0.5, 'AUMENTA', +25), (0.5, 2, 'MANTIENI', 0),
                                          (2, 5, 'RIDUCI', -25), (5, math.inf, 'SOSPENDI', -100)],
                                'nota': 'Cmin < 0,5: in alternativa accorciare a ogni 6–8 h; tra 0,5 e 0,7 '
                                        'sotto target, rivalutare. Cmin 2–5: −25/−50%.'},
    'aciclovir': {'tipo': 'range', 'analita': 'aciclovir', 'livello': 'Cmin', 'target': (1, 3),
                  'fasce': [(0, 0.5, 'AUMENTA', +25), (0.5, 3, 'MANTIENI', 0),
                            (3, 5, 'RIDUCI', -25), (5, math.inf, 'SOSPENDI', -100)],
                  'nota': 'Cmin < 0,5: in alternativa accorciare a ogni 6 h; tra 0,5 e 1 sotto target, '
                          'rivalutare. Cmin 3–5: −25/−50%.'},
    'teicoplanin': {'tipo': 'range', 'analita': 'teicoplanina', 'livello': 'Cmin', 'target': (15, 30),
                    'fasce': [(0, 10, 'AUMENTA', +25), (10, 15, 'AUMENTA', +25), (15, 30, 'MANTIENI', 0),
                              (30, 40, 'RIDUCI', -25), (40, math.inf, 'SOSPENDI', -100)],
                    'tox': 40, 'tox_txt': 'nefrotossicità',
                    'nota': 'Target 15–30 mg/L (MRSA grave), 20–40 mg/L (infezioni complicate); '
                            'Cmin > 30: sospendere, ridosare quando 10–30 mg/L. Nell\'insufficienza '
                            'renale preferire l\'allungamento dell\'intervallo.'},
    'colistin': {'tipo': 'range', 'analita': 'colistina', 'livello': 'Css', 'target': (2, 2.5),
                 'fasce': [(0, 2, 'AUMENTA', +25), (2, 2.5, 'MANTIENI', 0), (2.5, math.inf, 'RIDUCI', -25)],
                 'tox': 2.5, 'tox_txt': 'nefrotossicità',
                 'nota': 'Css,avg ≈ valle al giorno 3–4. Target pratico ~2 mg/L (AUCss,24 ~50 mg·h/L).'},
    'dalbavancin': {'tipo': 'range', 'analita': 'dalbavancina', 'livello': 'Cmin', 'target': (8.04, None),
                    'fasce': [(0, 8.04, 'AUMENTA', 0), (8.04, math.inf, 'MANTIENI', 0)],
                    'nota': 'Farmaco a lunga emivita: per la data della prossima dose usare il modulo '
                            'bayesiano (forecast). Sotto target: dose aggiuntiva da 1500 mg (adulto).'},
    # ---------------- AUC ----------------
    'vancomycin': {'tipo': 'auc', 'analita': 'vancomicina', 'auc_mic': (400, 600), 'tox_auc': 600,
                   'tox_valle': 20},
    'daptomycin': {'tipo': 'auc', 'analita': 'daptomicina', 'auc_mic': (666, 939), 'tox_valle': 24.3,
                   'nota': 'AUC24/MIC ≥ 666 (S. aureus); E. faecium fAUC/MIC 27–143. '
                           'Per la stima bayesiana usare il modulo dedicato.'},
    # ---------------- aminoglicosidi ----------------
    'gentamicin': {'tipo': 'amino', 'analita': 'gentamicina', 'picco_mic': 8, 'valle_max': 1.0,
                   'valle_tox': 2.0},
    'amikacin': {'tipo': 'amino', 'analita': 'amikacina', 'picco_mic': 8, 'valle_max': 2.5,
                 'valle_tox': 5.0},
    'tobramycin': {'tipo': 'amino', 'analita': 'tobramicina', 'picco_mic': 10, 'valle_max': 1.0,
                   'valle_tox': 2.0},
    # ---------------- fluorochinoloni ----------------
    'levofloxacin': {'tipo': 'fq', 'analita': 'levofloxacina', 'cmax_mic': 10, 'auc_mic': 125, 'valle_max': 3,
                     'nota': 'Valle > 2 mg/L: distanziare a ogni 36–48 h.'},
    'ciprofloxacin': {'tipo': 'fq', 'analita': 'ciprofloxacina', 'cmax_mic': 10, 'auc_mic': 125, 'valle_max': 2},
    # ---------------- altri ----------------
    'fosfomycin': {'tipo': 'fauc', 'analita': 'fosfomicina', 'fu': 1.0,
                   'target_per_gruppo': {'s_aureus': (20, 75), 'enterobacterales': (83.3, None),
                                         'p_aeruginosa': (40.8, None)}},
    'sulfamethoxazole-trimethoprim': {'tipo': 'smx', 'analita': 'cotrimossazolo',
                                      'tmp': (5, 15), 'smx': (100, 200)},
}


def classifica_fasce(valore, fasce):
    for lo, hi, azione, pct in fasce:
        if lo <= valore < hi:
            return azione, pct
    return fasce[-1][2], fasce[-1][3]


def arrotonda_dose(mg):
    """Arrotonda una dose giornaliera a un passo clinicamente leggibile."""
    if mg is None or mg <= 0:
        return mg
    if mg >= 4000:
        passo = 500
    elif mg >= 1000:
        passo = 250
    elif mg >= 200:
        passo = 25
    elif mg >= 50:
        passo = 5
    else:
        passo = 1 if mg >= 10 else 0.1
    return round(mg / passo) * passo


def auc_da_due_livelli(c1, t1, c2, t2, tau, tinf):
    """
    AUC su un intervallo (Sawchuk-Zaske, 1 cmt): C1 e C2 prelevati in fase di
    eliminazione ai tempi t1 < t2 (h dall'inizio dell'infusione).
    Restituisce (AUC24, Cmax a fine infusione, Cmin a fine intervallo, ke).
    """
    if not (c1 > 0 and c2 > 0 and t2 > t1 and tau > 0):
        return None
    ke = math.log(c1 / c2) / (t2 - t1)
    if ke <= 0:
        return None
    cmax = c1 * math.exp(ke * (t1 - tinf))
    cmin = c2 * math.exp(-ke * (tau - t2))
    auc_tau = (cmax + cmin) / 2 * tinf + (cmax - cmin) / ke
    return auc_tau * 24.0 / tau, cmax, cmin, ke


def descrivi_rene(egfr):
    if egfr is None:
        return None
    if egfr < 30:
        return "funzionalità renale gravemente ridotta"
    if egfr < 60:
        return "funzionalità renale ridotta"
    if egfr <= 130:
        return "funzionalità renale nella norma"
    return "clearance renale aumentata (ARC)"
