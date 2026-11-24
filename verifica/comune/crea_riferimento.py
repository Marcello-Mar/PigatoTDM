"""Raccoglie i casi verificati e i valori mapbayr in tests/dati/riferimento_mapbayr.json."""
import sys, json, csv
sys.path.insert(0, sys.argv[1])
from pigato_bayes import models
CIP = {m.split(' ')[0]: m for m in models.CIPRO_MODELLI}
DAL = {m.split(' ')[0]: m for m in models.DALBA_MODELLI}
casi = []
for farm, js, cs in [('Daptomicina', 'dapto/casi_dapto.json', 'dapto/mapbayr_dapto.csv'),
                     ('Ciprofloxacina', 'cipro/casi_cipro.json', 'cipro/mapbayr_cipro.csv'),
                     ('Dalbavancina', 'dalba_amino/casi_dalba.json', 'dalba_amino/mapbayr_dalba.csv'),
                     ('Aminoglicoside', 'dalba_amino/casi_amino.json', 'dalba_amino/mapbayr_amino.csv')]:
    R = {r['id']: r for r in csv.DictReader(open(cs))}
    for c in json.load(open(js)):
        r = R[c['id']]; cv = c['cov']; nome = c['modello']
        if farm == 'Daptomicina':
            chiamata = dict(funzione='prior_dapto', args=dict(sesso=cv['sesso'], gg=cv['gg'], peso=cv['peso'], crea=cv['crea'],
                            crrt=cv['crrt'], Sc=cv['Sc'], Qef=cv['Qef'], dialisi=cv['dialisi']))
        elif farm == 'Ciprofloxacina':
            mod = CIP[nome.split(' ')[0].replace(',', '')]
            chiamata = dict(funzione='prior_cipro', args=dict(sesso=cv['sesso'], gg=cv['gg'], peso=cv['peso'], crea=cv['crea'],
                            crrt=cv['crrt'], Sc=cv['Sc'], Qef=cv['Qef'], altezza=cv['altezza'], egfr=cv['egfr'], modello=mod))
        elif farm == 'Dalbavancina':
            mod = DAL['Mariani' if nome.startswith('Modello Mariani') else nome.split(' ')[0]]
            chiamata = dict(funzione='prior_dalba', args=dict(peso=cv['peso'], modello=mod))
        else:
            f = cv['farmaco']; aut = nome.split(' ')[0]
            mod = next(m for m in models.AMINO_MODELLI[f] if m.startswith(aut))
            chiamata = dict(funzione='prior_amino', args=dict(farmaco=f, modello=mod, sesso=cv['sesso'], gg=cv['gg'],
                            peso=cv['peso'], altezza=cv['altezza'], crea=cv['crea'], crrt=cv['crrt'], Sc=cv['Sc'],
                            Qef=cv['Qef'], eg_sett=cv['eg_sett']))
        PK = ('CL', 'V1', 'Q', 'V2')
        att = {k: float(r[k]) for k in PK if k in r}
        if 'V' in r: att['V1'] = float(r['V'])
        casi.append(dict(id=c['id'], descr=c['descr'], chiamata=chiamata, prior_atteso=c['prior'],
                         omega=c['omega'], storia=c['storia'], livelli=c['livelli'],
                         limiti_attivi_pigato=c['limiti_attivi'], finestra=c['finestra'],
                         atteso_mapbayr=dict(par=att, AUC24=float(r['AUC24']), Cmin=float(r['Cmin']),
                                             t_sotto_soglia=float(r['t_sotto_soglia']) if r.get('t_sotto_soglia') not in (None, '', 'NA') else None),
                         soglia=c.get('soglia'), t_val=c.get('t_val')))
meta = dict(descrizione="Valori di riferimento calcolati con R 4.3.3, mrgsolve 2.0.1 e mapbayr 0.10.2 (CRAN) il 24/11/2026, omega pubblicati senza riduzione; modelli scritti come ODE, AUC per integrazione. Script in verifica/ del rapporto.",
            casi=casi)
json.dump(meta, open(sys.argv[2], 'w'), indent=1, ensure_ascii=False)
print(len(casi), 'casi')
