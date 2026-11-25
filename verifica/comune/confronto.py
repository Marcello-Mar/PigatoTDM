"""Confronto PIGATO vs mapbayr (1 o 2 cmt). Uso: python confronto.py casi.json mapbayr.csv uscita.csv"""
import sys, json, csv
casi = json.load(open(sys.argv[1])); R = {r['id']: r for r in csv.DictReader(open(sys.argv[2]))}
rel = lambda a, b: 0.0 if a == b == 0 else 100 * (a - b) / b
righe, mx = [], {}
for c in casi:
    r = R[c['id']]; lim = bool(c['limiti_attivi'])
    p = c['map_senza_limiti'] if lim else c['map']
    st = c.get('stimati') or ['CL', 'V1']
    d = {k: rel(p[k], float(r[k if k in r else 'V'])) for k in st}
    d['AUC24'] = rel(c.get('auc24_senza_limiti') if lim and 'auc24_senza_limiti' in c else c['auc24'], float(r['AUC24']))
    d['Cmin'] = rel(c.get('cmin_senza_limiti') if lim and 'cmin_senza_limiti' in c else c['cmin'], float(r['Cmin']))
    tx = ''
    if c.get('soglia') is not None:
        a, b = c['t_sotto_soglia'], r.get('t_sotto_soglia')
        tx = f"{a - float(b):+.3f} h (PIGATO {a:.1f} h; ultimo dato {c['t_val']:.0f} h)" if b not in (None, '', 'NA') else f"PIGATO {a}, R NA"
        if b not in (None, '', 'NA'): mx['t_soglia_h'] = max(mx.get('t_soglia_h', 0), abs(a - float(b)))
    for k, v in d.items(): mx[k] = max(mx.get(k, 0), abs(v))
    righe.append([c['id'], c['modello'].split(',')[0], ';'.join(st)] + [f"{d.get(k, float('nan')):+.4f}" for k in ('CL', 'V1', 'Q', 'V2', 'AUC24', 'Cmin')]
                 + [','.join(c['limiti_attivi']), r['scarto_max_pop'], tx])
with open(sys.argv[3], 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['id', 'modello', 'stimati', 'dCL%', 'dV1%', 'dQ%', 'dV2%', 'dAUC24%', 'dCmin%', 'limiti_PIGATO_attivi', 'scarto_curva_pop', 'diff_tempo_sotto_soglia'])
    w.writerows(righe)
for x in righe: print(*x)
print('scarto massimo assoluto %:', {k: round(v, 4) for k, v in mx.items()})
