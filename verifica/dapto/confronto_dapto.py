"""Confronto PIGATO vs mapbayr, daptomicina. Uso: python confronto_dapto.py (nella cartella dei risultati)"""
import json, csv
casi = {c['id']: c for c in json.load(open('casi_dapto.json'))}
R = {r['id']: r for r in csv.DictReader(open('mapbayr_dapto.csv'))}
rel = lambda a, b: 100 * (a - b) / b
righe, mx = [], {}
for i, c in casi.items():
    r = R[i]; lim = bool(c['limiti_attivi'])
    p = c['map_senza_limiti'] if lim else c['map']
    d = dict(CL=rel(p['CL'], float(r['CL'])), V=rel(p['V1'], float(r['V'])),
             AUC24=rel(c['auc24_senza_limiti' if lim else 'auc24'], float(r['AUC24'])), Cmin=rel(c['cmin_senza_limiti' if lim else 'cmin'], float(r['Cmin'])))
    for k, v in d.items():
        mx[k] = max(mx.get(k, 0), abs(v))
    righe.append([i, c['modello'].split(' [')[0].split(',')[0] + (' [' + c['modello'].split('[')[1] if '[' in c['modello'] else ''),
                  f"{p['CL']:.4f}", r['CL'][:7], f"{d['CL']:+.4f}", f"{p['V1']:.3f}", r['V'][:6], f"{d['V']:+.4f}",
                  f"{c['auc24_senza_limiti' if lim else 'auc24']:.1f}", f"{float(r['AUC24']):.1f}", f"{d['AUC24']:+.3f}",
                  f"{c['cmin_senza_limiti' if lim else 'cmin']:.2f}", f"{float(r['Cmin']):.2f}", f"{d['Cmin']:+.4f}",
                  'senza limiti (PIGATO a ' + ','.join(c['limiti_attivi']) + ' al limite)' if lim else '',
                  f"{float(r['scarto_max_pop']):.1e}"])
hdr = ['id', 'modello', 'CL_PIG', 'CL_mapbayr', 'dCL%', 'V_PIG', 'V_mapbayr', 'dV%', 'AUC24_PIG', 'AUC24_mapbayr',
       'dAUC%', 'Cmin_PIG', 'Cmin_mapbayr', 'dCmin%', 'nota', 'scarto_max_curva_pop']
with open('confronto_dapto.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(righe)
for r in righe: print(r[0], r[4], r[7], r[10], r[13], r[14])
print('scarto massimo assoluto %:', {k: round(v, 4) for k, v in mx.items()})
d23 = casi['D23']
print('D23 PIGATO con limiti:', d23['map'], ' senza:', d23['map_senza_limiti'])
