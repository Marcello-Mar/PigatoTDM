"""Stile PIGATO (design system v7) e piccoli helper di formattazione."""

import html

import streamlit as st

EMBLEMA_SVG = """<svg viewBox="0 0 200 200" width="44" height="44" role="img" aria-label="PIGATO" fill="none" stroke="currentColor" stroke-width="4" stroke-linecap="round" stroke-linejoin="round" style="color:var(--pg-ink);flex:none">
<circle cx="100" cy="100" r="92"/><circle cx="100" cy="43" r="3.6" fill="currentColor" stroke="none"/>
<line x1="100" y1="46" x2="100" y2="92"/>
<path d="M99 53 C86 47 77 51 72 59 C81 55 90 56 98 60" stroke-width="3.2"/>
<path d="M101 53 C114 47 123 51 128 59 C119 55 110 56 102 60" stroke-width="3.2"/>
<g fill="var(--pg-wine)" stroke="none"><circle cx="100" cy="66" r="4.6"/><circle cx="92.5" cy="72.5" r="4.6"/><circle cx="107.5" cy="72.5" r="4.6"/><circle cx="96" cy="80" r="4.6"/><circle cx="104" cy="80" r="4.6"/><circle cx="100" cy="87.5" r="4.6"/></g>
<line x1="50" y1="151" x2="150" y2="151" stroke-width="3.4"/>
<path d="M65 108 C65 125 87 125 87 108Z" stroke-width="3.2"/><path d="M69 112 C70 120 82 120 83 112Z" fill="var(--pg-wine)" stroke="none"/>
<line x1="76" y1="123" x2="76" y2="147" stroke-width="3.2"/><path d="M68 149 Q76 144 84 149" stroke-width="3.2"/>
<rect x="101" y="90" width="7" height="18" rx="1.5" stroke-width="3.2"/>
<path d="M101 108 Q96 114 96 122 L96 150 L113 150 L113 122 Q113 114 108 108" stroke-width="3.2"/>
<rect x="99.5" y="127" width="10" height="12" rx="1.5" stroke-width="2.6"/>
<path d="M120 150 L120 143 Q120 140 123 140 L129 140 L133 134 L142 134 L145 140 L147 140 Q149 140 149 143 L149 150" stroke-width="3.2"/>
<circle cx="127" cy="150" r="4.6" stroke-width="3"/><circle cx="142" cy="150" r="4.6" stroke-width="3"/></svg>"""

CSS = """
<style>
/* Il font si imposta da .streamlit/config.toml (theme.font): qui niente regole
   globali su span/div, che sostituivano anche il font delle icone Material. */
:root{--pg-ink:#0E2C52;--pg-ink-soft:#41607F;--pg-ink-faint:#516A88;--pg-line:#D3DEED;
      --pg-wine:#8A6318;--pg-tint:#E7EFF9;--pg-surface:#FFFFFF;}
.pg-head{display:flex;flex-wrap:wrap;align-items:center;gap:10px 14px;margin:0 0 6px 0}
.pg-brand{display:flex;align-items:center;gap:12px;min-width:0}
.pg-word{font-family:"IBM Plex Serif",Georgia,serif;font-weight:600;font-size:30px;letter-spacing:.06em;
         color:var(--pg-ink);line-height:1}
.pg-tag{font-family:"IBM Plex Mono",monospace;font-size:11px;letter-spacing:.14em;text-transform:uppercase;
        color:var(--pg-wine);margin-top:4px}
.pg-inst{font-family:"IBM Plex Mono",monospace;font-size:12px;color:var(--pg-ink-faint);line-height:1.45;
         margin-left:auto;text-align:right}
@media (max-width:640px){.pg-inst{margin-left:0;text-align:left;flex-basis:100%}.pg-word{font-size:26px}}
.pg-step{font-family:"IBM Plex Mono",monospace;font-size:11px;letter-spacing:.14em;text-transform:uppercase;
         color:var(--pg-ink-faint);margin:18px 0 2px}
.pg-chip{display:inline-block;padding:3px 10px;border-radius:999px;font-size:13px;font-weight:600;white-space:nowrap}
.pg-g{background:#E1F0E6;color:#1F6E4A}.pg-y{background:#FBF0D6;color:#8A5A0E}.pg-o{background:#FBE5D3;color:#9A4A12}
.pg-r{background:#F7DCD8;color:#B23A2E}.pg-b{background:#E3ECF8;color:#1C4F8F}
/* schede numeriche: vanno a capo invece di troncare */
.pg-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:6px 0 12px}
.pg-card{border:1px solid var(--pg-line);border-radius:12px;padding:10px 12px;background:var(--pg-surface);min-width:0}
.pg-card.pg-sel{border-color:var(--pg-ink);box-shadow:inset 0 0 0 1px var(--pg-ink)}
.pg-k{font-size:12.5px;color:var(--pg-ink-soft);line-height:1.3;overflow-wrap:anywhere}
.pg-v{font-size:22px;font-weight:600;color:var(--pg-ink);line-height:1.2;margin-top:3px;
      font-variant-numeric:tabular-nums;overflow-wrap:anywhere}
.pg-u{font-size:13px;font-weight:500;color:var(--pg-ink-soft);margin-left:3px}
.pg-s{font-size:12px;color:var(--pg-ink-faint);margin-top:3px;line-height:1.3;overflow-wrap:anywhere}
.pg-s.ok{color:#1F6E4A}.pg-s.ko{color:#B23A2E}
/* scheda della proposta: il colore segue l'azione (verde mantieni, ambra cambia, rosso sospendi/tossicità) */
.pg-prop{border:1px solid var(--pg-p-line);border-radius:14px;background:var(--pg-p-bg);padding:14px 16px 12px;
         margin:8px 0 16px;box-shadow:0 1px 2px rgba(14,44,82,.06)}
.pg-prop-ok{--pg-p-bg:#EEF6F1;--pg-p-line:#B7D7C3;--pg-p-ink:#1A5E3F}
.pg-prop-cambia{--pg-p-bg:#FBF4E3;--pg-p-line:#E3CC93;--pg-p-ink:#74500B}
.pg-prop-stop{--pg-p-bg:#FBEDEB;--pg-p-line:#E5B3AC;--pg-p-ink:#952C22}
.pg-prop-riga{display:flex;flex-wrap:wrap;align-items:baseline;gap:4px 12px}
.pg-prop-az{font-size:21px;font-weight:700;color:var(--pg-p-ink);letter-spacing:.01em;white-space:nowrap}
.pg-prop-reg{font-size:18px;font-weight:600;color:var(--pg-ink);font-variant-numeric:tabular-nums;
             overflow-wrap:anywhere;min-width:0}
.pg-prop ul{margin:8px 0 0;padding-left:18px}
.pg-prop li{font-size:14.5px;color:var(--pg-ink);line-height:1.45;margin:1px 0;overflow-wrap:anywhere}
.pg-prop li.pg-prop-allerta{color:var(--pg-p-ink);font-weight:600}
@media (max-width:640px){.pg-prop-az{font-size:19px}.pg-prop-reg{font-size:16px}}
/* didascalie: il grigio di Streamlit (testo al 60%) scende sotto 4,5:1 sullo sfondo */
[data-testid="stCaptionContainer"]{color:var(--pg-ink-soft) !important;opacity:1 !important}
/* blocchi di codice del referto: mai oltre il riquadro */
[data-testid="stCode"] pre, [data-testid="stCode"] code{white-space:pre-wrap !important;overflow-wrap:anywhere}
</style>
"""


def applica_stile():
    st.markdown(CSS, unsafe_allow_html=True)


def intestazione(sottotitolo):
    st.markdown(
        f"""<div class="pg-head"><div class="pg-brand">{EMBLEMA_SVG}
        <div><div class="pg-word">PIGATO</div><div class="pg-tag">{sottotitolo}</div></div></div>
        <div class="pg-inst">IRCCS Istituto Giannina Gaslini – Genova<br>UOC Malattie Infettive</div></div>""",
        unsafe_allow_html=True)


def schede(voci):
    """
    Griglia di valori che va a capo sullo schermo stretto (al posto di st.metric,
    che tronca con "…"). voci: lista di dict {k, v, u?, s?, stato?: 'ok'|'ko', sel?: bool}.
    """
    from html import escape
    parti = []
    for x in voci:
        if x is None:
            continue
        sub = ""
        if x.get('s'):
            cls = {'ok': ' ok', 'ko': ' ko'}.get(x.get('stato'), '')
            sub = f'<div class="pg-s{cls}">{escape(str(x["s"]))}</div>'
        u = f'<span class="pg-u">{escape(x["u"])}</span>' if x.get('u') else ''
        parti.append(f'<div class="pg-card{" pg-sel" if x.get("sel") else ""}"><div class="pg-k">{escape(x["k"])}</div>'
                     f'<div class="pg-v">{escape(str(x["v"]))}{u}</div>{sub}</div>')
    st.markdown(f'<div class="pg-grid">{"".join(parti)}</div>', unsafe_allow_html=True)


def schede_egfr(eg):
    from .renal import ETICHETTE
    voci = []
    for kf, v in eg['tutte'].items():
        voci.append({'k': ETICHETTE[kf], 'v': n(v, 0), 'u': "mL/min" if kf == 'cg' else "mL/min/1,73 m²",
                     'sel': kf == eg['chiave'], 's': "formula raccomandata" if kf == eg['chiave'] else None})
    if eg.get('assoluta'):
        voci.append({'k': "Clearance assoluta", 'v': n(eg['assoluta'], 0), 'u': "mL/min",
                     's': "eGFR × BSA/1,73 (Mosteller)"})
    return voci


def passo(testo):
    st.markdown(f'<div class="pg-step">{testo}</div>', unsafe_allow_html=True)


def chip(testo, colore):
    return f'<span class="pg-chip pg-{colore}">{testo}</span>'


def n(x, dec=1):
    """Numero in formato italiano (virgola decimale), senza zeri inutili."""
    if x is None:
        return "N/D"
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    # niente separatore delle migliaia: "1.911" si leggerebbe come 1,9
    if dec == 0:
        s = f"{x:.0f}"
    else:
        s = f"{x:.{dec}f}".rstrip('0').rstrip('.')
    return s.replace('.', ',')


def n_semplice(x, dec=1):
    """Come n() ma senza separatore delle migliaia (per dosi in mg nel referto)."""
    if x is None:
        return "N/D"
    x = float(x)
    s = f"{x:.{dec}f}".rstrip('0').rstrip('.') if dec else f"{x:.0f}"
    return s.replace('.', ',')


AVVISO = ("Strumento a scopo educational, non validato né certificato come dispositivo medico: "
          "ogni proposta va verificata dal medico, che resta l'unico responsabile della terapia. "
          "Nessun dato viene conservato sul server; lo storico esiste solo nel file JSON che scarichi. "
          "Non inserire nomi, cognomi o date di nascita: usa un codice paziente non identificativo.")


def conferma(chiave, titolo, messaggio, esegui, cliccato):
    """
    Finestra di conferma per le azioni che cancellano i dati della sessione.
    Lo stato "aperta" vive in session_state, così la finestra resta visibile
    fra un rerun e l'altro e si chiude solo con Annulla o Conferma.
    """
    flag = f'_conferma_{chiave}'
    if cliccato:
        st.session_state[flag] = True
    if not st.session_state.get(flag):
        return

    def _finestra():
        st.markdown(messaggio)
        a, b = st.columns(2)
        if a.button("Annulla", width="stretch", key=f"{flag}_no"):
            st.session_state.pop(flag, None)
            st.rerun()
        if b.button("Conferma", type="primary", width="stretch", key=f"{flag}_si"):
            st.session_state.pop(flag, None)
            esegui()
            st.rerun()

    try:
        st.dialog(titolo, dismissible=False)(_finestra)()
    except TypeError:            # Streamlit senza il parametro dismissible
        st.dialog(titolo)(_finestra)()


def messaggio_perdita(sto, cosa):
    n = len(sto.get('monitoraggi', []))
    perdita = (f"Le **{n} valutazioni registrate** in questa sessione andranno perse se non hai già scaricato "
               "lo storico JSON." if n else "I dati inseriti e non registrati andranno persi.")
    return f"{cosa}\n\n{perdita}"


def tono_azione(azione, allerta=False):
    if azione == 'SOSPENDI' or (azione or '').startswith('NESSUN'):
        return 'stop'
    if allerta:
        return 'stop'
    return 'ok' if azione == 'MANTIENI' else 'cambia'


def proposta(azione, regime="", motivi=(), allerte=(), tono=None, dove=None):
    """
    Scheda della decisione: azione e regime in evidenza, motivi sotto.
    `allerte` sono righe che richiedono attenzione (tossicità, banda < 10%, target irraggiungibile).
    `dove` è il contenitore in cui scrivere (per metterla in cima ai risultati).
    """
    tono = tono or tono_azione(azione, bool(allerte))
    esc = html.escape
    righe = "".join(f'<li class="pg-prop-allerta">{esc(a)}</li>' for a in allerte if a)
    righe += "".join(f"<li>{esc(m)}</li>" for m in motivi if m)
    blocco = (f'<div class="pg-prop pg-prop-{tono}" role="status">'
              f'<div class="pg-prop-riga"><span class="pg-prop-az">{esc(azione)}</span>'
              + (f'<span class="pg-prop-reg">{esc(regime)}</span>' if regime else "")
              + "</div>" + (f"<ul>{righe}</ul>" if righe else "") + "</div>")
    (dove or st).markdown(blocco, unsafe_allow_html=True)

