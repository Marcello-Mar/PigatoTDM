"""Stile PIGATO (design system v7) e piccoli helper di formattazione."""

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
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Serif:wght@500;600&display=swap">
<style>
:root{--pg-ink:#0E2C52;--pg-ink-soft:#41607F;--pg-ink-faint:#6C86A5;--pg-line:#D3DEED;--pg-wine:#AE7E22;--pg-tint:#E7EFF9;}
@media (prefers-color-scheme:dark){:root{--pg-ink:#E8F0FB;--pg-ink-soft:#AEC4E0;--pg-ink-faint:#7C96B7;--pg-line:#20375A;--pg-wine:#D6B15A;--pg-tint:#122D51;}}
html, body, [class*="st-"], .stMarkdown, .stTextInput, .stNumberInput, .stSelectbox, button {font-family:"IBM Plex Sans",system-ui,sans-serif;}
h1, h2, h3 {font-family:"IBM Plex Serif",Georgia,serif !important; font-weight:600 !important; letter-spacing:.01em;}
code, pre, .pg-mono {font-family:"IBM Plex Mono",monospace !important;}
.pg-head{display:flex;align-items:center;gap:14px;margin:0 0 4px 0}
.pg-word{font-family:"IBM Plex Serif",Georgia,serif;font-weight:600;font-size:30px;letter-spacing:.06em;color:var(--pg-ink);line-height:1}
.pg-tag{font-family:"IBM Plex Mono",monospace;font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--pg-wine);margin-top:4px}
.pg-inst{font-family:"IBM Plex Mono",monospace;font-size:11px;color:var(--pg-ink-faint);line-height:1.45;margin-left:auto;text-align:right}
.pg-step{font-family:"IBM Plex Mono",monospace;font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--pg-ink-faint);margin:18px 0 2px}
.pg-chip{display:inline-block;padding:2px 9px;border-radius:999px;font-size:12px;font-weight:600;margin-right:6px}
.pg-g{background:#E1F0E6;color:#1F6E4A}.pg-y{background:#FBF0D6;color:#8A5A0E}.pg-o{background:#FBE5D3;color:#9A4A12}
.pg-r{background:#F7DCD8;color:#B23A2E}.pg-b{background:#E3ECF8;color:#1C4F8F}
</style>
"""


def applica_stile():
    st.markdown(CSS, unsafe_allow_html=True)


def intestazione(sottotitolo):
    st.markdown(
        f"""<div class="pg-head">{EMBLEMA_SVG}
        <div><div class="pg-word">PIGATO</div><div class="pg-tag">{sottotitolo}</div></div>
        <div class="pg-inst">IRCCS Istituto Giannina Gaslini – Genova<br>UOC Malattie Infettive</div></div>""",
        unsafe_allow_html=True)


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
    if dec == 0:
        s = f"{x:,.0f}"
    else:
        s = f"{x:,.{dec}f}".rstrip('0').rstrip('.')
    return s.replace(',', 'X').replace('.', ',').replace('X', '.')


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
