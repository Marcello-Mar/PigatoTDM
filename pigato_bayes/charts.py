"""Grafici Plotly con sfondo chiaro fisso (leggibili anche con browser in tema scuro)."""

import numpy as np
import plotly.graph_objects as go

from .pk_engine import conc

INK = '#0E2C52'
WINE = '#AE7E22'
FAINT = '#94A3B8'
DANGER = '#B23A2E'
OK = '#1F6E4A'

LAYOUT = dict(
    plot_bgcolor='#FFFFFF', paper_bgcolor='#FFFFFF',
    font=dict(color='#1F2937', size=13, family='IBM Plex Sans, sans-serif'),
    xaxis=dict(color='#1F2937', gridcolor='#E5E7EB', linecolor='#9CA3AF', automargin=True,
               title=dict(standoff=8)),
    yaxis=dict(color='#1F2937', gridcolor='#E5E7EB', linecolor='#9CA3AF', automargin=True,
               title=dict(standoff=10), rangemode='tozero'),
    hovermode="closest",
    # legenda sotto il grafico: in alto si sovrapponeva alla barra degli strumenti
    # e su telefono andava su quattro righe sopra la curva
    legend=dict(orientation="h", yanchor="top", y=-0.22, xanchor="left", x=0, font=dict(color='#1F2937', size=12)),
    margin=dict(l=60, r=16, t=24, b=24),
    height=440,
)


def _soglia(fig, val, etich, col):
    """Linea orizzontale con etichetta dentro l'area del grafico (a sinistra, sopra la linea)."""
    fig.add_hline(y=val, line_dash="dashdot", line_color=col)
    fig.add_annotation(x=0.01, xref="paper", y=val, yref="y", text=etich, showarrow=False,
                       xanchor="left", yanchor="bottom", font=dict(color=col, size=11),
                       bgcolor="rgba(255,255,255,0.8)")


def profilo(nome_farmaco, storia, livelli, par_ind, par_pop, t_inf, t_ultimo_storico=-1,
            soglie=None, t_fine=None):
    t_last = max([s['ore'] for s in storia] + [l['ore'] for l in livelli] + [0])
    t_fine = t_fine or (t_last + 24)
    tg = np.linspace(0, t_fine, 900)
    c_ind = conc(par_ind, storia, tg, t_inf)
    c_pop = conc(par_pop, storia, tg, t_inf)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=tg, y=c_pop, mode='lines', name='Popolazione',
                             line=dict(dash='dash', color=FAINT, width=2)))
    fig.add_trace(go.Scatter(x=tg, y=c_ind, mode='lines', name='Paziente (MAP)',
                             line=dict(color=INK, width=3)))
    ymax = float(max(c_ind.max(), c_pop.max(), max([l['conc'] for l in livelli] + [1])))
    if t_ultimo_storico >= 0:
        fig.add_vrect(x0=0, x1=t_ultimo_storico, fillcolor="#6EA8D8", opacity=0.08, layer="below", line_width=0)
    xs, ys = [], []
    for s in storia:
        xs += [s['ore'], s['ore'], None]
        ys += [0, ymax, None]
    fig.add_trace(go.Scatter(x=xs, y=ys, mode='lines', name='Dosi', hoverinfo='none',
                             line=dict(color='#475569', dash='dot', width=1)))
    if livelli:
        fig.add_trace(go.Scatter(
            x=[l['ore'] for l in livelli], y=[l['conc'] for l in livelli], mode='markers', name='Livelli misurati',
            marker=dict(symbol='diamond', size=12, color=WINE, line=dict(width=1, color='white')),
            text=[f"t={l['ore']:.1f} h<br>C={l['conc']:.3g} mg/L" + (f"<br>{l['nota']}" if l.get('nota') else "")
                  for l in livelli], hoverinfo='text'))
    for (val, etich, col) in (soglie or []):
        _soglia(fig, val, etich, col)
    fig.update_layout(xaxis_title_text="Ore dalla prima dose", yaxis_title_text=f"{nome_farmaco} (mg/L)", **LAYOUT)
    return fig


def forecast(nome_farmaco, tg, cg, soglia, storia, t_nuove=None, etichetta_soglia=None):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=tg, y=cg, mode='lines', name='Proiezione', line=dict(color=INK, width=2.5)))
    _soglia(fig, soglia, etichetta_soglia or f"Soglia {soglia:.2f} mg/L".replace('.', ','), DANGER)
    ymax = float(max(cg.max(), soglia) if len(cg) else soglia)
    xs, ys = [], []
    for s in storia:
        xs += [s['ore'], s['ore'], None]
        ys += [0, ymax, None]
    fig.add_trace(go.Scatter(x=xs, y=ys, mode='lines', name='Dosi', hoverinfo='none',
                             line=dict(color=FAINT, dash='dot', width=1)))
    for t in (t_nuove or []):
        fig.add_vline(x=t, line_dash="dash", line_color=OK)
    fig.update_layout(xaxis_title_text="Ore dalla prima dose", yaxis_title_text=f"{nome_farmaco} (mg/L)", **LAYOUT)
    return fig
