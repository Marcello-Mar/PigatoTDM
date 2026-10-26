"""Punto d'ingresso: stile, accesso con parola d'ordine, navigazione."""

import streamlit as st

from . import auth, bayes_page, referto_page, ui

VERSIONE = "7.0"


def main():
    st.set_page_config(page_title="PIGATO · TDM", layout="wide", page_icon=":material/science:")
    ui.applica_stile()
    auth.richiedi_accesso()
    pagine = st.navigation([
        st.Page(bayes_page.render, title="TDM bayesiana", url_path="bayes", default=True),
        st.Page(referto_page.render, title="Referto di consulenza", url_path="referto"),
    ], position="top")
    pagine.run()
    st.sidebar.caption(f"PIGATO TDM bayesiana v{VERSIONE}")
