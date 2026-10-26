"""
Accesso con parola d'ordine unica.

La parola d'ordine non sta nel codice: si legge da st.secrets["APP_PASSWORD"]
(Streamlit Community Cloud → App settings → Secrets) oppure, in locale, da
.streamlit/secrets.toml o dalla variabile d'ambiente PIGATO_PASSWORD.
Se non è configurata l'app resta chiusa.

È una barriera contro l'accesso casuale, non un'autenticazione per utente.
"""

import hmac
import os
import time

import streamlit as st


def _password_configurata():
    try:
        v = st.secrets.get("APP_PASSWORD")
    except Exception:
        v = None
    return v or os.environ.get("PIGATO_PASSWORD")


def richiedi_accesso():
    if st.session_state.get("_accesso_ok"):
        return True
    attesa = _password_configurata()
    if not attesa:
        st.error("Parola d'ordine non configurata. In Streamlit Community Cloud: App settings → "
                 "Secrets → aggiungi APP_PASSWORD = \"...\".")
        st.stop()
    with st.form("accesso", clear_on_submit=True):
        st.markdown("#### Accesso riservato")
        pw = st.text_input("Parola d'ordine", type="password")
        entra = st.form_submit_button("Entra", type="primary")
    if entra:
        tentativi = st.session_state.get("_tentativi", 0)
        if tentativi >= 5:
            time.sleep(min(2 * (tentativi - 4), 10))
        if hmac.compare_digest(pw.encode("utf-8"), str(attesa).encode("utf-8")):
            st.session_state["_accesso_ok"] = True
            st.session_state.pop("_tentativi", None)
            st.rerun()
        st.session_state["_tentativi"] = tentativi + 1
        st.error("Parola d'ordine errata.")
    st.stop()
