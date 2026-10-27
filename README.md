# PIGATO · TDM bayesiana

App Streamlit per la TDM degli antimicrobici, parte della famiglia PIGATO
(IRCCS Istituto Giannina Gaslini – UOC Malattie Infettive). Uso educational, non dispositivo medico.

## Moduli
- **TDM bayesiana** (MAP): daptomicina, ciprofloxacina, dalbavancina (forecast della prossima dose),
  gentamicina, amikacina, tobramicina (anche da un solo livello, es. a 10 h).
- **Referto di consulenza**: 32 farmaci di `tdm_data.json`. Funzione renale con la formula per età,
  penetrazione al sito (dataset PIGATO), MIC misurata o breakpoint EUCAST v16.1 in empirica,
  proposta riduci/mantieni/aumenta modificabile, testo del referto.

Nessun dato resta sul server: lo storico di ogni paziente × molecola × modulo è un file JSON che
l'operatore scarica e ricarica (compatibile con i file della versione 6).

## Avvio
```bash
pip install -r requirements.txt
streamlit run pigato.py
```

## Accesso
La parola d'ordine si imposta in Streamlit Community Cloud → App settings → Secrets:

```toml
APP_PASSWORD = "..."
```

In locale: copia `.streamlit/secrets.toml.example` in `.streamlit/secrets.toml` (escluso da git)
oppure usa la variabile d'ambiente `PIGATO_PASSWORD`.

## Struttura
- `pigato.py`: punto d'ingresso
- `pigato_bayes/models.py`: modelli di popolazione, con la nota di verifica sulle fonti (`VERIFICA`)
- `pigato_bayes/pk_engine.py`: 1/2 compartimenti, MAP, metriche
- `pigato_bayes/renal.py`: eGFR (stesse formule di PIGATO TDM), FFM, creatinina attesa
- `pigato_bayes/rules.py`: regole di interpretazione del referto
- `pigato_bayes/ui.py`, `charts.py`: stile PIGATO, schede, grafici
- `pigato_bayes/data/`: target TDM, penetrazione, breakpoint EUCAST
- `.streamlit/config.toml`: tema e font IBM Plex
