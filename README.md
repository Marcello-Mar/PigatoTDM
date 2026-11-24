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

## Versione nel browser (Cloudflare Pages)
Il codice Python gira nel browser con [stlite](https://github.com/whitphx/stlite) (Pyodide): nessun dato
passa da un server. Pensata per PC; il primo caricamento scarica circa 30–50 MB, poi resta in cache.

Impostazione su Cloudflare (Workers & Pages → Create → Pages → Connect to Git → questo repository):
- Build command: `python3 web/build.py`
- Build output directory: `web/dist`
- Settings → Variables and Secrets: `PIGATO_PASSWORD` (tipo Secret); facoltativo `PIGATO_UTENTE`

La parola d'ordine la chiede il browser prima di caricare qualsiasi file (`functions/_middleware.js`,
HTTP Basic). In locale: `python3 web/build.py && npx wrangler pages dev web/dist --binding PIGATO_PASSWORD=...`.

## Test
```bash
pip install pytest
python -m pytest tests -q
```
- `tests/test_calcoli.py`: eGFR, FFM, motore PK, MAP, prior di popolazione e regole contro valori calcolati a mano.
- `tests/test_referto.py`: decisioni per tipo di regola, referto completo sui 32 farmaci (bambino e adulto con IR), JSON andata e ritorno.
- `tests/test_bayes.py`: percorso completo dei 6 farmaci bayesiani, casi limite (neonato pretermine, CRRT, HD, IR).

