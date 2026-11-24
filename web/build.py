"""
Costruisce la versione nel browser (stlite) in web/dist, per Cloudflare Pages.

    python3 web/build.py

Cloudflare Pages: build command `python3 web/build.py`, output directory `web/dist`.
Il codice Python gira nel browser dell'utente (Pyodide): nessun dato passa da un server.
"""

import json
import shutil
import tomllib
from pathlib import Path

RADICE = Path(__file__).resolve().parents[1]
DIST = RADICE / "web" / "dist"
STLITE = "1.12.3"          # Streamlit 1.65, Pyodide 314 (Python 3.14)
REQUISITI = ["numpy", "scipy", "pandas", "plotly"]


def main():
    if DIST.exists():
        shutil.rmtree(DIST)
    DIST.mkdir(parents=True)
    file_app = ["pigato.py"] + sorted(
        str(p.relative_to(RADICE)) for p in (RADICE / "pigato_bayes").rglob("*")
        if p.is_file() and p.suffix in (".py", ".json") and "__pycache__" not in p.parts)
    for f in file_app:
        dest = DIST / "app" / f
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(RADICE / f, dest)

    # tema e opzioni da .streamlit/config.toml, in forma "sezione.chiave"
    cfg = tomllib.loads((RADICE / ".streamlit" / "config.toml").read_text())
    config = {f"{sez}.{k}": v for sez, voci in cfg.items() if sez in ("theme", "client") for k, v in voci.items()}
    # stlite 1.12 va in errore (HelmetDispatcher) con i font del tema dati come "Nome:URL":
    # al tema passa solo il nome, i fogli dei font li carica la pagina
    link_font = []
    for k in ("theme.font", "theme.headingFont", "theme.codeFont"):
        if isinstance(config.get(k), str) and ":http" in config[k]:
            nome, url = config[k].split(":", 1)
            config[k] = nome
            link_font.append(f'<link rel="stylesheet" href="{url}">')

    opzioni = {
        "entrypoint": "pigato.py",
        "requirements": REQUISITI,
        "files": {f: {"url": f"app/{f}"} for f in file_app},
        "streamlitConfig": config,
    }
    html = (RADICE / "web" / "index.template.html").read_text()
    html = html.replace("{{FONT}}", "\n  ".join(link_font)).replace("{{STLITE}}", STLITE).replace("{{OPZIONI}}", json.dumps(opzioni, indent=2, ensure_ascii=False))
    (DIST / "index.html").write_text(html)
    print(f"web/dist pronto: {len(file_app)} file Python/JSON, stlite {STLITE}")


if __name__ == "__main__":
    main()
