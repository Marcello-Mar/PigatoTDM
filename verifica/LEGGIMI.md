# Verifica PIGATO TDM contro mrgsolve/mapbayr (24/11/2026)
PIGATO: Marcello-Mar/PigatoTDM rev. bcad55f; numpy 2.5.3, scipy 1.18.1
Riferimento: R 4.3.3; mrgsolve 2.1.0.9000 (GitHub c529f2d); mapbayr 0.10.2.9002 (GitHub 0412b46)
  (CRAN non raggiungibile dal container: versioni di sviluppo da GitHub; da ripetere con le versioni CRAN)

comune/prior_indipendenti.R  prior, omega, sigma ricalcolati dalle formule delle fonti (senza codice PIGATO)
comune/mapbayr_2cmt.R        MAP a 2 compartimenti (ODE), errore log-additivo o combinato
comune/confronto.py          tabella degli scarti PIGATO vs mapbayr

Daptomicina (cartella dapto/):
  python casi_dapto.py <repo> .   ;  Rscript mapbayr_dapto.R casi_dapto.json mapbayr_dapto.csv
  python confronto_dapto.py       ;  Rscript ../comune/prior_indipendenti.R casi_dapto.json prior_dapto.csv
  python controlli_dapto.py <repo>
Ciprofloxacina (cartella cipro/):
  python casi_cipro.py <repo> .   ;  Rscript ../comune/mapbayr_2cmt.R casi_cipro.json mapbayr_cipro.csv
  python ../comune/confronto.py casi_cipro.json mapbayr_cipro.csv confronto_cipro.csv
  Rscript ../comune/prior_indipendenti.R casi_cipro.json prior_cipro.csv
  python controlli_cipro.py <repo> casi_cipro.json
Dalbavancina e aminoglicosidi (cartella dalba_amino/):
  python casi_dalba_amino.py <repo> .
  Rscript ../comune/mapbayr_2cmt.R casi_dalba.json mapbayr_dalba.csv   (idem per casi_amino.json)
  python ../comune/confronto.py casi_dalba.json mapbayr_dalba.csv confronto_dalba.csv
  Rscript ../comune/prior_indipendenti.R casi_dalba.json prior_dalba.csv

Riferimento per i test di regressione:
  python comune/crea_riferimento.py <repo> <repo>/tests/dati/riferimento_mapbayr.json
  (dalla cartella verifica/, dopo aver rigenerato i casi e i CSV di mapbayr)
  I test sono in tests/test_verifica_mapbayr.py.

Requisiti R: install.packages(c("mrgsolve", "mapbayr", "jsonlite")); serve un compilatore C++ (Rtools su Windows).
