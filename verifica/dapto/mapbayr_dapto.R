# Verifica PIGATO TDM, daptomicina: stima MAP indipendente con mrgsolve + mapbayr.
# Modello a 1 compartimento scritto come ODE (PIGATO usa la soluzione chiusa),
# infusione a velocità costante, ETA esponenziali su CL e V, errore residuo
# additivo su scala log (log_transformation di mapbayr = residuo log-normale di PIGATO).
# OMEGA e SIGMA come varianze (PIGATO usa SD: si elevano al quadrato).
# AUC calcolata per integrazione (compartimento AUC), non con i trapezi.
#
# Uso: Rscript mapbayr_dapto.R casi_dapto.json mapbayr_dapto.csv

.libPaths(c("~/Rlib", .libPaths()))  # libreria locale facoltativa
suppressPackageStartupMessages({
  library(jsonlite)
  library(mrgsolve)
  library(mapbayr)
})
args <- commandArgs(trailingOnly = TRUE)
casi <- fromJSON(args[1], simplifyVector = FALSE)

code <- '
$PROB Daptomicina 1 cmt, infusione, errore log-additivo
$PARAM @annotated @covariates
TVCL : 1  : CL tipica del prior (L/h)
TVV  : 10 : V tipico del prior (L)
$OMEGA 0.1 0.1
$SIGMA 0 0.01
$CMT @annotated
CENT : Compartimento centrale (mg) [ADM, OBS]
AUC  : AUC (mg*h/L) []
$MAIN
double CL = TVCL * exp(ETA(1));
double V  = TVV  * exp(ETA(2));
$ODE
dxdt_CENT = -CL / V * CENT;
dxdt_AUC  = CENT / V;
$TABLE
double CP = CENT / V;
double DV = CP * exp(EPS(2));
$CAPTURE DV CP
$SET rtol = 1e-10, atol = 1e-12, maxsteps = 1e6
'
mod0 <- mcode("dapto", code, quiet = TRUE)

dati <- function(cs, tvcl, tvv) {
  d <- do.call(rbind, lapply(cs$storia, function(s)
    data.frame(ID = 1, time = s$ore, evid = 1, cmt = 1, amt = s$dose,
               rate = s$dose / s$tinf, mdv = 1, DV = NA_real_)))
  o <- do.call(rbind, lapply(cs$livelli, function(l)
    data.frame(ID = 1, time = l$ore, evid = 0, cmt = 1, amt = 0, rate = 0, mdv = 0, DV = l$conc)))
  d <- rbind(d, o)
  d <- d[order(d$time, -d$evid), ]
  d$TVCL <- tvcl; d$TVV <- tvv
  d
}

# AUC e Cmin nella finestra dell'app (24 h dopo l'ultima dose, regime proseguito;
# per intervalli > 24 h: intero intervallo, AUC normalizzata a 24 h)
metriche <- function(mod, cl, v, fin) {
  ev <- do.call(rbind, lapply(fin$dosi, function(s)
    data.frame(ID = 1, time = s$ore, evid = 1, cmt = 1, amt = s$dose, rate = s$dose / s$tinf)))
  m <- param(mod, TVCL = cl, TVV = v) %>% zero_re()
  t0 <- fin$start; t1 <- fin$start + fin$len
  out <- mrgsim_df(m, data = ev, tgrid = tgrid(t0, t1, fin$len / 4000), obsonly = TRUE)
  out <- out[out$time >= t0 - 1e-9 & out$time <= t1 + 1e-9, ]
  c(auc24 = (out$AUC[nrow(out)] - out$AUC[1]) * fin$fattore, cmin = min(out$CP), cmax = max(out$CP))
}

ris <- list()
for (cs in casi) {
  om <- unlist(cs$omega)[c("CL", "V1")]
  mod <- mod0 %>% omat(dmat(om[["CL"]]^2, om[["V1"]]^2)) %>% smat(dmat(0, cs$sigma_prop^2))
  tvcl <- cs$prior$CL; tvv <- cs$prior$V1
  d <- dati(cs, tvcl, tvv)

  # curva di popolazione (ETA = 0) sulla griglia di PIGATO
  pop <- mrgsim_df(zero_re(mod), data = d[d$evid == 1, ], tgrid = unlist(cs$tgrid), obsonly = TRUE)
  pop <- pop[match(round(unlist(cs$tgrid), 4), round(pop$time, 4)), ]
  sc_pop <- max(abs(pop$CP - unlist(cs$pop_curve)) / pmax(unlist(cs$pop_curve), 1e-3))

  stima <- function(qb) {
    est <- mapbayest(mod, data = d, quantile_bound = qb, verbose = FALSE, progress = FALSE)
    eta <- get_eta(est)
    c(CL = tvcl * exp(eta[["ETA1"]]), V = tvv * exp(eta[["ETA2"]]))
  }
  p_def <- stima(0.001)   # limiti predefiniti mapbayr: ETA entro ±3,09 omega
  p_lar <- stima(1e-300)  # limiti di fatto assenti
  mt <- metriche(mod, p_lar[["CL"]], p_lar[["V"]], cs$finestra)
  ris[[length(ris) + 1]] <- data.frame(
    id = cs$id, CL = p_lar[["CL"]], V = p_lar[["V"]], AUC24 = mt[["auc24"]], Cmin = mt[["cmin"]],
    CL_limdef = p_def[["CL"]], V_limdef = p_def[["V"]], scarto_max_pop = sc_pop)
}
res <- do.call(rbind, ris)
write.csv(res, args[2], row.names = FALSE)
cat(sprintf("mrgsolve %s, mapbayr %s, R %s\n", packageVersion("mrgsolve"),
            packageVersion("mapbayr"), getRversion()))
print(res, digits = 5)
