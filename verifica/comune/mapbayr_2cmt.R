# Verifica PIGATO TDM: stima MAP indipendente con mrgsolve + mapbayr, modello a 2 compartimenti
# (ODE; PIGATO usa la soluzione chiusa in alfa/beta). Usabile per ciprofloxacina, dalbavancina,
# aminoglicosidi. Per caso legge dal JSON: prior (CL, V1, Q, V2), omega (SD log; 0 = fisso),
# sigma_prop, sigma_add, storia delle dosi, livelli, finestra per AUC24/Cmin.
# Errore residuo:
#   sigma_add = 0 -> additivo su scala log (DV = CP*exp(EPS(2)); log_transformation di mapbayr)
#   sigma_add > 0 -> combinato su scala lineare, var = (sp*pred)^2 + sa^2 (come PIGATO)
# OMEGA e SIGMA come varianze. AUC per integrazione (compartimento AUC).
#
# Uso: Rscript mapbayr_2cmt.R casi.json risultati.csv

.libPaths(c("~/Rlib", .libPaths()))  # libreria locale facoltativa
suppressPackageStartupMessages({
  library(jsonlite)
  library(mrgsolve)
  library(mapbayr)
})
args <- commandArgs(trailingOnly = TRUE)
casi <- fromJSON(args[1], simplifyVector = FALSE)

base <- '
$PARAM @annotated @covariates
TVCL : 1  : CL (L/h)
TVV1 : 10 : V1 (L)
TVQ  : 1  : Q (L/h)
TVV2 : 10 : V2 (L)
$OMEGA 0.1 0.1 0.1 0.1
$SIGMA 0.01 0.01
$CMT @annotated
CENT : Centrale (mg) [ADM, OBS]
PERI : Periferico (mg) []
AUC  : AUC (mg*h/L) []
$MAIN
double CL = TVCL * exp(ETA(1));
double V1 = TVV1 * exp(ETA(2));
double Q  = TVQ  * exp(ETA(3));
double V2 = TVV2 * exp(ETA(4));
$ODE
dxdt_CENT = -(CL + Q) / V1 * CENT + Q / V2 * PERI;
dxdt_PERI =  Q / V1 * CENT - Q / V2 * PERI;
dxdt_AUC  = CENT / V1;
$TABLE
double CP = CENT / V1;
ERRORE
$CAPTURE DV CP
$SET rtol = 1e-10, atol = 1e-12, maxsteps = 1e7
'
mod_log <- mcode("pk2_log", sub("ERRORE", "double DV = CP * exp(EPS(2));", base), quiet = TRUE)
mod_cmb <- mcode("pk2_cmb", sub("ERRORE", "double DV = CP * (1 + EPS(1)) + EPS(2);", base), quiet = TRUE)

ev_df <- function(dosi) do.call(rbind, lapply(dosi, function(s)
  data.frame(ID = 1, time = s$ore, evid = 1, cmt = 1, amt = s$dose, rate = s$dose / s$tinf,
             mdv = 1, DV = NA_real_)))

ris <- list()
for (cs in casi) {
  pr <- cs$prior
  tv <- c(CL = pr$CL, V1 = pr$V1, Q = pr$Q, V2 = pr$V2)
  om <- sapply(names(tv), function(k) if (is.null(cs$omega[[k]])) 0 else cs$omega[[k]])
  if (cs$sigma_add > 0) {
    mod <- mod_cmb %>% smat(dmat(cs$sigma_prop^2, cs$sigma_add^2))
  } else {
    mod <- mod_log %>% smat(dmat(0, cs$sigma_prop^2))
  }
  mod <- mod %>% omat(dmat(om[["CL"]]^2, om[["V1"]]^2, om[["Q"]]^2, om[["V2"]]^2))

  d <- ev_df(cs$storia)
  o <- do.call(rbind, lapply(cs$livelli, function(l)
    data.frame(ID = 1, time = l$ore, evid = 0, cmt = 1, amt = 0, rate = 0, mdv = 0, DV = l$conc)))
  d <- rbind(d, o); d <- d[order(d$time, -d$evid), ]
  d$TVCL <- tv[["CL"]]; d$TVV1 <- tv[["V1"]]; d$TVQ <- tv[["Q"]]; d$TVV2 <- tv[["V2"]]

  # curva di popolazione sulla griglia di PIGATO
  pop <- mrgsim_df(zero_re(mod), data = d[d$evid == 1, ], tgrid = unlist(cs$tgrid), obsonly = TRUE)
  pop <- pop[match(round(unlist(cs$tgrid), 4), round(pop$time, 4)), ]
  pc <- unlist(cs$pop_curve)
  sc_pop <- max(abs(pop$CP - pc) / pmax(pc, 1e-3 * max(pc)))

  est <- mapbayest(mod, data = d, quantile_bound = 1e-300, verbose = FALSE, progress = FALSE)
  eta <- get_eta(est)
  sel <- which(om > 0)
  p <- tv
  for (i in sel) p[i] <- tv[i] * exp(eta[[paste0("ETA", i)]])

  # metriche nella finestra dell'app
  fin <- cs$finestra
  m <- param(zero_re(mod), TVCL = p[["CL"]], TVV1 = p[["V1"]], TVQ = p[["Q"]], TVV2 = p[["V2"]])
  t0 <- fin$start; t1 <- fin$start + fin$len
  out <- mrgsim_df(m, data = ev_df(fin$dosi)[, 1:6], tgrid = tgrid(t0, t1, fin$len / 4000), obsonly = TRUE)
  out <- out[out$time >= t0 - 1e-9 & out$time <= t1 + 1e-9, ]
  # dalbavancina: primo tempo dopo t_val con C < soglia (griglia 0,05 h, poi uniroot)
  t_x <- NA
  if (!is.null(cs$soglia)) {
    tv0 <- cs$t_val
    sim <- mrgsim_df(m, data = ev_df(cs$storia)[, 1:6], tgrid = tgrid(tv0, tv0 + 2400, 0.05), obsonly = TRUE)
    sim <- sim[sim$time >= tv0 - 1e-9, ]
    i <- which(sim$CP < cs$soglia)[1]
    if (!is.na(i)) {
      if (i == 1) t_x <- sim$time[1] else {
        f <- function(tt) mrgsim_df(m, data = ev_df(cs$storia)[, 1:6], tgrid = tt, obsonly = TRUE)$CP[1] - cs$soglia
        t_x <- uniroot(f, c(sim$time[i - 1], sim$time[i]), tol = 1e-8)$root
      }
    }
  }
  ris[[length(ris) + 1]] <- data.frame(
    id = cs$id, CL = p[["CL"]], V1 = p[["V1"]], Q = p[["Q"]], V2 = p[["V2"]],
    AUC24 = (out$AUC[nrow(out)] - out$AUC[1]) * fin$fattore, Cmin = min(out$CP),
    t_sotto_soglia = t_x, scarto_max_pop = sc_pop)
}
res <- do.call(rbind, ris)
write.csv(res, args[2], row.names = FALSE)
cat(sprintf("mrgsolve %s, mapbayr %s, R %s\n", packageVersion("mrgsolve"),
            packageVersion("mapbayr"), getRversion()))
print(res, digits = 5)
