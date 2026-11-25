# Verifica PIGATO TDM: ricalcolo indipendente dei valori di popolazione (prior), di omega e sigma
# dalle formule delle pubblicazioni, a partire dalle covariate del paziente. Non usa codice PIGATO.
# Confronta con i valori che PIGATO ha passato al MAP (campi prior, omega_pubblicati/omega, sigma).
#
# Uso: Rscript prior_indipendenti.R casi.json confronto_prior.csv

suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
casi <- fromJSON(args[1], simplifyVector = FALSE)
w_cv <- function(cv) sqrt(log(1 + (cv / 100)^2))          # CV% log-normale -> SD log
nz <- function(x, d) if (is.null(x)) d else x

# ---- funzione renale (formule originali)
cockcroft_gault <- function(eta, peso, scr, sesso) {
  if (is.null(scr) || eta < 18) return(NA)
  v <- (140 - eta) * peso / (72 * scr); if (sesso == "Femmina") v * 0.85 else v }
ckd_epi_2021 <- function(eta, scr, sesso) {   # Inker LA, et al. N Engl J Med 2021;385:1737-1749
  if (is.null(scr) || eta < 18) return(NA)
  f <- sesso == "Femmina"; k <- if (f) 0.7 else 0.9; a <- if (f) -0.241 else -0.302
  142 * min(scr / k, 1)^a * max(scr / k, 1)^-1.2 * 0.9938^eta * (if (f) 1.012 else 1) }
mosteller <- function(peso, h) if (is.null(h)) NA else sqrt(peso * h / 3600)

# ---- daptomicina
WEI <- list(`3-6m` = c(19.7, 5.5, 0.17), `7-12m` = c(19.6, 5.8, 0.17), `13-24m` = c(21.8, 3.0, 0.11),
            `2-6y` = c(20.0, 6.1, 0.11), `7-11y` = c(17.0, 5.7, 0.11), `12-17y` = c(11.0, 2.1, 0.11))
YE <- rbind(`12-17y` = c(0.98, 1.08, 1.24, 1.33), `7-11y` = c(1.00, 1.19, 1.49, 1.73),
            `2-6y` = c(1.02, 1.28, 1.74, 2.20), `1-2y` = c(1.04, 1.29, 1.72, 2.13))
colnames(YE) <- c("Lieve", "Moderata", "Grave", "ESRD")
prior_dapto <- function(cv) {
  gg <- cv$gg; peso <- cv$peso; eta <- gg / 365.25; dia <- nz(cv$dialisi, "Nessuna")
  if (gg >= 4380 && gg < 6570 && grepl("HD|CCPD", dia)) {          # Lim 2021, scalato a 40 kg
    s <- peso / 40
    r <- list(CL = (0.0226 + 0.158) * s, V1 = 5.90 * s, wCL = 0.0878, wV = 0.112, sig = 0.085, tinf = 0.5)
  } else if (gg < 6570) {                                           # Wei 2020 (+ Ye 2022)
    g <- max(gg, 91)
    b <- if (g < 213) "3-6m" else if (g < 365) "7-12m" else if (g < 730) "13-24m" else
         if (g < 2190) "2-6y" else if (g < 4380) "7-11y" else "12-17y"
    p <- WEI[[b]]
    cl <- p[1] * peso / 1000
    cat_ye <- sub(" .*", "", dia)
    if (cat_ye %in% colnames(YE)) {
      by <- if (b %in% c("3-6m", "7-12m", "13-24m")) "1-2y" else b
      cl <- cl / YE[by, cat_ye]
    }
    r <- list(CL = cl, V1 = p[3] * peso, wCL = w_cv(100 * p[2] / p[1]), wV = 0.25, sig = 0.12,
              tinf = if (gg < 2190) 1 else 0.5)
  } else {                                                          # Soraluce 2018
    clcr <- cockcroft_gault(eta, peso, cv$crea, cv$sesso); if (is.na(clcr)) clcr <- 60
    r <- list(CL = 0.160 + 0.367 * clcr / 49, V1 = 12.30, wCL = w_cv(36.7), wV = w_cv(27.8), sig = 0.123, tinf = 0.5)
  }
  if (isTRUE(cv$crrt)) r$CL <- r$CL + cv$Sc * cv$Qef
  list(par = c(CL = r$CL, V1 = r$V1), om = c(CL = r$wCL, V1 = r$wV), sp = r$sig, sa = 0, tinf = r$tinf)
}

# ---- ciprofloxacina
prior_cipro <- function(cv, modello) {
  peso <- cv$peso; eta <- cv$gg / 365.25; eg <- cv$egfr
  if (grepl("^Hirt", modello)) {
    fw <- peso / 20; e <- nz(eg, 199)
    par <- c(CL = 9.07 * fw^0.75 * (e / 199)^0.337 * (max(eta, 0.01) / 7.3)^0.0828,
             V1 = 26.5 * fw, Q = 3.9 * fw^0.75, V2 = 17.4 * fw)
    om <- c(CL = 0.347, V1 = 0, Q = 0, V2 = 0); sp <- 0.423; sa <- 0
  } else if (grepl("^Rajagopalan", modello)) {
    fw <- peso / 70
    cl <- 30.3 * fw^0.75 * (1 + 0.045 * (eta - 2.5))
    if (!is.null(eg) && eg < 100) { f <- 0.145 / (0.145 + 0.167); cl <- cl * ((1 - f) + f * eg / 100) }
    par <- c(CL = cl, V1 = 56.7 * fw, Q = 37.5 * fw^0.75, V2 = 89.8 * fw)
    om <- c(CL = 0.40, V1 = 0.40, Q = 0.30, V2 = 0.30); sp <- 0.27; sa <- 0.04
  } else {                                                          # Forrest 1993
    clcr <- nz(eg, ckd_epi_2021(eta, cv$crea, cv$sesso)); if (is.na(clcr)) clcr <- 90
    bsa <- mosteller(peso, cv$altezza); if (is.na(bsa)) bsa <- 1.73
    par <- c(CL = (0.00145 * clcr + 0.167) * peso, V1 = 0.69 * peso, Q = 38 * bsa / 1.73, V2 = 0.51 * peso)
    om <- c(CL = w_cv(44), V1 = w_cv(26), Q = w_cv(24), V2 = w_cv(33)); sp <- 0.15; sa <- 0
  }
  if (isTRUE(cv$crrt)) par[["CL"]] <- par[["CL"]] + cv$Sc * cv$Qef
  list(par = par, om = om, sp = sp, sa = sa, tinf = 1)
}

# ---- dalbavancina (allometria a 70 kg: CL e Q esponente 0,75, volumi lineari)
prior_dalba <- function(cv, modello) {
  if (grepl("Mariani", modello)) { tv <- c(CL = 0.0305, V1 = 4.81, Q = 0.0081, V2 = 4.81); allo <- TRUE
    om <- c(CL = 0.198, V1 = 0, Q = 0, V2 = 0); sp <- 0.375
  } else if (grepl("^Cojutti", modello)) { tv <- c(CL = 0.106, V1 = 17.4, Q = 0.103, V2 = 15.1); allo <- FALSE
    om <- c(CL = w_cv(36.21), V1 = w_cv(44.27), Q = 0, V2 = w_cv(62.34)); sp <- 0.18
  } else { tv <- c(CL = 0.0273, V1 = 3.6, Q = 0.0225, V2 = 6.4); allo <- TRUE
    om <- c(CL = w_cv(22), V1 = w_cv(17.3), Q = w_cv(55.9), V2 = w_cv(30.1)); sp <- 0.144 }
  if (allo) { f <- cv$peso / 70; tv <- tv * c(f^0.75, f, f^0.75, f) }
  list(par = tv, om = om, sp = sp, sa = 0, tinf = 0.5)
}

# ---- aminoglicosidi
ffm_kg <- function(peso, h, eta, sesso) {      # Al-Sallami 2015 (< 18 anni), Janmahasatian 2005
  bmi <- peso / (h / 100)^2
  if (sesso == "Femmina") { j <- 9270 * peso / (8780 + 244 * bmi)
    if (eta < 18) j <- j * (1.11 + (1 - 1.11) / (1 + (max(eta, 0.01) / 7.1)^-1.1))
  } else { j <- 9270 * peso / (6680 + 216 * bmi)
    if (eta < 18) j <- j * (0.88 + (1 - 0.88) / (1 + (max(eta, 0.01) / 13.4)^-12.7)) }
  j }
scr_attesa <- function(eta, sesso) {            # Ceriotti 2008, forma clinPK (umol/L)
  a <- max(eta, 0.05)
  if (a <= 15) return(-2.37330 - 12.91367 * log(a) + 23.93581 * sqrt(a))
  if (a < 18) return(if (sesso == "Maschio") 9.5471 * a - 87.847 else 4.7137 * a - 15.347)
  if (sesso == "Maschio") 84 else 69.5 }
prior_amino <- function(cv, modello) {
  eta <- cv$gg / 365.25; peso <- cv$peso
  if (grepl("^Llanos", modello)) {
    ffm <- ffm_kg(peso, cv$altezza, eta, cv$sesso)
    pma <- nz(cv$eg_sett, 40) + cv$gg / 7
    gfr <- (ffm / 70)^0.75 * pma^3.33 / (55.4^3.33 + pma^3.33) * 112
    cl <- 5.77 * gfr / 100 * (scr_attesa(eta, cv$sesso) / (cv$crea * 88.4))^0.55
    par <- c(CL = cl, V1 = 21.6 * ffm / 70, Q = 0.62 * (ffm / 70)^0.75, V2 = 13.8 * ffm / 70)
    om <- c(CL = w_cv(sqrt(16^2 + 20.7^2)), V1 = w_cv(21.5), Q = 0, V2 = w_cv(62.4)); sp <- 0.275; sa <- 0.04
  } else if (grepl("^Alhadab", modello)) {
    f <- peso / 70
    par <- c(CL = 11.1 * f^0.75, V1 = 30.2 * f, Q = 4.26 * f^0.75, V2 = 14.9 * f)
    om <- c(CL = w_cv(33), V1 = w_cv(63.1), Q = 0, V2 = 0); sp <- 0.25; sa <- 0.5
  } else {
    f <- peso / 70
    par <- c(CL = 6.37 * f^0.75, V1 = 18.7 * f, Q = 0.39, V2 = 1.32)
    om <- c(CL = w_cv(sqrt(11.7^2 + 6.47^2)), V1 = w_cv(11.66), Q = 0, V2 = 0); sp <- 0.19; sa <- 0
  }
  if (isTRUE(cv$crrt)) par[["CL"]] <- par[["CL"]] + cv$Sc * cv$Qef
  list(par = par, om = om, sp = sp, sa = sa, tinf = 0.5)
}

rel <- function(a, b) ifelse(b == 0 & a == 0, 0, 100 * (a - b) / b)
ris <- list()
for (cs in casi) {
  dapto <- grepl("^(Wei|Lim|Soraluce)", cs$modello)
  ind <- if (dapto) prior_dapto(cs$cov) else if (grepl("Mariani|^Cojutti|^Baiardi", cs$modello))
    prior_dalba(cs$cov, cs$modello) else if (grepl("^(Llanos|Alhadab|Hennig)", cs$modello))
    prior_amino(cs$cov, cs$modello) else prior_cipro(cs$cov, cs$modello)
  pig_om <- nz(cs$omega_pubblicati, cs$omega)
  k <- names(ind$par)
  dpar <- sapply(k, function(x) rel(cs$prior[[x]], ind$par[[x]]))
  dom <- sapply(k, function(x) rel(nz(pig_om[[x]], 0), ind$om[[x]]))
  ris[[length(ris) + 1]] <- data.frame(
    id = cs$id, modello = sub(",.*", "", cs$modello),
    max_scarto_prior_pct = max(abs(dpar)), max_scarto_omega_pct = max(abs(dom)),
    scarto_sigma_prop_pct = rel(cs$sigma_prop, ind$sp), scarto_sigma_add_pct = rel(cs$sigma_add, ind$sa),
    tinf_uguale = isTRUE(all.equal(cs$t_inf, ind$tinf)))
}
res <- do.call(rbind, ris)
write.csv(res, args[2], row.names = FALSE)
print(res, digits = 3)
