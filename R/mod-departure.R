#' @rdname moduleset-gmc19
#' @export
deaths_covid_gmc19 <- function(dat, at) {

  ## Attributes
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  ## Parameters
  mort.rates <- get_param(dat, "mort.rates")
  mort.dis.mult <- get_param(dat, "mort.dis.mult", override.null.error = TRUE)
  if (is.null(mort.dis.mult)) mort.dis.mult <- 1   # 1 = no excess disease mortality
  age.breaks <- get_param(dat, "age.breaks")
  nAgeGrp <- length(age.breaks) - 1

  idsElig <- which(as.logical(active))

  nDeaths <- 0L
  nDisDeaths <- 0L
  age_grp_disdep <- integer(0)

  if (length(idsElig) > 0) {
    # Age-indexed background mortality (ages >= 86y use index 86).
    age_idx_elig <- pmin(ceiling(age[idsElig]), 86)
    death_rates <- mort.rates[age_idx_elig]

    # Excess disease mortality for severe (hospitalized) cases: daily death
    # probability is the background age rate scaled by mort.dis.mult (capped at 1).
    # This is the infection -> hospitalization -> death pathway; the mort.dis.mult
    # magnitude is a calibration target (HFR/IFR by age).
    severe <- status[idsElig] == "h"
    death_rates[severe] <- pmin(1, death_rates[severe] * mort.dis.mult)

    dep_local <- which(runif(length(death_rates)) < death_rates)
    idsDep <- idsElig[dep_local]

    if (length(idsDep) > 0) {
      dis_dep <- severe[dep_local]          # which departures were disease deaths
      nDeaths <- length(idsDep)
      nDisDeaths <- sum(dis_dep)
      if (nDisDeaths > 0) {
        age_grp_disdep <- cut(age[idsDep[dis_dep]], age.breaks,
                              labels = FALSE, right = FALSE)
      }
      dat <- set_attr(dat, "active", 0, posit_ids = idsDep)
      dat <- depart_nodes(dat, departures = idsDep)
      attr.length <- unique(vapply(get_attr_list(dat), length, numeric(1)))
      if (attr.length != attributes(dat$run$el[[1]])[["n"]]) {
        stop("mismatch between el and attr length in departures mod")
      }
    }
  }

  ## Summary output: total departures, disease deaths, disease deaths by age group
  dat <- set_epi(dat, "d.flow", at, nDeaths)
  dat <- set_epi(dat, "d.dis.flow", at, nDisDeaths)
  for (g in seq_len(nAgeGrp)) {
    dat <- set_epi(dat, paste0("d.dis.flow.age", g), at,
                   sum(age_grp_disdep == g, na.rm = TRUE))
  }

  return(dat)
}
