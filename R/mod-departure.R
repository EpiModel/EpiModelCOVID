#' @rdname moduleset-gmc19
#' @export
deaths_covid_gmc19 <- function(dat, at) {

  ## Attributes
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  ## Parameters
  mort.rates <- get_param(dat, "mort.rates")
  age.breaks <- get_param(dat, "age.breaks")
  nAgeGrp <- length(age.breaks) - 1
  # Direct disease-death hazard for symptomatic cases (clinical "ic" or hospitalized
  # "h"), by age group, set per-pathogen via the scenario. Decoupled from
  # hospitalization so the IFR is hit directly while prop.hospit stays a literal
  # hospitalization rate; this captures out-of-hospital disease mortality (which
  # dominates for flu/RSV). The daily hazard is calibrated to age-specific IFR
  # targets (dis.death.rate_1..6 in scenarios_pathogen.csv). Absent -> no excess.
  dis.death.rate <- get_param(dat, "dis.death.rate", override.null.error = TRUE)
  if (is.null(dis.death.rate)) dis.death.rate <- rep(0, nAgeGrp)

  idsElig <- which(as.logical(active))

  nDeaths <- 0L
  nDisDeaths <- 0L
  age_grp_disdep <- integer(0)

  if (length(idsElig) > 0) {
    # Age-indexed background mortality (ages >= 86y use index 86).
    age_idx_elig <- pmin(ceiling(age[idsElig]), 86)
    bg_rates <- mort.rates[age_idx_elig]

    # Disease-death hazard applies to symptomatic cases only (ic or h), indexed by
    # model age group (not hospitalization status).
    age_grp_elig <- cut(age[idsElig], age.breaks, labels = FALSE, right = FALSE)
    sympt <- status[idsElig] %in% c("ic", "h")
    dis_rates <- numeric(length(idsElig))
    dis_rates[sympt] <- pmin(1, dis.death.rate[age_grp_elig[sympt]])

    # Background and disease deaths are independent competing draws; a node that
    # draws a disease death is counted as a disease death regardless of background.
    bg_dep  <- runif(length(idsElig)) < bg_rates
    dis_dep <- runif(length(idsElig)) < dis_rates
    dep_mask <- bg_dep | dis_dep
    idsDep <- idsElig[dep_mask]

    if (length(idsDep) > 0) {
      is_dis <- dis_dep[dep_mask]           # which departures were disease deaths
      nDeaths <- length(idsDep)
      nDisDeaths <- sum(is_dis)
      if (nDisDeaths > 0) {
        age_grp_disdep <- cut(age[idsDep[is_dis]], age.breaks,
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
