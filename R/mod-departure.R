#' @rdname moduleset-gmc19
#' @export
deaths_covid_gmc19 <- function(dat, at) {

  ## Attributes
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")
  vax <- get_attr(dat, "vax")

  ## Parameters
  mort.rates <- get_param(dat, "mort.rates")
  age.breaks <- get_param(dat, "age.breaks")
  nAgeGrp <- length(age.breaks) - 1
  # Disease deaths (GlobalMix India #71). Whether an episode is fatal is drawn
  # once, at symptom onset, in progress_general(), with probability ifr /
  # prop.clinical by severity band and the vaccine's and prior infection's
  # protection against death applied there, so the IFR per infection is the
  # input whatever the stage durations. A fatal case dies when it would
  # otherwise have recovered: at the end of its symptomatic stay (icr.rate) if
  # not admitted, or of its hospital stay (hr.rate). This replaces a daily
  # hazard on every day in ic or h, which made the IFR per infection 2.2 to
  # 4.4 times the daily input, and moved with any change to the stage
  # durations, prop.hospit or hospitalization VE.
  fatal <- get_attr(dat, "fatal")
  hospit <- get_attr(dat, "hospit")
  statusTime <- get_attr(dat, "statusTime")
  inf.count <- get_attr(dat, "inf.count")
  is_infant <- get_attr(dat, "is_infant")

  idsElig <- which(as.logical(active))

  nDeaths <- 0L
  nDisDeaths <- 0L
  age_grp_disdep <- integer(0)
  age_disdep <- numeric(0)
  vax_disdep <- logical(0)
  infant_disdep <- logical(0)

  if (length(idsElig) > 0) {
    # Background mortality by single year of age: the daily rate at index
    # ceiling(age), the last element for every age beyond the vector. The floor
    # at 1 keeps a node of age exactly 0 (a newborn, if a module order ever put
    # departures before aging) on the first year's rate rather than dropping it
    # from the index (GlobalMix India #64).
    age_idx_elig <- pmin(pmax(ceiling(age[idsElig]), 1), length(mort.rates))
    bg_rates <- mort.rates[age_idx_elig]
    if (length(bg_rates) != length(idsElig) || anyNA(bg_rates)) {
      stop("background mortality has no rate for some living node; mort.rates has ",
           length(mort.rates), " elements", call. = FALSE)
    }

    # A fatal episode ends in death at the rate it would have ended in recovery
    st <- status[idsElig]
    due <- fatal[idsElig] %in% 1 & statusTime[idsElig] < at
    dis_rates <- numeric(length(idsElig))
    if (any(due)) {
      dis_rates[due & st == "ic" & hospit[idsElig] %in% 0] <- get_param(dat, "icr.rate")
      dis_rates[due & st == "h"] <- get_param(dat, "hr.rate")
    }

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
        dd <- idsDep[is_dis]
        age_grp_disdep <- cut(age[dd], age.breaks, labels = FALSE, right = FALSE)
        vax_disdep <- vax[dd] >= 1   # vaccinated at death?
        if (!is.null(is_infant)) infant_disdep <- is_infant[dd]
        age_disdep <- age[dd]   # issue #46: for the infant 0-6 / 6-12mo split
        if (sev_diagnostics(dat)) dat <- sev_flows(dat, at, "death", dd, inf.count[dd] == 1)
      }
      # Mark the dead only. nwupdate.FUN removes them from every layer after
      # the arrivals module has run (EpiModel's convention), so a module that
      # runs between here and nwupdate.FUN must skip active == 0.
      dat <- set_attr(dat, "active", 0, posit_ids = idsDep)
      dat <- set_attr(dat, "exitTime", at, posit_ids = idsDep)
    }
  }

  ## Summary output: total departures, disease deaths, disease deaths by age group,
  ## and by vaccination status (incl age x vax cross-tab). The deaths-by-vax split
  ## is what lets the strategy comparison distinguish "vaccinated cases protected"
  ## from "fewer cases overall" and read direct vs indirect mortality protection;
  ## the age marginal alone cannot.
  dat <- set_epi(dat, "d.flow", nDeaths, at = at)
  dat <- set_epi(dat, "d.dis.flow", nDisDeaths, at = at)
  dat <- set_epi(dat, "d.dis.flow.vax", sum(vax_disdep, na.rm = TRUE), at = at)
  dat <- set_epi(dat, "d.dis.flow.unvax", sum(!vax_disdep, na.rm = TRUE), at = at)
  # Infant disease deaths, reported separately so the infant burden (RSV) is not
  # diluted in the coarse youngest age band (d.dis.flow.age1).
  dat <- set_epi(dat, "d.dis.flow.infant", sum(infant_disdep, na.rm = TRUE), at = at)
  # Issue #46: infant deaths split at 6 months, matching the graded severity
  # inputs (the 0-6 and 6-12 month severity bands, #71). Pooling under 1y hid the fact
  # that 0-6mo carries most of the RSV burden, which is exactly the quantity the
  # direct-infant-product arms are meant to move.
  dat <- set_epi(dat, "d.dis.flow.infant.young",
                 sum(infant_disdep & age_disdep <  0.5, na.rm = TRUE), at = at)
  dat <- set_epi(dat, "d.dis.flow.infant.old",
                 sum(infant_disdep & age_disdep >= 0.5, na.rm = TRUE), at = at)
  for (g in seq_len(nAgeGrp)) {
    dat <- set_epi(dat, paste0("d.dis.flow.age", g),
                   sum(age_grp_disdep == g, na.rm = TRUE), at = at)
    dat <- set_epi(dat, paste0("d.dis.flow.unvax.age", g),
                   sum(age_grp_disdep == g & !vax_disdep, na.rm = TRUE), at = at)
  }

  return(dat)
}
