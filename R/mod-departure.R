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
  # Direct disease-death hazard for symptomatic cases (clinical "ic" or hospitalized
  # "h"), by age group, set per-pathogen via the scenario. Decoupled from
  # hospitalization so the IFR is hit directly while prop.hospit stays a literal
  # hospitalization rate; this captures out-of-hospital disease mortality (which
  # dominates for flu/RSV). The daily hazard is calibrated to age-specific IFR
  # targets (dis.death.rate_1..6 in scenarios_pathogen.csv). Absent -> no excess.
  dis.death.rate <- get_param(dat, "dis.death.rate", override.null.error = TRUE)
  if (is.null(dis.death.rate)) dis.death.rate <- rep(0, nAgeGrp)
  # Infant-specific disease-death hazard (issue #23), graded by sub-age because
  # RSV severity is sharply front-loaded: 0-6mo carries ~64-80% of infant burden
  # (dis.death.rate.infant), 6-12mo is intermediate (dis.death.rate.infant6).
  # Both set per pathogen in scenarios_pathogen.csv; absent (covid/flu) -> infants
  # use the age-band rate.
  dis.death.rate.infant  <- get_param(dat, "dis.death.rate.infant", override.null.error = TRUE)
  dis.death.rate.infant6 <- get_param(dat, "dis.death.rate.infant6", override.null.error = TRUE)
  is_infant <- get_attr(dat, "is_infant")

  idsElig <- which(as.logical(active))

  nDeaths <- 0L
  nDisDeaths <- 0L
  age_grp_disdep <- integer(0)
  vax_disdep <- logical(0)
  infant_disdep <- logical(0)

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
    # Infant override (before death-VE), graded by sub-age: 0-6mo (age < 0.5) take
    # the higher dis.death.rate.infant, 6-12mo take dis.death.rate.infant6. Applied
    # to the base rate so a (rare) vaccinated infant still gets the death-VE
    # reduction below.
    .ok <- function(x) !is.null(x) && !is.na(x) && x > 0
    if (!is.null(is_infant) && (.ok(dis.death.rate.infant) || .ok(dis.death.rate.infant6))) {
      inf <- is_infant[idsElig]
      young <- sympt & inf & age[idsElig] <  0.5    # 0-6 months
      older <- sympt & inf & age[idsElig] >= 0.5    # 6-12 months
      if (.ok(dis.death.rate.infant))  dis_rates[young] <- pmin(1, dis.death.rate.infant)
      if (.ok(dis.death.rate.infant6)) dis_rates[older] <- pmin(1, dis.death.rate.infant6)
    }

    # Direct severity/mortality protection: a vaccinated symptomatic case dies at
    # the VE-reduced rate (death-VE), the dominant arm of the real vaccines and the
    # mechanism by which directly vaccinating high-IFR groups averts their deaths.
    if (any(sympt)) {
      ve_death <- compute_ve(at = at, ids = idsElig[sympt], vax = vax,
        vax.age.group = vax_age_group_for(dat),
        last.dose.time = get_attr(dat, "last.dose.time"),
        vax.schedule = build_vax_schedule(dat, get_pathogen(dat)),
        outcome = "death")
      dis_rates[sympt] <- dis_rates[sympt] * ve_death$rr
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
      }
      dat <- set_attr(dat, "active", 0, posit_ids = idsDep)
      dat <- depart_nodes(dat, departures = idsDep)
      attr.length <- unique(vapply(get_attr_list(dat), length, numeric(1)))
      if (attr.length != attributes(dat$run$el[[1]])[["n"]]) {
        stop("mismatch between el and attr length in departures mod")
      }
    }
  }

  ## Summary output: total departures, disease deaths, disease deaths by age group,
  ## and by vaccination status (incl age x vax cross-tab). The deaths-by-vax split
  ## is what lets the strategy comparison distinguish "vaccinated cases protected"
  ## from "fewer cases overall" and read direct vs indirect mortality protection;
  ## the age marginal alone cannot.
  dat <- set_epi(dat, "d.flow", at, nDeaths)
  dat <- set_epi(dat, "d.dis.flow", at, nDisDeaths)
  dat <- set_epi(dat, "d.dis.flow.vax", at, sum(vax_disdep, na.rm = TRUE))
  dat <- set_epi(dat, "d.dis.flow.unvax", at, sum(!vax_disdep, na.rm = TRUE))
  # Infant disease deaths, reported separately so the infant burden (RSV) is not
  # diluted in the coarse youngest age band (d.dis.flow.age1).
  dat <- set_epi(dat, "d.dis.flow.infant", at, sum(infant_disdep, na.rm = TRUE))
  for (g in seq_len(nAgeGrp)) {
    dat <- set_epi(dat, paste0("d.dis.flow.age", g), at,
                   sum(age_grp_disdep == g, na.rm = TRUE))
    dat <- set_epi(dat, paste0("d.dis.flow.unvax.age", g), at,
                   sum(age_grp_disdep == g & !vax_disdep, na.rm = TRUE))
  }

  return(dat)
}
