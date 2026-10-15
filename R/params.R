# Parameters the GMC19 model reads (GlobalMix India #67)
#
# One list of every parameter name the GMC19 modules read. The project's
# parameter table used to carry a third of its rows that nothing read
# (CorporateMix vaccine rows, half.life, ...), which is how the inherited
# modifiers removed in #66 went unnoticed. init_gmc19() warns on any name
# outside this list, and tests/testthat/test-gmc19-params.R scans the package
# source and fails if a module reads a name missing here or a name here is
# never read.
#
# Names are as the modules read them, after param.net() reassembles a table's
# name_1, ..., name_n rows into one vector.

# Read by the GMC19 modules on both designs
gmc19_params_common <- c(
  # demography and aging
  "time.unit", "age.breaks", "age.grps", "arrival.age", "birth.prop.male", "infant.age.max",
  "mort.rates", "birth.rate",
  # layers and transmission
  "layer.names", "inf.prob", "exposure.rate", "inf.prob.a.rr", "sympt.contact.rr", "susc.age.rr",
  "susc.infant.rr", "susc.infant.age",
  "contact.weight", "contact.weight.equal", "seasonal.amp", "seasonal.phase", "import.rate",
  "bg.hazard", "bg.mix",
  # natural history and severity
  "ea.rate", "ar.rate", "eip.rate", "ipic.rate", "ich.rate", "icr.rate", "hr.rate", "rs.rate",
  "prop.clinical", "prop.hospit", "prop.hospit.infant", "prop.hospit.infant6",
  "dis.death.rate", "dis.death.rate.infant", "dis.death.rate.infant6",
  # vaccination and targeting; the per-dose schedule is in vax_dose_fields
  "pathogen", "vax.strategy", "vax.supply.rate", "vax.supply.total", "vax.age.breaks",
  "vax.infant.eligible", "vax.infant.product", "vax.season.length", "vax.target.uptake",
  "vax.hybrid.frac", "vax.degree.score", "compute.degree"
)

# Read only by the relationship model (#63): the activation step, the traits and
# the new contacts, and the birth placement (#65). No name is read by the legacy
# four-layer design alone.
gmc19_params_relationship <- c(
  "act.o", "act.kappa", "act.kappa.mult", "act.q", "act.q.mult", "act.p11", "act.s",
  "act.pi", "act.theta.F", "act.theta.I", "act.rho", "act.class.mult", "act.theta.N",
  "act.pi.N", "act.episodes", "act.window", "act.schedule", "act.diagnostics",
  "act.seed.newborns", "act.reset.cross", "act.occ.age", "new.rate", "new.rate.mult", "new.mix", "new.share.norel",
  # birth placement (#65), with a population that carries sex
  "birth.age.weight", "birth.hh.coef", "birth.interval.min"
)

# Set by EpiModel itself: param.net() adds act.rate (1) and vital, init_nets()
# sets groups (which the edges correction reads), and use_scenario() sets
# .param.updater.list and .scenario.id
epimodel_params <- c("groups", "vital", "act.rate", ".param.updater.list", ".scenario.id")

# The per-dose vaccine schedule, read by build_vax_schedule() as
# vax.<pathogen>.<field>.dose<i> for doses 1 to vax.<pathogen>.n.doses. start
# and rate hold one value per vaccine age group, the others one per dose. The
# death VE (ve.*.death) scales the disease-death hazard of a vaccinated
# symptomatic case, the dominant and most durable effect of the real vaccines.
vax_dose_fields <- c(
  "start", "rate", "interval", "annual",
  "ve.peak.infect", "ve.halflife.infect", "ve.floor.infect",
  "ve.peak.clinical", "ve.halflife.clinical", "ve.floor.clinical",
  "ve.peak.hosp", "ve.halflife.hosp", "ve.floor.hosp",
  "ve.peak.death", "ve.halflife.death", "ve.floor.death",
  "ve.delay"
)
vax_dose_by_age <- c("start", "rate")

# Read only by the CorporateMix modules still in the package
# (infect_covid_corporate(), progress_covid_corporate(), vax_covid_corporate(),
# dx_covid_corporate(), isolate_covid_corporate()), which are not part of the
# GMC19 model, so a GMC19 run that sets one warns. The parameters of the
# removed modifiers they also read are in retired_params, and stop the run.
corporate_params <- c(
  "half.life",
  paste0("vax", 1:4, ".rate"), paste0("vax", 1:4, ".boost"), paste0("vax", 1:4, ".boost.start"),
  paste0("vax", 1:4, ".rr.infect"), paste0("vax", 1:4, ".rr.clinical"),
  paste0("vax", 1:4, ".rr.hosp"),
  "vax1.start", "vax3.start", "vax4.start", "vax2.interval", "vax3.interval", "vax4.interval"
)

#' Parameters the GMC19 model reads
#'
#' @description
#' Every parameter name that the GMC19 modules, or EpiModel's netsim internals
#' under them, read: the names as [EpiModel::param.net()] stores them, after a
#' parameter table's `name_1, ..., name_n` rows are reassembled into one vector.
#' [init_gmc19()] warns, once, on any other name in the parameters or their
#' updaters (GlobalMix India #67).
#'
#' @details
#' The relationship design (GlobalMix India #63) reads every name the legacy
#' four-layer design reads, plus the activation, trait and new-contact
#' parameters (`act.*`, `new.*`) and the birth placement (`birth.age.weight`,
#' `birth.hh.coef`, `birth.interval.min`; #65); no name is read by the legacy
#' design alone.
#' `groups`, `vital`, `act.rate`, `.param.updater.list` and `.scenario.id` are
#' set by EpiModel itself, and none is a row of a parameter table. The per-dose
#' vaccine schedule is read as `vax.<pathogen>.<field>.dose<i>`, for doses 1 to
#' `vax.<pathogen>.n.doses`.
#'
#' Not listed: `hh.pairs` and the parameters of the removed CorporateMix
#' modifiers, which are read only to stop a run that sets them, and the
#' parameters read only by the CorporateMix modules still in the package
#' (`vax1.rate`, `half.life`, ...), which are not part of the GMC19 model.
#'
#' @param design `"all"` (the default), `"relationship"` or `"legacy"`: the
#'   names read on that design.
#' @param n.doses Number of doses of each pathogen's vaccine schedule, named by
#'   pathogen (the `vax.<pathogen>.n.doses` parameters). The default is the
#'   GlobalMix India study's: three COVID-19 doses and one each for influenza
#'   and RSV. [init_gmc19()] passes the run's own values.
#' @return A sorted character vector.
#' @export
gmc19_known_params <- function(design = c("all", "relationship", "legacy"),
                               n.doses = c(covid = 3, flu = 1, rsv = 1)) {
  design <- match.arg(design)
  if (length(n.doses) > 0 && (is.null(names(n.doses)) || any(names(n.doses) == "") ||
                              anyNA(n.doses) || any(n.doses < 0))) {
    stop("`n.doses` must be non-negative dose counts named by pathogen", call. = FALSE)
  }
  schedule <- unlist(lapply(names(n.doses), function(p) {
    c(paste0("vax.", p, ".n.doses"),
      as.vector(outer(vax_dose_fields, seq_len(n.doses[[p]]),
                      function(f, i) paste0("vax.", p, ".", f, ".dose", i))))
  }))
  sort(unique(c(gmc19_params_common, if (design != "legacy") gmc19_params_relationship,
                epimodel_params, schedule)))
}

# The model's rates and probabilities are per day, and EpiModel no longer sets
# time.unit (EpiModel #1097). inf.prob is set by the pathogen scenario; its
# default in the parameter table is NA, so a run that skips the scenario stops
# here rather than running a stale value with no disease mortality (#67).
check_required_params <- function(dat) {
  time.unit <- get_param(dat, "time.unit", override.null.error = TRUE)
  if (is.null(time.unit)) {
    stop("the parameters have no time.unit; add time.unit = 1 (days per time step), ",
         "which EpiModel no longer sets", call. = FALSE)
  }
  if (!isTRUE(length(time.unit) == 1 && time.unit == 1)) {
    stop("time.unit is ", toString(time.unit), ", but the model's rates and probabilities ",
         "are per day: set time.unit = 1", call. = FALSE)
  }
  inf.prob <- get_param(dat, "inf.prob", override.null.error = TRUE)
  if (is.null(inf.prob) || anyNA(inf.prob)) {
    stop("inf.prob is set by the pathogen scenario, scenarios_pathogen.csv; the registry ",
         "default is NA. Apply the pathogen's scenario with use_scenario() before netsim().",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Births and background deaths (GlobalMix India #64). Births are at the crude
# birth rate, independent of deaths, so the relationship design must set
# birth.rate (0 for the closed cohort); births used to replace each step's
# deaths, which only the legacy four-layer model still does. mort.rates is the
# daily death probability by single year of age.
check_demography <- function(dat) {
  birth.rate <- get_param(dat, "birth.rate", override.null.error = TRUE)
  if (is.null(birth.rate)) {
    if (activation_on(dat)) {
      stop("the relationship design needs birth.rate, the crude birth rate per person-day ",
           "(0 for a closed population); births no longer replace deaths (GlobalMix India #64)",
           call. = FALSE)
    }
  } else if (!(is.numeric(birth.rate) && length(birth.rate) == 1 && is.finite(birth.rate) &&
               birth.rate >= 0 && birth.rate < 1)) {
    stop("birth.rate must be one daily probability in [0, 1), not ", toString(birth.rate),
         call. = FALSE)
  }
  mort.rates <- get_param(dat, "mort.rates", override.null.error = TRUE)
  if (!is.null(mort.rates) && !(is.numeric(mort.rates) && length(mort.rates) > 0 &&
                                !anyNA(mort.rates) && all(mort.rates >= 0 & mort.rates <= 1))) {
    stop("mort.rates must be daily death probabilities by single year of age, each in [0, 1]",
         call. = FALSE)
  }
  invisible(TRUE)
}

# One warning listing every parameter, set in the table or by an updater, that
# nothing on this run's design reads
check_known_params <- function(dat) {
  updated <- unlist(lapply(dat$param$.param.updater.list, function(u) names(u$param)))
  dose_names <- grep("^vax\\..+\\.n\\.doses$", names(dat$param), value = TRUE)
  n.doses <- vapply(dose_names, function(p) as.numeric(dat$param[[p]])[1], numeric(1))
  names(n.doses) <- sub("^vax\\.(.+)\\.n\\.doses$", "\\1", dose_names)
  n.doses[is.na(n.doses)] <- 0
  known <- gmc19_known_params(if (activation_on(dat)) "relationship" else "legacy", n.doses)
  unknown <- sort(unique(setdiff(c(names(dat$param), updated), known)))
  if (length(unknown) > 0) {
    warning("no module reads the parameter(s) ", toString(unknown),
            "; remove them from the parameters (GlobalMix India #67)", call. = FALSE)
  }
  invisible(unknown)
}
