# The RSV infant products and the cocoon (GlobalMix India #76, decided October
# 25, 2026) -----------------------------------------------------------------------
#
# The infant products are a dose line of their own, beside the adult budget of
# vax_general() and whatever arm orders it (#58 items 10 and 13): every arm of a
# cell can carry the same infant program, and a strategy only reorders the adult
# doses. They used to be strategies drawing on the adult budget, so the direct
# product against cocooning comparison turned on the rollout (#76).
# vax.infant.program names the product:
#   nirsevimab  the long-acting monoclonal, given to the infant: a catch-up dose
#               on the day the campaign opens to the infants under
#               vax.infant.catchup.months (8), and from that day a dose at birth
#   maternal    the maternal vaccine, given in pregnancy, which protects the
#               newborn from birth by passive transfer. The program is taken as
#               already running: the infants of the starting population were
#               born to a vaccinated mother with the coverage probability and
#               are protected from their birth, as is every covered newborn of
#               the run. With vax.maternal.mother, the mother on record is also
#               vaccinated, with dose 1 of the run's RSV product,
#               vax.maternal.lead days (35) before the birth.
# Coverage is vax.infant.coverage. An infant is covered when its draw iprog.u,
# made once at initialization or birth on a side stream (with_side_rng()), is
# below it: the infants covered at 50% are among those covered at 80%, the two
# products reach the same infants, and every arm of a cell covers the same
# infants while their random streams agree. The adult product's willingness does
# not apply. A covered infant's iprog.time is the day of its dose (nirsevimab)
# or of its birth (maternal), and its protection follows passive_ve().

infant_program <- function(dat) param_or(dat, "vax.infant.program", "none")

# The product's efficacy profile: the peak marginal VE by outcome, as
# vax_ve_outcomes, and the onset delay, the end of the plateau and the day the
# protection reaches 0, in days from the dose (or the birth)
infant_prog_pars <- function(dat, prog = infant_program(dat)) {
  f <- function(x) get_param(dat, paste0("vax.", prog, ".", x))
  list(peak = vapply(stats::setNames(nm = vax_ve_outcomes), function(o) f(paste0("ve.peak.", o)), numeric(1)),
       delay = f("ve.delay"), plateau = f("ve.plateau"), duration = f("ve.duration"))
}

# What compute_ve() needs to protect the infants who hold the product; NULL on a
# run with no infant program
infant_prog_state <- function(dat) {
  prog <- infant_program(dat)
  if (identical(prog, "none")) return(NULL)
  list(time = get_attr(dat, "iprog.time"), pars = infant_prog_pars(dat, prog))
}

# Protection of a passive antibody product t days after the dose, or after the
# birth for the maternal vaccine: 0 before ve.delay, the peak to day ve.plateau,
# then a straight line to 0 at day ve.duration, and 0 after. One shape for every
# outcome, so the stages stay ordered on every day. The exponential waning of
# the active products never reaches 0, and for nirsevimab it cannot draw the
# plateau the trials show: efficacy is flat through 150 to 180 days and gone by
# the second season (the project's rsv_vaccine_efficacy_2026-10-20.md,
# Section 2c).
passive_ve <- function(t, pars, outcome) {
  v <- pars$peak[[outcome]] * pmin(1, pmax(0, (pars$duration - t) / (pars$duration - pars$plateau)))
  v[t < pars$delay] <- 0
  v
}

# Checked at initialization
check_infant_program <- function(dat) {
  prog <- param_or(dat, "vax.infant.program", "none")
  if (length(prog) != 1 || !prog %in% c("none", infant_products)) {
    stop("vax.infant.program must be 'none' or one of ", toString(infant_products), " (GlobalMix India #76)",
         call. = FALSE)
  }
  if (prog == "none") return(invisible(TRUE))
  pa <- get_pathogen(dat)
  if (pa != "rsv") {
    stop("vax.infant.program '", prog, "' is an RSV product, and the run's pathogen is ", pa,
         " (GlobalMix India #76)", call. = FALSE)
  }
  one <- function(name, lo, hi) {
    v <- get_param(dat, name, override.null.error = TRUE)
    if (is.null(v) || !is.numeric(v) || length(v) != 1 || is.na(v) || v < lo || v > hi) {
      stop(name, " must be one number in [", lo, ", ", hi, "] (GlobalMix India #76)", call. = FALSE)
    }
    v
  }
  one("vax.infant.coverage", 0, 1)
  one("vax.infant.catchup.months", 0, 12)
  for (o in vax_ve_outcomes) one(paste0("vax.", prog, ".ve.peak.", o), 0, 1)
  delay <- one(paste0("vax.", prog, ".ve.delay"), 0, Inf)
  plateau <- one(paste0("vax.", prog, ".ve.plateau"), 0, Inf)
  duration <- one(paste0("vax.", prog, ".ve.duration"), 0, Inf)
  if (plateau < delay || duration <= plateau) {
    stop("the ", prog, " profile needs ve.delay <= ve.plateau < ve.duration (GlobalMix India #76)", call. = FALSE)
  }
  p <- infant_prog_pars(dat, prog)$peak
  if (p[["infect"]] > p[["clinical"]] || p[["clinical"]] > p[["hosp"]] || p[["clinical"]] > p[["death"]]) {
    stop("the ", prog, " efficacy must not fall along the stages (infect <= clinical <= hosp, and ",
         "clinical <= death) (GlobalMix India #76)", call. = FALSE)
  }
  if (prog == "maternal") {
    m <- as.logical(param_or(dat, "vax.maternal.mother", FALSE))
    if (length(m) != 1 || is.na(m)) stop("vax.maternal.mother must be TRUE or FALSE (GlobalMix India #76)", call. = FALSE)
    one("vax.maternal.lead", 0, 280)
    if (m && is.null(get_param(dat, paste0("vax.", pa, ".n.doses"), override.null.error = TRUE))) {
      stop("vax.maternal.mother vaccinates the mother with the run's RSV product, and the run has none ",
           "(GlobalMix India #76)", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# The maternal program at initialization: the covered infants of the starting
# population hold the product from their birth, 1 less their age in days (as
# init_last_birth() dates it), and their mothers on record are vaccinated
init_infant_program <- function(dat) {
  if (infant_program(dat) != "maternal") return(dat)
  is_infant <- get_attr(dat, "is_infant")
  ids <- which(get_attr(dat, "active") == 1 & is_infant %in% TRUE &
                 get_attr(dat, "iprog.u") < get_param(dat, "vax.infant.coverage"))
  birth <- 1 - round(get_attr(dat, "age")[ids] * days_per_year)
  time <- get_attr(dat, "iprog.time")
  time[ids] <- birth
  dat <- set_attr(dat, "iprog.time", time)
  cm <- credit_mothers(dat, ids, birth)
  dat <- set_epi(cm$dat, "iprog.init", length(ids), at = 1)
  set_epi(dat, "iprog.init.mother", cm$n, at = 1)
}

# The nirsevimab catch-up, once, on the day the campaign opens: every covered
# living infant under vax.infant.catchup.months who holds no dose yet, whatever
# their infection state (no one is tested). Returns dat and the doses given.
infant_program_catchup <- function(dat, at) {
  if (infant_program(dat) != "nirsevimab" || is.null(dat$run$vax.campaign) || !is.null(dat$run$iprog.catchup)) {
    return(list(dat = dat, n = 0L))
  }
  time <- get_attr(dat, "iprog.time")
  ids <- which(get_attr(dat, "active") == 1 & get_attr(dat, "age") * 12 < get_param(dat, "vax.infant.catchup.months") &
                 is.na(time) & get_attr(dat, "iprog.u") < get_param(dat, "vax.infant.coverage"))
  time[ids] <- at
  dat <- set_attr(dat, "iprog.time", time)
  dat$run$iprog.catchup <- at
  list(dat = dat, n = length(ids))
}

# Birth doses (arrival(), after the birth placement): every covered newborn,
# from the day the campaign opens for nirsevimab and on every day for the
# maternal program, which also vaccinates the mother. The step's birth doses
# and mothers vaccinated are nVax.iprog.birth and nVax.iprog.mother.
infant_program_births <- function(dat, at, new_ids) {
  prog <- infant_program(dat)
  n_birth <- 0L
  n_mother <- 0L
  if (prog != "none" && length(new_ids) > 0 && (prog == "maternal" || at >= get_param(dat, "vax.start"))) {
    cov <- new_ids[get_attr(dat, "iprog.u")[new_ids] < get_param(dat, "vax.infant.coverage")]
    if (length(cov) > 0) {
      time <- get_attr(dat, "iprog.time")
      time[cov] <- at
      dat <- set_attr(dat, "iprog.time", time)
      n_birth <- length(cov)
      if (prog == "maternal") {
        cm <- credit_mothers(dat, cov, at)
        dat <- cm$dat
        n_mother <- cm$n
      }
    }
  }
  dat <- set_epi(dat, "nVax.iprog.birth", n_birth, at = at)
  set_epi(dat, "nVax.iprog.mother", n_mother, at = at)
}

# Vaccinates the living mothers on record (mother.uid) of the covered children
# child_ids, born on the days birth, when vax.maternal.mother: dose 1 of the
# run's RSV product, given vax.maternal.lead days before the birth. A mother
# already dosed keeps her dose. Returns dat and the mothers vaccinated.
credit_mothers <- function(dat, child_ids, birth) {
  if (length(child_ids) == 0 || !isTRUE(as.logical(param_or(dat, "vax.maternal.mother", FALSE)))) {
    return(list(dat = dat, n = 0L))
  }
  birth <- rep_len(birth, length(child_ids))
  m <- match(get_attr(dat, "mother.uid")[child_ids], get_attr(dat, "unique_id"))
  vax <- get_attr(dat, "vax")
  ok <- !is.na(m)
  ok[ok] <- get_attr(dat, "active")[m[ok]] == 1 & vax[m[ok]] == 0
  sel <- which(ok & !duplicated(m))
  if (length(sel) == 0) return(list(dat = dat, n = 0L))
  ldt <- get_attr(dat, "last.dose.time")
  vax[m[sel]] <- 1
  ldt[m[sel]] <- birth[sel] - param_or(dat, "vax.maternal.lead", 35)
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "last.dose.time", ldt)
  list(dat = dat, n = length(sel))
}

# The infant line's outputs, every step (vax_general()): today's catch-up doses,
# and today's doses on the line, the catch-up plus the birth doses (a maternal
# dose is given to the mother and counted once); and the living infants who hold
# the product
infant_program_outputs <- function(dat, at, n_catchup) {
  nb <- get_epi(dat, "nVax.iprog.birth", at = at, override.null.error = TRUE)
  nb <- if (is.null(nb) || is.na(nb)) 0 else nb
  dat <- set_epi(dat, "nVax.iprog.catchup", n_catchup, at = at)
  dat <- set_epi(dat, "nVax.iprog", n_catchup + nb, at = at)
  time <- get_attr(dat, "iprog.time", override.null.error = TRUE)
  is_infant <- get_attr(dat, "is_infant", override.null.error = TRUE)
  n <- if (is.null(time) || is.null(is_infant)) 0L else
    sum(get_attr(dat, "active") == 1 & is_infant %in% TRUE & !is.na(time))
  set_epi(dat, "n.iprog.infant", n, at = at)
}

# The cocoon tier: the living adults vax.cocoon.age.min (18) to under
# vax.cocoon.age.max (50) of the households with a living infant under
# vax.cocoon.infant.months (6) today, 1 for a member and 0 otherwise. Read on
# every arm, so the doses an arm gives the tier are known (nVax.cocoon).
cocoon_tier <- function(dat) {
  active <- get_attr(dat, "active")
  out <- integer(length(active))
  hh <- get_attr(dat, "hh.ids", override.null.error = TRUE)
  if (is.null(hh) || all(is.na(hh))) return(out)
  age <- get_attr(dat, "age")
  live <- active == 1 & !is.na(hh)
  young <- live & age * 12 < param_or(dat, "vax.cocoon.infant.months", 6)
  if (!any(young)) return(out)
  as.integer(live & hh %in% unique(hh[young]) & age >= param_or(dat, "vax.cocoon.age.min", 18) &
               age < param_or(dat, "vax.cocoon.age.max", 50))
}
