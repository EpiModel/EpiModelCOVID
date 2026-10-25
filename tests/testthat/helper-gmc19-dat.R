# Fixtures for the GMC19 tests (episodes, immunity, severity, vaccination):
# a netsim_dat object with no contacts, and mocked vaccine effects

ages12 <- c(0.2, 0.7, 3, 7, 15, 25, 35, 45, 55, 65, 75, 85)  # one age in each severity band

# A netsim_dat object with N people and no contacts: one tie type with an empty
# edgelist, so only imports and the background hazard infect. Everyone is
# susceptible, unvaccinated and 30 unless set otherwise.
nh_dat <- function(N, age = 30, ...) {
  param <- do.call(param.net, c(list(
    time.unit = 1, inf.prob = 0.1, exposure.rate = 1, inf.prob.a.rr = 0.5, layer.names = "community",
    age.breaks = c(0, 10, 20, 30, 40, 60, Inf), arrival.age = 0, mort.rates = rep(0, 101),
    ea.rate = 0.25, ar.rate = 0.2, eip.rate = 0.25, ipic.rate = 0.5, ich.rate = 0.3,
    icr.rate = 0.25, hr.rate = 0.15), sev_test_params(), list(...)))
  control <- control.net(type = NULL, nsteps = 400, nsims = 1, verbose = FALSE, initialize.FUN = init_gmc19)
  dat <- create_dat_object(param, init.net(), control)
  dat <- append_core_attr(dat, 1, N)
  na <- rep(NA, N)
  attrs <- list(status = rep("s", N), statusTime = na, infTime = na, clinical = na, hospit = na,
                fatal = na, recTime = na, inf.count = rep(0, N), ever.inf = rep(0, N),
                ever.sympt = rep(0, N), prior.immune = rep(0, N), vax = rep(0, N), vax.willing = rep(1L, N),
                last.dose.time = na, age = rep_len(age, N), deg_work = rep(0, N),
                degree_quartile = rep(1, N), is_bridge = rep(FALSE, N), hh.ids = na)
  for (a in names(attrs)) dat <- set_attr(dat, a, attrs[[a]])
  dat <- set_attr(dat, "is_infant", get_attr(dat, "age") < 1)
  dat$run$el <- list(community = matrix(integer(0), ncol = 2))
  dat$num.nw <- 1
  dat$nwparam <- list(list())
  dat
}

# Vaccine effects mocked: a vaccinated person (vax >= 1) has relative risk
# ve_rr[[outcome]], everyone else 1
mock_ve <- function(ve_rr = c(infect = 1, clinical = 1, hosp = 1, death = 1, infectious = 1), env = parent.frame()) {
  local_mocked_bindings(
    build_vax_schedule = function(dat, disease) NULL,
    compute_ve = function(at, ids, vax, ..., outcome) {
      r <- if (outcome %in% names(ve_rr)) ve_rr[[outcome]] else 1
      list(rr = ifelse(vax[ids] >= 1, r, 1))
    },
    vax_age_group_for = function(dat) rep(NA_integer_, length(get_attr(dat, "active"))),
    .package = "EpiModelCOVID", .env = env)
}


# One pathogen's schedule (the vaccine and allocation tests) of n.doses doses with marginal VE by outcome; a value
# may be one per vaccine age group (5)
ve_params <- function(infect = 0.3, clinical = 0.5, hosp = 0.8, death = 0.9, infectious = 0.5,
                      halflife = 200, floor = 0, delay = 14, n.doses = 1, ...) {
  p <- list(vax.covid.n.doses = n.doses, vax.age.breaks = c(0, 5, 18, 50, 65, Inf), pathogen = "covid",
            # no campaign unless a test sets one (#74)
            vax.budget = 0, vax.delivery.days = 1, vax.start = 1, vax.willing = 1,
            vax.covid.age.min = 0, vax.covid.age.max = Inf)
  for (i in seq_len(n.doses)) {
    d <- paste0(".dose", i)
    p[[paste0("vax.covid.interval", d)]] <- if (i == 1) NA else 21
    p[[paste0("vax.covid.annual", d)]] <- FALSE
    p[[paste0("vax.covid.ve.delay", d)]] <- delay
    for (o in c("infect", "clinical", "hosp", "death", "infectious")) {
      p[[paste0("vax.covid.ve.peak.", o, d)]] <- get(o)
      p[[paste0("vax.covid.ve.halflife.", o, d)]] <- halflife
      p[[paste0("vax.covid.ve.floor.", o, d)]] <- floor
    }
  }
  utils::modifyList(p, list(...))
}


# The RSV infant products (GlobalMix India #76) at the registry's values: the
# program, its coverage and catch-up age, the maternal vaccine's mother, and
# both products' efficacy profiles
infant_params <- function(program = "nirsevimab", coverage = 0.8, ...) {
  p <- list(vax.infant.program = program, vax.infant.coverage = coverage, vax.infant.catchup.months = 8,
            vax.maternal.mother = TRUE, vax.maternal.lead = 35)
  peaks <- list(nirsevimab = c(0.10, 0.50, 0.80, 0.85, 0), maternal = c(0.30, 0.70, 0.90, 0.90, 0))
  shape <- list(nirsevimab = c(1, 180, 365), maternal = c(0, 0, 270))
  for (pr in names(peaks)) {
    for (k in seq_along(peaks[[pr]])) {
      p[[paste0("vax.", pr, ".ve.peak.", c("infect", "clinical", "hosp", "death", "infectious")[k])]] <- peaks[[pr]][k]
    }
    p[paste0("vax.", pr, ".ve.", c("delay", "plateau", "duration"))] <- as.list(shape[[pr]])
  }
  utils::modifyList(p, list(...))
}

# A vaccine parameter list moved from COVID-19 to RSV
as_rsv <- function(p) {
  names(p) <- sub("^vax\\.covid\\.", "vax.rsv.", names(p))
  p$pathogen <- "rsv"
  p
}
