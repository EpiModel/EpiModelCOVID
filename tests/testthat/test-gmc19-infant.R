# The RSV infant products, their dose line and their efficacy (GlobalMix India
# #76). The cocoon ordering is tested with the other strategies
# (test-gmc19-allocation.R), and the products on a full small run with the
# instrumentation (test-gmc19-instrumentation.R).

# N people on an RSV run with the infant program as given, no adult doses
# unless a budget is set, and coverage draws u (spread over (0, 1) by default)
infant_dat <- function(N, age, program = "nirsevimab", coverage = 0.8, start = 1, u = NULL, budget = 0, ...) {
  vp <- as_rsv(ve_params(vax.budget = budget, vax.delivery.days = 1, vax.start = start, vax.strategy = "age_oldest"))
  vp$vax.rsv.age.min <- 60
  dat <- do.call(nh_dat, c(list(N = N, age = age), vp, infant_params(program, coverage, ...)))
  for (a in c("iprog.time", "mother.uid")) dat <- set_attr(dat, a, rep(NA_real_, N))
  dat <- set_attr(dat, "iprog.u", if (is.null(u)) (seq_len(N) * 0.6180339887) %% 1 else u)
  dat <- set_attr(dat, "degree_total", rep(1, N))
  set_attr(dat, "n_layers_active", rep(1L, N))
}
run_vax <- function(dat, steps) {
  for (at in steps) dat <- vax_general(dat, at)
  dat
}
held <- function(dat) which(!is.na(get_attr(dat, "iprog.time")))

test_that("the passive profile is flat to the plateau, falls on a straight line and is 0 from the stated duration", {
  dat <- infant_dat(5, 0.5)
  nir <- EpiModelCOVID:::infant_prog_pars(dat, "nirsevimab")
  mat <- EpiModelCOVID:::infant_prog_pars(dat, "maternal")
  pv <- EpiModelCOVID:::passive_ve
  # nirsevimab: nothing on the day of the dose, 0.80 against admission from the
  # next day to day 180, halfway down at day 272.5, 0 from day 365
  expect_equal(pv(c(0, 1, 90, 180, 272.5, 365, 500), nir, "hosp"), c(0, 0.8, 0.8, 0.8, 0.4, 0, 0))
  # maternal: the peak at birth, falling to 0 at day 270
  expect_equal(pv(c(0, 135, 270, 400), mat, "death"), c(0.9, 0.45, 0, 0))
  t <- 0:500
  for (p in list(nir, mat)) {
    v <- sapply(c("infect", "clinical", "hosp", "death", "infectious"), function(o) pv(t, p, o))
    expect_true(all(v[, "infect"] <= v[, "clinical"] & v[, "clinical"] <= v[, "hosp"] &
                      v[, "clinical"] <= v[, "death"]))
    expect_true(all(v[t >= p$duration, ] == 0))
    expect_true(all(v[, "infectious"] == 0))
  }
})

test_that("compute_ve() protects the infants who hold the product, at every stage", {
  N <- 4
  dat <- infant_dat(N, age = 0.3)
  dat <- set_attr(dat, "iprog.time", c(10, 10, NA, -400))   # dosed on day 10; none; dosed long ago
  st <- EpiModelCOVID:::infant_prog_state(dat)
  rr <- function(o) compute_ve(at = 20, ids = 1:N, vax = rep(0, N), vax.age.group = rep(1L, N),
                               last.dose.time = rep(NA, N), vax.schedule = NULL, outcome = o, iprog = st)$rr
  expect_equal(rr("infect"), c(0.9, 0.9, 1, 1))
  expect_equal(rr("clinical"), c(0.5 / 0.9, 0.5 / 0.9, 1, 1))
  expect_equal(rr("hosp"), c(0.4, 0.4, 1, 1))
  expect_equal(rr("death"), c(0.3, 0.3, 1, 1))
  expect_equal(rr("infectious"), rep(1, N))
  # the susceptibility every contact, import and the background hazard take
  expect_equal(EpiModelCOVID:::susceptibility_rr(dat, 1:N, 20), c(0.9, 0.9, 1, 1))
  # and no protection on a run without the program
  expect_equal(EpiModelCOVID:::susceptibility_rr(set_param(dat, "vax.infant.program", "none"), 1:N, 20), rep(1, N))
})

test_that("the nirsevimab catch-up goes once, on the day the campaign opens, to the covered infants under 8 months", {
  N <- 300
  a <- c(seq(0.01, 0.99, length.out = 100), rep(30, 200))
  u <- (seq_len(N) * 0.6180339887) %% 1
  dat <- run_vax(infant_dat(N, a, start = 5, u = u), 1:10)
  want <- which(a * 12 < 8 & u < 0.8)
  expect_identical(held(dat), want)
  expect_true(all(get_attr(dat, "iprog.time")[want] == 5))
  expect_equal(get_epi(dat, "nVax.iprog.catchup")[1:10], replace(rep(0, 10), 5, length(want)))
  expect_equal(get_epi(dat, "nVax.iprog")[5], length(want))
  expect_equal(get_epi(dat, "n.iprog.infant")[c(4, 10)], c(0, length(want)))
  # at 50% the infants covered are among those covered at 80%
  d50 <- run_vax(infant_dat(N, a, start = 5, u = u, coverage = 0.5), 1:6)
  expect_identical(held(d50), which(a * 12 < 8 & u < 0.5))
  expect_true(all(held(d50) %in% want))
  # the adult product is untouched, and the maternal program has no catch-up
  expect_true(all(get_attr(dat, "vax") == 0))
  dm <- run_vax(infant_dat(N, a, program = "maternal", start = 5, u = u), 1:6)
  expect_length(held(dm), 0)
})

test_that("birth doses: nirsevimab from the open day, the maternal vaccine every day with the mother on record", {
  N <- 50
  a <- c(rep(0, 10), rep(28, 20), rep(70, 20))     # 1-10 newborns, 11-30 women of 28
  u <- rep(c(0.1, 0.9), 25)                          # the odd ids are covered
  nb <- 1:10
  covered <- c(1L, 3L, 5L, 7L, 9L)
  births <- EpiModelCOVID:::infant_program_births
  dat <- infant_dat(N, a, start = 5, u = u)
  d <- births(dat, 3, nb)
  expect_length(held(d), 0)
  expect_equal(get_epi(d, "nVax.iprog.birth", at = 3), 0)
  d <- births(dat, 5, nb)
  expect_identical(held(d), covered)
  expect_equal(get_epi(d, "nVax.iprog.birth", at = 5), 5)
  expect_equal(get_epi(d, "nVax.iprog.mother", at = 5), 0)
  # maternal: from the first day; child k's mother is node 10 + k, and of the
  # mothers of the covered, node 15 is already dosed and node 17 has died
  dat <- infant_dat(N, a, program = "maternal", start = 50, u = u)
  mu <- rep(NA_real_, N)
  mu[nb] <- get_attr(dat, "unique_id")[10 + nb]
  dat <- set_attr(dat, "mother.uid", mu)
  dat <- set_attr(dat, "vax", replace(rep(0, N), 15, 1))
  dat <- set_attr(dat, "active", replace(rep(1, N), 17, 0))
  d <- births(dat, 3, nb)
  expect_identical(held(d), covered)
  expect_true(all(get_attr(d, "iprog.time")[covered] == 3))
  moms <- c(11, 13, 19)
  expect_identical(which(get_attr(d, "vax") == 1), c(11L, 13L, 15L, 19L))
  expect_equal(get_attr(d, "last.dose.time")[moms], rep(3 - 35, 3))
  expect_equal(get_epi(d, "nVax.iprog.birth", at = 3), 5)
  expect_equal(get_epi(d, "nVax.iprog.mother", at = 3), 3)
  d <- births(set_param(dat, "vax.maternal.mother", FALSE), 3, nb)
  expect_identical(which(get_attr(d, "vax") == 1), 15L)
})

test_that("the maternal program protects the covered starting infants from their birth", {
  N <- 40
  a <- c(seq(0.05, 0.95, length.out = 10), rep(28, 10), rep(60, 20))
  u <- rep(c(0.1, 0.9), 20)
  dat <- infant_dat(N, a, program = "maternal", u = u)
  mu <- rep(NA_real_, N)
  mu[1:10] <- get_attr(dat, "unique_id")[11:20]
  d <- EpiModelCOVID:::init_infant_program(set_attr(dat, "mother.uid", mu))
  cov <- c(1L, 3L, 5L, 7L, 9L)
  birth <- 1 - round(a[cov] * 365)
  expect_identical(held(d), cov)
  expect_equal(get_attr(d, "iprog.time")[cov], birth)
  expect_equal(get_attr(d, "last.dose.time")[10 + cov], birth - 35)
  expect_equal(get_epi(d, "iprog.init", at = 1), 5)
  expect_equal(get_epi(d, "iprog.init.mother", at = 1), 5)
  # on day 2 the protection follows the infant's age
  st <- EpiModelCOVID:::infant_prog_state(d)
  ve <- EpiModelCOVID:::vax_protection(2, cov, rep(0, N), rep(NA, N), NULL, rep(1L, N), "death", st)
  expect_equal(ve, 0.9 * pmax(0, (270 - (2 - birth)) / 270))
  # nirsevimab gives nothing at initialization
  expect_length(held(EpiModelCOVID:::init_infant_program(infant_dat(N, a, u = u))), 0)
})

test_that("check_infant_program() stops a bad program", {
  check <- EpiModelCOVID:::check_infant_program
  ok <- infant_dat(5, 0.5)
  expect_true(check(ok))
  expect_true(check(set_param(ok, "vax.infant.program", "none")))
  expect_error(check(set_param(ok, "vax.infant.program", "palivizumab")), "vax.infant.program must be")
  expect_error(check(set_param(ok, "pathogen", "covid")), "is an RSV product")
  expect_error(check(set_param(ok, "vax.infant.coverage", 1.5)), "vax.infant.coverage must be")
  expect_error(check(set_param(ok, "vax.infant.catchup.months", NA)), "vax.infant.catchup.months must be")
  expect_error(check(set_param(ok, "vax.nirsevimab.ve.peak.clinical", 0.9)), "must not fall along the stages")
  expect_error(check(set_param(ok, "vax.nirsevimab.ve.duration", 100)), "ve.delay <= ve.plateau < ve.duration")
  m <- set_param(ok, "vax.infant.program", "maternal")
  expect_true(check(m))
  expect_error(check(set_param(m, "vax.maternal.lead", -1)), "vax.maternal.lead must be")
  expect_error(check(set_param(m, "vax.maternal.mother", NA)), "vax.maternal.mother must be")
})
