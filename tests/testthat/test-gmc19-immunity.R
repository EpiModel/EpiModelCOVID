# Episodes, post-infection immunity (GlobalMix India #69), imports and prior
# immunity (#70), and severity per infection with mortality per episode (#71)

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
                ever.sympt = rep(0, N), prior.immune = rep(0, N), vax = rep(0, N),
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
mock_ve <- function(ve_rr = c(infect = 1, clinical = 1, hosp = 1, death = 1), env = parent.frame()) {
  local_mocked_bindings(
    build_vax_schedule = function(dat, disease) NULL,
    compute_ve = function(at, ids, vax, ..., outcome) list(rr = ifelse(vax[ids] >= 1, ve_rr[[outcome]], 1)),
    vax_age_group_for = function(dat) rep(NA_integer_, length(get_attr(dat, "active"))),
    .package = "EpiModelCOVID", .env = env)
}

infected <- function(dat) get_attr(dat, "status") == "e"

test_that("each episode draws its own severity, with the vaccine's effect on a reinfection", {
  mock_ve(c(infect = 1, clinical = 0.5, hosp = 1, death = 1))
  set.seed(1)
  N <- 20000
  dat <- nh_dat(N, prop.clinical = rep(0.6, 12), immune.full = 5)
  dat <- start_episode(dat, seq_len(N), 1)
  dat <- progress_general(dat, 1)
  expect_equal(mean(get_attr(dat, "clinical")), 0.6, tolerance = 0.02)

  # recovered on day 1, back in s on day 6, vaccinated, reinfected on day 7
  dat <- set_attr(dat, "status", rep("r", N))
  dat <- set_attr(dat, "recTime", rep(1, N))
  dat <- set_attr(dat, "hospit", rep(1, N))
  dat <- set_attr(dat, "fatal", rep(0, N))
  dat <- progress_general(dat, 5)
  expect_true(all(get_attr(dat, "status") == "r"))
  dat <- progress_general(dat, 6)
  expect_true(all(get_attr(dat, "status") == "s"))
  dat <- set_attr(dat, "vax", rep(1, N))
  dat <- start_episode(dat, seq_len(N), 7)
  expect_true(all(is.na(get_attr(dat, "clinical"))))
  expect_true(all(is.na(get_attr(dat, "hospit"))))
  expect_true(all(is.na(get_attr(dat, "fatal"))))
  expect_true(all(get_attr(dat, "inf.count") == 2))
  dat <- progress_general(dat, 7)
  # the first episode's draw is not replayed: clinical VE reaches the reinfection
  expect_equal(mean(get_attr(dat, "clinical")), 0.3, tolerance = 0.02)
})

test_that("admission and death are drawn once at onset, per infection by severity band", {
  mock_ve()
  set.seed(2)
  N <- 60000
  pc <- seq(0.3, 0.8, length.out = 12)
  ihr <- pc * seq(0.02, 0.4, length.out = 12)
  ifr <- pc * seq(0.001, 0.1, length.out = 12)
  dat <- nh_dat(N, age = ages12, prop.clinical = pc, ihr = ihr, ifr = ifr)
  dat <- set_attr(dat, "status", rep("ic", N))
  dat <- set_attr(dat, "statusTime", rep(1, N))
  dat <- set_attr(dat, "clinical", rep(1, N))
  dat <- set_attr(dat, "infTime", rep(1, N))
  dat <- progress_general(dat, 1)
  band <- rep_len(1:12, N)
  expect_equal(as.numeric(tapply(get_attr(dat, "hospit"), band, mean)), ihr / pc, tolerance = 0.05)
  expect_equal(as.numeric(tapply(get_attr(dat, "fatal"), band, mean)), ifr / pc, tolerance = 0.1)
  # drawn once: a later step leaves both draws alone
  h <- get_attr(dat, "hospit"); f <- get_attr(dat, "fatal")
  dat <- progress_general(dat, 2)
  expect_identical(get_attr(dat, "fatal"), f)
  expect_identical(get_attr(dat, "hospit"), h)
})

# A cohort infected on day 1, run through progression and deaths with no
# transmission: deaths per infection by band
cohort_ifr <- function(seed, ...) {
  set.seed(seed)
  N <- 48000
  pc <- rep(c(0.4, 0.7), each = 6)
  ifr <- pc * rep(c(0.01, 0.2), each = 6)
  dat <- nh_dat(N, age = ages12, prop.clinical = pc, ihr = pc * 0.3, ifr = ifr, ...)
  dat <- start_episode(dat, seq_len(N), 1)
  reached_r <- logical(N)
  for (at in 1:300) {
    dat <- deaths_covid_gmc19(dat, at)
    dat <- progress_general(dat, at)
    reached_r <- reached_r | get_attr(dat, "status") == "r"
  }
  band <- rep_len(1:12, N)
  list(ifr = as.numeric(tapply(get_attr(dat, "active") == 0, band, mean)), target = ifr,
       fatal_recovered = sum(reached_r & get_attr(dat, "fatal") %in% 1),
       open = sum(get_attr(dat, "active") == 1 & get_attr(dat, "status") %in% c("e", "a", "ip", "ic", "h")))
}

test_that("a fatal episode ends in death, so the IFR per infection is the input whatever the durations", {
  mock_ve()
  short <- cohort_ifr(3)
  long <- cohort_ifr(4, icr.rate = 0.1, hr.rate = 0.05, ich.rate = 0.1)
  expect_identical(short$open, 0L)
  expect_identical(long$open, 0L)
  expect_identical(short$fatal_recovered, 0L)
  expect_equal(short$ifr, short$target, tolerance = 0.12)
  # longer symptomatic and hospital stays used to raise the IFR in proportion
  expect_equal(long$ifr, long$target, tolerance = 0.12)
  expect_equal(mean(long$ifr[7:12]) / mean(short$ifr[7:12]), 1, tolerance = 0.05)
})

test_that("the vaccine's protection against death applies to the episode's death draw", {
  mock_ve(c(infect = 1, clinical = 1, hosp = 1, death = 0.25))
  set.seed(5)
  N <- 40000
  dat <- nh_dat(N, prop.clinical = rep(0.5, 12), ihr = rep(0.1, 12), ifr = rep(0.1, 12))
  dat <- set_attr(dat, "vax", rep(0:1, N / 2))
  dat <- set_attr(dat, "status", rep("ic", N))
  dat <- set_attr(dat, "statusTime", rep(1, N))
  dat <- set_attr(dat, "clinical", rep(1, N))
  dat <- progress_general(dat, 1)
  f <- tapply(get_attr(dat, "fatal"), get_attr(dat, "vax"), mean)
  expect_equal(f[["0"]], 0.2, tolerance = 0.05)
  expect_equal(f[["1"]], 0.05, tolerance = 0.1)
})

test_that("death.hosp.rr moves deaths into hospital and keeps their total", {
  p_death <- rep(0.1, 4); p_hosp <- rep(0.2, 4)
  d <- EpiModelCOVID:::death_given_admission(p_death, p_hosp, c(0, 1, 0, 1), 3)
  expect_equal(d[2] / d[1], 3)
  expect_equal(0.2 * d[2] + 0.8 * d[1], 0.1)
  expect_identical(EpiModelCOVID:::death_given_admission(p_death, p_hosp, c(0, 1, 0, 1), 1), p_death)
})

test_that("full protection lasts immune.full days after recovery, Inf keeps people immune", {
  mock_ve()
  dat <- nh_dat(10, immune.full = 30)
  dat <- set_attr(dat, "status", rep("r", 10))
  dat <- set_attr(dat, "recTime", 1:10)
  back <- integer(10)
  for (at in 2:60) {
    before <- get_attr(dat, "status") == "r"
    dat <- progress_general(dat, at)
    back[before & get_attr(dat, "status") == "s"] <- at
  }
  expect_identical(back - 1:10, rep(30L, 10))

  dat <- nh_dat(10, immune.full = Inf)
  dat <- set_attr(dat, "status", rep("r", 10))
  dat <- set_attr(dat, "recTime", rep(1, 10))
  for (at in 2:400) dat <- progress_general(dat, at)
  expect_true(all(get_attr(dat, "status") == "r"))
})

imm_pars <- list(immune.full = 30,
                 immune.peak.infect = 0.8, immune.halflife.infect = 200, immune.floor.infect = 0.2,
                 immune.peak.clinical = 0.85, immune.halflife.clinical = 200, immune.floor.clinical = 0.3,
                 immune.peak.hosp = 0.95, immune.halflife.hosp = Inf, immune.floor.hosp = 0.95,
                 immune.peak.death = 0.95, immune.halflife.death = Inf, immune.floor.death = 0.95)

test_that("partial protection is marginal, wanes as specified and gives conditional stage risks", {
  dat <- do.call(nh_dat, c(list(N = 3), imm_pars))
  dat <- set_attr(dat, "recTime", c(NA, -99, -399))
  dat <- set_attr(dat, "infTime", c(1, 1, 1))
  d <- c(NA, 100, 400)
  P <- function(peak, hl, fl) ifelse(is.na(d), 0, fl + (peak - fl) * 0.5^(d / hl))
  Pi <- P(0.8, 200, 0.2); Pc <- P(0.85, 200, 0.3); Ph <- P(0.95, Inf, 0.95)
  expect_equal(immune_rr(dat, 1:3, "infect", at = 1), 1 - Pi)
  expect_equal(immune_rr(dat, 1:3, "clinical"), (1 - Pc) / (1 - Pi))
  expect_equal(immune_rr(dat, 1:3, "hosp"), (1 - Ph) / (1 - Pc))
  expect_equal(immune_rr(dat, 1:3, "death"), (1 - Ph) / (1 - Pc))
})

test_that("time from recovery to reinfection follows the protection curve under a constant hazard", {
  mock_ve()
  set.seed(6)
  N <- 30000
  h <- 0.02
  dat <- do.call(nh_dat, c(list(N = N, import.rate = h), imm_pars))
  dat <- set_attr(dat, "status", rep("r", N))
  dat <- set_attr(dat, "recTime", rep(1, N))
  dat <- set_attr(dat, "inf.count", rep(1, N))
  reinf <- rep(NA_integer_, N)
  for (at in 2:150) {
    dat <- infect_general(dat, at)
    new <- is.na(reinf) & infected(dat)
    reinf[new] <- at
    # the reinfected leave the risk set: hold them in e
    dat <- set_attr(dat, "clinical", ifelse(infected(dat), 0, get_attr(dat, "clinical")))
    dat <- set_attr(dat, "statusTime", ifelse(infected(dat), 1e6, get_attr(dat, "statusTime")))
    dat <- progress_general(dat, at)
  }
  # no reinfection within immune.full days, then hazard h (1 - P_infect(d))
  expect_true(all(is.na(reinf) | reinf - 1 > 30))
  d <- 31:149
  haz <- h * (1 - (0.2 + 0.6 * 0.5^(d / 200)))
  surv <- cumprod(1 - haz)
  for (t in c(60, 100, 149)) {
    expect_equal(mean(is.na(reinf) | reinf > t + 1), surv[d == t], tolerance = 0.02)
  }
  expect_equal(sum(get_epi(dat, "se.flow.first"), na.rm = TRUE), 0)
  expect_equal(sum(get_epi(dat, "se.flow.reinf"), na.rm = TRUE), sum(!is.na(reinf)))
})

test_that("imports take the vaccine's and the age band's reduction in susceptibility (#70)", {
  mock_ve(c(infect = 0.4, clinical = 1, hosp = 1, death = 1))
  set.seed(7)
  N <- 80000
  dat <- nh_dat(N, age = c(5, 25), import.rate = 0.1, susc.age.rr = c(0.5, 1, 1, 1, 1, 1))
  vax <- rep(c(0, 0, 1, 1), N / 4)
  dat <- set_attr(dat, "vax", vax)
  dat <- infect_general(dat, 1)
  inf <- infected(dat)
  age <- get_attr(dat, "age")
  rate <- function(v, a) mean(inf[vax == v & age == a])
  expect_equal(rate(0, 25), 0.1, tolerance = 0.04)
  expect_equal(rate(1, 25) / rate(0, 25), 0.4, tolerance = 0.06)
  expect_equal(rate(0, 5) / rate(0, 25), 0.5, tolerance = 0.06)
  expect_equal(rate(1, 5) / rate(0, 25), 0.2, tolerance = 0.08)
  expect_equal(get_epi(dat, "se.import.flow", at = 1), sum(inf))
  expect_equal(get_epi(dat, "se.flow.first", at = 1), sum(inf))
  expect_true(all(get_attr(dat, "inf.count")[inf] == 1))
})

test_that("the introduction pulse brings import.pulse expected introductions in its window", {
  mock_ve()
  set.seed(8)
  dat <- nh_dat(100000, import.pulse = 1000, import.pulse.start = 3, import.pulse.days = 5)
  for (at in 1:10) dat <- infect_general(dat, at)
  imp <- get_epi(dat, "se.import.flow")[1:10]
  expect_equal(imp[c(1:2, 8:10)], rep(0, 5))
  expect_equal(sum(imp), 1000, tolerance = 0.1)
})

test_that("prior immunity holds each band's target share and follows the propensity", {
  set.seed(9)
  N <- 120000
  target <- c(0, 0, 0.2, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6)
  dat <- nh_dat(N, age = ages12, init.immune = target, init.immune.attr = "immune.w")
  w <- rgamma(N, shape = 1)
  dat <- set_attr(dat, "immune.w", w)
  status <- get_attr(dat, "status"); status[1:24] <- "e"
  dat <- set_attr(dat, "status", status)
  dat <- init_immunity(dat)
  band <- rep_len(1:12, N)
  prior <- get_attr(dat, "prior.immune")
  expect_equal(as.numeric(tapply(prior, band, mean)), target, tolerance = 0.03)
  expect_true(all(get_attr(dat, "status")[prior == 1] == "r"))
  expect_true(all(is.na(get_attr(dat, "recTime"))))
  # more likely immune the higher the propensity
  b <- band == 6
  expect_gt(mean(prior[b & w > median(w[b])]), mean(prior[b & w <= median(w[b])]) + 0.2)
  # the seeds are the run's first infections, not immune
  expect_identical(get_attr(dat, "ever.inf")[1:24], rep(1, 24))
  expect_identical(get_attr(dat, "inf.count")[1:24], rep(1, 24))
  expect_identical(prior[1:24], rep(0, 24))
  # weighted inclusion probabilities hold the share exactly and are capped at 1
  p <- EpiModelCOVID:::weighted_inclusion(c(1, 1, 1, 100), 0.5)
  expect_equal(mean(p), 0.5)
  expect_identical(p[4], 1)
  expect_error(EpiModelCOVID:::weighted_inclusion(c(0, 0, 0, 1), 0.5), "positive propensity")
})

test_that("newborns get a default for every propensity the population carries", {
  # The population is drawn once and carries one propensity per pathogen, so a
  # run of one pathogen holds the others too, and the arrivals module needs a
  # default for each: without them a two-pathogen population died at the first
  # step with a birth rather than at initialization
  dat <- nh_dat(5, init.immune.attr = "immune.w.flu")
  dat <- set_attr(dat, "immune.w.flu", rep(1, 5))
  dat <- set_attr(dat, "immune.w.rsv", rep(1, 5))
  expect_setequal(EpiModelCOVID:::immune_weight_attrs(dat), c("immune.w.flu", "immune.w.rsv"))
  d <- EpiModelCOVID:::get_default_attrs(dat)
  expect_identical(d$immune.w.flu, NA_real_)
  expect_identical(d$immune.w.rsv, NA_real_)
  # a run naming none of them still carries what the population holds
  no_attr <- set_attr(nh_dat(5), "immune.w.rsv", rep(1, 5))
  expect_identical(EpiModelCOVID:::immune_weight_attrs(no_attr), "immune.w.rsv")
  expect_identical(EpiModelCOVID:::get_default_attrs(no_attr)$immune.w.rsv, NA_real_)
})

test_that("the severity, immunity, import and prior-immunity parameters are checked", {
  ok <- nh_dat(5)
  expect_true(EpiModelCOVID:::check_severity_params(ok))
  expect_true(EpiModelCOVID:::check_immunity_params(ok))
  expect_true(EpiModelCOVID:::check_import_params(ok))
  expect_true(EpiModelCOVID:::check_prior_immunity(ok))
  expect_error(EpiModelCOVID:::check_severity_params(nh_dat(5, ihr = rep(0.6, 12))),
               "cannot exceed prop.clinical")
  expect_error(EpiModelCOVID:::check_severity_params(nh_dat(5, ifr = rep(0.1, 6))),
               "one probability per severity band")
  expect_error(EpiModelCOVID:::check_severity_params(nh_dat(5, sev.age.breaks = c(0, 5, 60))),
               "closed by Inf")
  expect_error(EpiModelCOVID:::check_immunity_params(nh_dat(5, immune.full = 0)), "at least 1")
  expect_error(EpiModelCOVID:::check_immunity_params(
    nh_dat(5, immune.peak.infect = 0.9, immune.peak.clinical = 0.5)), "must not fall")
  expect_error(EpiModelCOVID:::check_import_params(nh_dat(5, import.pulse = -1)), "non-negative")
  expect_error(EpiModelCOVID:::check_prior_immunity(nh_dat(5, init.immune = rep(0.5, 6))),
               "one share in \\[0, 1\\] per severity band")
  expect_error(EpiModelCOVID:::check_prior_immunity(nh_dat(5, init.immune = rep(0.5, 12),
                                                           init.immune.attr = "nope")),
               "does not carry")
})

test_that("the parameters the redo replaced stop the run", {
  for (p in c("rs.rate", "prop.hospit", "dis.death.rate", "dis.death.rate.infant")) {
    dat <- set_param(nh_dat(5), p, 0.1)
    expect_error(EpiModelCOVID:::check_retired_params(dat), paste0("no longer reads ", p))
  }
})
