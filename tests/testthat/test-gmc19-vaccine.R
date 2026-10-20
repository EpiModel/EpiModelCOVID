# Vaccine efficacy (GlobalMix India #73): marginal VE inputs converted to the
# conditional stage relative risks, the onset delay, age-specific VE, protection
# against infectiousness, and the schedule check

# N people, the first half dosed `since` days before day 1
dosed_dat <- function(N, since = 30, age = 30, ...) {
  dat <- do.call(nh_dat, c(list(N = N, age = age), ve_params(...)))
  dat <- set_attr(dat, "vax", rep(1:0, each = N / 2))
  set_attr(dat, "last.dose.time", rep(c(1 - since, NA), each = N / 2))
}

ve_at <- function(dat, at, outcome, ids = seq_len(length(get_attr(dat, "active")))) {
  compute_ve(at = at, ids = ids, vax = get_attr(dat, "vax"), last.dose.time = get_attr(dat, "last.dose.time"),
             vax.schedule = build_vax_schedule(dat, "covid"), outcome = outcome,
             vax.age.group = vax_age_group_for(dat))
}

test_that("a VE field holds one value per dose or one per vaccine age group, and the schedule has one row per dose", {
  dat <- do.call(nh_dat, c(list(N = 5, age = c(2, 10, 30, 55, 70)),
                           ve_params(vax.covid.ve.peak.infect.dose1 = c(0.6, 0.5, 0.5, 0.4, 0.3))))
  vs <- build_vax_schedule(dat, "covid")
  expect_identical(nrow(vs), 1L)
  expect_identical(vs$ve.peak.infect[[1]], c(0.6, 0.5, 0.5, 0.4, 0.3))
  expect_identical(vs$ve.peak.clinical[[1]], 0.5)
  ag <- vax_age_group_for(dat)
  expect_identical(ag, 1:5)
  expect_identical(EpiModelCOVID:::get_schedule_value(vs, "ve.peak.infect", 1, ag), c(0.6, 0.5, 0.5, 0.4, 0.3))
  expect_identical(EpiModelCOVID:::get_schedule_value(vs, "ve.peak.clinical", 1, ag), rep(0.5, 5))
  # two doses, two rows, each with its own vector
  dat2 <- do.call(nh_dat, c(list(N = 5), ve_params(n.doses = 2, vax.covid.ve.peak.death.dose2 = rep(0.95, 5))))
  vs2 <- build_vax_schedule(dat2, "covid")
  expect_identical(nrow(vs2), 2L)
  expect_identical(vs2$ve.peak.death[[2]], rep(0.95, 5))
  expect_identical(vs2$interval[[2]], 21)
  # the wrong length stops the run
  bad <- do.call(nh_dat, c(list(N = 5), ve_params(vax.covid.ve.peak.infect.dose1 = c(0.6, 0.5, 0.4))))
  expect_error(build_vax_schedule(bad, "covid"), "one, or one per vaccine age group, 5")
  bad <- do.call(nh_dat, c(list(N = 5), ve_params(vax.covid.ve.delay.dose1 = c(14, 14))))
  expect_error(build_vax_schedule(bad, "covid"), "needs one")
})

test_that("VE is 0 before the onset delay, the peak at it, and wanes toward the floor from there", {
  dat <- dosed_dat(4, since = 0, delay = 14, halflife = 100, floor = 0.1, infect = 0.6)
  ve <- function(at) ve_at(dat, at, "infect")$ve
  expect_identical(ve(1)[1], 0)              # dosed today
  expect_identical(ve(14)[1], 0)             # the day before the delay ends
  expect_equal(ve(15)[1], 0.6)               # 14 days after the dose
  expect_equal(ve(15 + 100)[1], 0.1 + 0.5 * 0.5)
  expect_equal(ve(15 + 300)[1], 0.1 + 0.5 * 0.125)
  expect_identical(ve(15)[3], 0)             # unvaccinated
  expect_equal(ve_at(dat, 15, "infect")$time.since.last.dose, c(14, 14, 0, 0))
  # a passive product (delay 0) protects from its dose
  dat0 <- dosed_dat(2, since = 0, delay = 0, infect = 0.8)
  expect_equal(ve_at(dat0, 1, "infect")$ve[1], 0.8)
})

test_that("the stage relative risks reproduce the marginal VE of the schedule, at any time and by age group", {
  dat <- do.call(nh_dat, c(list(N = 10, age = rep(c(2, 10, 30, 55, 70), 2)),
                           ve_params(halflife = 150, floor = 0.05,
                                     vax.covid.ve.peak.infect.dose1 = c(0.5, 0.4, 0.3, 0.3, 0.2),
                                     vax.covid.ve.peak.clinical.dose1 = c(0.6, 0.6, 0.5, 0.5, 0.4),
                                     vax.covid.ve.peak.hosp.dose1 = c(0.9, 0.9, 0.8, 0.8, 0.7),
                                     vax.covid.ve.peak.death.dose1 = c(0.95, 0.95, 0.9, 0.9, 0.85))))
  dat <- set_attr(dat, "vax", rep(1:0, each = 5))
  dat <- set_attr(dat, "last.dose.time", rep(c(-13, NA), each = 5))   # day 1 is 14 days after the dose
  for (at in c(1, 60, 200, 800)) {
    V <- function(o) ve_at(dat, at, o)$ve
    R <- function(o) ve_at(dat, at, o)$rr
    # the marginal VE of each stage is recovered from the chain of stage risks
    expect_equal((1 - V("infect")) * R("clinical"), 1 - V("clinical"), tolerance = 1e-12)
    expect_equal((1 - V("clinical")) * R("hosp"), 1 - V("hosp"), tolerance = 1e-12)
    expect_equal((1 - V("clinical")) * R("death"), 1 - V("death"), tolerance = 1e-12)
    expect_equal(R("infect"), 1 - V("infect"))
    expect_equal(R("infectious"), 1 - V("infectious"))
    expect_true(all(R("clinical") <= 1 & R("hosp") <= 1 & R("death") <= 1))
    # the unvaccinated carry no protection and a relative risk of 1
    expect_identical(V("death")[6:10], rep(0, 5))
    expect_identical(R("death")[6:10], rep(1, 5))
  }
  # at the peak, by age group
  V1 <- ve_at(dat, 1, "infect")$ve
  expect_equal(V1[1:5], c(0.5, 0.4, 0.3, 0.3, 0.2))
  expect_equal(ve_at(dat, 1, "clinical")$rr[1:5], (1 - c(0.6, 0.6, 0.5, 0.5, 0.4)) / (1 - c(0.5, 0.4, 0.3, 0.3, 0.2)))
})

test_that("the schedule check stops VE that falls along the stages, and bad profiles", {
  expect_true(EpiModelCOVID:::check_vax_schedule(do.call(nh_dat, c(list(N = 2), ve_params()))))
  expect_error(EpiModelCOVID:::check_vax_schedule(do.call(nh_dat, c(list(N = 2), ve_params(infect = 0.6, clinical = 0.5)))),
               "must not fall along the stages")
  expect_error(EpiModelCOVID:::check_vax_schedule(do.call(nh_dat, c(list(N = 2), ve_params(
    vax.covid.ve.peak.hosp.dose1 = c(0.9, 0.9, 0.4, 0.9, 0.9))))), "vaccine age group\\(s\\) 3")
  expect_error(EpiModelCOVID:::check_vax_schedule(do.call(nh_dat, c(list(N = 2), ve_params(
    vax.covid.ve.floor.death.dose1 = 0.2, vax.covid.ve.floor.clinical.dose1 = 0.3)))), "ve.floor")
  expect_error(EpiModelCOVID:::check_vax_schedule(do.call(nh_dat, c(list(N = 2), ve_params(halflife = 0)))), "positive")
  expect_error(EpiModelCOVID:::check_vax_schedule(do.call(nh_dat, c(list(N = 2), ve_params(death = 1.2)))), "in \\[0, 1\\]")
  expect_error(EpiModelCOVID:::check_vax_schedule(do.call(nh_dat, c(list(N = 2), ve_params(delay = -1)))), "non-negative")
  # a run with no schedule for its pathogen has nothing to check
  expect_true(EpiModelCOVID:::check_vax_schedule(nh_dat(2)))
})

test_that("realized protection reproduces the marginal VE: infection under a constant hazard, then each stage among the infected", {
  set.seed(21)
  N <- 60000
  h <- 0.01
  dat <- dosed_dat(N, since = 30, halflife = Inf, import.rate = h, prop.clinical = rep(0.8, 12), ihr = rep(0.4, 12),
                   ifr = rep(0.24, 12))
  vax <- get_attr(dat, "vax")
  # infection: 20 days of a constant import hazard h, h (1 - VE_infect) for the dosed
  for (at in 1:20) {
    dat <- infect_general(dat, at)
    dat <- set_attr(dat, "statusTime", ifelse(get_attr(dat, "status") == "e", 1e6, get_attr(dat, "statusTime")))
  }
  inf <- get_attr(dat, "status") == "e"
  p <- tapply(inf, vax, mean)
  expect_equal(1 - log(1 - p[["1"]]) / log(1 - p[["0"]]), 0.3, tolerance = 0.05)

  # the stages: everyone infected on day 1, progression and deaths with no transmission
  dat <- dosed_dat(N, since = 30, halflife = Inf, prop.clinical = rep(0.8, 12), ihr = rep(0.4, 12), ifr = rep(0.24, 12),
                   vax.diagnostics = TRUE)
  dat <- start_episode(dat, seq_len(N), 1)
  for (at in 1:200) {
    dat <- deaths_covid_gmc19(dat, at)
    dat <- progress_general(dat, at)
  }
  # the vaccine diagnostics count the same events by vaccine age group (3, ages 30) and status
  vd <- function(what, v) sum(get_epi(dat, paste0("vd.", what, ".3.", v)), na.rm = TRUE)
  expect_identical(vd("sympt", 1), sum(get_attr(dat, "clinical")[vax == 1]))
  expect_identical(vd("sympt", 0), sum(get_attr(dat, "clinical")[vax == 0]))
  expect_identical(vd("hosp", 1), sum(get_attr(dat, "hospit")[vax == 1] %in% 1))
  expect_identical(vd("death", 1) + vd("death", 0), sum(get_attr(dat, "active") == 0))
  sympt <- tapply(get_attr(dat, "clinical"), vax, mean)
  hosp <- tapply(get_attr(dat, "hospit") %in% 1, vax, mean)
  died <- tapply(get_attr(dat, "active") == 0, vax, mean)
  # among the infected, the conditional protection; the marginal VE follows by
  # multiplying the infection stage back in
  expect_equal(sympt[["0"]], 0.8, tolerance = 0.02)
  expect_equal(1 - (1 - 0.3) * sympt[["1"]] / sympt[["0"]], 0.5, tolerance = 0.02)
  expect_equal(hosp[["0"]], 0.4, tolerance = 0.03)
  expect_equal(1 - (1 - 0.3) * hosp[["1"]] / hosp[["0"]], 0.8, tolerance = 0.02)
  expect_equal(died[["0"]], 0.24, tolerance = 0.04)
  expect_equal(1 - (1 - 0.3) * died[["1"]] / died[["0"]], 0.9, tolerance = 0.02)
})

test_that("a vaccinated infector transmits at 1 - VE_infectious, on contacts and in the background hazard", {
  set.seed(22)
  M <- 40000
  dat <- do.call(nh_dat, c(list(N = 2 * M), ve_params(infectious = 0.5, infect = 0, halflife = Inf, vax.diagnostics = TRUE)))
  dat <- set_param(dat, "inf.prob", 0.2)
  # infectors 1..M (the first half dosed), each tied to one susceptible
  dat <- set_attr(dat, "vax", c(rep(1:0, each = M / 2), rep(0, M)))
  dat <- set_attr(dat, "last.dose.time", c(rep(-30, M / 2), rep(NA, M + M / 2)))
  dat <- set_attr(dat, "status", c(rep("ip", M), rep("s", M)))
  dat <- set_attr(dat, "clinical", c(rep(1, M), rep(NA, M)))
  dat <- set_attr(dat, "statusTime", c(rep(1, M), rep(NA, M)))
  dat <- set_attr(dat, "infTime", c(rep(1, M), rep(NA, M)))
  dat$run$el <- list(community = cbind(seq_len(M), M + seq_len(M)))
  dat$control$tergmLite <- TRUE   # get_edgelist() reads dat$run$el
  dat <- infect_general(dat, 2)
  inf <- get_attr(dat, "status")[M + seq_len(M)] == "e"
  expect_identical(get_epi(dat, "vd.inf.3.0", at = 2), sum(inf))
  expect_identical(get_epi(dat, "vd.inf.3.1", at = 2), 0L)
  expect_equal(mean(inf[(M / 2 + 1):M]), 0.2, tolerance = 0.03)
  # the exposure accounting: every susceptible has one discordant contact, and
  # the infections they imply are the pair probabilities (0.2 from an undosed
  # infector, 0.1 from a dosed one), with no vaccine of their own to take out
  expect_identical(get_epi(dat, "vd.exp.3.0", at = 2), M)
  expect_equal(get_epi(dat, "vd.einf.3.0", at = 2), 0.1 * M / 2 + 0.2 * M / 2)
  expect_equal(get_epi(dat, "vd.einf0.3.0", at = 2), 0.1 * M / 2 + 0.2 * M / 2)
  # with the susceptibles dosed instead (VE against infection 0.5), the implied
  # infections halve and einf0 gives them back
  dat2 <- do.call(nh_dat, c(list(N = 2 * M), ve_params(infectious = 0, infect = 0.5, halflife = Inf, vax.diagnostics = TRUE)))
  dat2 <- set_param(dat2, "inf.prob", 0.2)
  dat2 <- set_attr(dat2, "vax", c(rep(0, M), rep(1, M)))
  dat2 <- set_attr(dat2, "last.dose.time", c(rep(NA, M), rep(-30, M)))
  dat2 <- set_attr(dat2, "status", c(rep("ip", M), rep("s", M)))
  dat2 <- set_attr(dat2, "clinical", c(rep(1, M), rep(NA, M)))
  dat2 <- set_attr(dat2, "statusTime", c(rep(1, M), rep(NA, M)))
  dat2 <- set_attr(dat2, "infTime", c(rep(1, M), rep(NA, M)))
  dat2$run$el <- list(community = cbind(seq_len(M), M + seq_len(M)))
  dat2$control$tergmLite <- TRUE
  dat2 <- infect_general(dat2, 2)
  expect_identical(get_epi(dat2, "vd.exp.3.1", at = 2), M)
  expect_equal(get_epi(dat2, "vd.einf.3.1", at = 2), 0.1 * M)
  expect_equal(get_epi(dat2, "vd.einf0.3.1", at = 2), 0.2 * M)
  expect_equal(get_epi(dat2, "vd.inf.3.1", at = 2) / M, 0.1, tolerance = 0.05)
  expect_equal(mean(inf[1:(M / 2)]) / mean(inf[(M / 2 + 1):M]), 0.5, tolerance = 0.05)

  # the background hazard weights an infector's prevalence the same way
  dat <- do.call(nh_dat, c(list(N = 2 * M), ve_params(infectious = 0.5, infect = 0, halflife = Inf)))
  dat <- set_param(dat, "inf.prob", 0)
  dat <- set_param(dat, "bg.hazard", 0.5)
  dat <- set_param(dat, "bg.mix", matrix(1, 6, 6))
  dat <- set_param(dat, "age.grps", 1:6)
  dat <- set_attr(dat, "age.grp", rep(3L, 2 * M))
  dat$control$tergmLite <- TRUE
  dat <- set_attr(dat, "status", c(rep("ip", M), rep("s", M)))
  dat <- set_attr(dat, "clinical", c(rep(1, M), rep(NA, M)))
  dat <- set_attr(dat, "statusTime", c(rep(1, M), rep(NA, M)))
  dat <- set_attr(dat, "infTime", c(rep(1, M), rep(NA, M)))
  one <- function(dosed) {
    d <- set_attr(dat, "vax", c(rep(dosed, M), rep(0, M)))
    d <- set_attr(d, "last.dose.time", c(rep(if (dosed) -30 else NA, M), rep(NA, M)))
    d <- infect_general(d, 2)
    mean(get_attr(d, "status")[M + seq_len(M)] == "e")
  }
  # everyone is in band 3, where half are infectious at weight 1 (prevalence
  # 0.5) and the within-band mixing weight is 2/7 (the diagonal is doubled
  # before row normalization); dosed infectors count at half the weight
  expect_equal(-log(1 - one(0)), 0.5 * 2 / 7 * 0.5, tolerance = 0.05)
  expect_equal(-log(1 - one(1)) / -log(1 - one(0)), 0.5, tolerance = 0.1)
})
