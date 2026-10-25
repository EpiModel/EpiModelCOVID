# Instrumentation for the redo estimands (GlobalMix India #77): decedent-inclusive
# person counts and the day-1 cohort, the household SAR on the start-of-step
# state, dose accounting by band and state, years of life lost from exact age at
# death, eta-squared of degree on age, the baseline degree quartile, and the
# partition identities of the new flows on a full small run

counts <- EpiModelCOVID:::departed_counts

test_that("the person counts do not fall when an ever-infected person dies", {
  mock_ve()
  N <- 240
  ages <- rep(c(5, 30, 70), length.out = N)   # bands 1, 4, 6
  dat <- nh_dat(N, age = ages)
  # everyone 61 and over dies on the next departures step
  mr <- rep(0, 101); mr[61:101] <- 1
  dat <- set_param(dat, "mort.rates", mr)
  inf <- c(which(ages == 5)[1:20], which(ages == 30)[1:20], which(ages == 70)[1:20])
  sympt <- inf[seq(1, 60, by = 2)]
  dat <- start_episode(dat, inf, 1)
  dat <- set_attr(dat, "ever.inf", 1, inf)
  dat <- set_attr(dat, "ever.sympt", 1, sympt)
  dat <- set_attr(dat, "status", "r", inf)   # recovered, so they can die of background mortality only

  dat <- prevalence_covid_corporate(dat, 2)
  expect_equal(get_epi(dat, "ever.inf.pers")[2], 60)
  expect_equal(get_epi(dat, "ever.inf.dead")[2], 0)
  expect_equal(get_epi(dat, "ever.sympt.pers")[2], 30)
  expect_equal(get_epi(dat, "n.cohort")[2], N)
  expect_equal(get_epi(dat, "ever.inf.cohort.pers")[2], 60)

  set.seed(1)
  dat <- deaths_covid_gmc19(dat, 3)
  expect_equal(get_epi(dat, "d.flow")[3], sum(ages == 70))
  expect_identical(counts(dat)$inf, 20)
  dat <- prevalence_covid_corporate(dat, 3)
  # the living-only count falls by the 20 infected 70-year-olds; the person count does not
  expect_equal(get_epi(dat, "ever.inf.num")[3], 40)
  expect_equal(get_epi(dat, "ever.inf.dead")[3], 20)
  expect_equal(get_epi(dat, "ever.inf.pers")[3], 60)
  expect_equal(get_epi(dat, "ever.inf.pers.age6")[3], 20)
  expect_equal(get_epi(dat, "ever.inf.dead.age6")[3], 20)
  expect_equal(get_epi(dat, "ever.inf.pers.age1")[3], 20)
  expect_equal(get_epi(dat, "ever.inf.age6")[3], 0)
  expect_equal(get_epi(dat, "ever.sympt.num")[3], 20)
  expect_equal(get_epi(dat, "ever.sympt.dead")[3], 10)
  expect_equal(get_epi(dat, "ever.sympt.pers")[3], 30)
  expect_equal(get_epi(dat, "ever.sympt.pers.age6")[3], 10)
  # the cohort: everyone was present on day 1
  expect_equal(get_epi(dat, "n.cohort")[3], N - sum(ages == 70))
  expect_equal(get_epi(dat, "d.cohort")[3], sum(ages == 70))
  expect_equal(get_epi(dat, "ever.inf.cohort")[3], 40)
  expect_equal(get_epi(dat, "ever.inf.cohort.pers")[3], 60)
  # a later step with no deaths leaves the departed counts where they are
  dat <- set_param(dat, "mort.rates", rep(0, 101))
  dat <- deaths_covid_gmc19(dat, 4)
  dat <- prevalence_covid_corporate(dat, 4)
  expect_equal(get_epi(dat, "ever.inf.pers")[4], 60)
  expect_equal(get_epi(dat, "ever.inf.dead")[4], 20)
})

test_that("person-time is emitted by severity band and vaccine age group", {
  mock_ve()
  N <- 120
  ages <- rep(c(0.2, 0.7, 3, 7, 15, 25, 35, 45, 55, 65, 75, 85), each = 10)
  dat <- nh_dat(N, age = ages, vax.age.breaks = c(0, 5, 18, 50, 65, Inf))
  dat <- set_attr(dat, "active", 0, posit_ids = 1:5)
  dat <- prevalence_covid_corporate(dat, 2)
  sev <- vapply(1:12, function(k) get_epi(dat, paste0("n.sev", k))[2], numeric(1))
  expect_equal(sev, c(5, rep(10, 11)))
  vax <- vapply(1:5, function(k) get_epi(dat, paste0("n.vaxage", k))[2], numeric(1))
  expect_equal(vax, c(25, 20, 30, 10, 30))
  expect_equal(sum(sev), get_epi(dat, "num")[2])
})

# A tergmLite netsim_dat object at t = 1 on a community layer and a household
# clique layer (15 households of 4), as the modifiers tests build it, with the
# statuses of the first eight people set for the household SAR hand computation:
#   household 1: 1 symptomatic, 2 and 3 susceptible, 4 admitted
#   household 2: 5 admitted, 6 susceptible
#   household 3: 7 susceptible, 8 recovered
# Everyone else is recovered, so nothing else is exposed or infectious.
hh_dat <- function(inf.prob = c(0, 1)) {
  N <- 60
  hh <- rep(1:15, each = 4)
  nw <- network_initialize(n = N)
  est_com <- suppressMessages(netest(nw, formation = ~edges, target.stats = 40,
                                     coef.diss = dissolution_coefs(~offset(edges), 20), verbose = FALSE))
  nw_hh <- network_initialize(n = N)
  nw_hh <- set_vertex_attribute(nw_hh, "hh.ids", hh)
  x <- list(community = est_com, home = netclique(nw_hh, "hh.ids"))
  param <- param.net(inf.prob = 0.1, act.rate = 1)
  init <- init.net(i.num = 5)
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE)
  EpiModel:::crosscheck.net(x, param, init, control)
  control <- EpiModel:::netsim_validate_control(control)
  dat <- initialize.net(x, param, init, control, s = 1)
  dat <- set_param(dat, "layer.names", names(x))
  dat <- set_param(dat, "inf.prob", inf.prob)
  dat <- set_param(dat, "exposure.rate", c(1, 1))
  dat <- set_param(dat, "inf.prob.a.rr", 0.5)
  dat <- set_param(dat, "age.breaks", c(0, 10, 20, 30, 40, 60, Inf))
  dat <- set_param(dat, "arrival.age", 0)
  status <- rep("r", N)
  status[c(1, 2, 3, 4, 5, 6, 7, 8)] <- c("ic", "s", "s", "h", "h", "s", "s", "r")
  dat <- set_attr(dat, "status", status)
  for (a in c("vax", "ever.inf", "inf.count")) dat <- set_attr(dat, a, rep(0, N))
  dat <- set_attr(dat, "last.dose.time", rep(NA, N))
  dat <- set_attr(dat, "age", rep(30, N))
  dat <- set_attr(dat, "is_infant", rep(FALSE, N))
  dat <- set_attr(dat, "degree_quartile", rep(1, N))
  dat <- set_attr(dat, "is_bridge", rep(FALSE, N))
  EpiModelCOVID:::check_layer_params(dat)
}

test_that("the household SAR counts the exposed susceptibles at the start of the step, with the transmitting infectious set", {
  mock_ve()
  set.seed(3)
  # every exposed susceptible is infected through the household (inf.prob 1)
  dat <- infect_general(hh_dat(c(0, 1)), 2)
  expect_equal(get_epi(dat, "hh.exposed.n")[2], 1)     # household 1 only: the admitted do not transmit
  expect_equal(get_epi(dat, "hh.sar.den")[2], 2)       # persons 2 and 3, counted before they were infected
  expect_equal(get_epi(dat, "hh.sar.num")[2], 2)
  expect_equal(get_epi(dat, "se.flow.home")[2], 2)
  expect_identical(get_attr(dat, "status")[2:3], c("e", "e"))
  # no transmission in the household: the same denominator, no numerator
  dat <- infect_general(hh_dat(c(0, 0)), 2)
  expect_equal(get_epi(dat, "hh.sar.den")[2], 2)
  expect_equal(get_epi(dat, "hh.sar.num")[2], 0)
  # with the symptomatic member admitted, no household is exposed
  d0 <- hh_dat(c(0, 1))
  d0 <- set_attr(d0, "status", "h", posit_ids = 1)
  dat <- infect_general(d0, 2)
  expect_equal(get_epi(dat, "hh.exposed.n")[2], 0)
  expect_equal(get_epi(dat, "hh.sar.den")[2], 0)
  expect_equal(get_epi(dat, "hh.sar.num")[2], 0)
})

test_that("layer flows by age partition the layer flows, and the baseline quartile flows sum to se.flow", {
  mock_ve()
  set.seed(4)
  dat <- hh_dat(c(0, 1))
  dat <- set_attr(dat, "age", c(5, 15, 70, rep(30, 57)))
  dat <- set_attr(dat, "degree_quartile_base", c(rep(1:4, 14), NA, NA, NA, NA))
  dat <- infect_general(dat, 2)
  by_age <- vapply(1:6, function(g) get_epi(dat, paste0("se.flow.home.age", g))[2], numeric(1))
  expect_equal(by_age, c(0, 1, 0, 0, 0, 1))          # persons 2 (15y) and 3 (70y)
  expect_equal(sum(by_age), get_epi(dat, "se.flow.home")[2])
  expect_equal(sum(vapply(1:6, function(g) get_epi(dat, paste0("se.flow.community.age", g))[2], numeric(1))), 0)
  qb <- vapply(0:4, function(q) get_epi(dat, paste0("se.flow.degQb", q))[2], numeric(1))
  expect_equal(qb, c(0, 0, 1, 1, 0))                 # persons 2 and 3 hold quartiles 2 and 3
  expect_equal(sum(qb), get_epi(dat, "se.flow")[2])
})

# N people under a campaign, as the allocation tests build them
dose_dat <- function(N, age, ...) {
  dat <- do.call(nh_dat, c(list(N = N, age = age), ve_params(...)))
  dat <- set_attr(dat, "degree_total", rep(1, N))
  dat <- set_attr(dat, "degree_score", rep(1, N))
  dat <- set_attr(dat, "n_layers_active", rep(1L, N))
  dat
}

test_that("doses by band and by recipient state match the dose flows", {
  set.seed(74)
  N <- 1000
  ages <- sample(c(0.5, 3, 10, 25, 40, 55, 62, 70, 80), N, replace = TRUE)
  dat <- dose_dat(N, ages, vax.budget = 0.3, vax.delivery.days = 5, vax.start = 3, vax.strategy = "random")
  # a fifth recovered, half of them immune from before the run; a few waned susceptibles
  rec <- seq(1, N, by = 5)
  dat <- set_attr(dat, "status", "r", posit_ids = rec)
  dat <- set_attr(dat, "prior.immune", 1, posit_ids = rec[seq_along(rec) %% 2 == 0])
  dat <- set_attr(dat, "ever.inf", 1, posit_ids = c(rec, seq(2, N, by = 5)))
  for (at in 1:10) dat <- vax_general(dat, at)
  e <- function(k) get_epi(dat, k)[1:10]
  total <- e("daily_vax_total")
  expect_equal(total, e("nVax1"))
  expect_equal(total, c(0, 0, rep(60, 5), 0, 0, 0))
  expect_equal(Reduce(`+`, lapply(1:6, function(g) e(paste0("nVax.age", g)))), total)
  expect_equal(Reduce(`+`, lapply(1:5, function(k) e(paste0("nVax.vaxage", k)))), total)
  expect_equal(e("nVax.sus") + e("nVax.prior") + e("nVax.rec") + e("nVax.inf"), total)
  expect_equal(e("nVax.inf"), rep(0, 10))
  expect_true(all(e("nVax.immune") <= e("nVax.rec")))
  expect_gt(sum(e("nVax.immune")), 0)
  expect_gt(sum(e("nVax.prior")), 0)
  expect_equal(cumsum(total), e("vax.doses"))
  expect_equal(sum(total), get_epi(dat, "vax.budget.doses")[10])
  # infants dosed are the recipients under 1
  vax <- get_attr(dat, "vax")
  expect_equal(sum(e("nVax.infant")), sum(vax >= 1 & ages < 1))
  expect_equal(sum(e("nVax.age1")), sum(vax >= 1 & ages < 10))
  # the campaign's opening day and its pools
  expect_equal(get_epi(dat, "alloc.day")[3], 3)
  expect_true(all(is.na(get_epi(dat, "alloc.day")[c(1, 2, 4:10)])))
  expect_equal(get_epi(dat, "n.elig.alloc")[3], N)
  expect_equal(get_epi(dat, "n.elig.willing.alloc")[3], N)
  # with vax.start 1 the campaign opens on step 2, the first step the modules run
  d1 <- dose_dat(100, 30, vax.budget = 0.1, vax.delivery.days = 1, vax.start = 1, vax.strategy = "random")
  for (at in 1:3) d1 <- vax_general(d1, at)
  expect_equal(get_epi(d1, "alloc.day")[2], 2)
})

test_that("years of life lost come from the expectancy at exact age, interpolated and discounted", {
  le <- c(70, 71, 70, seq(69, 0, length.out = 83))   # ages 0 to 85
  expect_equal(EpiModelCOVID:::yll_remaining(c(0, 1, 2), le), c(70, 71, 70))
  expect_equal(EpiModelCOVID:::yll_remaining(0.5, le), 70.5)
  expect_equal(EpiModelCOVID:::yll_remaining(1.25, le), 71 - 0.25)
  expect_equal(EpiModelCOVID:::yll_remaining(c(85, 90, 120), le), rep(le[86], 3))
  expect_equal(EpiModelCOVID:::yll_remaining(numeric(0), le), numeric(0))
  expect_equal(EpiModelCOVID:::yll_discounted(10, 0.03), (1 - exp(-0.3)) / 0.03)
  expect_equal(EpiModelCOVID:::yll_discounted(c(10, 0), 0), c(10, 0))
  expect_lt(EpiModelCOVID:::yll_discounted(40, 0.03), 40)

  # the departures module emits the flows from the exact age of each disease death
  mock_ve()
  N <- 30
  ages <- c(rep(0.5, 10), rep(64.5, 10), rep(90, 10))
  dat <- nh_dat(N, age = ages, yll.le = le, yll.discount = 0.03)
  # every one of them is a fatal symptomatic case due to die this step
  dat <- set_attr(dat, "status", rep("ic", N))
  dat <- set_attr(dat, "statusTime", rep(1, N))
  dat <- set_attr(dat, "fatal", rep(1, N))
  dat <- set_attr(dat, "hospit", rep(0, N))
  dat <- set_param(dat, "icr.rate", 1)
  dat <- deaths_covid_gmc19(dat, 2)
  expect_equal(get_epi(dat, "d.dis.flow")[2], N)
  L <- EpiModelCOVID:::yll_remaining(ages, le)
  expect_equal(get_epi(dat, "yll.flow")[2], sum(L))
  expect_equal(get_epi(dat, "yll.disc.flow")[2], sum((1 - exp(-0.03 * L)) / 0.03))
  expect_equal(get_epi(dat, "yll.flow.age1")[2], 10 * 70.5)
  expect_equal(get_epi(dat, "yll.flow.age6")[2], sum(L[11:30]))
  expect_equal(sum(vapply(1:6, function(g) get_epi(dat, paste0("yll.flow.age", g))[2], numeric(1))), sum(L))
  # disease deaths by severity band: 6-11 months, 60-69 and 80+
  sev <- vapply(1:12, function(k) get_epi(dat, paste0("d.dis.flow.sev", k))[2], numeric(1))
  expect_equal(sev, c(0, 10, rep(0, 7), 10, 0, 10))
  # without yll.le no flow is emitted; a malformed one or a bad rate stops at init
  d0 <- nh_dat(5)
  d0 <- deaths_covid_gmc19(d0, 2)
  expect_false("yll.flow" %in% names(d0$epi))
  expect_error(EpiModelCOVID:::check_yll_params(nh_dat(5, yll.le = c(70, NA))), "yll.le")
  expect_error(EpiModelCOVID:::check_yll_params(nh_dat(5, yll.le = 70)), "yll.le")
  expect_error(EpiModelCOVID:::check_yll_params(nh_dat(5, yll.discount = 1.5)), "yll.discount")
  expect_true(EpiModelCOVID:::check_yll_params(nh_dat(5, yll.le = le, yll.discount = 0)))
})

test_that("eta-squared is the between-group share of the variance", {
  eta <- EpiModelCOVID:::eta_squared
  expect_equal(eta(c(1, 2, 3, 4), c(1, 1, 2, 2)), 0.8)
  expect_equal(eta(c(1, 2, 3, 4), c(1, 2, 1, 2)), 0.2)
  expect_equal(eta(c(1, 1, 2, 2), c("a", "a", "b", "b")), 1)
  expect_equal(eta(c(1, 2, 1, 2), c(1, 1, 2, 2)), 0)
  expect_true(is.na(eta(c(3, 3, 3), c(1, 1, 2))))
  expect_true(is.na(eta(5, 1)))
  # against the ANOVA decomposition
  set.seed(9)
  x <- rnorm(500) + rep(1:5, 100); g <- rep(1:5, 100)
  a <- anova(lm(x ~ factor(g)))
  expect_equal(eta(x, g), a[["Sum Sq"]][1] / sum(a[["Sum Sq"]]))
  # NA values and groups are left out
  expect_equal(eta(c(1, 2, 3, 4, NA), c(1, 1, 2, 2, 2)), 0.8)
})

test_that("the baseline quartile breaks ties at random without touching the run's random stream", {
  q4 <- EpiModelCOVID:::score_quartile
  side <- EpiModelCOVID:::with_side_rng
  set.seed(11)
  before <- .Random.seed
  q <- side(q4(rep(1, 100), rep(TRUE, 100)))
  expect_identical(.Random.seed, before)
  expect_equal(as.numeric(table(q)), rep(25, 4))
  # the draws after it are the draws without it
  set.seed(11); a <- runif(3)
  set.seed(11); invisible(side(runif(5))); b <- runif(3)
  expect_identical(a, b)
  # different run states break ties differently; the same state the same way
  set.seed(1); q1 <- side(q4(rep(1, 100), rep(TRUE, 100)))
  set.seed(2); q2 <- side(q4(rep(1, 100), rep(TRUE, 100)))
  set.seed(1); q3 <- side(q4(rep(1, 100), rep(TRUE, 100)))
  expect_false(identical(q1, q2))
  expect_identical(q1, q3)
  # distinct scores rank as they should; the dead have no quartile
  q <- q4(c(10, 1, 5, 7, 2, 8, 3, 9), c(rep(TRUE, 7), FALSE))
  expect_equal(q, c(4L, 1L, 3L, 3L, 2L, 4L, 2L, NA))
  # under L'Ecuyer-CMRG, the stream the forked checks use
  old_kind <- RNGkind("L'Ecuyer-CMRG")[1]
  on.exit(RNGkind(old_kind), add = TRUE)
  set.seed(5); before <- .Random.seed
  invisible(side(runif(10)))
  expect_identical(.Random.seed, before)
})

# A small relationship model with sex and occupation class, as the population
# tests build it, run end to end with the GMC19 module set (the control of
# make_control_gmc19() in GlobalMix-India-NetVax). The test fits carry no
# age-pair coefficients, so the networks advance by EpiModel's resim_nets()
# rather than the exact step; the identities tested do not depend on that.
bands_i <- c("0-9y", "10-19y", "20-29y", "30-39y", "40-59y", "60+y")
act_pars_i <- list(o = c(0.9, 0.95, 0.9, 0.85, 0.9, 0.8), kappa = 5, q = rep(0.75, 6),
                   p11 = 0.5, s = 0.05, pi = c(0.1, 0, 0, 0.2, 0, 0.3), theta_F = rep(2, 6),
                   theta_I = 0.3, rho = -0.5, new_rate = rep(0.1, 6), new_share = 0.5,
                   new_mix = matrix(1, 6, 6) + diag(4, 6))
run_small <- function(nsteps = 40, seed = 77, young = 0, rsv = FALSE, ...) {
  N <- 600
  set.seed(5)
  band <- rep(1:6, length.out = N)
  sex <- rep(c("F", "M"), each = N / 2)
  occ <- c("child", "student", "worker", "home")[c(1, 2, 3, 3, 4, 4)][band]
  tr <- draw_relationship_traits(band, act_pars_i, seed = 6)
  nw <- network_initialize(n = N)
  age <- c(5, 15, 25, 35, 50, 70)[band] + runif(N, -4, 4)
  # `young` infants in the youngest band, for the infant products (#76)
  if (young > 0) age[which(band == 1)[seq_len(young)]] <- seq(0.02, 0.98, length.out = young)
  nw <- set_vertex_attribute(nw, "age", age)
  nw <- set_vertex_attribute(nw, "age.grp", bands_i[band])
  nw <- set_vertex_attribute(nw, "never", tr$never)
  nw <- set_vertex_attribute(nw, "log_aF", tr$log_aF)
  nw <- set_vertex_attribute(nw, "log_aI", tr$log_aI)
  nw <- set_vertex_attribute(nw, "sex", sex)
  nw <- set_vertex_attribute(nw, "occ", occ)
  fit <- function(attr, edges) {
    f <- as.formula(paste0("~edges + offset(nodecov(\"", attr, "\")) + offset(nodefactor(\"never\", levels = I(1)))"))
    est <- suppressMessages(suppressWarnings(netest(nw, formation = f, target.stats = edges, coef.form = c(1, -Inf),
                                                    coef.diss = dissolution_coefs(~offset(edges), 100), verbose = FALSE)))
    est$cell_theta <- matrix(est$coef.form.crude[["edges"]], 6, 6)
    est
  }
  nw_hh <- network_initialize(n = N)
  nw_hh <- set_vertex_attribute(nw_hh, "hh.ids", rep(seq_len(N / 4), each = 4))
  x <- list(freq = fit("log_aF", 900), infreq = fit("log_aI", 600),
            home = netclique(nw_hh, "hh.ids", arrivals = "join", arrivals.FUN = hh_join_birth))
  vp <- ve_params(vax.budget = 0.3, vax.delivery.days = 5, vax.start = 5, vax.strategy = "degree")
  if (rsv) vp <- as_rsv(vp)
  param <- do.call(param.net, c(list(
    time.unit = 1, inf.prob = rep(0.15, 4), exposure.rate = rep(1, 4), inf.prob.a.rr = 0.5,
    layer.names = c("freq", "infreq", "new", "home"),
    age.breaks = c(0, 10, 20, 30, 40, 60, Inf), age.grps = bands_i, arrival.age = 0,
    birth.prop.male = 0.52, birth.rate = 40 / 1000 / 365, mort.rates = rep(0.003, 86),
    ea.rate = 0.25, ar.rate = 0.2, eip.rate = 0.25, ipic.rate = 0.5, ich.rate = 0.3, icr.rate = 0.25, hr.rate = 0.15,
    import.pulse = 20, import.pulse.start = 2, import.pulse.days = 3,
    act.o = act_pars_i$o, act.kappa = act_pars_i$kappa, act.q = act_pars_i$q, act.p11 = act_pars_i$p11,
    act.s = act_pars_i$s, act.pi = act_pars_i$pi, act.theta.F = act_pars_i$theta_F,
    act.theta.I = act_pars_i$theta_I, act.rho = act_pars_i$rho, new.rate = act_pars_i$new_rate,
    new.mix = act_pars_i$new_mix, new.share.norel = act_pars_i$new_share,
    birth.age.weight = c(5.9, 84.7, 94.9, 48.8, 12.0, 3.5, 0.8),
    birth.hh.coef = setNames(rep(0, 7), EpiModelCOVID:::birth_hh_terms), birth.interval.min = 304,
    yll.le = seq(80, 5, length.out = 86), yll.discount = 0.03),
    sev_test_params(prop.clinical = rep(0.5, 12), ihr = rep(0.1, 12), ifr = rep(0.1, 12)),
    utils::modifyList(vp, list(...))))
  control <- control.net(
    type = NULL, nsteps = nsteps, nsims = 1, ncores = 1, verbose = FALSE,
    initialize.FUN = init_gmc19, aging.FUN = aging, departures.FUN = deaths_covid_gmc19,
    arrivals.FUN = arrival, nwupdate.FUN = nwupdate.net, resim_nets.FUN = resim_nets,
    summary_nets.FUN = summary_nets, activation.FUN = activate_contacts, degree.FUN = netdegree,
    infection.FUN = infect_general, progress.FUN = progress_general, vax.FUN = vax_general,
    prevalence.FUN = prevalence_covid_corporate,
    module.order = c("aging.FUN", "departures.FUN", "arrivals.FUN", "nwupdate.FUN", "resim_nets.FUN",
                     "summary_nets.FUN", "activation.FUN", "degree.FUN", "infection.FUN", "progress.FUN",
                     "vax.FUN", "prevalence.FUN"),
    resimulate.network = TRUE, tergmLite = TRUE, save.nwstats = FALSE, save.network = FALSE,
    save.transmat = FALSE, save.other = "attr", edges.correct.attr = "act.eligible")
  set.seed(seed)
  suppressWarnings(netsim(x, param, init.net(e.num = 0), control))
}

test_that("the new flows partition on a full small run, and persons ever infected are the seeds plus the first infections", {
  sim <- run_small()
  d <- as.data.frame(sim)
  d <- d[d$time > 1, ]
  rs <- function(cols) Reduce(`+`, lapply(cols, function(k) d[[k]]))
  types <- c("freq", "infreq", "new", "home")
  # an epidemic, deaths among the infected, doses and births all happened
  expect_gt(sum(d$se.flow), 100)
  expect_gt(tail(d$ever.inf.dead, 1), 0)
  expect_gt(sum(d$d.dis.flow), 0)
  expect_gt(sum(d$nNew), 0)
  # the budget is 30% of the living on the day the campaign opens, delivered in full
  expect_equal(sum(d$nVax1), tail(d$vax.budget.doses, 1))
  expect_gt(sum(d$nVax1), 170)
  # persons: the decedent-inclusive count is the living plus the departed, sums over
  # the bands, never falls, and equals the seeds (none) plus the cumulative first infections
  expect_equal(d$ever.inf.pers, d$ever.inf.num + d$ever.inf.dead)
  expect_equal(rs(paste0("ever.inf.pers.age", 1:6)), d$ever.inf.pers)
  expect_equal(rs(paste0("ever.inf.dead.age", 1:6)), d$ever.inf.dead)
  expect_true(all(diff(d$ever.inf.pers) >= 0))
  expect_equal(d$ever.inf.pers, cumsum(d$se.flow.first))
  expect_equal(d$ever.sympt.pers, d$ever.sympt.num + d$ever.sympt.dead)
  expect_equal(rs(paste0("ever.sympt.pers.age", 1:6)), d$ever.sympt.pers)
  expect_equal(d$ever.sympt.pers, cumsum(d$ipic.flow.first))
  expect_equal(d$se.flow.first + d$se.flow.reinf, d$se.flow)
  expect_equal(rs(paste0("se.flow.first.age", 1:6)), d$se.flow.first)
  # the day-1 cohort: living plus dead is the day-1 population
  n1 <- as.data.frame(sim)$num[1]
  expect_equal(d$n.cohort + d$d.cohort, rep(n1, nrow(d)))
  expect_true(all(d$ever.inf.cohort.pers <= d$ever.inf.pers))
  expect_true(all(diff(d$ever.inf.cohort.pers) >= 0))
  # layers: by tie type, and by tie type and age
  expect_equal(rs(paste0("se.flow.", types)) + d$se.import.flow, d$se.flow)
  for (k in types) expect_equal(rs(paste0("se.flow.", k, ".age", 1:6)), d[[paste0("se.flow.", k)]])
  expect_equal(rs(paste0("se.flow.degQb", 0:4)), d$se.flow)
  expect_gt(sum(d$se.flow.degQb0), 0)   # infections among the people born in the run
  # doses
  expect_equal(rs(paste0("nVax.age", 1:6)), d$nVax1)
  expect_equal(rs(paste0("nVax.vaxage", 1:5)), d$nVax1)
  expect_equal(d$nVax.sus + d$nVax.prior + d$nVax.rec + d$nVax.inf, d$nVax1)
  expect_equal(d$nVax.inf, rep(0, nrow(d)))
  expect_equal(cumsum(d$nVax1), d$vax.doses)
  expect_equal(d$alloc.day[d$time == 5], 5)
  # deaths and years of life lost
  expect_equal(rs(paste0("d.dis.flow.sev", 1:12)), d$d.dis.flow)
  expect_equal(rs(paste0("d.dis.flow.age", 1:6)), d$d.dis.flow)
  expect_equal(rs(paste0("yll.flow.age", 1:6)), d$yll.flow)
  expect_equal(rs(paste0("yll.disc.flow.age", 1:6)), d$yll.disc.flow)
  expect_true(all(d$yll.disc.flow <= d$yll.flow))
  expect_gt(sum(d$yll.flow), 0)
  # person-time
  expect_equal(rs(paste0("n.sev", 1:12)), d$num)
  expect_equal(rs(paste0("n.vaxage", 1:5)), d$num)
  expect_equal(rs(paste0("n.age", 1:6)), d$num)
  # the household SAR
  expect_true(all(d$hh.sar.num <= d$hh.sar.den))
  # fixed-day scalars: day 2 and the allocation day only
  for (s in c("total", "xlayer", "hazard")) {
    v <- d[[paste0("eta2.", s, ".d2")]]
    expect_true(!is.na(v[d$time == 2]) && v[d$time == 2] >= 0 && v[d$time == 2] <= 1)
    expect_true(all(is.na(v[d$time != 2])))
    expect_true(!is.na(d[[paste0("eta2.", s, ".alloc")]][d$time == 5]))
  }
  expect_equal(d$mean_degree_freq.d2[d$time == 2], d$mean_degree_freq[d$time == 2])
  expect_equal(d$mean_degree_hh.alloc[d$time == 5], d$mean_degree_hh[d$time == 5])
  # the baseline quartile: every survivor of day 1 holds one, nobody born since does
  attr <- sim$attr[[1]]
  born <- attr$entrTime > 1
  expect_true(any(born))
  expect_true(all(is.na(attr$degree_quartile_base[born])))
  expect_true(all(attr$degree_quartile_base[!born] %in% 1:4))
})

test_that("every #75 arm runs end to end and delivers the same doses on the same days", {
  coef <- c("(Intercept)" = 1, "band10-19y" = 0.2, sexM = 0.05, occstudent = 0.6, occworker = 0.2, hh_size = -0.1)
  hx <- list(vax.degree.score = "hazard_xlayer")
  arms <- list(
    age_oldest = list(vax.strategy = "age_oldest"),
    age_transmission = c(list(vax.strategy = "age_transmission"), hx),
    random = list(vax.strategy = "random"),
    degree = c(list(vax.strategy = "degree", vax.degree.frozen = TRUE), hx),
    precision = c(list(vax.strategy = "degree", vax.degree.frozen = TRUE, vax.rank.spearman = 0.5), hx),
    acquaintance = list(vax.strategy = "acquaintance"),
    proxy_obs = list(vax.strategy = "proxy_obs", vax.proxy.coef = coef),
    hybrid = c(list(vax.strategy = "hybrid", vax.hybrid.f = 0.5), hx))
  runs <- lapply(arms, function(a) as.data.frame(do.call(run_small, c(list(nsteps = 20), a))))
  doses <- sapply(runs, function(d) d$nVax1[-1])
  for (k in seq_len(ncol(doses))[-1]) expect_identical(doses[, k], doses[, 1], info = names(arms)[k])
  expect_gt(sum(doses[, 1]), 170)
  # the precision arm's realized rank correlation on its opening day, and the hybrid's reserve
  expect_equal(max(runs$precision$vax.rank.spearman.alloc, na.rm = TRUE), 0.5, tolerance = 0.01)
  expect_gt(max(runs$hybrid$vax.reserve.doses, na.rm = TRUE), 0)
  expect_equal(sum(runs$hybrid$nVax.reserve[-1]), max(runs$hybrid$vax.reserve.doses, na.rm = TRUE))
  expect_equal(sum(runs$age_oldest$nVax.reserve[-1]), 0)
  # the non-household hazard is part of the hazard, and age carries less of it
  expect_true(is.finite(runs$degree$eta2.hazard_xlayer.d2[2]))
})

test_that("the infant products run end to end beside the adult arms, and only the cocoon moves the adult doses (#76)", {
  # 40 infants in the youngest band, on the all-ages RSV product (the cocoon's)
  base <- list(nsteps = 20, young = 40, rsv = TRUE, vax.rsv.age.min = 5, vax.budget = 0.1)
  arms <- list(
    direct = c(infant_params("nirsevimab"), vax.strategy = "age_oldest"),
    combined = c(infant_params("nirsevimab"), vax.strategy = "cocoon"),
    mat_direct = c(infant_params("maternal"), vax.strategy = "age_oldest"),
    mat_combined = c(infant_params("maternal"), vax.strategy = "cocoon"))
  runs <- lapply(arms, function(a) as.data.frame(do.call(run_small, c(base, a))))
  # the adult line: the same doses on the same days under each program
  expect_identical(runs$direct$nVax1, runs$combined$nVax1)
  expect_identical(runs$mat_direct$nVax1, runs$mat_combined$nVax1)
  expect_gt(sum(runs$direct$nVax1, na.rm = TRUE), 50)
  # the infant line: the catch-up on the open day, the same in both arms
  cu <- runs$direct$nVax.iprog.catchup
  expect_identical(cu, runs$combined$nVax.iprog.catchup)
  expect_identical(which(cu > 0), 5L)
  expect_equal(runs$direct$n.iprog.infant[c(4, 6)] > 0, c(FALSE, TRUE))
  # the cocoon arm doses the tier first; the oldest-first arm reaches no one in it here
  d <- runs$combined[runs$combined$time == 5, ]
  expect_gt(d$n.cocoon.elig, 0)
  expect_equal(d$nVax.cocoon, min(d$nVax1, d$n.cocoon.elig))
  expect_equal(sum(runs$direct$nVax.cocoon, na.rm = TRUE), 0)
  # the maternal program holds from the first day, with mothers on record vaccinated
  m <- runs$mat_direct
  expect_gt(max(m$iprog.init, na.rm = TRUE), 0)
  expect_gt(max(m$iprog.init.mother, na.rm = TRUE), 0)
  expect_gt(m$n.iprog.infant[2], 0)
  expect_equal(sum(m$nVax.iprog.catchup, na.rm = TRUE), 0)
})
