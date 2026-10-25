# The allocation model (GlobalMix India #74): one dose budget released on a
# calendar schedule and never refilled, a willingness flag, eligibility by
# product, and the strategies as orderings of the willing eligibles

# N people under a campaign: ages, degrees and willingness as given
alloc_dat <- function(N, age = 30, budget = 0.2, days = 1, start = 1, willing = 1, strategy = "random",
                      age.min = 0, age.max = Inf, degree = NULL, ...) {
  dat <- do.call(nh_dat, c(list(N = N, age = age),
                           ve_params(vax.budget = budget, vax.delivery.days = days, vax.start = start,
                                     vax.willing = willing, vax.strategy = strategy,
                                     vax.covid.age.min = age.min, vax.covid.age.max = age.max, ...)))
  deg <- if (is.null(degree)) rep(1, N) else degree
  dat <- set_attr(dat, "degree_total", deg)
  dat <- set_attr(dat, "degree_score", deg)
  dat <- set_attr(dat, "n_layers_active", rep(1L, N))
  dat
}
run_vax <- function(dat, steps) {
  for (at in steps) dat <- vax_general(dat, at)
  dat
}
dosed <- function(dat) which(get_attr(dat, "vax") >= 1)

set.seed(74)
ages <- sample(c(0.5, 3, 10, 25, 40, 55, 62, 70, 80), 2000, replace = TRUE)
degs <- rpois(2000, 5)

arms_nh <- c("age_oldest", "age_transmission", "random", "degree", "hybrid")

test_that("cumulative doses are identical across arms in a cell, on every day", {
  given <- lapply(arms_nh, function(s) {
    dat <- alloc_dat(2000, age = ages, degree = degs, budget = 0.2, days = 10, start = 3, strategy = s,
                     vax.hybrid.f = 0.5)
    dat <- run_vax(dat, 1:20)
    list(n = length(dosed(dat)), daily = get_epi(dat, "nVax1")[1:20], cum = get_epi(dat, "vax.doses")[1:20])
  })
  for (g in given) expect_identical(g$n, 400L)
  # the throughput schedule: nothing before the start, 40 a day for 10 days
  expect_equal(given[[1]]$daily, c(0, 0, rep(40, 10), rep(0, 8)))
  for (g in given[-1]) expect_identical(given[[1]]$daily, g$daily)
  expect_equal(given[[1]]$cum, c(0, 0, 40 * 1:10, rep(400, 8)))
  # when the willing pool is short of the budget, every arm doses the whole pool
  short <- vapply(arms_nh, function(s) {
    dat <- alloc_dat(2000, age = ages, degree = degs, budget = 0.5, days = 1, strategy = s, vax.hybrid.f = 0.5)
    dat <- set_attr(dat, "vax.willing", rep(c(1L, 0L), c(300, 1700)))
    length(dosed(run_vax(dat, 1:3)))
  }, integer(1))
  expect_identical(unname(short), rep(300L, length(arms_nh)))
})

test_that("with willingness 1 and unlimited supply, every arm vaccinates the same set", {
  sets <- lapply(arms_nh, function(s) {
    dat <- alloc_dat(500, age = ages[1:500], degree = degs[1:500], budget = 1, days = 1, strategy = s, age.min = 18,
                     vax.hybrid.f = 1)
    dosed(run_vax(dat, 1:2))
  })
  elig <- which(ages[1:500] >= 18)
  for (s in sets) expect_identical(s, elig)
})

test_that("no dose goes outside the eligible set: the product's age window, the state and willingness", {
  N <- 600
  dat <- alloc_dat(N, age = ages[1:N], degree = degs[1:N], budget = 1, days = 1, strategy = "degree",
                   age.min = 18, age.max = 65)
  status <- rep(c("s", "e", "a", "ip", "ic", "h", "r"), length.out = N)
  dat <- set_attr(dat, "status", status)
  willing <- rep(c(1L, 1L, 0L), length.out = N)
  dat <- set_attr(dat, "vax.willing", willing)
  dat <- run_vax(dat, 1:2)
  elig <- which(ages[1:N] >= 18 & ages[1:N] < 65 & status %in% c("s", "r") & willing == 1)
  expect_identical(dosed(dat), elig)
  expect_equal(get_epi(dat, "n.elig")[1], as.numeric(sum(ages[1:N] >= 18 & ages[1:N] < 65 & status %in% c("s", "r"))))
  expect_equal(get_epi(dat, "n.elig.willing")[1], as.numeric(length(elig)))
  # recipients by state at the dose
  expect_equal(get_epi(dat, "nVax.sus")[1] + get_epi(dat, "nVax.rec")[1], as.numeric(length(elig)))
  expect_equal(get_epi(dat, "nVax.rec")[1], as.numeric(sum(status[elig] == "r")))
})

test_that("the eligibility vector never contains NA, and an NA stops the run", {
  dat <- alloc_dat(6, age = c(30, 30, 30, 30, 30, NA))
  e <- EpiModelCOVID:::vax_eligible(dat, "covid", active = c(1, 1, 1, 1, 1, 0), status = c("s", "r", "e", "ic", "s", "s"),
                                    age = c(30, 30, 30, 30, 30, NA))
  expect_false(anyNA(e))
  expect_identical(e, c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_error(EpiModelCOVID:::vax_eligible(dat, "covid", active = rep(1, 6), status = rep("s", 6),
                                            age = c(30, 30, 30, 30, 30, NA)), "eligibility is NA for 1 node")
})

test_that("a waned, untested susceptible stays eligible and is dosed", {
  dat <- alloc_dat(4, budget = 1)
  dat <- set_attr(dat, "status", c("s", "s", "r", "s"))
  dat <- set_attr(dat, "ever.inf", c(0, 1, 1, 1))
  dat <- set_attr(dat, "recTime", c(NA, -40, -5, -100))
  dat <- run_vax(dat, 1)
  expect_identical(dosed(dat), 1:4)
  expect_equal(get_epi(dat, "nVax.prior")[1], 2)
})

test_that("the budget is fixed from the living when the campaign opens and is not refilled by deaths", {
  dat <- alloc_dat(1000, budget = 0.1, days = 1, start = 2)
  dat <- set_attr(dat, "active", rep(c(1, 0), c(900, 100)))
  dat <- run_vax(dat, 1:2)
  expect_equal(get_epi(dat, "vax.budget.doses")[2], 90)
  expect_identical(length(dosed(dat)), 90L)
  # fifty recipients die: no dose comes back
  active <- get_attr(dat, "active"); active[dosed(dat)[1:50]] <- 0
  dat <- set_attr(dat, "active", active)
  dat <- run_vax(dat, 3:5)
  expect_equal(get_epi(dat, "nVax1")[3:5], c(0, 0, 0))
  expect_equal(get_epi(dat, "vax.doses")[5], 90)
  expect_identical(sum(get_attr(dat, "vax")[active == 1] >= 1), 40L)
})

test_that("the throughput schedule carries a day's shortfall forward", {
  dat <- alloc_dat(1000, budget = 0.1, days = 10)        # 10 doses a day
  dat <- set_attr(dat, "vax.willing", rep(c(1L, 0L), c(15, 985)))
  dat <- run_vax(dat, 1:2)
  expect_equal(get_epi(dat, "nVax1")[1:2], c(10, 5))
  willing <- get_attr(dat, "vax.willing"); willing[16:200] <- 1L
  dat <- set_attr(dat, "vax.willing", willing)
  dat <- run_vax(dat, 3:12)
  expect_equal(get_epi(dat, "nVax1")[3], 15)        # 30 allowed by day 3, 15 given
  expect_equal(get_epi(dat, "vax.doses")[12], 100)
  expect_equal(get_epi(dat, "vax.allowed")[1:12], c(10 * 1:10, 100, 100))
})

test_that("willingness is drawn by vaccine age group, once per person", {
  dat <- alloc_dat(5, age = c(2, 10, 30, 55, 70), willing = c(1, 0, 1, 0, 1))
  expect_identical(EpiModelCOVID:::draw_vax_willing(dat, c(2, 10, 30, 55, 70)), c(1L, 0L, 1L, 0L, 1L))
  dat0 <- alloc_dat(5, willing = 0)
  expect_identical(EpiModelCOVID:::draw_vax_willing(dat0, rep(30, 5)), rep(0L, 5))
  set.seed(1)
  p <- mean(EpiModelCOVID:::draw_vax_willing(alloc_dat(5, willing = 0.3), rep(30, 10000)))
  expect_equal(p, 0.3, tolerance = 0.05)
  bad <- alloc_dat(5, willing = c(0.5, 0.5))
  expect_error(EpiModelCOVID:::draw_vax_willing(bad, rep(30, 5)), "one per vaccine age group")
  expect_error(check_vax_campaign(bad), "one per vaccine age group")
})

test_that("the strategies are orderings of the willing eligibles", {
  N <- 1000
  a <- ages[1:N]; d <- degs[1:N]
  # age_oldest: oldest first
  dat <- run_vax(alloc_dat(N, age = a, degree = d, budget = 0.1, strategy = "age_oldest"), 1)
  expect_true(min(a[dosed(dat)]) >= max(a[-dosed(dat)]))
  # degree: highest score first
  dat <- run_vax(alloc_dat(N, age = a, degree = d, budget = 0.1, strategy = "degree"), 1)
  expect_true(min(d[dosed(dat)]) >= max(d[-dosed(dat)]))
  # random: a uniform draw, so the dosed share is the same across age bands
  dat <- run_vax(alloc_dat(N, age = a, degree = d, budget = 0.5, strategy = "random"), 1)
  share <- tapply(seq_len(N) %in% dosed(dat), a >= 55, mean)
  expect_equal(unname(share[1]), unname(share[2]), tolerance = 0.1)
  # one eligible person: dosed (sample(x, n) used to draw from 1:x)
  dat <- alloc_dat(3, budget = 1, strategy = "random")
  dat <- set_attr(dat, "vax.willing", c(0L, 1L, 0L))
  expect_identical(dosed(run_vax(dat, 1)), 2L)
})

test_that("cocoon doses the 18-49 of young-infant households first, at random, then the oldest (#76)", {
  # 100 households of 4: in households 1-20 the first member is under 6 months,
  # in 21-30 between 6 and 12 months (no cocoon), elsewhere an adult
  set.seed(76)
  N <- 400
  hh <- rep(1:100, each = 4)
  a <- runif(N, 18, 80)
  first <- seq(1, N, by = 4)
  a[first[1:20]] <- 0.3
  a[first[21:30]] <- 0.8
  tier <- hh <= 20 & a >= 18 & a < 50
  cocoon_dat <- function(budget, days = 1) {
    d <- alloc_dat(N, age = a, budget = budget, days = days, strategy = "cocoon", age.min = 18)
    set_attr(d, "hh.ids", hh)
  }
  expect_identical(EpiModelCOVID:::cocoon_tier(cocoon_dat(0)), as.integer(tier))
  dat <- run_vax(cocoon_dat(0.25), 1)                  # 100 doses
  got <- dosed(dat)
  expect_true(all(which(tier) %in% got))
  expect_equal(get_epi(dat, "nVax.cocoon")[1], sum(tier))
  expect_equal(get_epi(dat, "n.cocoon.tier")[1], sum(tier))
  expect_equal(get_epi(dat, "n.cocoon.elig")[1], sum(tier))
  # the rest oldest first, with the infants (outside the product's window) left out
  rest <- setdiff(got, which(tier))
  others <- setdiff(which(!tier & a >= 18), got)
  expect_length(rest, 100 - sum(tier))
  expect_gte(min(a[rest]), max(a[others]))
  # at random within the tier: with 8 doses a day, every member is dosed in some draw
  hits <- integer(N)
  for (r in 1:50) {
    g <- dosed(run_vax(cocoon_dat(0.02), 1))
    expect_true(all(tier[g]))
    hits[g] <- hits[g] + 1L
  }
  expect_true(all(hits[tier] > 0))
  # the tier is read each day: a newborn in household 50 brings its adults in
  dat <- run_vax(cocoon_dat(0.25, days = 2), 1)
  a2 <- a
  a2[first[50]] <- 0.01
  dat <- set_attr(dat, "age", a2)
  new_tier <- hh == 50 & a2 >= 18 & a2 < 50 & !seq_len(N) %in% dosed(dat)
  dat <- run_vax(dat, 2)
  expect_true(all(which(new_tier) %in% dosed(dat)))
  # every arm reads the tier, and the oldest-first arm gives it nothing here
  dat <- alloc_dat(N, age = a, budget = 0.1, strategy = "age_oldest", age.min = 18)
  dat <- run_vax(set_attr(dat, "hh.ids", hh), 1)
  expect_equal(get_epi(dat, "n.cocoon.tier")[1], sum(tier))
  expect_equal(get_epi(dat, "nVax.cocoon")[1], sum(tier[dosed(dat)]))
  # a product that does not cover the tier stops the run; the knobs are checked
  expect_error(check_vax_campaign(alloc_dat(5, strategy = "cocoon", age.min = 60)), "does not cover")
  expect_true(check_vax_campaign(alloc_dat(5, strategy = "cocoon", age.min = 5)))
  expect_error(check_vax_campaign(alloc_dat(5, vax.cocoon.age.min = 50)), "vax.cocoon.age.min")
  expect_error(check_vax_campaign(alloc_dat(5, vax.cocoon.infant.months = 0)), "vax.cocoon.infant.months")
})

test_that("a second dose counts against the same budget, and an annual dose reopens each season", {
  dat <- alloc_dat(100, budget = 1.5, n.doses = 2)
  dat <- run_vax(dat, 1:30)
  expect_identical(sum(get_attr(dat, "vax") >= 1), 100L)
  expect_identical(sum(get_attr(dat, "vax") == 2), 50L)
  expect_equal(get_epi(dat, "nVax2")[21], 0)
  expect_equal(get_epi(dat, "nVax2")[22], 50)        # interval 21; the budget's last 50 doses
  expect_equal(get_epi(dat, "vax.doses")[30], 150)
  # annual: one dose a season, re-up from the season start, from one budget
  dat <- alloc_dat(100, budget = 1, vax.covid.annual.dose1 = TRUE, vax.season.length = 50)
  dat <- set_attr(dat, "vax.willing", rep(c(1L, 0L), 50))
  dat <- run_vax(dat, 1:60)
  expect_equal(get_epi(dat, "nVax1")[c(1, 2, 51, 52)], c(50, 0, 50, 0))
  expect_equal(get_epi(dat, "vax.doses")[60], 100)
})

test_that("check_vax_campaign() stops a bad campaign", {
  expect_true(check_vax_campaign(alloc_dat(5)))
  expect_error(check_vax_campaign(alloc_dat(5, budget = -1)), "vax.budget must be one number")
  expect_error(check_vax_campaign(alloc_dat(5, days = 0)), "vax.delivery.days must be")
  expect_error(check_vax_campaign(alloc_dat(5, start = 0)), "vax.start must be")
  expect_error(check_vax_campaign(alloc_dat(5, willing = 1.5)), "vax.willing must be")
  expect_error(check_vax_campaign(alloc_dat(5, age.min = 65, age.max = 60)), "age.min must be below")
  # no campaign parameters at all: the run does not start
  dat <- do.call(nh_dat, c(list(N = 5), ve_params()))
  dat$param$vax.budget <- NULL
  expect_error(check_vax_campaign(dat), "vax.budget must be")
})

test_that("the retired supply caps, uptake rates and infant gate stop the run", {
  dat <- alloc_dat(5)
  for (nm in c("vax.supply.total", "vax.supply.rate", "vax.infant.eligible", "vax.covid.rate.dose1", "vax.flu.start.dose1",
               "vax.infant.product")) {
    expect_error(EpiModelCOVID:::check_retired_params(set_param(dat, nm, 1)), "replaced by", info = nm)
  }
})

# The strategy set (GlobalMix India #75) -------------------------------------

test_that("retired and unknown strategies stop the run, and the knobs are checked", {
  for (s in c("age", "bridge", "age_then_degree", "proxy_hh", "proxy_fit", "age_plus_proxy", "age_plus_degree",
              "infant_hh", "infant_direct", "infant_direct_cocoon")) {
    expect_error(check_vax_campaign(alloc_dat(5, strategy = s)), "was retired", info = s)
  }
  expect_error(check_vax_campaign(alloc_dat(5, strategy = "oldest")), "unknown vax.strategy")
  expect_error(check_vax_campaign(alloc_dat(5, strategy = "degree", vax.rank.spearman = 0)), "vax.rank.spearman")
  expect_error(check_vax_campaign(alloc_dat(5, strategy = "hybrid", vax.hybrid.f = 1.2)), "vax.hybrid.f")
  expect_error(check_vax_campaign(alloc_dat(5, strategy = "acquaintance")), "needs the relationship model")
  expect_error(check_vax_campaign(alloc_dat(5, strategy = "proxy_obs")), "vax.proxy.coef")
  for (nm in c("vax.target.uptake", "vax.hybrid.frac")) {
    expect_error(EpiModelCOVID:::check_retired_params(set_param(alloc_dat(5), nm, 1)), "replaced by", info = nm)
  }
})

test_that("age_transmission ranks five-year bands by their mean score on the opening day, at random within a band", {
  N <- 2000
  set.seed(11)
  a <- runif(N, 18, 80)
  b1 <- a >= 30 & a < 35; b2 <- a >= 45 & a < 50
  # the 30-34 band has the highest score, then 45-49; everyone else 1
  d <- ifelse(b1, 10, ifelse(b2, 5, 1))
  dat <- alloc_dat(N, age = a, degree = d, budget = 0.5, days = 5, strategy = "age_transmission")
  dat <- run_vax(dat, 1)
  got <- dosed(dat)                                   # 200 doses on day 1
  expect_true(all(which(b1) %in% got))
  expect_true(all(b2[setdiff(got, which(b1))]))
  # within a band the order is random, not by age
  expect_lt(abs(mean(a[setdiff(got, which(b1))]) - 47.5), 0.6)
  # the band order was fixed on the opening day: a later change of score moves nothing
  dat <- set_attr(dat, "degree_score", ifelse(a >= 60, 100, 1))
  dat <- run_vax(dat, 2)
  got <- dosed(dat)
  expect_true(all(which(b2) %in% got))
  rest <- setdiff(got, which(b1 | b2))
  expect_lt(mean(a[rest] >= 60), 0.6)                 # about 20 / 52 of the other ages, not all 60+
  expect_equal(EpiModelCOVID:::age_band5(c(0, 4.9, 5, 17.99, 74.9, 75, 99)), c(0, 0, 5, 15, 70, 75, 75))
})

test_that("a frozen degree arm ranks once on the opening day; an unfrozen one follows today's score", {
  N <- 1000
  d0 <- rep(c(10, 1), c(100, 900))
  d1 <- rep(c(1, 10), c(100, 900))                  # the scores swap after the first day
  run <- function(frozen) {
    dat <- alloc_dat(N, degree = d0, budget = 0.2, days = 4, strategy = "degree", vax.degree.frozen = frozen)
    dat <- run_vax(dat, 1)
    dat <- set_attr(dat, "degree_score", d1)
    run_vax(dat, 2)
  }
  frozen <- dosed(run(TRUE))
  live <- dosed(run(FALSE))
  expect_true(all(1:100 %in% frozen))               # 50 a day: all of the first hundred by day 2
  expect_identical(sum(frozen <= 100), 100L)
  expect_identical(sum(live <= 100), 50L)           # day 1 from them, day 2 from the new top
  # someone born after the opening day ranks last on a frozen arm
  dat <- alloc_dat(10, degree = rep(1, 10), budget = 0.5, days = 5, strategy = "degree", vax.degree.frozen = TRUE)
  dat <- run_vax(dat, 1)
  camp <- dat$run$vax.campaign
  s <- EpiModelCOVID:::strategy_score(dat, camp, "degree", rep(1, 11), rep(30, 11), c(get_attr(dat, "unique_id"), 999))
  expect_identical(s[11], -Inf)
})

test_that("the precision arm's noise is fixed per person and calibrated to the target Spearman correlation", {
  N <- 4000
  set.seed(75)
  d <- rpois(N, 4)                                  # a tied integer score, as degree is
  ns <- EpiModelCOVID:::normal_score(d, sort(d))
  expect_equal(rank(ns), rank(d))                   # the normal score keeps the order and the ties
  for (rho in c(0.9, 0.7, 0.5, 0.3)) {
    dat <- alloc_dat(N, degree = d, budget = 0.25, days = 10, start = 2, strategy = "degree", vax.rank.spearman = rho)
    set.seed(1)
    dat <- run_vax(dat, 1:2)
    camp <- dat$run$vax.campaign
    expect_equal(camp$spearman, rho, tolerance = 0.002)
    expect_equal(get_epi(dat, "vax.rank.spearman.alloc")[2], camp$spearman)
    # the dosed are the top of the noisy order, not of the score
    noisy <- ns + camp$sigma * camp$z
    expect_true(min(noisy[dosed(dat)]) >= max(noisy[-dosed(dat)]) - 1e-12)
    if (rho == 0.9) z9 <- camp$z else expect_identical(camp$z, z9)   # one noise for every level
  }
  # an exact arm draws no noise
  dat <- run_vax(alloc_dat(N, degree = d, budget = 0.25, strategy = "degree"), 1)
  expect_identical(dat$run$vax.campaign$sigma, 0)
  expect_true(min(d[dosed(dat)]) >= max(d[-dosed(dat)]))
  # the noise is drawn on a side stream: the main stream is where it was
  set.seed(9); r0 <- runif(1)
  set.seed(9); invisible(EpiModelCOVID:::with_side_rng(rnorm(10))); expect_identical(runif(1), r0)
})

test_that("the hybrid reserves f of the willing 60+ on the opening day, released with the budget, oldest first", {
  N <- 3000
  a <- rep(c(25, 45, 65, 75), length.out = N)
  d <- rep(c(9, 5, 1, 1), length.out = N)
  dat <- alloc_dat(N, age = a, degree = d, budget = 0.3, days = 10, start = 2, strategy = "hybrid", vax.hybrid.f = 0.5)
  dat <- run_vax(dat, 1:13)
  camp <- dat$run$vax.campaign
  n60 <- sum(a >= 60)
  expect_identical(camp$reserve, as.integer(round(0.5 * n60)))      # 750 of the 900-dose budget
  expect_equal(get_epi(dat, "vax.reserve.doses")[2], 750)
  res <- get_epi(dat, "nVax.reserve")[1:13]
  expect_equal(sum(res), 750)
  expect_equal(res[2:11], diff(c(0, floor(750 * (90 * 1:10) / 900))))   # in step with the budget
  got <- dosed(dat)
  expect_identical(length(got), 900L)
  expect_identical(sum(a[got] >= 60), 750L)
  expect_identical(sum(a[got] == 75), 750L)                          # oldest first
  expect_identical(sum(d[got] == 9), 150L)                           # the rest by degree
  # a reserve larger than the budget takes the whole budget
  dat <- run_vax(alloc_dat(N, age = a, degree = d, budget = 0.1, strategy = "hybrid", vax.hybrid.f = 1), 1)
  expect_identical(dat$run$vax.campaign$reserve, 300L)
  expect_true(all(a[dosed(dat)] >= 60))
  # f = 0 is the degree arm
  dat <- run_vax(alloc_dat(N, age = a, degree = d, budget = 0.1, strategy = "hybrid", vax.hybrid.f = 0), 1)
  expect_true(all(d[dosed(dat)] == 9))
  # under a younger eligibility window the 60+ are absent and the queue is degree
  dat <- run_vax(alloc_dat(N, age = a, degree = d, budget = 0.1, strategy = "hybrid", vax.hybrid.f = 1,
                           age.max = 60), 1)
  expect_identical(dat$run$vax.campaign$reserve, 0L)
})

test_that("acquaintance nominees come in proportion to the friendship-paradox weight", {
  # a star of 1 hub and 20 leaves, and 79 people with no frequent relationship
  N <- 100
  dat <- alloc_dat(N, strategy = "acquaintance", act.o = 1, budget = 0.3)
  dat$run$el <- list(freq = cbind(1L, 2:21))
  w <- EpiModelCOVID:::acquaintance_weight(dat)
  expect_equal(w[1], 20)                     # each leaf names the hub
  expect_equal(w[2:21], rep(1 / 20, 20))     # the hub names one leaf in 20
  expect_equal(w[22:N], rep(0, N - 21))
  # the named come first, the hub almost always first of them (20 / 21), and
  # then the unnamed at random
  ords <- replicate(2000, EpiModelCOVID:::allocation_strategy(1:N, 30, "acquaintance", rep(30, N), score = w))
  expect_true(all(apply(ords, 2, function(o) setequal(o[1:21], 1:21))))
  expect_equal(mean(ords[1, ] == 1), 20 / 21, tolerance = 0.015)
  dat <- run_vax(dat, 1)
  expect_true(all(1:21 %in% dosed(dat)))
  # a weighted draw: a person named twice as often is taken first two times in three
  x <- replicate(20000, EpiModelCOVID:::allocation_strategy(1:2, 1, "acquaintance", c(30, 30), score = c(2, 1)))
  expect_equal(mean(x == 1), 2 / 3, tolerance = 0.02)
})

test_that("proxy_obs ranks on the frozen linear predictor of the recorded traits", {
  N <- 8
  dat <- alloc_dat(N, age = c(5, 15, 25, 35, 45, 65, 25, 25), strategy = "proxy_obs", budget = 0.25,
                   age.grps = c("0-9y", "10-19y", "20-29y", "30-39y", "40-59y", "60+y"),
                   vax.proxy.coef = c("(Intercept)" = 1, "band10-19y" = 0.5, "band60+y" = -1, sexM = 0.1,
                                      occstudent = 0.6, occworker = 0.2, hh_size = -0.1))
  dat <- set_attr(dat, "age.grp", c("0-9y", "10-19y", "20-29y", "30-39y", "40-59y", "60+y", "20-29y", "20-29y"))
  dat <- set_attr(dat, "sex", c("F", "M", "F", "M", "F", "M", "F", "F"))
  dat <- set_attr(dat, "occ", c("child", "student", "worker", "worker", "home", "home", "worker", "worker"))
  dat <- set_attr(dat, "hh.ids", c(1, 1, 2, 2, 2, 3, 4, NA))
  s <- EpiModelCOVID:::proxy_obs_score(dat)
  expect_equal(s, c(1 - 0.2, 1 + 0.5 + 0.1 + 0.6 - 0.2, 1 + 0.2 - 0.3, 1 + 0.1 + 0.2 - 0.3, 1 - 0.3,
                    1 - 1 + 0.1 - 0.1, 1 + 0.2 - 0.1, 1 + 0.2 - 0.1))
  dat <- run_vax(dat, 1)
  got <- dosed(dat)                                  # 2 doses: the student, then a worker living alone
  expect_true(2L %in% got)
  expect_true(setdiff(got, 2L) %in% c(7L, 8L))
  # an unknown term stops
  bad <- set_param(dat, "vax.proxy.coef", c("(Intercept)" = 1, hh_size = 0, educ = 1))
  expect_error(EpiModelCOVID:::check_proxy_coef(bad, get_param(bad, "vax.proxy.coef")), "must be named")
})

test_that("the noise calibration reaches its target and is monotone", {
  set.seed(3)
  ns <- qnorm(ppoints(5000)); z <- rnorm(5000)
  for (rho in c(0.95, 0.6, 0.2)) {
    s <- EpiModelCOVID:::calibrate_rank_noise(ns, z, rho)
    expect_equal(EpiModelCOVID:::spearman(ns, ns + s * z), rho, tolerance = 0.002)
  }
  expect_identical(EpiModelCOVID:::calibrate_rank_noise(ns, z, 1), 0)
})
