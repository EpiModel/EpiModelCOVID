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

test_that("cumulative doses are identical across arms in a cell, on every day", {
  given <- lapply(c("age", "random", "degree"), function(s) {
    dat <- alloc_dat(2000, age = ages, degree = degs, budget = 0.2, days = 10, start = 3, strategy = s)
    dat <- run_vax(dat, 1:20)
    list(n = length(dosed(dat)), daily = get_epi(dat, "nVax1")[1:20], cum = get_epi(dat, "vax.doses")[1:20])
  })
  for (g in given) expect_identical(g$n, 400L)
  # the throughput schedule: nothing before the start, 40 a day for 10 days
  expect_equal(given[[1]]$daily, c(0, 0, rep(40, 10), rep(0, 8)))
  expect_identical(given[[1]]$daily, given[[2]]$daily)
  expect_identical(given[[1]]$daily, given[[3]]$daily)
  expect_equal(given[[1]]$cum, c(0, 0, 40 * 1:10, rep(400, 8)))
  # when the willing pool is short of the budget, every arm doses the whole pool
  short <- vapply(c("age", "random", "degree"), function(s) {
    dat <- alloc_dat(2000, age = ages, degree = degs, budget = 0.5, days = 1, strategy = s)
    dat <- set_attr(dat, "vax.willing", rep(c(1L, 0L), c(300, 1700)))
    length(dosed(run_vax(dat, 1:3)))
  }, integer(1))
  expect_identical(unname(short), rep(300L, 3))
})

test_that("with willingness 1 and unlimited supply, every arm vaccinates the same set", {
  sets <- lapply(c("age", "random", "degree", "age_then_degree"), function(s) {
    dat <- alloc_dat(500, age = ages[1:500], degree = degs[1:500], budget = 1, days = 1, strategy = s, age.min = 18)
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
  # age: oldest first
  dat <- run_vax(alloc_dat(N, age = a, degree = d, budget = 0.1, strategy = "age"), 1)
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
  # the parallel hybrid: half of today's doses to the 65 and over, oldest
  # first, the rest by degree over everyone not taken
  dat <- run_vax(alloc_dat(N, age = a, degree = d, budget = 0.1, strategy = "age_plus_degree",
                           vax.hybrid.frac = 0.5), 1)
  old <- dosed(dat)[a[dosed(dat)] >= 65]
  expect_true(length(old) >= 50)
  expect_true(min(a[old]) >= 70 || sum(a >= 65) <= 50)
  # a degraded network arm: half of today's doses to the top-ranked, the rest at random
  dat <- run_vax(alloc_dat(N, age = a, degree = d, budget = 0.1, strategy = "degree", vax.target.uptake = 0.5), 1)
  d50 <- sort(d, decreasing = TRUE)[50]                 # the 50 top-ranked share the boundary degree
  expect_true(sum(d[dosed(dat)] >= d50) >= 50)
  expect_true(sum(d[dosed(dat)] < d50) >= 25)           # the random half reaches below it
  expect_identical(length(dosed(dat)), 100L)
})

test_that("the infant products dose infants only, never more than the infant pool", {
  N <- 500
  a <- ages[1:N]
  dat <- alloc_dat(N, age = a, budget = 0.5, strategy = "infant_direct")
  dat <- set_attr(dat, "hh.ids", rep(1:100, length.out = N))
  dat <- run_vax(dat, 1)
  infants <- which(a < 1)
  expect_identical(dosed(dat), infants)
  expect_equal(get_epi(dat, "nVax1")[1], as.numeric(length(infants)))
  # cocooning: infants first, then their households' adults, then the rest
  dat <- alloc_dat(N, age = a, budget = 0.2, strategy = "infant_direct_cocoon")
  hh <- rep(1:100, length.out = N)
  dat <- set_attr(dat, "hh.ids", hh)
  dat <- run_vax(dat, 1)
  got <- dosed(dat)
  expect_true(all(infants %in% got))
  rest <- setdiff(got, infants)
  expect_true(all(hh[rest] %in% hh[infants]))
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
  for (nm in c("vax.supply.total", "vax.supply.rate", "vax.infant.eligible", "vax.covid.rate.dose1", "vax.flu.start.dose1")) {
    expect_error(EpiModelCOVID:::check_retired_params(set_param(dat, nm, 1)), "replaced by", info = nm)
  }
})
