# Per-layer parameters and the removed CorporateMix modifiers (GlobalMix India #66)

N <- 60
set.seed(2)
hh <- rep(1:15, each = 4)
nw <- network_initialize(n = N)
est_com <- suppressMessages(netest(nw, formation = ~edges, target.stats = 40,
                                   coef.diss = dissolution_coefs(~offset(edges), 20),
                                   verbose = FALSE))
nw_hh <- network_initialize(n = N)
nw_hh <- set_vertex_attribute(nw_hh, "hh.ids", hh)
est_hh <- netclique(nw_hh, "hh.ids")

# A tergmLite netsim_dat object at t = 1 on a community layer and a household
# clique layer, with the attributes and parameters infect_general() reads. One
# member of each household is symptomatic, so both layers have discordant ties.
layer_dat <- function(x = list(community = est_com, home = est_hh),
                      layer.names = names(x), inf.prob = c(0.1, 0.1),
                      exposure.rate = c(1, 1)) {
  param <- param.net(inf.prob = 0.1, act.rate = 1)
  init <- init.net(i.num = 5)
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE)
  EpiModel:::crosscheck.net(x, param, init, control)
  control <- EpiModel:::netsim_validate_control(control)
  dat <- initialize.net(x, param, init, control, s = 1)
  if (!is.null(layer.names)) dat <- set_param(dat, "layer.names", layer.names)
  dat <- set_param(dat, "inf.prob", inf.prob)
  dat <- set_param(dat, "exposure.rate", exposure.rate)
  dat <- set_param(dat, "inf.prob.a.rr", 0.5)
  dat <- set_param(dat, "age.breaks", c(0, 10, 20, 30, 40, 60, Inf))
  dat <- set_param(dat, "arrival.age", 0)
  status <- rep("s", N)
  status[seq(1, N, by = 4)] <- "ic"
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "vax", rep(0, N))
  dat <- set_attr(dat, "ever.inf", rep(0, N))
  dat <- set_attr(dat, "last.dose.time", rep(NA, N))
  dat <- set_attr(dat, "age", rep(30, N))
  dat <- set_attr(dat, "is_infant", rep(FALSE, N))
  dat <- set_attr(dat, "degree_quartile", rep(1, N))
  dat <- set_attr(dat, "is_bridge", rep(FALSE, N))
  dat
}

# Run infect_general() once and return the per-tie daily transmission
# probabilities it passed to rbinom(), one vector per layer with discordant
# ties. Vaccine effects are mocked away (rr = 1) and no transmission happens.
tie_probs <- function(dat, at) {
  probs <- list()
  local_mocked_bindings(
    rbinom = function(n, size, prob) {
      probs[[length(probs) + 1]] <<- prob
      integer(n)
    },
    build_vax_schedule = function(dat, disease) NULL,
    compute_ve = function(at, ids, ...) list(rr = rep(1, length(ids))),
    vax_age_group_for = function(dat) rep(NA_integer_, length(get_attr(dat, "active"))),
    .package = "EpiModelCOVID"
  )
  infect_general(dat, at)
  probs
}

test_that("every tie carries the same daily hazard on every layer and every day", {
  dat <- EpiModelCOVID:::check_layer_params(layer_dat())
  # day 20 is after the old day-15 cut on every layer but the first, day 60 is
  # inside the old inf.add window, day 140 inside the old inf.sub window
  for (at in c(2, 20, 60, 140, 300)) {
    probs <- tie_probs(dat, at)
    expect_length(probs, 2)
    expect_equal(unlist(probs), rep(0.1, length(unlist(probs))))
  }
})

test_that("exposure.rate and inf.prob act on their own layer", {
  dat <- EpiModelCOVID:::check_layer_params(
    layer_dat(inf.prob = c(0.1, 0.2), exposure.rate = c(1, 2)))
  probs <- tie_probs(dat, 20)
  expect_equal(unique(probs[[1]]), 0.1)
  expect_equal(unique(probs[[2]]), 1 - 0.8^2)
})

test_that("a transmission probability that is NA or outside [0, 1] stops the run", {
  dat <- EpiModelCOVID:::check_layer_params(layer_dat(inf.prob = c(0.1, NA)))
  expect_error(tie_probs(dat, 2), "probabilities on layer `home` at step 2 are NA")
  dat <- EpiModelCOVID:::check_layer_params(layer_dat(inf.prob = c(1.5, 0.1)))
  expect_error(tie_probs(dat, 2), "on layer `community`.*outside \\[0, 1\\]")
})

test_that("per-layer parameters are named by layer and must match the layers", {
  dat <- EpiModelCOVID:::check_layer_params(layer_dat())
  expect_named(get_param(dat, "inf.prob"), c("community", "home"))
  expect_named(get_param(dat, "exposure.rate"), c("community", "home"))

  expect_error(EpiModelCOVID:::check_layer_params(layer_dat(layer.names = NULL)),
               "no layer.names")
  expect_error(EpiModelCOVID:::check_layer_params(layer_dat(layer.names = c("home", "community"))),
               "are not the layers passed to netsim\\(\\) \\(community, home\\)")
  expect_error(EpiModelCOVID:::check_layer_params(layer_dat(x = list(est_com, est_hh),
                                                            layer.names = c("community", "home"))),
               "name every layer")
  expect_error(EpiModelCOVID:::check_layer_params(layer_dat(inf.prob = c(0.1, 0.1, 0.1))),
               "`inf.prob` has 3 values for 2 layers")
  expect_error(EpiModelCOVID:::check_layer_params(
    layer_dat(exposure.rate = c(home = 1, community = 5))),
    "`exposure.rate` is named home, community")

  # after initialization, as an updater would change it: an unnamed vector is
  # read in layer order, a misnamed one stops
  dat <- set_param(dat, "inf.prob", c(0.2, 0.3))
  expect_equal(EpiModelCOVID:::layer_param(dat, "inf.prob"), c(community = 0.2, home = 0.3))
  dat <- set_param(dat, "inf.prob", c(home = 0.2, community = 0.3))
  expect_error(tie_probs(dat, 2), "`inf.prob` is named home, community")
})

test_that("a retired CorporateMix parameter stops the run", {
  dat <- layer_dat()
  expect_true(EpiModelCOVID:::check_retired_params(dat))
  expect_error(EpiModelCOVID:::check_retired_params(set_param(dat, "hosp.boost.mult", 5)),
               "no longer reads hosp.boost.mult")
  # a neutral value is refused too: nothing would read it
  expect_error(EpiModelCOVID:::check_retired_params(set_param(dat, "inf.add", 0)),
               "no longer reads inf.add")
  dat$param$.param.updater.list <- list(list(at = 100, param = list(dx.rate.sympt = 0.1)))
  expect_error(EpiModelCOVID:::check_retired_params(dat), "no longer reads dx.rate.sympt")
})

test_that("the model carries no testing or isolation state", {
  attrs <- names(EpiModelCOVID:::get_default_attrs(layer_dat()))
  expect_false(any(c("dxStatus", "dxTime", "isolate", "isoTime") %in% attrs))
})

test_that("people returned from R to S stay eligible for vaccination", {
  dose_row <- data.frame(dose = 1, interval = NA, annual = FALSE)
  dose_row$start <- I(list(rep(20, 5)))
  # an R to S return used to leave dxStatus NA, which made the node ineligible
  elig <- EpiModelCOVID:::get_ids_eligible_for_dose(
    dose_row = dose_row, dose_i = 1, active = rep(1, 4),
    status = c("s", "r", "ic", "h"), vax = rep(0, 4), last.dose.time = rep(NA, 4),
    vax.age.group = rep(3, 4), at = 30)
  expect_equal(elig, 1:2)
})
