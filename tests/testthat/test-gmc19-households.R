
# Network membership and the household clique layer (GlobalMix India #85)

N <- 60
set.seed(1)
hh <- rep(1:15, each = 4)
age.grp <- rep(c("0-9y", "20-29y", "30-39y", "60+y"), 15)
# households 1-3 have no 0-9y member, so newborns can never join them
age.grp[hh %in% 1:3] <- "40-59y"

nw <- network_initialize(n = N)
nw <- set_vertex_attribute(nw, "age.grp", age.grp)
est_com <- suppressMessages(netest(nw, formation = ~edges, target.stats = 40,
                                   coef.diss = dissolution_coefs(~offset(edges), 20),
                                   verbose = FALSE))
nw_hh <- network_initialize(n = N)
nw_hh <- set_vertex_attribute(nw_hh, "hh.ids", hh)
est_hh <- netclique(nw_hh, "hh.ids", arrivals = "join",
                    arrivals.FUN = hh_join_young)

# A tergmLite netsim_dat object at t = 1, school-style layer first and the
# household clique layer second
dat_at_start <- function(x = list(est_com, est_hh)) {
  param <- param.net(inf.prob = 0.1, act.rate = 1)
  init <- init.net(i.num = 5)
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE)
  EpiModel:::crosscheck.net(x, param, init, control)
  control <- EpiModel:::netsim_validate_control(control)
  initialize.net(x, param, init, control, s = 1)
}

# One step of the engine's vital dynamics without the disease modules: mark
# `dead` (departures module), append `n_new` newborns (arrivals module), then
# nwupdate_gmc19()
vital_step <- function(dat, at, dead, n_new) {
  dat <- set_current_timestep(dat, at)
  dat <- set_attr(dat, "active", 0, posit_ids = dead)
  dat <- set_attr(dat, "exitTime", at, posit_ids = dead)
  if (n_new > 0) {
    dat <- append_core_attr(dat, at, n_new)
    n <- length(get_attr(dat, "active"))
    for (a in names(get_attr_list(dat))) {
      v <- get_attr(dat, a)
      if (length(v) < n) {
        val <- switch(a, age.grp = "0-9y", status = "s", hh.ids = NA, NA)
        dat <- append_attr(dat, a, val, n - length(v))
      }
    }
  }
  nwupdate_gmc19(dat, at)
}

pairs_of <- function(el) {
  el <- matrix(as.integer(el), ncol = 2)
  sort(paste(pmin(el[, 1], el[, 2]), pmax(el[, 1], el[, 2])))
}

test_that("hh_layer finds the household layer by its grouping attribute", {
  expect_equal(EpiModelCOVID:::hh_layer(dat_at_start()), 2L)
  expect_equal(EpiModelCOVID:::hh_layer(dat_at_start(list(est_hh, est_com))), 1L)
  expect_true(is.na(EpiModelCOVID:::hh_layer(dat_at_start(list(est_com)))))
})

test_that("hh_join_young draws only living, existing 0-9y members' households", {
  dat <- dat_at_start()
  dat <- set_current_timestep(dat, 2)
  # kill every 0-9y member except household 7's
  young <- which(age.grp == "0-9y")
  dat <- set_attr(dat, "active", 0, posit_ids = setdiff(young, which(hh == 7)))
  dat <- append_core_attr(dat, 2, 3)
  dat <- append_attr(dat, "age.grp", "0-9y", 3)
  dat <- append_attr(dat, "hh.ids", NA, 3)
  new_ids <- N + 1:3
  # a length-1 pool must not trigger sample()'s 1:value behavior
  expect_equal(hh_join_young(dat, 2, new_ids, 2), rep(7, 3))

  # an empty pool leaves the newborns in no household
  dat <- set_attr(dat, "active", 0, posit_ids = which(hh == 7))
  expect_equal(hh_join_young(dat, 2, new_ids, 2), rep(NA, 3))
})

test_that("nwupdate_gmc19 removes the dead and wires newborns into full cliques", {
  dat <- dat_at_start()
  # Every 0-9y member dies except household 5's (hh5[1]), so household 5 is the
  # only one newborns can join and all three are same-day siblings there. One
  # adult of household 5 and one member of household 2 die too.
  young <- which(age.grp == "0-9y")
  hh5 <- which(hh == 5)
  dead <- c(setdiff(young, hh5[1]), hh5[2], which(hh == 2)[1])
  dat <- vital_step(dat, at = 2, dead = dead, n_new = 3)

  n1 <- N - length(dead) + 3
  hh1 <- get_attr(dat, "hh.ids")
  expect_length(hh1, n1)
  expect_true(all(get_attr(dat, "active") == 1))
  for (network in 1:2) {
    expect_equal(attr(dat$run$el[[network]], "n"), n1)
    expect_equal(dat$run$net_attr[[network]][["n"]], n1)
  }
  # the newborns are the last three nodes and all joined household 5
  expect_equal(hh1[n1 - 2:0], rep(5, 3))
  expect_equal(sum(hh1 == 5), 3 + 3)
  expect_equal(sum(hh1 == 2), 3)
  # the household layer is exactly the cliques of hh.ids: the newborns are
  # connected to every living member and to each other, and nothing is left of
  # the dead
  expect_equal(pairs_of(dat$run$el[[2]]),
               pairs_of(EpiModel:::group_edgelist(hh1)))
})

test_that("nwupdate_gmc19 stops when an attribute was not set for newborns", {
  dat <- dat_at_start()
  dat <- set_current_timestep(dat, 2)
  dat <- append_core_attr(dat, 2, 1)
  expect_error(nwupdate_gmc19(dat, 2), "unequal lengths")
})
