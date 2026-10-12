# Relationships, daily activation and new contacts (GlobalMix India #63)

bands <- c("0-9y", "10-19y", "20-29y", "30-39y", "40-59y", "60+y")
act_pars <- list(o = c(0.9, 0.95, 0.9, 0.85, 0.9, 0.8), kappa = 5, q = rep(0.75, 6),
                 p11 = 0.5, s = 0.05, pi = c(0.1, 0, 0, 0.2, 0, 0.3), theta_F = rep(2, 6),
                 theta_I = 0.3, rho = -0.5, new_rate = rep(0.1, 6), new_share = 0.5,
                 new_mix = matrix(1, 6, 6) + diag(4, 6), episodes = TRUE, window = 30)

# A small relationship model: frequent and infrequent TERGMs with the activity and
# never-active offsets, and a household clique layer
N <- 600
set.seed(3)
band <- rep(1:6, length.out = N)
tr <- draw_relationship_traits(band, act_pars, seed = 4)
rel_nw <- function() {
  nw <- network_initialize(n = N)
  nw <- set_vertex_attribute(nw, "age", c(5, 15, 25, 35, 50, 70)[band] + 0.5)
  nw <- set_vertex_attribute(nw, "age.grp", bands[band])
  nw <- set_vertex_attribute(nw, "never", tr$never)
  nw <- set_vertex_attribute(nw, "log_aF", tr$log_aF)
  nw <- set_vertex_attribute(nw, "log_aI", tr$log_aI)
  nw
}
rel_fit <- function(attr, edges) {
  f <- as.formula(paste0("~edges + offset(nodecov(\"", attr, "\")) + offset(nodefactor(\"never\", levels = I(1)))"))
  suppressMessages(suppressWarnings(netest(rel_nw(), formation = f, target.stats = edges,
                                           coef.form = c(1, -Inf),
                                           coef.diss = dissolution_coefs(~offset(edges), 100),
                                           verbose = FALSE)))
}
est_F <- rel_fit("log_aF", 900)
est_I <- rel_fit("log_aI", 600)
est_F$cell_theta <- matrix(est_F$coef.form.crude[["edges"]], 6, 6)
est_I$cell_theta <- matrix(est_I$coef.form.crude[["edges"]], 6, 6)
nw_hh <- network_initialize(n = N)
nw_hh <- set_vertex_attribute(nw_hh, "hh.ids", rep(seq_len(N / 4), each = 4))
est_H <- netclique(nw_hh, "hh.ids")

engine_params <- function(dat) {
  dat <- set_param(dat, "act.o", act_pars$o)
  dat <- set_param(dat, "act.kappa", act_pars$kappa)
  dat <- set_param(dat, "act.q", act_pars$q)
  dat <- set_param(dat, "act.p11", act_pars$p11)
  dat <- set_param(dat, "act.s", act_pars$s)
  dat <- set_param(dat, "act.pi", act_pars$pi)
  dat <- set_param(dat, "act.theta.F", act_pars$theta_F)
  dat <- set_param(dat, "act.theta.I", act_pars$theta_I)
  dat <- set_param(dat, "act.rho", act_pars$rho)
  dat <- set_param(dat, "new.rate", act_pars$new_rate)
  dat <- set_param(dat, "new.mix", act_pars$new_mix)
  dat <- set_param(dat, "new.share.norel", act_pars$new_share)
  dat
}

# A tergmLite netsim_dat at t = 1 with the engine's attributes and the activation
# parameters, as init_gmc19() leaves it (without the disease modules)
rel_dat <- function(layer.names = c("freq", "infreq", "new", "home"), activation = TRUE,
                    x = list(freq = est_F, infreq = est_I, home = est_H)) {
  param <- param.net(inf.prob = 0.1, act.rate = 1)
  init <- init.net(i.num = 5, e.num = 0)
  control <- control.net(type = "SI", nsteps = 30, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE)
  EpiModel:::crosscheck.net(x, param, init, control)
  control <- EpiModel:::netsim_validate_control(control)
  dat <- initialize.net(x, param, init, control, s = 1)
  dat <- set_param(dat, "layer.names", layer.names)
  dat <- set_param(dat, "inf.prob", rep(0.1, length(layer.names)))
  dat <- set_param(dat, "exposure.rate", rep(1, length(layer.names)))
  dat <- set_param(dat, "inf.prob.a.rr", 0.5)
  dat <- set_param(dat, "age.breaks", c(0, 10, 20, 30, 40, 60, Inf))
  dat <- set_param(dat, "age.grps", bands)
  dat <- set_param(dat, "arrival.age", 0)
  if (activation) dat <- engine_params(dat)
  dat <- EpiModelCOVID:::check_layer_params(dat)
  dat <- EpiModelCOVID:::init_attrs(dat)
  dat <- set_attr(dat, "status", rep("s", N))
  dat
}

test_that("relationship traits follow their parameters", {
  b <- rep(1:6, each = 20000)
  x <- draw_relationship_traits(b, act_pars, seed = 1)
  expect_equal(as.numeric(tapply(x$never, b, mean)), act_pars$pi, tolerance = 0.02)
  expect_equal(mean(exp(x$log_aF)), 1, tolerance = 0.03)
  expect_equal(mean(exp(x$log_aI)), 1, tolerance = 0.1)
  expect_lt(stats::cor(x$log_aF, x$log_aI, method = "spearman"), -0.3)
  o <- draw_out_propensity(b, act_pars$o, 5)
  expect_equal(as.numeric(tapply(o, b, mean)), act_pars$o, tolerance = 0.01)
  expect_identical(draw_out_propensity(b, act_pars$o, 1e4), act_pars$o[b])
})

test_that("episodes are stationary at s and persist with probability p11", {
  key <- as.numeric(1:20000)
  st <- EpiModelCOVID:::episode_update(key, NULL, 0.6, 0.1)
  on <- matrix(FALSE, length(key), 40)
  for (d in 1:40) {
    st <- EpiModelCOVID:::episode_update(key, list(key = key, on = st), 0.6, 0.1)
    on[, d] <- st
  }
  expect_equal(mean(on), 0.1, tolerance = 0.03)
  expect_equal(sum(on[, -1] & on[, -40]) / sum(on[, -40]), 0.6, tolerance = 0.02)
  # A tie not seen yesterday starts at the stationary probability, and with
  # episodes off every day is independent
  expect_equal(mean(EpiModelCOVID:::episode_update(key, list(key = -1, on = TRUE), 0.6, 0.1)), 0.1,
               tolerance = 0.01)
  flat <- EpiModelCOVID:::episode_update(key, list(key = key, on = rep(TRUE, length(key))), 0.6, 0.1,
                                         episodes = FALSE)
  expect_equal(mean(flat), 0.1, tolerance = 0.01)
})

test_that("new contacts have the expected number of edge ends by band", {
  N_b <- c(1000, 2000, 1500, 1500, 3000, 1000)
  E <- new_contact_edges(act_pars$new_mix, act_pars$new_rate, N_b)
  expect_equal(rowSums(E) + diag(E), act_pars$new_rate * N_b, tolerance = 1e-8)
  b <- rep(1:6, N_b)
  w <- rep(1, length(b)); w[seq(1, length(b), by = 2)] <- 0
  ends <- 0
  for (r in 1:50) {
    el <- draw_new_contacts(E, b, w)
    expect_true(all(el[, 1] < el[, 2]))
    expect_true(all(w[el] > 0))
    ends <- ends + tabulate(b[el], 6)
  }
  expect_equal(ends / 50, act_pars$new_rate * N_b, tolerance = 0.05)
})

test_that("the day's contacts are relationships met by two people out, plus new contacts", {
  dat <- rel_dat()
  dat <- set_current_timestep(dat, 2)
  dat <- activate_contacts(dat, 2)
  ct <- dat$run$contacts
  expect_identical(names(ct), c("freq", "infreq", "new", "home"))
  expect_identical(dat$run$contacts_at, 2)
  key <- function(el) paste(el[, 1], el[, 2])
  expect_true(all(key(ct$freq) %in% key(dat$run$el$freq)))
  expect_true(all(key(ct$infreq) %in% key(dat$run$el$infreq)))
  expect_identical(ct$home, dat$run$el$home)
  never <- get_attr(dat, "never")
  expect_false(any(never[c(ct$freq, ct$infreq)] == 1))
  # half of the new-contact ends go to people with no relationship met today
  dat2 <- dat
  ends <- 0; norel <- 0
  for (t in 3:12) {
    dat2 <- set_current_timestep(dat2, t)
    dat2 <- activate_contacts(dat2, t)
    ends <- ends + get_epi(dat2, "act.new.ends", t); norel <- norel + get_epi(dat2, "act.new.norel", t)
  }
  expect_equal(norel / ends, 0.5, tolerance = 0.1)
  # the episode state is carried for exactly today's infrequent ties
  expect_length(dat$run$episode$on, nrow(dat$run$el$infreq))
  expect_equal(get_attr(dat, "act.contact.days"),
               tabulate(c(ct$freq, ct$infreq, ct$new), nbins = N))
  dat <- set_current_timestep(dat, 3)
  dat <- activate_contacts(dat, 3)
  expect_false(is.null(get_epi(dat, "act.both.F.den", 3)))
})

test_that("frequent relationships are met at q when both people are out", {
  dat <- rel_dat()
  dat <- set_attr(dat, "o_i", rep(1, N))
  met <- 0; ties <- 0
  for (t in 2:21) {
    dat <- set_current_timestep(dat, t)
    dat <- activate_contacts(dat, t)
    eF <- dat$run$el$freq
    ok <- get_attr(dat, "never")[eF[, 1]] == 0 & get_attr(dat, "never")[eF[, 2]] == 0
    met <- met + nrow(dat$run$contacts$freq); ties <- ties + sum(ok)
  }
  expect_equal(met / ties, 0.75, tolerance = 0.03)
})

test_that("infection runs on the day's contacts, by tie type, with the age-pair weights", {
  dat <- rel_dat()
  dat <- set_current_timestep(dat, 2)
  dat <- activate_contacts(dat, 2)
  status <- rep("s", N); status[seq(1, N, by = 3)] <- "ic"
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "vax", rep(0, N))
  dat <- set_attr(dat, "ever.inf", rep(0, N))
  dat <- set_attr(dat, "last.dose.time", rep(NA, N))
  dat <- set_attr(dat, "is_infant", rep(FALSE, N))
  dat <- set_param(dat, "inf.prob", c(0.1, 0.2, 0.3, 0.4))
  W <- matrix(2, 6, 6)
  dat <- set_param(dat, "contact.weight", list(freq = W, infreq = W, new = W, home = W / 2))
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
  dat <- infect_general(dat, 2)
  # freq 0.1 and infreq 0.2 with weight 2 (two acts), new 0.3 with weight 2,
  # home 0.4 with weight 1
  p <- unique(round(unlist(probs), 10))
  expect_true(all(p %in% round(c(1 - 0.9^2, 1 - 0.8^2, 1 - 0.7^2, 0.4), 10)))
  expect_true(all(round(c(1 - 0.9^2, 0.4), 10) %in% p))
  expect_false(is.null(get_epi(dat, "se.flow.new", 2)))
  # the equal-intensity sensitivity removes the weights
  dat <- set_param(dat, "contact.weight.equal", TRUE)
  expect_null(EpiModelCOVID:::contact_weights(dat))
})

test_that("contacts must come from this step's activation", {
  local_mocked_bindings(
    build_vax_schedule = function(dat, disease) NULL,
    vax_age_group_for = function(dat) rep(NA_integer_, length(get_attr(dat, "active"))),
    .package = "EpiModelCOVID"
  )
  dat <- rel_dat()
  dat <- set_current_timestep(dat, 2)
  dat <- activate_contacts(dat, 2)
  expect_error(infect_general(dat, 3), "no contacts for step 3")
})

test_that("relationship networks without the activation step stop the run", {
  expect_error(rel_dat(activation = FALSE), "need the activation")
  expect_error(rel_dat(layer.names = c("freq", "infreq", "home")), "new-contact")
})

test_that("degree scores rank on the named score", {
  dat <- rel_dat()
  dat <- set_current_timestep(dat, 2)
  dat <- activate_contacts(dat, 2)
  dat <- set_param(dat, "vax.degree.score", "relationship")
  dat <- netdegree(dat, 2)
  expect_equal(get_attr(dat, "degree_score"), tabulate(dat$run$el$freq, nbins = N))
  dat <- set_param(dat, "vax.degree.score", "active")
  dat <- netdegree(dat, 2)
  expect_true(all(get_attr(dat, "degree_score") <= get_attr(dat, "degree_freq") +
                    get_attr(dat, "degree_infreq") + 1e-12))
  dat <- set_param(dat, "vax.degree.score", "nonsense")
  expect_error(netdegree(dat, 2), "unknown vax.degree.score")
})

test_that("newborns get relationships at arrival and the edgelists stay sorted", {
  dat <- rel_dat()
  dat <- set_current_timestep(dat, 2)
  n_new <- 20
  dat <- append_core_attr(dat, 2, n_new)
  for (a in names(get_attr_list(dat))) {
    v <- get_attr(dat, a)
    if (length(v) < N + n_new) {
      val <- switch(a, age.grp = "0-9y", status = "s", never = 0, log_aF = 0.5, log_aI = 0.5, NA)
      dat <- append_attr(dat, a, val, N + n_new - length(v))
    }
  }
  dat <- EpiModel::nwupdate.net(dat, 2)
  before <- nrow(dat$run$el$freq)
  dat <- EpiModelCOVID:::seed_newborn_ties(dat, 2)
  el <- dat$run$el$freq
  expect_gt(nrow(el), before)
  expect_true(all(el[, 1] < el[, 2]))
  expect_false(is.unsorted(el[, 1] * (N + n_new + 1) + el[, 2]))
  expect_false(anyDuplicated(paste(el[, 1], el[, 2])) > 0)
  expect_gt(sum(tabulate(el, nbins = N + n_new)[N + seq_len(n_new)]), 0)
  expect_identical(attr(el, "n"), attr(dat$run$el$infreq, "n"))
})
