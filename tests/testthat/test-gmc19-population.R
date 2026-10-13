# Aging, sex and occupation class, and the class multipliers of the relationship
# activities (GlobalMix India #61)

bands <- c("0-9y", "10-19y", "20-29y", "30-39y", "40-59y", "60+y")
act_pars <- list(o = c(0.9, 0.95, 0.9, 0.85, 0.9, 0.8), kappa = 5, q = rep(0.75, 6),
                 p11 = 0.5, s = 0.05, pi = c(0.1, 0, 0, 0.2, 0, 0.3), theta_F = rep(2, 6),
                 theta_I = 0.3, rho = -0.5, new_rate = rep(0.1, 6), new_share = 0.5,
                 new_mix = matrix(1, 6, 6) + diag(4, 6))

# Multipliers for every class in the population below and for newborns
classes <- expand.grid(band = 1:6, sex = c("F", "M"), occ = c("child", "student", "worker", "home"),
                       stringsAsFactors = FALSE)
classes$mult_F <- ifelse(classes$occ == "worker", 1.5, 0.8)
classes$mult_I <- ifelse(classes$sex == "M", 1.2, 0.9)

N <- 600
set.seed(5)
band <- rep(1:6, length.out = N)
sex <- rep(c("F", "M"), each = N / 2)
occ <- c("child", "student", "worker", "home")[c(1, 2, 3, 3, 4, 4)][band]
tr <- draw_relationship_traits(band, act_pars, seed = 6)

# A small relationship model: the two relationship networks and a household
# clique layer, with or without the population's sex and occupation class
pop_nets <- function(with_class = TRUE) {
  nw <- network_initialize(n = N)
  nw <- set_vertex_attribute(nw, "age", c(5, 15, 25, 35, 50, 70)[band] + 0.5)
  nw <- set_vertex_attribute(nw, "age.grp", bands[band])
  nw <- set_vertex_attribute(nw, "never", tr$never)
  nw <- set_vertex_attribute(nw, "log_aF", tr$log_aF)
  nw <- set_vertex_attribute(nw, "log_aI", tr$log_aI)
  if (with_class) {
    nw <- set_vertex_attribute(nw, "sex", sex)
    nw <- set_vertex_attribute(nw, "occ", occ)
  }
  fit <- function(attr, edges) {
    f <- as.formula(paste0("~edges + offset(nodecov(\"", attr, "\")) + ",
                           "offset(nodefactor(\"never\", levels = I(1)))"))
    est <- suppressMessages(suppressWarnings(netest(nw, formation = f, target.stats = edges,
                                                    coef.form = c(1, -Inf),
                                                    coef.diss = dissolution_coefs(~offset(edges), 100),
                                                    verbose = FALSE)))
    est$cell_theta <- matrix(est$coef.form.crude[["edges"]], 6, 6)
    est
  }
  nw_hh <- network_initialize(n = N)
  nw_hh <- set_vertex_attribute(nw_hh, "hh.ids", rep(seq_len(N / 4), each = 4))
  list(freq = fit("log_aF", 900), infreq = fit("log_aI", 600),
       home = netclique(nw_hh, "hh.ids", arrivals = "join", arrivals.FUN = hh_join_young))
}
nets_class <- pop_nets(TRUE)
nets_plain <- pop_nets(FALSE)

pop_param <- function(...) {
  param.net(
    time.unit = 1, inf.prob = rep(0.1, 4), exposure.rate = rep(1, 4), inf.prob.a.rr = 0.5,
    layer.names = c("freq", "infreq", "new", "home"),
    age.breaks = c(0, 10, 20, 30, 40, 60, Inf), age.grps = bands, arrival.age = 0,
    birth.prop.male = 0.52,
    act.o = act_pars$o, act.kappa = act_pars$kappa, act.q = act_pars$q, act.p11 = act_pars$p11,
    act.s = act_pars$s, act.pi = act_pars$pi, act.theta.F = act_pars$theta_F,
    act.theta.I = act_pars$theta_I, act.rho = act_pars$rho, new.rate = act_pars$new_rate,
    new.mix = act_pars$new_mix, new.share.norel = act_pars$new_share, ...)
}

# init_gmc19() on the small model, as netsim() calls it
pop_init <- function(x = nets_class, param = pop_param()) {
  init <- init.net(e.num = 0)
  control <- control.net(type = NULL, nsteps = 5, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE,
                         initialize.FUN = init_gmc19, aging.FUN = aging)
  EpiModel:::crosscheck.net(x, param, init, control)
  control <- EpiModel:::netsim_validate_control(control)
  init_gmc19(x, param, init, control, s = 1)
}

# Newborns as the arrivals module adds them
add_newborns <- function(dat, n_new, at = 2) {
  dat <- set_current_timestep(dat, at)
  EpiModelCOVID:::init_new_nodes_attrs(dat, n_new)
}

test_that("aging adds one day, 1/365 of a year, per step", {
  dat <- pop_init()
  age <- get_attr(dat, "age")
  dat <- aging(dat, 2)
  expect_equal(get_attr(dat, "age") - age, rep(1 / 365, N))
  # an age a day short of a band's top crosses it on the step
  dat <- set_attr(dat, "age", rep(10 - 0.5 / 365, N))
  dat <- aging(dat, 3)
  expect_true(all(get_attr(dat, "age.grp") == "10-19y"))
})

test_that("the population's sex and occupation class reach the attributes", {
  dat <- pop_init()
  expect_identical(get_attr(dat, "sex"), sex)
  expect_identical(get_attr(dat, "occ"), occ)
  # and are absent when the population carries neither
  dat <- pop_init(nets_plain)
  expect_false(any(c("sex", "occ") %in% names(get_attr_list(dat))))
})

test_that("newborns are children, male with probability birth.prop.male", {
  dat <- add_newborns(pop_init(), 40000)
  new <- N + seq_len(40000)
  expect_equal(mean(get_attr(dat, "sex")[new] == "M"), 0.52, tolerance = 0.02)
  expect_true(all(get_attr(dat, "sex")[new] %in% c("F", "M")))
  expect_true(all(get_attr(dat, "occ")[new] == "child"))
  expect_identical(get_attr(dat, "sex")[seq_len(N)], sex)
})

test_that("a population without sex and occ draws nothing for them", {
  dat <- pop_init(nets_plain)
  dat <- set_current_timestep(dat, 2)
  set.seed(21)
  out <- EpiModelCOVID:::init_new_nodes_attrs(dat, 50)
  expect_false(any(c("sex", "occ") %in% names(get_attr_list(out))))
  # the newborns' draws are the traits, the propensity to go out and the
  # new-contact propensity, in that order, as before #61
  pars <- EpiModelCOVID:::activation_params(dat)
  set.seed(21)
  t1 <- draw_relationship_traits(rep(1L, 50), pars)
  o1 <- draw_out_propensity(rep(1L, 50), pars$o, pars$kappa)
  new <- N + seq_len(50)
  expect_identical(get_attr(out, "log_aF")[new], t1$log_aF)
  expect_identical(get_attr(out, "log_aI")[new], t1$log_aI)
  expect_identical(get_attr(out, "o_i")[new], o1)
})

test_that("the population's class attributes are checked at initialization", {
  expect_error(pop_init(param = pop_param(birth.prop.male = NULL)), "birth.prop.male")
  expect_error(pop_init(param = pop_param(birth.prop.male = 1.2)), "birth.prop.male")
  bad <- nets_class
  bad$freq$newnetwork <- set_vertex_attribute(bad$freq$newnetwork, "sex", rep("X", N))
  expect_error(pop_init(bad), "sex must be F or M; the population has X")
  one <- nets_plain
  one$freq$newnetwork <- set_vertex_attribute(one$freq$newnetwork, "occ", occ)
  expect_error(pop_init(one), "carries occ but not sex")
  expect_error(pop_init(nets_plain, pop_param(act.class.mult = classes)), "needs each person's sex and occ")
})

test_that("class multipliers scale each activity's mean, shape unchanged", {
  b <- rep(3L, 80000)
  s <- rep(c("F", "M"), 40000)
  o <- rep(c("worker", "worker", "home", "home"), 20000)
  pars <- act_pars
  pars$class_mult <- EpiModelCOVID:::normalize_class_mult(classes, bands)
  x <- draw_relationship_traits(b, pars, seed = 2, sex = s, occ = o)
  y <- draw_relationship_traits(b, act_pars, seed = 2)
  m <- EpiModelCOVID:::class_multipliers(b, s, o, pars$class_mult)
  # the same draws, times the class's multiplier (away from the floor on a
  # variate near zero)
  expect_identical(x$never, y$never)
  expect_equal(x$log_aF - y$log_aF, log(m$mult_F))
  ok <- y$log_aI > log(1e-6)
  expect_equal((x$log_aI - y$log_aI)[ok], log(m$mult_I)[ok])
  # so each class's mean activity is its multiplier
  aF <- tapply(exp(x$log_aF), paste(s, o), mean)
  expect_equal(as.numeric(aF[c("F home", "M worker")]), c(0.8, 1.5), tolerance = 0.03)
  aI <- tapply(exp(x$log_aI), s, mean)
  expect_equal(as.numeric(aI), c(0.9, 1.2), tolerance = 0.1)
})

test_that("without class multipliers the traits are those drawn before #61", {
  b <- rep(1:6, 500)
  s <- rep(c("F", "M"), 1500)
  o <- rep("worker", 3000)
  ref <- draw_relationship_traits(b, act_pars, seed = 3)
  expect_identical(draw_relationship_traits(b, act_pars, seed = 3, sex = s, occ = o), ref)
  set.seed(3)
  expect_identical(draw_relationship_traits(b, act_pars), ref)
})

test_that("a class missing from the multipliers is an error naming it", {
  pars <- act_pars
  pars$class_mult <- EpiModelCOVID:::normalize_class_mult(
    classes[!(classes$band == 2 & classes$sex == "M" & classes$occ == "home"), ], bands)
  expect_error(draw_relationship_traits(c(1L, 2L), pars, sex = c("F", "M"), occ = c("child", "home")),
               "no row for band 2, sex M, occ home")
  expect_error(draw_relationship_traits(1:2, pars), "need each node's sex and occ")
  # newborns' classes are checked at initialization
  no_child <- classes[!(classes$band == 1 & classes$occ == "child"), ]
  expect_error(pop_init(param = pop_param(act.class.mult = no_child)),
               "no row for band 1, sex F, occ child; band 1, sex M, occ child")
  dup <- rbind(classes, classes[1, ])
  expect_error(pop_init(param = pop_param(act.class.mult = dup)), "more than one row")
  neg <- classes; neg$mult_F[3] <- 0
  expect_error(pop_init(param = pop_param(act.class.mult = neg)), "must be positive")
})

test_that("newborns draw their traits with their own class", {
  dat <- pop_init(param = pop_param(act.class.mult = classes))
  n_new <- 40000
  dat <- add_newborns(dat, n_new)
  new <- N + seq_len(n_new)
  s <- get_attr(dat, "sex")[new]
  aI <- tapply(exp(get_attr(dat, "log_aI")[new]), s, mean)
  # newborns are children, so mult_F is 0.8 for both sexes and mult_I is by sex
  expect_equal(mean(exp(get_attr(dat, "log_aF")[new])), 0.8, tolerance = 0.03)
  expect_equal(as.numeric(aI), c(0.9, 1.2), tolerance = 0.1)
})

test_that("class multipliers must be those the networks were fitted with", {
  fitted <- nets_class
  fitted$freq$trait_pars <- list(pi = act_pars$pi, theta_F = act_pars$theta_F,
                                 theta_I = act_pars$theta_I, rho = act_pars$rho,
                                 class_mult = classes)
  # the same table, in another row order and with band labels, passes
  shuffled <- classes[rev(seq_len(nrow(classes))), ]
  shuffled$band <- bands[shuffled$band]
  expect_silent(pop_init(fitted, pop_param(act.class.mult = shuffled)))
  other <- classes; other$mult_F[1] <- 2
  expect_error(pop_init(fitted, pop_param(act.class.mult = other)), "\\(class_mult\\) differ")
  expect_error(pop_init(fitted), "\\(class_mult\\) differ")
  # and none on the run against none on the fit
  fitted$freq$trait_pars$class_mult <- NULL
  expect_error(pop_init(fitted, pop_param(act.class.mult = classes)), "\\(class_mult\\) differ")
  expect_silent(pop_init(fitted))
})
