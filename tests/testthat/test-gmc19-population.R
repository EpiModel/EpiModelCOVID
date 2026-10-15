# Aging, sex and occupation class, and the class multipliers of the relationship
# activities (GlobalMix India #61); births and background deaths (#64)

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
    birth.prop.male = 0.52, birth.rate = 12 / 1000 / 365, mort.rates = rep(1e-4, 86),
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

# Births and background deaths (GlobalMix India #64) ------------------------------

test_that("births are binomial on the living at birth.rate, whatever the deaths", {
  dat <- pop_init(param = pop_param(birth.rate = 0.2))
  dat <- set_current_timestep(dat, 2)
  # this step's dead are marked but still in the attributes until nwupdate.FUN
  dat <- set_attr(dat, "active", 0, posit_ids = 1:100)
  dat <- set_epi(dat, "d.flow", 100, at = 2)
  set.seed(31)
  out <- arrival(dat, 2)
  set.seed(31)
  n_new <- rbinom(1, N - 100, 0.2)
  expect_equal(get_epi(out, "nNew", at = 2), n_new)
  expect_length(get_attr(out, "active"), N + n_new)
  # the same draw with no deaths marked: births do not follow d.flow
  dat <- set_epi(dat, "d.flow", 0, at = 2)
  set.seed(31)
  expect_equal(get_epi(arrival(dat, 2), "nNew", at = 2), n_new)
  # a closed population has none
  dat <- set_param(dat, "birth.rate", 0)
  expect_equal(get_epi(arrival(dat, 2), "nNew", at = 2), 0L)
})

test_that("without birth.rate births replace deaths, as on the legacy model", {
  dat <- pop_init()
  dat <- set_current_timestep(dat, 2)
  dat$param$birth.rate <- NULL
  dat <- set_epi(dat, "d.flow", 7, at = 2)
  expect_equal(get_epi(arrival(dat, 2), "nNew", at = 2), 7L)
})

test_that("the relationship design needs birth.rate, and both rates are checked", {
  expect_error(pop_init(param = pop_param(birth.rate = NULL)), "needs birth.rate")
  for (bad in list(-0.1, 1, c(0.1, 0.2), NA_real_, "0.1")) {
    expect_error(pop_init(param = pop_param(birth.rate = bad)), "birth.rate must be one daily probability")
  }
  for (bad in list(c(1e-4, NA), c(1e-4, -1e-4), 2, numeric(0))) {
    expect_error(pop_init(param = pop_param(mort.rates = bad)), "mort.rates must be daily death probabilities")
  }
  expect_silent(pop_init(param = pop_param(birth.rate = 0, mort.rates = rep(0, 86))))
})

test_that("background mortality is indexed by single year of age, floored at 1 and capped", {
  dat <- pop_init(param = pop_param(mort.rates = c(1, rep(0, 84), 1)))
  dat <- set_current_timestep(dat, 2)
  age <- rep(5.5, N)
  age[1] <- 0          # index 1, not dropped from the index
  age[2] <- 0.4        # index 1
  age[3] <- 85.2       # index 86
  age[4] <- 120        # beyond the vector: its last element
  dat <- set_attr(dat, "age", age)
  out <- deaths_covid_gmc19(dat, 2)
  expect_identical(which(get_attr(out, "active") == 0), 1:4)
  expect_identical(get_attr(out, "exitTime")[1:4], rep(2, 4))
  expect_equal(get_epi(out, "d.flow", at = 2), 4)
  # a vector with no rate for an age stops the run rather than dropping the node
  dat <- set_param(dat, "mort.rates", c(1e-4, NA))
  expect_error(deaths_covid_gmc19(dat, 2), "no rate for some living node")
})

test_that("aging keeps the infant flag on age", {
  dat <- pop_init()
  age <- get_attr(dat, "age")
  age[1:3] <- c(0.2, 1 - 0.5 / 365, 1.5)
  dat <- set_attr(dat, "age", age)
  dat <- set_attr(dat, "is_infant", age < 1 | seq_len(N) == 3)     # a stale flag on node 3
  dat <- set_attr(dat, "active", 0, posit_ids = 4)
  dat <- set_attr(dat, "is_infant", TRUE, posit_ids = 4)
  dat <- aging(dat, 2)
  inf <- get_attr(dat, "is_infant")
  expect_identical(inf[1:3], c(TRUE, FALSE, FALSE))
  alive <- get_attr(dat, "active") == 1
  expect_identical(inf[alive], get_attr(dat, "age")[alive] < 1)
  # the dead keep the flag they died with
  expect_true(inf[4])
})

# Band tops (GlobalMix India #64) ----------------------------------------------------

# Trait parameters that differ by band, so a crossing has something to map
cross_pars <- function() {
  p <- act_pars
  p$theta_F <- c(2, 4, 1, 3, 2, 1.5)
  p$pi <- c(0.1, 0, 0, 0.2, 0, 0.3)
  p$o <- c(0.9, 0.95, 0.9, 0.7, 0.9, 0.8)
  cm <- classes
  cm$mult_F <- cm$mult_F * c(1, 1.2, 0.9, 1.1, 1, 0.8)[cm$band]
  p$class_mult <- EpiModelCOVID:::normalize_class_mult(cm, bands)
  p
}

test_that("traits keep their values within a band and their rank across a band top", {
  p <- cross_pars()
  n <- 60000
  s <- rep(c("F", "M"), n / 2)
  o <- rep(c("worker", "home", "student"), n / 3)
  draw <- function(b, seed) {
    tr <- draw_relationship_traits(rep(b, n), p, seed = seed, sex = s, occ = o)
    tr$o_i <- draw_out_propensity(rep(b, n), p$o, p$kappa)
    tr
  }
  x <- draw(3, 4)
  # no crossing, no change
  expect_identical(remap_relationship_traits(rep(3L, n), rep(3L, n), x, p, s, o), x)

  # 3 -> 4: the mapped traits have band 4's distributions...
  set.seed(8)
  y <- remap_relationship_traits(rep(3L, n), rep(4L, n), x, p, s, o)
  z <- draw(4, 5)
  expect_equal(mean(y$never), p$pi[4], tolerance = 0.02)
  expect_true(all(y$never[x$never == 1] == 1))          # pi rises, so nobody becomes active
  for (k in unique(paste(s, o))) {
    j <- paste(s, o) == k
    expect_equal(mean(exp(y$log_aF[j])), mean(exp(z$log_aF[j])), tolerance = 0.03, info = k)
    expect_equal(var(exp(y$log_aF[j])) / mean(exp(y$log_aF[j]))^2, 1 / p$theta_F[4], tolerance = 0.06, info = k)
    # ...and each person keeps their rank within their class
    expect_equal(cor(y$log_aF[j], x$log_aF[j], method = "spearman"), 1, tolerance = 1e-6, info = k)
  }
  expect_identical(y$log_aI, x$log_aI)                  # theta_I and mult_I do not change
  expect_equal(mean(y$o_i), p$o[4], tolerance = 0.005)
  expect_equal(cor(y$o_i, x$o_i, method = "spearman"), 1, tolerance = 1e-6)

  # 4 -> 5: pi falls to 0, so every never-active person becomes active
  y5 <- remap_relationship_traits(rep(4L, n), rep(5L, n), y, p, s, o)
  expect_identical(sum(y5$never), 0L)
  # with no variation in going out, the new band's mean
  p$kappa <- 1e4
  expect_identical(remap_relationship_traits(rep(4L, n), rep(5L, n), y, p, s, o)$o_i, rep(p$o[5], n))
})

# Ties of a set of nodes in a relationship network's edgelist
ties_of <- function(el, ids) el[el[, 1] %in% ids | el[, 2] %in% ids, , drop = FALSE]
key_of <- function(el) paste(el[, 1], el[, 2])

test_that("a person crossing a band top takes the new band's traits and relationships", {
  dat <- pop_init(param = pop_param(act.class.mult = classes))
  age <- get_attr(dat, "age")
  cross <- which(age < 10)[1:20]
  age[cross] <- 10 - 0.5 / 365
  dat <- set_attr(dat, "age", age)
  el0 <- dat$run$el
  tr0 <- lapply(c(never = "never", log_aF = "log_aF", o_i = "o_i"), function(a) get_attr(dat, a))
  dat <- set_current_timestep(dat, 2)
  set.seed(12)
  dat1 <- aging(dat, 2)
  expect_identical(which(get_attr(dat1, "rel.cross") == 2), cross)
  expect_true(all(get_attr(dat1, "age.grp")[cross] == "10-19y"))
  expect_equal(get_epi(dat1, "n.band.cross", at = 2), 20)
  # only the crossers' traits move; band 1 and 2 share theta_F and the class
  # multipliers here, so the activity is mapped only through o and pi
  others <- setdiff(seq_len(N), cross)
  for (a in names(tr0)) expect_identical(get_attr(dat1, a)[others], tr0[[a]][others], info = a)
  expect_true(all(get_attr(dat1, "o_i")[cross] != tr0$o_i[cross]))
  expect_identical(get_attr(dat1, "act.eligible"), 1 - get_attr(dat1, "never"))
  # no one else's ties move, and the edgelists stay sorted and unique
  for (L in c("freq", "infreq")) {
    el <- dat1$run$el[[L]]
    expect_identical(key_of(el[!(el[, 1] %in% cross | el[, 2] %in% cross), , drop = FALSE]),
                     key_of(el0[[L]][!(el0[[L]][, 1] %in% cross | el0[[L]][, 2] %in% cross), , drop = FALSE]))
    expect_true(all(el[, 1] < el[, 2]))
    expect_false(anyDuplicated(key_of(el)) > 0)
    expect_identical(order(el[, 1], el[, 2]), seq_len(nrow(el)))
  }
  # the newborn draw leaves the crossers' ties alone
  expect_identical(EpiModelCOVID:::seed_relationship_ties(dat1, 2)$run$el, dat1$run$el)
})

test_that("a mover's ties change only as far as the new cross-section needs", {
  dat <- pop_init()
  dat <- set_current_timestep(dat, 2)
  band <- EpiModelCOVID:::band_index(dat)
  ids <- which(band == 1 & get_attr(dat, "never") == 0)[1:15]
  old <- list(never = get_attr(dat, "never")[ids], log_aF = get_attr(dat, "log_aF")[ids],
              log_aI = get_attr(dat, "log_aI")[ids])
  # no change in band or traits: nothing moves
  expect_identical(EpiModelCOVID:::move_relationship_ties(dat, ids, band[ids], old)$run$el, dat$run$el)

  # the movers enter band 3 with a new activity; the expected ties of each
  # mover after the move, given the ties before it, against the simulated mean
  grp <- get_attr(dat, "age.grp"); grp[ids] <- bands[3]
  dat <- set_attr(dat, "age.grp", grp)
  dat <- set_attr(dat, "log_aF", get_attr(dat, "log_aF")[ids] + log(1.5), posit_ids = ids)
  band1 <- EpiModelCOVID:::band_index(dat)
  la0 <- get_attr(dat, "log_aF"); la0[ids] <- old$log_aF
  la1 <- get_attr(dat, "log_aF")
  elig <- which(get_attr(dat, "never") == 0)
  cell <- unname(nets_class$freq$cell_theta)
  el0 <- dat$run$el$freq
  expected <- vapply(ids, function(i) {
    j <- setdiff(elig, i)
    p0 <- plogis(cell[cbind(rep(1, length(j)), band1[j])] + la0[i] + la0[j])
    p0[j %in% ids] <- plogis(cell[1, 1] + la0[i] + la0[j[j %in% ids]])
    p1 <- plogis(cell[cbind(rep(3, length(j)), band1[j])] + la1[i] + la1[j])
    tied <- j %in% c(el0[el0[, 1] == i, 2], el0[el0[, 2] == i, 1])
    sum(ifelse(tied, pmin(1, p1 / p0), pmax(0, p1 - p0) / (1 - p0)))
  }, numeric(1))
  sims <- replicate(300, {
    out <- EpiModelCOVID:::move_relationship_ties(dat, ids, rep(1L, length(ids)), old)$run$el$freq
    tabulate(c(out[out[, 1] %in% ids, 1], out[out[, 2] %in% ids, 2]), nbins = N)[ids]
  })
  # pairs of two movers count for both, so compare the total
  expect_equal(sum(rowMeans(sims)), sum(expected), tolerance = 0.03)
})

test_that("with act.reset.cross FALSE a crosser keeps their traits and ties", {
  dat <- pop_init(param = pop_param(act.reset.cross = FALSE))
  age <- get_attr(dat, "age")
  cross <- which(age < 10)[1:20]
  age[cross] <- 10 - 0.5 / 365
  dat <- set_attr(dat, "age", age)
  el0 <- dat$run$el
  o0 <- get_attr(dat, "o_i")
  dat <- set_current_timestep(dat, 2)
  dat <- aging(dat, 2)
  expect_true(all(get_attr(dat, "age.grp")[cross] == "10-19y"))
  expect_true(all(is.na(get_attr(dat, "rel.cross"))))
  expect_identical(get_attr(dat, "o_i"), o0)
  dat <- EpiModelCOVID:::seed_relationship_ties(dat, 2)
  expect_identical(dat$run$el, el0)
})

# Occupation class by age (GlobalMix India #64) --------------------------------------

# Class multipliers in which a student is more active than a child in band 1
cross_mult <- function() {
  cm <- classes
  cm$mult_F[cm$occ == "student"] <- 1.3
  cm$mult_F[cm$occ == "child"] <- 0.7
  cm
}


# Class shares by age on a 0.5-year grid to 100: children to 5, then students,
# workers from 20 and at home from 60
occ_age_tab <- function() {
  ages <- seq(0, 100, by = 0.5)
  p <- function(a) {
    if (a < 5) c(1, 0, 0, 0) else if (a < 20) c(0, 1, 0, 0) else if (a < 60) c(0, 0.2, 0.6, 0.2)
    else c(0, 0, 0.3, 0.7)
  }
  m <- t(vapply(ages, p, numeric(4)))
  do.call(rbind, lapply(c("F", "M"), function(x) {
    data.frame(sex = x, age = ages, child = m[, 1], student = m[, 2], worker = m[, 3], home = m[, 4])
  }))
}

test_that("class transitions carry the class shares from one age to the next", {
  n <- 100000
  P <- matrix(c(0.5, 0.3, 0.2, 0), n, 4, byrow = TRUE)
  Q <- matrix(c(0.2, 0.3, 0.3, 0.2), n, 4, byrow = TRUE)
  lv <- EpiModelCOVID:::occ_levels
  occ <- lv[sample.int(4, n, replace = TRUE, prob = P[1, ])]
  set.seed(3)
  out <- EpiModelCOVID:::next_occ(occ, P, Q)
  expect_equal(as.numeric(table(factor(out, lv)) / n), Q[1, ], tolerance = 0.01)
  # only the excess moves: nobody leaves student or worker, whose shares do not fall
  expect_true(all(out[occ %in% c("student", "worker")] == occ[occ %in% c("student", "worker")]))
  expect_equal(mean(out != occ), 0.3, tolerance = 0.01)
  # no change in the distribution, no change of class
  expect_identical(EpiModelCOVID:::next_occ(occ, P, P), occ)
})

test_that("act.occ.age is checked and looked up by age step", {
  oa <- EpiModelCOVID:::normalize_occ_age(occ_age_tab())
  expect_equal(oa$step, 0.5)
  pr <- EpiModelCOVID:::occ_age_probs(c(4.99, 5, 25, 150), c("F", "M", "F", "M"), oa)
  expect_equal(pr[, 1], c(1, 0, 0, 0))
  expect_equal(pr[, 2], c(0, 1, 0.2, 0))
  bad <- occ_age_tab(); bad$home[3] <- 0.5
  expect_error(EpiModelCOVID:::normalize_occ_age(bad), "summing to 1")
  expect_error(EpiModelCOVID:::normalize_occ_age(occ_age_tab()[-2, ]), "same ages for each sex")
  expect_error(pop_init(param = pop_param(act.occ.age = bad)), "summing to 1")
})

test_that("a child who reaches school age enrolls and takes a student's traits and ties", {
  dat <- pop_init(param = pop_param(act.class.mult = cross_mult(), act.occ.age = occ_age_tab()))
  age <- get_attr(dat, "age")
  occ <- get_attr(dat, "occ")
  kids <- which(age < 10)[1:20]
  age[kids] <- 5 - 0.5 / 365
  occ[kids] <- "child"
  dat <- set_attr(dat, "age", age)
  dat <- set_attr(dat, "occ", occ)
  aF0 <- get_attr(dat, "log_aF")
  dat <- set_current_timestep(dat, 2)
  dat <- aging(dat, 2)
  expect_true(all(get_attr(dat, "occ")[kids] == "student"))
  expect_equal(get_epi(dat, "n.class.change", at = 2), 20)
  expect_identical(which(get_attr(dat, "rel.cross") == 2), kids)
  # same band, so the activity moves by the class multiplier ratio
  cm <- cross_mult()
  r <- cm$mult_F[cm$band == 1 & cm$occ == "student"][1] / cm$mult_F[cm$band == 1 & cm$occ == "child"][1]
  expect_equal(get_attr(dat, "log_aF")[kids] - aF0[kids], rep(log(r), 20), tolerance = 1e-8)
  # nobody else changes class: their distribution is flat over the step
  others <- setdiff(seq_len(N), kids)
  expect_identical(get_attr(dat, "occ")[others], occ[others])
})
