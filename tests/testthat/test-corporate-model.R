
# Model Setup -------------------------------------------------------------

library("EpiModelCOVID")

mortality_rate <- c(588.45, 24.8, 11.7, 14.55, 47.85, 88.2, 105.65, 127.2,
                    154.3, 206.5, 309.3, 495.1, 736.85, 1051.15, 1483.45,
                    2294.15, 3642.95, 6139.4, 13938.3)

mr_pp_pd <- mortality_rate / 1e5 / 365

age_spans <- c(1, 4, rep(5, 16), 1)
mr_vec <- rep(mr_pp_pd, times = age_spans)

# Set network
n <- 1000
nw <- network_initialize(n)

age <- round(rnorm(n, 40, 21), 1)
age <- pmax(age, 0)
age <- pmin(age, 99)

age.breaks <- seq(0, 100, 10)
age.grp <- cut(age, age.breaks, labels = FALSE, right = FALSE)

nw <- set_vertex_attribute(nw, "age", age)
nw <- set_vertex_attribute(nw, "age.grp", age.grp)

## Within household network
md.hh <- 2.5
target.stats.hh <- md.hh * n/2

formation.hh <- ~edges

coef.diss.hh <- dissolution_coefs(dissolution = ~offset(edges), duration = 1e5)
suppressMessages(
  est.hh <- netest(nw, formation.hh, target.stats.hh, coef.diss.hh,
                   set.control.ergm = control.ergm(MCMLE.maxit = 500))
)

## Within office network
md.oo <- 5
target.stats.oo <- md.oo * n/2

formation.oo <- ~edges

coef.diss.oo <- dissolution_coefs(dissolution = ~offset(edges), duration = 1e5)

suppressMessages(
  est.oo <- netest(nw, formation.oo, target.stats.oo, coef.diss.oo,
                   set.control.ergm = control.ergm(MCMLE.maxit = 500))
)

## Community network
md.cc <- 2
target.stats.cc <- md.cc * n/2

formation.cc <- ~edges

coef.diss.cc <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)

suppressMessages(
  est.cc <- netest(nw, formation.cc, target.stats.cc, coef.diss.cc,
                   set.control.ergm = control.ergm(MCMLE.maxit = 500))
)

est <- list(est.hh, est.oo, est.cc)


# Model Testing -----------------------------------------------------------

test_that("base corporate model parameterization", {
  param <- param.net(inf.prob = c(0.11, 0.11, 0.11),
                     act.rate = c(5, 1, 1),
                     inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                     inf.prob.inter.time = c(Inf, 15, 15),
                     act.rate.inter.rr = c(1, 1, 1),
                     act.rate.inter.time = c(Inf, Inf, Inf),
                     inf.prob.a.rr = 0.5,

                     act.rate.dx.inter.rr = 0.1,
                     act.rate.dx.inter.time = 1,
                     act.rate.sympt.inter.rr = 0.5,
                     act.rate.sympt.inter.time = 1,

                     prop.clinical = c(0.40, 0.25, 0.37, 0.42, 0.51, 0.59, 0.72, 0.76),
                     prop.hospit = c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),

                     ea.rate = 1/4.0,
                     ar.rate = 1/5.0,
                     eip.rate = 1/4.0,
                     ipic.rate = 1/1.5,
                     icr.rate = 1/3.5,
                     ich.rate = 1/3.5,
                     hr.rate = 1/7,

                     pcr.sens = 0.8,
                     dx.rate.sympt = 0.1,
                     dx.rate.other = 0.01,
                     allow.rescreen = TRUE,

                     vax.start = Inf,
                     vax1.rate = 0.01,
                     vax2.interval = 21,
                     vax1.rr.infect = 0.75,
                     vax2.rr.infect = 0.25,
                     vax.rr.clinical = 0.05,
                     vax1.immune = 7,
                     vax2.immune = 14,

                     a.rate = mean(mr_vec),
                     arrival.age = 0,
                     mort.rates = mr_vec,
                     mort.dis.mult = 180)
  init <- init.net(e.num = 100)
  control <- control.net(nsteps = 100,
                         nsims = 1,
                         ncores = 1,
                         initialize.FUN = init_covid_corporate,
                         aging.FUN = aging_covid,
                         departures.FUN = deaths_covid_corporate,
                         arrivals.FUN = arrival_covid_corporate,
                         edges_correct.FUN = NULL,
                         resim_nets.FUN = resim_nets_covid_corporate,
                         infection.FUN = infect_covid_corporate,
                         recovery.FUN = progress_covid,
                         dx.FUN = dx_covid,
                         vax.FUN = vax_covid,
                         prevalence.FUN = prevalence_covid_corporate,
                         module.order = c("aging.FUN",
                                          "departures.FUN",
                                          "arrivals.FUN",
                                          "resim_nets.FUN",
                                          "infection.FUN",
                                          "recovery.FUN",
                                          "dx.FUN",
                                          "vax.FUN",
                                          "prevalence.FUN"),
                         resimulate.network = TRUE,
                         skip.check = TRUE,
                         tergmLite = TRUE,
                         verbose = FALSE)

  sim <- netsim(est, param, init, control)

  expect_s3_class(sim, "netsim")

  df <- as.data.frame(sim)
  expect_equal(nrow(df), 100)
  expect_gt(sum(df$se.flow, na.rm = TRUE), 0)
  expect_equal(sum(df$v1.num + df$v2.num, na.rm = TRUE), 0)

})
