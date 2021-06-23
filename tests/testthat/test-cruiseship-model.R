

# Model Setup -------------------------------------------------------------

mortality_rate <- c(588.45, 24.8, 11.7, 14.55, 47.85, 88.2, 105.65, 127.2,
                    154.3, 206.5, 309.3, 495.1, 736.85, 1051.15, 1483.45,
                    2294.15, 3642.95, 6139.4, 13938.3)

mr_pp_pd <- mortality_rate / 1e5 / 365

age_spans <- c(1, 4, rep(5, 16), 1)
mr_vec <- rep(mr_pp_pd, times = age_spans)

n.sectors <- 10
prop.within.sector.post <- 0.98
n.sectors.pass <- NULL  # If NULL, only within-cabin pass/pass mixing

n.crew <- 1045
n.pass <- 2666
n <- n.crew + n.pass

n.rooms <- 1337
n.pass.per.room <- n.pass/n.rooms
n.pass.per.room

pass.ids <- 1:n.pass

room.ids <- 1:n.rooms
room.ids.pass <- apportion_lr(n.pass, room.ids, rep(1/n.rooms, n.rooms))

type.attr <- rep(c("p", "c"), times = c(n.pass, n.crew))

crew.ids <- which(type.attr == "c")

room.sectors <- apportion_lr(length(room.ids), 1:n.sectors,
                             rep(1/n.sectors, n.sectors))
df.match <- data.frame(room.ids, room.sectors)

sector <- rep(NA, n)
for (id in pass.ids) {
  sector[id] <- df.match[df.match$room.ids == room.ids.pass[id], "room.sectors"]
}

room.sectors.c <- apportion_lr(length(crew.ids), 1:n.sectors,
                               rep(1/n.sectors, n.sectors), shuffled = TRUE)
sector[crew.ids] <- room.sectors.c

ages.pass <- round(rnorm(n.pass, 69, 6), 1)
ages.crew <- pmax(round(rnorm(n.crew, 36, 10), 1), 18)

age <- rep(0, n)
age[type.attr == "p"] <- ages.pass
age[type.attr == "c"] <- ages.crew

nw <- network_initialize(n)
nw <- set_vertex_attribute(nw, "type", type.attr)
nw <- set_vertex_attribute(nw, "pass.room", room.ids.pass, pass.ids)
nw <- set_vertex_attribute(nw, "pass.room", 0, crew.ids)
nw <- set_vertex_attribute(nw, "age", age)
nw <- set_vertex_attribute(nw, "sector", sector)

## PP network
md.pre <- 5
edges.pre <- n.pass * md.pre/2
target.stats1.pre <- edges.pre

formation1.pre <- ~edges + offset(nodefactor("type", levels = -2))

coef.diss1 <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
est1 <- suppressMessages(
  netest(nw, formation1.pre, target.stats1.pre, coef.diss1, coef.form = -Inf,
         set.control.ergm = control.ergm(MCMLE.maxit = 500))
)

## CC network
md.pre <- 10
edges.pre <- n.crew * md.pre/2
prop.within.sector.pre <- 0.50
target.stats2.pre <- edges.pre

formation2.pre <- ~edges + offset(nodefactor("type", levels = -1))
coef.diss2 <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
est2 <- suppressMessages(
  netest(nw, formation2.pre, target.stats2.pre, coef.diss2, coef.form = -Inf,
         set.control.ergm = control.ergm(MCMLE.maxit = 500))
)

## CP network
md.pre <- 8
edges.pre <- md.pre * n.rooms
target.stats3.pre <- c(edges.pre, 0)

formation3.pre <- ~edges + nodematch("type")
coef.diss3 <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
est3 <- suppressMessages(
  netest(nw, formation3.pre, target.stats3.pre, coef.diss3,
         set.control.ergm = control.ergm(MCMLE.maxit = 500))
)

est <- list(est1, est2, est3, est1, est2, est3)


# Model Testing -----------------------------------------------------------

test_that("base ship model parameterization", {
  param <- param.net(inf.prob.pp = 0.10455227,
                     inf.prob.pp.inter.rr = 0.6,
                     inf.prob.pp.inter.time = Inf,
                     act.rate.pp = 5,
                     act.rate.pp.inter.rr = 1,
                     act.rate.pp.inter.time = Inf,
                     inf.prob.pc = 0.10455227,
                     inf.prob.pc.inter.rr = 0.6,
                     inf.prob.pc.inter.time = 15,
                     act.rate.pc = 1,
                     act.rate.pc.inter.rr = 1,
                     act.rate.pc.inter.time = Inf,
                     inf.prob.cc = 0.10455227,
                     inf.prob.cc.inter.rr = 0.6,
                     inf.prob.cc.inter.time = 15,
                     act.rate.cc = 1,
                     act.rate.cc.inter.rr = 1,
                     act.rate.cc.inter.time = Inf,
                     inf.prob.a.rr = 0.5,
                     prop.clinical = c(0.40, 0.25, 0.37, 0.42, 0.51, 0.59, 0.72, 0.76),
                     act.rate.dx.inter.rr = 0.1,
                     act.rate.dx.inter.time = 15,
                     act.rate.sympt.inter.rr = 0.1,
                     act.rate.sympt.inter.time = 15,
                     network.lockdown.time = 15,
                     ea.rate = 1/4.0,
                     ar.rate = 1/5.0,
                     eip.rate = 1/4.0,
                     ipic.rate = 1/1.5,
                     icr.rate = 1/3.5,
                     pcr.sens = 0.8,
                     dx.rate.sympt = c(rep(0, 15), rep(0.21430869, 5), rep(0.30194363, 5), rep(0.80075530, 100)),
                     dx.rate.other = c(rep(0, 15), rep(0, 5), rep(0.09680617, 5), rep(0.21919574, 100)),
                     allow.rescreen = FALSE,
                     mort.rates = mr_vec,
                     mort.dis.mult = 180,
                     exit.rate.pass = 0,
                     exit.rate.crew = 0,
                     exit.elig.status = c("ip", "ic"),
                     exit.require.dx = FALSE)
  init <- init.net(e.num.pass = 8,
                   e.num.crew = 0)
  control <- control.net(nsteps = 31,
                         nsims = 1,
                         ncores = 1,
                         initialize.FUN = init_covid_ship,
                         aging.FUN = aging_covid,
                         departures.FUN = deaths_covid_ship,
                         arrivals.FUN = NULL,
                         edges_correct.FUN = NULL,
                         resim_nets.FUN = resim_nets_covid_ship,
                         infection.FUN = infect_covid_ship,
                         recovery.FUN = progress_covid_ship,
                         dx.FUN = dx_covid,
                         prevalence.FUN = prevalence_covid_ship,
                         nwupdate.FUN = NULL,
                         module.order = c("aging.FUN",
                                          "departures.FUN",
                                          "resim_nets.FUN",
                                          "infection.FUN",
                                          "recovery.FUN",
                                          "dx.FUN",
                                          "prevalence.FUN"),
                         resimulate.network = TRUE,
                         skip.check = TRUE,
                         tergmLite = TRUE,
                         mcmc.control.ergm.1 = control.simulate.formula(),
                         mcmc.control.ergm.2 = control.simulate.formula(),
                         mcmc.control.ergm.3 = control.simulate.formula(),
                         mcmc.control.ergm.4 = control.simulate.formula(),
                         mcmc.control.ergm.5 = control.simulate.formula(),
                         mcmc.control.ergm.6 = control.simulate.formula(),
                         verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
})

