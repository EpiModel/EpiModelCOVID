
##
## COVID-19 Corporate Network Model
## Epidemic Model
##
## Authors: Samuel M. Jenness
## Date: February 2021
##

library("EpiModelCOVID")


# Read in fitted network models
est <- readRDS("data/input/est.simple.rds")

# Model parameters
source("R/01-epi-params.R")
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
                   vax2.interval = 21 + 14,
                   vax1.rr.infect = 0.75,
                   vax2.rr.infect = 0.25,
                   vax.rr.clinical = 0.05,

                   a.rate = mean(mr_vec),
                   arrival.age = 0,
                   mort.rates = mr_vec,
                   mort.dis.mult = 180)

# Initial conditions
init <- init.net(e.num = 100)

# Control settings
pkgload::load_all("~/git/EpiModelCOVID")
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
                       resimulate.networks = TRUE,
                       skip.check = TRUE,
                       tergmLite = TRUE)

sim <- netsim(est, param, init, control)
# print(sim)


plot(sim, y = c("s.num", "e.num", "a.num"))

sim <- mutate_epi(sim, se.cuml = cumsum(se.flow),
                       dx.cuml = cumsum(nDx),
                       dx.pos.cuml = cumsum(nDx.pos),
                       totI = e.num + a.num + ip.num + ic.num)
df <- as.data.frame(sim, out = "mean")
names(df)
names(df)
