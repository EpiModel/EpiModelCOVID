# The parameters the GMC19 model reads, and the checks init_gmc19() runs on
# them (GlobalMix India #67)

# The GMC19 model's modules (make_control_gmc19() in GlobalMix-India-NetVax) and
# the household layer's arrival rule. Every function they call is part of the
# model; the package's other functions are the CorporateMix modules.
gmc19_modules <- c("init_gmc19", "aging", "deaths_covid_gmc19", "arrival", "hh_join_birth",
                   "resim_relationships", "activate_contacts", "netdegree", "infect_general",
                   "progress_general", "vax_general", "prevalence_covid_corporate")

ns <- asNamespace("EpiModelCOVID")
pkg_funs <- Filter(function(f) is.function(get(f, ns)), ls(ns, all.names = TRUE))

# Parameter names a function reads by name: get_param(dat, "x"), the engine's
# param_or(dat, "x", default) and layer_param(dat, "x"), and dat$param$x. The
# per-dose vaccine schedule, read under names built with paste0(), is tested
# below on build_vax_schedule() itself.
param_reads <- function(fun) {
  out <- character(0)
  walk <- function(e) {
    if (!is.call(e)) return()
    f <- e[[1]]
    if (is.name(f) && as.character(f) %in% c("get_param", "param_or", "layer_param")) {
      args <- as.list(e)[-1]
      item <- if ("item" %in% names(args)) args$item else if ("name" %in% names(args)) args$name else args[[2]]
      if (is.character(item)) out <<- c(out, item)
    }
    if (identical(f, as.name("$")) && identical(e[[2]], quote(dat$param))) {
      out <<- c(out, as.character(e[[3]]))
    }
    for (a in as.list(e)[-1]) if (!missing(a)) walk(a)
  }
  walk(body(fun))
  unique(out)
}
reads <- lapply(setNames(nm = pkg_funs), function(f) param_reads(get(f, ns)))

# The functions the modules reach, through any chain of calls
reachable <- character(0)
todo <- gmc19_modules
while (length(todo) > 0) {
  f <- todo[1]; todo <- todo[-1]
  if (f %in% reachable) next
  reachable <- c(reachable, f)
  todo <- c(todo, intersect(codetools::findGlobals(get(f, ns)), pkg_funs))
}
gmc19_reads <- unique(unlist(reads[reachable]))
corporate_reads <- unique(unlist(reads[setdiff(pkg_funs, reachable)]))
# read only to stop a run that sets them
refused <- c(EpiModelCOVID:::retired_params, names(EpiModelCOVID:::replaced_params), "hh.pairs")

test_that("every parameter the modules read is listed, and every listed one is read", {
  expect_true(all(gmc19_modules %in% pkg_funs))
  known <- gmc19_known_params(n.doses = numeric(0))
  expect_identical(setdiff(gmc19_reads, c(known, refused)), character(0))
  listed <- c(EpiModelCOVID:::gmc19_params_common, EpiModelCOVID:::gmc19_params_relationship)
  expect_identical(setdiff(listed, gmc19_reads), character(0))
  expect_identical(anyDuplicated(listed), 0L)
  expect_identical(intersect(listed, refused), character(0))
})

test_that("the relationship design reads the legacy design's names plus act.*, new.* and birth.*", {
  rel <- gmc19_known_params("relationship")
  legacy <- gmc19_known_params("legacy")
  expect_identical(gmc19_known_params(), rel)
  expect_true(all(legacy %in% rel))
  expect_true(all(grepl("^(act|new|birth)\\.", setdiff(rel, legacy))))
})

test_that("the CorporateMix modules' own parameters are listed apart", {
  corp <- EpiModelCOVID:::corporate_params
  expect_setequal(setdiff(corporate_reads, c(gmc19_known_params(), refused)), corp)
  expect_identical(intersect(corp, gmc19_known_params()), character(0))
})

test_that("the names EpiModel sets are listed", {
  p <- param.net(inf.prob = 0.1)
  expect_true(all(c("act.rate", "vital") %in% names(p)))
  nw <- network_initialize(n = 20)
  est <- suppressMessages(netest(nw, formation = ~edges, target.stats = 10,
                                 coef.diss = dissolution_coefs(~offset(edges), 10), verbose = FALSE))
  dat <- EpiModel::init_nets(create_dat_object(p, init.net(i.num = 1), control.net(type = "SI", nsteps = 2)),
                             est)
  expect_true("groups" %in% names(dat$param))
  expect_true(all(EpiModelCOVID:::epimodel_params %in% gmc19_known_params()))
})

test_that("build_vax_schedule() reads exactly the listed per-dose names", {
  known <- gmc19_known_params(n.doses = c(covid = 2))
  # the product's age window (#74) is read by check_vax_campaign(), the rest
  # by build_vax_schedule()
  expect_true(all(c("vax.covid.age.min", "vax.covid.age.max") %in% known))
  sched <- grep("^vax\\.covid\\.(n\\.doses|.+\\.dose[0-9]+)$", known, value = TRUE)
  expect_length(sched, 1 + 2 * length(EpiModelCOVID:::vax_dose_fields))
  value <- function(nm) {
    if (grepl("annual", nm)) FALSE
    else if (nm == "vax.covid.n.doses") 2
    else 0.5
  }
  param <- setNames(lapply(sched, value), sched)
  dat <- create_dat_object(param = param)
  vs <- EpiModelCOVID:::build_vax_schedule(dat, "covid")
  expect_identical(nrow(vs), 2L)
  expect_setequal(names(vs), c("dose", EpiModelCOVID:::vax_dose_fields))
  expect_equal(vs$interval[[2]], 0.5)
  expect_equal(unlist(vs$ve.peak.death), c(0.5, 0.5))
  # every listed name is read: without it the schedule cannot be built
  for (nm in sched) {
    expect_error(EpiModelCOVID:::build_vax_schedule(create_dat_object(param = param[names(param) != nm]),
                                                    "covid"), info = nm)
  }
})

test_that("infant_prog_pars() reads exactly the listed infant product names", {
  known <- gmc19_known_params(n.doses = numeric(0))
  for (pr in EpiModelCOVID:::infant_products) {
    nm <- grep(paste0("^vax\\.", pr, "\\.ve\\."), known, value = TRUE)
    expect_length(nm, length(EpiModelCOVID:::infant_prog_fields))
    param <- setNames(as.list(seq_along(nm) / 10), nm)
    expect_length(EpiModelCOVID:::infant_prog_pars(create_dat_object(param = param), pr)$peak, 5)
    for (x in nm) {
      expect_error(EpiModelCOVID:::infant_prog_pars(create_dat_object(param = param[names(param) != x]), pr),
                   info = x)
    }
  }
})

# init_gmc19() up to its first checks; nothing after them needs the networks
init_param <- function(...) {
  param.net(inf.prob = rep(0.1, 2), layer.names = c("community", "home"), ...)
}
run_init <- function(param) {
  init_gmc19(NULL, param, init.net(e.num = 0), control.net(type = "SI", nsteps = 2), s = 1)
}

test_that("time.unit must be set, to 1", {
  expect_error(run_init(init_param()), "no time.unit; add time.unit = 1")
  expect_error(run_init(init_param(time.unit = 7)), "time.unit is 7, but the model's rates")
  expect_error(run_init(init_param(time.unit = NA)), "time.unit is NA")
})

test_that("an unset inf.prob stops the run", {
  msg <- "inf.prob is set by the pathogen scenario, scenarios_pathogen.csv; the registry default is NA"
  p <- init_param(time.unit = 1)
  p$inf.prob <- c(NA, NA)
  expect_error(run_init(p), msg, fixed = TRUE)
  p$inf.prob <- c(0.1, NA)
  expect_error(run_init(p), msg, fixed = TRUE)
  p$inf.prob <- NULL
  expect_error(run_init(p), msg, fixed = TRUE)
})

test_that("one warning names every parameter nothing reads", {
  dat <- create_dat_object(param = init_param(time.unit = 1, vax1.rate = 0.1, half.life = 5,
                                              act.q.mult = 1, vax.covid.n.doses = 1,
                                              vax.covid.ve.delay.dose1 = 14,
                                              vax.covid.ve.delay.dose2 = 14))
  dat$param$.param.updater.list <- list(list(at = 20, param = list(sympt.contact.rr = 1, zzz = 2)))
  # act.q.mult is read only by the relationship model, and dose 2 is past
  # vax.covid.n.doses; no warning names a parameter that is read
  expect_warning(
    unknown <- EpiModelCOVID:::check_known_params(dat),
    "reads the parameter\\(s\\) act.q.mult, half.life, vax.covid.ve.delay.dose2, vax1.rate, zzz;")
  expect_identical(unknown, c("act.q.mult", "half.life", "vax.covid.ve.delay.dose2", "vax1.rate", "zzz"))
  dat$param$act.o <- rep(0.9, 6)
  expect_warning(EpiModelCOVID:::check_known_params(dat), "half.life, vax.covid.ve.delay.dose2, vax1.rate, zzz;")
  dat$param[c("vax1.rate", "half.life", "vax.covid.ve.delay.dose2")] <- NULL
  dat$param$.param.updater.list <- NULL
  expect_silent(EpiModelCOVID:::check_known_params(dat))
})
