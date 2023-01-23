
#' @rdname moduleset-ship
#' @export
resim_nets_covid_ship <- function(dat, at) {
  ## Edges correction
  dat <- edges_correct_covid(dat, at)

  if (at < dat$param$network.lockdown.time) {
    nets <- 1:3
  } else {
    nets <- 4:6
  }

  # Network Resimulation
  # controls
  set.control.tergm <- get_control(dat, "set.control.tergm")
  set.control.ergm <- get_control(dat, "set.control.ergm")
  save.nwstats <- get_control(dat, "save.nwstats")
  nwstats.formulas <- get_control(dat, "nwstats.formulas")

  ## Edges correction
  dat <- edges_correct_covid(dat, at)

  # Network Resimulation
  for (i in nets) {
    nwparam <- EpiModel::get_nwparam(dat, network = i)
    isTERGM <- nwparam$coef.diss$duration > 1

    nwL <- networkLite(dat[["el"]][[i]], dat[["attr"]])

    if (get_control(dat, "tergmLite.track.duration")) {
      nwL %n% "time" <- dat[["nw"]][[i]] %n% "time"
      nwL %n% "lasttoggle" <- dat[["nw"]][[i]] %n% "lasttoggle"
    }

    if (isTERGM) {
      dat[["nw"]][[i]] <- simulate(
        nwL ~ Form(nwparam[["formation"]]) +
              Persist(nwparam[["coef.diss"]][["dissolution"]]),
        coef = c(nwparam[["coef.form"]], nwparam[["coef.diss"]][["coef.adj"]]),
        constraints = nwparam[["constraints"]],
        time.start = at - 1,
        time.slices = 1,
        time.offset = 1,
        monitor = nwstats.formulas[[i]],
        control = set.control.tergm,
        output = "final",
        dynamic = TRUE
      )
    } else {
      modified_controls <- set.control.tergm
      modified_controls$MCMC.prop.args <- list(discordance_fraction = 0)
      dat[["nw"]][[i]] <- simulate(
        basis = nwL,
        object = nwparam[["formation"]],
        coef = nwparam[["coef.form"]],
        constraints = nwparam[["constraints"]],
        monitor = nwstats.formulas[[i]],
        control = modified_controls,
        time.start = at - 1,
        time.slices = 1,
        time.offset = 1,
        dynamic = TRUE,
        output = "final"
      )
      rm(modified_controls)
    }

    dat[["el"]][[i]] <- as.edgelist(dat[["nw"]][[i]])
  }

  if (save.nwstats) {
    dat <- update_nwstats(dat)
  }

  return(dat)
}


#' @rdname moduleset-corporate
#' @export
resim_nets_covid_corporate <- function(dat, at) {

  # controls
  tergmLite.track.duration <- get_control(dat, "tergmLite.track.duration")
  set.control.tergm <- get_control(dat, "set.control.tergm")
  set.control.ergm <- get_control(dat, "set.control.ergm")
  save.nwstats <- get_control(dat, "save.nwstats")
  nwstats.formulas <- get_control(dat, "nwstats.formulas")

  ## Edges correction
  dat <- edges_correct_covid(dat, at)

  # Network Resimulation
  for (i in seq_along(dat[["el"]])) {
    nwparam <- get_nwparam(dat, network = i)
    isTERGM <- nwparam[["isTERGM"]]

    nwL <- networkLite(dat[["el"]][[i]], dat[["attr"]])

    if (isTERGM) {
      dat[["nw"]][[i]] <- simulate(
        nwL ~ Form(nwparam[["formation"]]) +
              Persist(nwparam[["coef.diss"]][["dissolution"]]),
        coef = c(nwparam[["coef.form"]], nwparam[["coef.diss"]][["coef.adj"]]),
        constraints = nwparam[["constraints"]],
        time.start = at - 1,
        time.slices = 1,
        time.offset = 1,
        monitor = nwstats.formulas[[i]],
        control = set.control.tergm,
        output = "final",
        dynamic = TRUE
      )
    } else {
      modified_controls <- set.control.tergm
      modified_controls$MCMC.prop.args <- list(discordance_fraction = 0)
      dat[["nw"]][[i]] <- simulate(
        basis = nwL,
        object = nwparam[["formation"]],
        coef = nwparam[["coef.form"]],
        constraints = nwparam[["constraints"]],
        monitor = nwstats.formulas[[i]],
        control = modified_controls,
        time.start = at - 1,
        time.slices = 1,
        time.offset = 1,
        dynamic = TRUE,
        output = "final"
      )
      rm(modified_controls)
    }

    dat[["el"]][[i]] <- as.edgelist(dat[["nw"]][[i]])
  }

  if (get_control(dat, "cumulative.edgelist")) {
    for (n_network in seq_along(dat$nwparam)) {
      dat <- update_cumulative_edgelist(dat, n_network,
                                        get_control(dat, "truncate.el.cuml"))
    }
  }

  if (save.nwstats) {
    dat <- update_nwstats(dat)
  }


  return(dat)
}


edges_correct_covid <- function(dat, at) {

  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)
  adjust <- log(old.num) - log(new.num)

  for (i in seq_along(dat$nwparam)) {
    coef.form1 <- get_nwparam(dat, network = i)$coef.form
    coef.form1[1] <- coef.form1[1] + adjust
    dat$nwparam[[i]]$coef.form <- coef.form1
  }

  return(dat)
}

update_nwstats <- function(dat) {
  for (i in seq_along(dat[["nwparam"]])) {
    new.nwstats <- tail(attributes(dat$nw[[i]])$stats, 1)
    keep.cols <- which(!duplicated(colnames(new.nwstats)))
    new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
    dat$stats$nwstats[[i]] <- rbind(dat$stats$nwstats[[i]], new.nwstats)
  }
  return(dat)
}
