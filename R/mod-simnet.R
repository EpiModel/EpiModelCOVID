
#' @rdname moduleset-ship
#' @export
resim_nets_covid_ship <- function(dat, at) {

  # controls
  set.control.stergm <- get_control(dat, "set.control.stergm")
  set.control.ergm <- get_control(dat, "set.control.ergm")

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
        control = set.control.tergm,
        output = "final",
        dynamic = TRUE
      )
    } else {
      dat[["nw"]][[i]] <- simulate(
        basis = nwL,
        object = nwparam[["formation"]],
        coef = nwparam[["coef.form"]],
        constraints = nwparam[["constraints"]],
        control = set.control.ergm,
        dynamic = FALSE,
        nsim = 1,
        output = "network"
      )
    }

    dat[["el"]][[i]] <- as.edgelist(dat[["nw"]][[i]])

    if (save.nwstats) {
      if (isTERGM) {
        term.options <- set.control.tergm$term.options
      } else {
        term.options <- set.control.ergm$term.options
      }
      dat$stats$nwstats[[i]] <- rbind(
        dat$stats$nwstats[[i]],
        summary(
          dat$control$nwstats.formulas[[i]],
          basis = nwL,
          term.options = term.options,
          dynamic = isTERGM
        )
      )
    }
  }

  return(dat)
}


#' @rdname moduleset-corporate
#' @export
resim_nets_covid_corporate <- function(dat, at) {

  # controls
  set.control.tergm <- get_control(dat, "set.control.tergm")
  set.control.ergm <- get_control(dat, "set.control.ergm")
  save.nwstats <- get_control(dat, "save.nwstats")

  ## Edges correction
  dat <- edges_correct_covid(dat, at)

  # Network Resimulation
  for (i in 1:length(dat$el)) {
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
        control = set.control.tergm,
        output = "final",
        dynamic = TRUE
      )
    } else {
      dat[["nw"]][[i]] <- simulate(
        basis = nwL,
        object = nwparam[["formation"]],
        coef = nwparam[["coef.form"]],
        constraints = nwparam[["constraints"]],
        control = set.control.ergm,
        dynamic = FALSE,
        nsim = 1,
        output = "network"
      )
    }

    dat[["el"]][[i]] <- as.edgelist(dat[["nw"]][[i]])

    if (save.nwstats) {
      if (isTERGM) {
        term.options <- set.control.tergm$term.options
      } else {
        term.options <- set.control.ergm$term.options
      }
      dat$stats$nwstats[[i]] <- rbind(
        dat$stats$nwstats[[i]],
        summary(
          dat$control$nwstats.formulas[[i]],
          basis = nwL,
          term.options = term.options,
          dynamic = isTERGM
        )
      )
    }

  }

  return(dat)
}


edges_correct_covid <- function(dat, at) {

  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)
  adjust <- log(old.num) - log(new.num)

  for (i in 1:length(dat$nwparam)) {
    coef.form1 <- get_nwparam(dat, network = i)$coef.form
    coef.form1[1] <- coef.form1[1] + adjust
    dat$nwparam[[i]]$coef.form <- coef.form1
  }

  return(dat)
}

#' @rdname moduleset-contacttrace
#' @export
resim_nets_covid_contacttrace <- function(dat, at) {
  # Control
  tergmLite.track.duration <- get_control(dat, "tergmLite.track.duration")
  set.control.stergm <- get_control(dat, "set.control.stergm")
  set.control.ergm <- get_control(dat, "set.control.ergm")

  ## Edges correction
  dat <- edges_correct_covid(dat, at)

  # Network Resimulation
  for (i in 1:length(dat$el)) {
    nwparam <- EpiModel::get_nwparam(dat, network = i)
    isTERGM <- ifelse(nwparam$coef.diss$duration > 1, TRUE, FALSE)

    nwL <- networkLite(dat[["el"]][[i]], dat[["attr"]])

    if (tergmLite.track.duration == TRUE) {
      nwL %n% "time" <- dat[["nw"]][[i]] %n% "time"
      nwL %n% "lasttoggle" <- dat[["nw"]][[i]] %n% "lasttoggle"
    }

    if (isTERGM == TRUE) {
      dat[["nw"]][[i]] <- simulate(
        nwL,
        formation = nwparam[["formation"]],
        dissolution = nwparam[["coef.diss"]][["dissolution"]],
        coef.form = nwparam[["coef.form"]],
        coef.diss = nwparam[["coef.diss"]][["coef.adj"]],
        constraints = nwparam[["constraints"]],
        time.start = at - 1,
        time.slices = 1,
        time.offset = 1,
        control = set.control.stergm,
        output = "final"
      )
    } else {
      dat[["nw"]][[i]] <- simulate(
        basis = nwL,
        object = nwparam[["formation"]],
        coef = nwparam[["coef.form"]],
        constraints = nwparam[["constraints"]],
        control = set.control.ergm,
        dynamic = FALSE,
        nsim = 1,
        output = "network"
      )
    }

    dat[["el"]][[i]] <- as.edgelist(dat[["nw"]][[i]])
  }

  if (dat$control$save.nwstats == TRUE) {
    if (get_control(dat, "save.nwstats") == TRUE) {
      term.options <- if (isTERGM == TRUE) {
        set.control.stergm$term.options
      } else {
        set.control.ergm$term.options
      }
      dat$stats$nwstats[[i]] <- rbind(
        dat$stats$nwstats[[i]],
        summary(
          dat$control$nwstats.formulas[[i]],
          basis = nwL,
          term.options = term.options,
          dynamic = isTERGM
        )
      )
    }
  }

  for (n_network in seq_along(dat[["nw"]])) {
    dat <- update_cumulative_edgelist(dat, n_network)
  }

  return(dat)
}
