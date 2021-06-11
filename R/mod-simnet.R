
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
  for (i in nets) {
    nwparam <- EpiModel::get_nwparam(dat, network = i)
    isTERGM <- all(nwparam$coef.diss$duration > 1)
    dat <- tergmLite::updateModelTermInputs(dat, network = i)
    if (isTERGM == TRUE) {
      rv <- tergmLite::simulate_network(state = dat$p[[i]]$state,
                                        coef = c(nwparam$coef.form,
                                                 nwparam$coef.diss$coef.adj),
                                        control = dat$control$mcmc.control[[i]],
                                        save.changes = TRUE)
      dat$el[[i]] <- rv$el
    } else {
      rv <- tergmLite::simulate_ergm(state = dat$p[[i]]$state,
                                     coef = nwparam$coef.form,
                                     control = dat$control$mcmc.control[[i]])

      dat$el[[i]] <- rv$el
    }

    if (get_control(dat, "save.nwstats") == TRUE) {
      nwL <- tergmLite::networkLite(dat$el[[i]], dat$attr)
      dat$stats$nwstats[[i]] <- rbind(dat$stats$nwstats[[i]],
                                      summary(dat$control$nwstats.formulas[[i]],
                                              basis = nwL,
                                              term.options = dat$control$mcmc.control[[i]]$term.options,
                                              dynamic = isTERGM))
    }
  }

  return(dat)
}


#' @rdname moduleset-corporate
#' @export
resim_nets_covid_corporate <- function(dat, at) {

  ## Edges correction
  dat <- edges_correct_covid(dat, at)

  # Network Resimulation
  for (i in 1:length(dat$el)) {
    nwparam <- EpiModel::get_nwparam(dat, network = i)
    isTERGM <- all(nwparam$coef.diss$duration > 1)
    dat <- tergmLite::updateModelTermInputs(dat, network = i)
    if (isTERGM == TRUE) {
      rv <- tergmLite::simulate_network(state = dat$p[[i]]$state,
                                        coef = c(nwparam$coef.form,
                                                 nwparam$coef.diss$coef.adj),
                                        control = dat$control$mcmc.control[[i]],
                                        save.changes = TRUE)
      dat$el[[i]] <- rv$el
    } else {
      rv <- tergmLite::simulate_ergm(state = dat$p[[i]]$state,
                                     coef = nwparam$coef.form,
                                     control = dat$control$mcmc.control[[i]])

      dat$el[[i]] <- rv$el
    }

    if (get_control(dat, "save.nwstats") == TRUE) {
      nwL <- tergmLite::networkLite(dat$el[[i]], dat$attr)
      dat$stats$nwstats[[i]] <- rbind(dat$stats$nwstats[[i]],
                                      summary(dat$control$nwstats.formulas[[i]],
                                              basis = nwL,
                                              term.options = dat$control$mcmc.control[[i]]$term.options,
                                              dynamic = isTERGM))
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

