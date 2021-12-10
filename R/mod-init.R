
#' @rdname moduleset-ship
#' @export
init_covid_ship <- function(x, param, init, control, s) {

  # Master Data List
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  dat$attr <- list()
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()

  ## Network Setup ##
  # Initial network simulations
  dat$nw <- list()
  for (i in 1:length(x)) {
    dat$nw[[i]] <- simulate(x[[i]]$fit, basis = x[[i]]$fit$newnetwork,
                            dynamic = FALSE)
  }
  nw <- dat$nw

  # Pull Network parameters
  dat$nwparam <- list()
  for (i in 1:length(x)) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }

  ## Nodal Attributes Setup ##
  num <- network.size(nw[[1]])
  dat$attr$active <- rep(1, num)
  dat$attr$arrival.time <- rep(1, num)
  dat$attr$uid <- 1:num

  # Pull in attributes on network
  nwattr.all <- names(nw[[1]][["val"]][[1]])
  nwattr.use <- nwattr.all[!nwattr.all %in% c("na", "vertex.names")]
  for (i in seq_along(nwattr.use)) {
    dat$attr[[nwattr.use[i]]] <- get.vertex.attribute(nw[[1]], nwattr.use[i])
  }

  # Convert to tergmLite method
  dat <- init_tergmLite(dat)

  ## Infection Status and Time Modules
  dat <- init_status_covid_ship(dat)

  ## Get initial prevalence
  dat <- prevalence_covid_ship(dat, at = 1)

  # Network stats
  if (get_control(dat, "save.nwstats") == TRUE) {
    for (i in 1:length(x)) {
      nwL <- networkLite(dat$el[[i]], dat$attr)
      nwstats <- summary(dat$control$nwstats.formulas[[i]],
                         basis = nwL,
                         term.options = dat$control$mcmc.control[[i]]$term.options,
                         dynamic = i < 3)

      dat$stats$nwstats[[i]] <- matrix(nwstats, nrow = 1,
                                       ncol = length(nwstats),
                                       dimnames = list(NULL, names(nwstats)))

      dat$stats$nwstats[[i]] <- as.data.frame(dat$stats$nwstats[[i]])
    }
  }

  return(dat)
}


init_status_covid_ship <- function(dat) {

  e.num.pass <- dat$init$e.num.pass
  e.num.crew <- dat$init$e.num.crew

  type <- dat$attr$type
  active <- dat$attr$active
  num <- sum(dat$attr$active)

  ## Disease status
  status <- rep("s", num)
  if (e.num.pass > 0) {
    status[sample(which(active == 1 & type == "p"), size = e.num.pass)] <- "e"
  }
  if (e.num.crew > 0) {
    status[sample(which(active == 1 & type == "c"), size = e.num.crew)] <- "e"
  }

  dat$attr$status <- status
  dat$attr$active <- rep(1, length(status))
  dat$attr$entrTime <- rep(1, length(status))
  dat$attr$exitTime <- rep(NA, length(status))

  # Infection Time
  idsInf <- which(status == "e")
  infTime <- rep(NA, length(status))
  clinical <- rep(NA, length(status))
  statusTime <- rep(NA, length(status))
  statusTime[idsInf] <- 1
  dxStatus <- rep(0, length(status))
  transmissions <- rep(0, length(status))

  dat$attr$statusTime <- statusTime
  dat$attr$infTime <- infTime
  dat$attr$clinical <- clinical
  dat$attr$dxStatus <- dxStatus
  dat$attr$transmissions <- transmissions

  return(dat)
}


#' @rdname moduleset-corporate
#' @export
init_covid_corporate <- function(x, param, init, control, s) {

  # Master Data List
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  dat$attr <- list()
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()

  ## Network Setup ##
  # Initial network simulations
  dat$nw <- list()
  for (i in 1:length(x)) {
    dat$nw[[i]] <- simulate(x[[i]]$fit, basis = x[[i]]$fit$newnetwork,
                            dynamic = FALSE)
  }
  nw <- dat$nw

  # Pull Network parameters
  dat$nwparam <- list()
  for (i in 1:length(x)) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }

  ## Nodal Attributes Setup ##
  num <- network.size(nw[[1]])
  dat <- append_core_attr(dat, 1, num)

  # Pull in attributes on network
  nwattr.all <- names(nw[[1]][["val"]][[1]])
  nwattr.use <- nwattr.all[!nwattr.all %in% c("na", "vertex.names")]
  for (i in seq_along(nwattr.use)) {
    dat$attr[[nwattr.use[i]]] <- get.vertex.attribute(nw[[1]], nwattr.use[i])
  }

  # Convert to tergmLite method
  dat <- init_tergmLite(dat)

  ## Infection Status and Time Modules
  dat <- init_status_covid_corporate(dat)

  ## Get initial prevalence
  dat <- prevalence_covid_corporate(dat, at = 1)

  # Network stats
  if (get_control(dat, "save.nwstats") == TRUE) {
    for (i in 1:length(x)) {
      nwL <- networkLite(dat$el[[i]], dat$attr)
      nwstats <- summary(dat$control$nwstats.formulas[[i]],
                         basis = nwL,
                         term.options = dat$control$mcmc.control[[i]]$term.options,
                         dynamic = i < 3)

      dat$stats$nwstats[[i]] <- matrix(nwstats, nrow = 1,
                                       ncol = length(nwstats),
                                       dimnames = list(NULL, names(nwstats)))

      dat$stats$nwstats[[i]] <- as.data.frame(dat$stats$nwstats[[i]])
    }
  }

  return(dat)
}


init_status_covid_corporate <- function(dat) {

  e.num <- get_init(dat, "e.num")

  active <- get_attr(dat, "active")
  num <- sum(active)

  ## Disease status
  status <- rep("s", num)
  if (e.num > 0) {
    status[sample(which(active == 1), size = e.num)] <- "e"
  }

  dat <- set_attr(dat, "status", status)

  # Infection Time and related attributes
  idsInf <- which(status == "e")
  infTime <- rep(NA, num)
  clinical <- rep(NA, num)
  hospit <- rep(NA, num)
  statusTime <- rep(NA, num)
  statusTime[idsInf] <- 1
  dxStatus <- rep(0, num)
  vax <- rep(0, num)
  vax1Time <- rep(NA, num)
  vax2Time <- rep(NA, num)

  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "infTime", infTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "hospit", hospit)
  dat <- set_attr(dat, "dxStatus", dxStatus)
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vax1Time", vax1Time)
  dat <- set_attr(dat, "vax2Time", vax2Time)

  return(dat)
}



#' @rdname moduleset-contacttrace
#' @export
init_covid_contacttrace <- function(x, param, init, control, s) {
  dat <- create_dat_object(param, init, control)

  # Master Data List
  dat$stats$nwstats <- list()

  ## Network Setup ##
  # Initial network simulations
  dat$nw <- list()
  for (i in 1:length(x)) {
    dat$nw[[i]] <- simulate(x[[i]]$fit, basis = x[[i]]$fit$newnetwork)
  }
  nw <- dat$nw

  # Pull Network parameters
  dat$nwparam <- list()
  for (i in 1:length(x)) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }

  ## Nodal Attributes Setup ##
  num <- network.size(nw[[1]])
  dat <- append_core_attr(dat, 1, num)

  # Pull in attributes on network
  nwattr.all <- names(nw[[1]][["val"]][[1]])
  nwattr.use <- nwattr.all[!nwattr.all %in% c("na", "vertex.names")]
  for (i in seq_along(nwattr.use)) {
    dat$attr[[nwattr.use[i]]] <- get.vertex.attribute(nw[[1]], nwattr.use[i])
  }

  # Convert to tergmLite method
  dat <- init_tergmLite(dat)

  ## Infection Status and Time Modules
  dat <- init_status_covid_contacttrace(dat)

  ## Get initial prevalence
  dat <- prevalence_covid_contacttrace(dat, at = 1)

  # Network stats
  if (get_control(dat, "save.nwstats")) {
    for (i in seq_along(x)) {
      nwL <- networkLite(dat$el[[i]], dat$attr)
      nwstats <- summary(
        dat$control$nwstats.formulas[[i]],
        basis = nwL,
        term.options = dat$control$mcmc.control[[i]]$term.options,
        dynamic = i < 3
      )

      dat$stats$nwstats[[i]] <- matrix(
        nwstats, 
        nrow = 1, ncol = length(nwstats),
        dimnames = list(NULL, names(nwstats))
      )

      dat$stats$nwstats[[i]] <- as.data.frame(dat$stats$nwstats[[i]])
    }
  }

  for (n_network in seq_along(dat[["nw"]])) {
    dat <- update_cumulative_edgelist(dat, n_network)
  }

  return(dat)
}


init_status_covid_contacttrace <- function(dat) {

  e.num <- get_init(dat, "e.num")

  active <- get_attr(dat, "active")
  num <- sum(active)

  ## Disease status
  status <- rep("s", num)
  if (e.num > 0) {
    status[sample(which(active == 1), size = e.num)] <- "e"
  }

  dat <- set_attr(dat, "status", status)

  # Infection Time and related attributes
  idsInf <- which(status == "e")
  infTime <- rep(NA, num)
  eligible.case <- rep(NA, num)
  clinical <- rep(NA, num)
  # hospit <- rep(NA, num)
  intensive <- rep(NA, num)
  branch <- rep(NA, num)
  statusTime.Ic <- rep(NA, num)
  statusTime <- rep(NA, num)
  statusTime[idsInf] <- 1
  dxStatus <- rep(0, num)
  dxTime <- rep(NA, num)
  # vax <- rep(0, num)
  # vax1Time <- rep(NA, num)

  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "statusTime.Ic", statusTime.Ic)
  dat <- set_attr(dat, "infTime", infTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "eligible.case", eligible.case)
  # dat <- set_attr(dat, "hospit", hospit)
  dat <- set_attr(dat, "intensive", intensive)
  dat <- set_attr(dat, "branch", branch)
  dat <- set_attr(dat, "dxStatus", dxStatus)
  dat <- set_attr(dat, "dxTime", dxTime)
  # dat <- set_attr(dat, "vax", vax)
  # dat <- set_attr(dat, "vax1Time", vax1Time)

  return(dat)
}
