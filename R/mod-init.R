
#' @rdname moduleset-ship
#' @export
init_covid_ship <- function(x, param, init, control, s) {

  # Master Data List
  dat <- create_dat_object(param, init, control)

  ## Network Setup ##
  # Initial network simulations
  dat[["nw"]] <- list()
  for (i in seq_along(x)) {
    dat[["nw"]][[i]] <- simulate(
      x[[i]][["formula"]],
      coef = x[[i]][["coef.form.crude"]],
      basis = x[[i]][["newnetwork"]],
      constraints = x[[i]][["constraints"]],
      control = get_control(dat, "set.control.ergm"),
      dynamic = FALSE
    )
  }
  nw <- dat[["nw"]]

  # Pull Network parameters
  dat[["nwparam"]] <- list()
  for (i in seq_along(x)) {
    dat[["nwparam"]][i] <- list(x[[i]][!(names(x[[i]]) %in% c("fit", "newnetwork"))])
    dat[["nwparam"]][[i]]["isTERGM"] <- all(x[[i]][["coef.diss"]][["duration"]] > 1)
  }

  ## Nodal Attributes Setup ##
  num <- network.size(nw[[1]])
  dat <- append_core_attr(dat, 1, num)

  # Pull in attributes on network
  nwattr.all <- list.vertex.attributes(nw[[1]])
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
  if (get_control(dat, "save.nwstats")) {
    dat <- initialize_nwstats(dat)
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

  ## Master Data List Setup ##
  dat <- create_dat_object(param, init, control)

  dat$stats$nwstats <- list()
  ## Network Setup ##
  # Initial network simulations
  dat[["nw"]] <- list()
  for (i in 1:3) {
    dat[["nw"]][[i]] <- simulate(
      x[[i]][["formula"]],
      coef = x[[i]][["coef.form.crude"]],
      basis = x[[i]][["newnetwork"]],
      constraints = x[[i]][["constraints"]],
      control = get_control(dat, "set.control.ergm"),
      dynamic = FALSE
    )
  }
  nw <- dat[["nw"]]

  # Pull Network parameters
  dat[["nwparam"]] <- list()
  for (i in seq_along(x)) {
    dat[["nwparam"]][i] <- list(x[[i]][!(names(x[[i]]) %in% c("fit", "newnetwork"))])
    dat[["nwparam"]][[i]]["isTERGM"] <- all(x[[i]][["coef.diss"]][["duration"]] > 1)
  }

  ## Nodal Attributes Setup ##
  num <- network.size(nw[[1]])
  dat <- append_core_attr(dat, 1, num)

  # Pull in attributes on network
  nwattr.all <- list.vertex.attributes(nw[[1]])
  nwattr.use <- nwattr.all[!nwattr.all %in% c("na", "vertex.names")]
  for (i in seq_along(nwattr.use)) {
    dat[["attr"]][[nwattr.use[i]]] <- get.vertex.attribute(nw[[1]], nwattr.use[i])
  }

  # Convert to tergmLite method
  dat <- init_tergmLite(dat)

  ## Infection Status and Time Modules
  dat <- init_status_covid_corporate(dat)

  ## Get initial prevalence
  dat <- prevalence_covid_corporate(dat, at = 1)

  # Network stats
  if (get_control(dat, "save.nwstats")) {
    dat <- initialize_nwstats(dat)
  }

  class(dat) <- "dat"
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
  eligible.case <- rep(NA, num)
  traced.cc <- rep(NA, num)
  quar <- rep(NA, num)
  clinical <- rep(NA, num)
  eligible.cc <- rep(NA, num)
  hospit <- rep(NA, num)
  symendTime <- rep(NA, num)
  statusTime.Ic <- rep(NA, num)
  statusTime <- rep(NA, num)
  statusTime[idsInf] <- 1
  dxStatus <- rep(0, num)
  dxTime <- rep(NA, num)
  iso.end <- rep(NA, num)
  tracedTime <- rep(NA, num)
  quarEnd <- rep(NA, num)
  vax <- rep(0, num)
  vax1Time <- rep(NA, num)
  vax2Time <- rep(NA, num)

  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "statusTime.Ic", statusTime.Ic)
  dat <- set_attr(dat, "symendTime", symendTime)
  dat <- set_attr(dat, "infTime", infTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "eligible.case", eligible.case)
  dat <- set_attr(dat, "eligible.cc", eligible.cc)
  dat <- set_attr(dat, "traced.cc", traced.cc)
  dat <- set_attr(dat, "quar", quar)
  dat <- set_attr(dat, "hospit", hospit)
  dat <- set_attr(dat, "dxStatus", dxStatus)
  dat <- set_attr(dat, "dxTime", dxTime)
  dat <- set_attr(dat, "iso.end", iso.end)
  dat <- set_attr(dat, "tracedTime", tracedTime)
  dat <- set_attr(dat, "quarEnd", quarEnd)
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vax1Time", vax1Time)
  dat <- set_attr(dat, "vax2Time", vax2Time)

  return(dat)
}

initialize_nwstats <- function(dat) {
  dat[["stats"]][["nwstats"]] <- list()
  for (i in seq_along(dat[["nwparam"]])) {
    new.nwstats <- attributes(dat$nw[[i]])$stats
    keep.cols <- which(!duplicated(colnames(new.nwstats)))
    new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
    dat$stats$nwstats[[i]] <- new.nwstats
  }
  return(dat)
}
