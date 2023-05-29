#' @rdname moduleset-vaxDecisions
#' @export
init_covid_vax_decisions <- function(x, param, init, control, s) {

  ## Master Data List Setup ##
  dat <- create_dat_object(param, init, control)

  ## Network Setup ##
  # Initial network simulations
  dat[["nw"]] <- list()
  for (i in 1:length(x)) {
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
  
  #Add household network edgelist
  dat$el[[length(x) + 1]] <- as.matrix(dat$param$hh.pairs)
  attr(dat$el[[length(x) + 1]], 'n') <- num

  ## Infection Status and Time Modules
  dat <- init_status_covid_vax_decisions(dat)

  ## Get initial prevalence
  dat <- prevalence_covid_vax_decisions(dat, at = 1)

  # Network stats
  if (get_control(dat, "save.nwstats")) {
    dat <- initialize_nwstats(dat)
  }

  class(dat) <- "dat"
  return(dat)
}


init_status_covid_vax_decisions <- function(dat) {

  e.num <- get_init(dat, "e.num")

  active <- get_attr(dat, "active")
  num <- sum(active)

  ## Disease status
  status <- rep("s", num)
  if (e.num > 0) {
    status[sample(which(active == 1), size = e.num)] <- "e"
  }

  dat <- set_attr(dat, "status", status)
  
  # Age group for vaccination processes
  age <- get_attr(dat, "age")
  
  vax.age.group <- rep(NA, length(age))
  vax.age.group[age < 5] <- 1
  vax.age.group[age >= 5 & age < 18] <- 2
  vax.age.group[age >= 18 & age < 50] <- 3
  vax.age.group[age >= 50 & age < 65] <- 4
  vax.age.group[age >= 65] <- 5
  
  dat <- set_attr(dat, "vax.age.group", vax.age.group)
  
  ## Vaccine willingness vs. resistance
  vax.willing.prob <- get_param(dat, "vax.willing.prob")
  vaxType <- rep(NA, num)
  willing.time <- rep(NA, num)
  idsAdults <- which(vax.age.group >= 3)
  vaxType.adults <- rbinom(length(idsAdults), 1, 
                           vax.willing.prob[vax.age.group[idsAdults] - 2])
  vaxType[idsAdults] <- vaxType.adults
  willing.time[which(vaxType == 1)] <- 1
  dat <- set_attr(dat, "vaxType", vaxType)
  dat <- set_attr(dat, "willing.time", willing.time)

  # Infection Time and related attributes
  idsInf <- which(status == "e")
  infTime <- rep(NA, num)
  symptStartTime <- rep(NA, num)
  clinical <- rep(NA, num)
  hospit <- rep(NA, num)
  statusTime <- rep(NA, num)
  statusTime[idsInf] <- 1
  dxStatus <- rep(0, num)
  dxTime <- rep(NA, num)
  vax <- rep(0, num)
  vax1Time <- rep(NA, num)
  vax2Time <- rep(NA, num)
  vax3Time <- rep(NA, num)
  vax4Time <- rep(NA, num)
  vaxSE <- rep(NA, num)

  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "infTime", infTime)
  dat <- set_attr(dat, "symptStartTime", symptStartTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "hospit", hospit)
  dat <- set_attr(dat, "dxStatus", dxStatus)
  dat <- set_attr(dat, "dxTime", dxTime)
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vax1Time", vax1Time)
  dat <- set_attr(dat, "vax2Time", vax2Time)
  dat <- set_attr(dat, "vax3Time", vax3Time)
  dat <- set_attr(dat, "vax4Time", vax4Time)
  dat <- set_attr(dat, "vaxSE", vaxSE)

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
