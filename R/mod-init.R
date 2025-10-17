
#' @rdname moduleset-ship
#' @export
init_covid_ship <- function(x, param, init, control, s) {

  # Master Data List
  dat <- create_dat_object(param, init, control)

  dat <- init_nets(dat, x)

  # simulate first time step
  dat <- sim_nets_t1(dat)
  dat <- summary_nets(dat, at = 1L)

  ## Infection Status and Time Modules
  dat <- init_status_covid_ship(dat)

  ## Get initial prevalence
  dat <- prevalence_covid_ship(dat, at = 1)

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

  dat <- init_nets(dat, x)

  # simulate first time step
  dat <- sim_nets_t1(dat)
  dat <- summary_nets(dat, at = 1L)

  ## Infection Status and Time Modules
  dat <- init_status_covid_corporate(dat)

  ## Get initial prevalence
  dat <- prevalence_covid_corporate(dat, at = 1)

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


#' @rdname moduleset-gmc19
#' @export
init_gmc19 <- function(x, param, init, control, s) {
  # talk to Adrien on the home layer. more stream line way to pull the attribute. have attribute defined in dat.
  ## Master Data List Setup ##
  dat <- create_dat_object(param, init, control)
  dat <- init_nets(dat, x)
  
  
  # define initial nodal attributes
  node_age_grp <- as.character(init$attr$node.age.grp)           
  node_age     <- as.numeric  (init$attr$node.age)               
  deg_work   <- as.integer  (init$attr$contact_attribute_Work) 
  deg_school <- as.integer  (init$attr$contact_attribute_School)
  deg_nonhome   <- as.integer  (init$attr$contact_attribute_Nonhome)
  
  no_contact <- 1L - deg_nonhome
  

  ## set attributes
  dat <- set_attr(dat, "age.grp", node_age_grp)
  dat <- set_attr(dat, "age",     node_age)
  dat <- set_attr(dat, "deg_school", deg_school)
  dat <- set_attr(dat, "deg_work",   deg_work)
  dat <- set_attr(dat, "no.contact",         no_contact)
  
  # simulate first time step
  dat <- sim_nets_t1(dat)
  dat <- summary_nets(dat, at = 1L)
  
  ## add home layer edgelist
  dat$num.nw <- dat$num.nw + 1
  
  dat$nw[["home"]] <- dat$param$nw_home   
  dat$num.nw <- length(dat$nw)
  
  # addition to keep a place holder for home layer
  # keep HOME last (helps with any "num.nw - 1" tricks)
  
  nm <- c("school",  "work"  ,  "nonhome" ,"home"   )
  
  # helper: ensure a list has one slot per layer, with names
  ensure_slots <- function(lst, n, nm) {
    if (is.null(lst)) lst <- vector("list", n)
    if (length(lst) < n) length(lst) <- n
    names(lst) <- nm
    lst
  }
  
  # grow/align per-layer bookkeeping so [[home]] exists
  dat$run$el   <- ensure_slots(dat$run$el,   dat$num.nw, nm)  # cumulative edgelist
  dat$control$nwstats.formula[[4]] <- NA
  
  # end addition
  
  ## Infection Status and Time Modules
  dat <- init_status_covid_corporate(dat)
  
  ## Get initial prevalence
  dat <- prevalence_covid_corporate(dat, at = 1)
  
  return(dat)
}
