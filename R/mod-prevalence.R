
#' @rdname moduleset-ship
#' @export
prevalence_covid_ship <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  type <- dat$attr$type
  dxStatus <- dat$attr$dxStatus

  nsteps <- dat$control$nsteps

  # Initialize Outputs
  var.names <- c("num", "s.num", "e.num", "a.num", "ip.num", "ic.num", "r.num",
                 "i.pass.num", "i.crew.num",
                 "dx.pos.num",
                 "se.flow", "ea.flow", "ar.flow", "Rt",
                 "eip.flow", "ipic.flow", "icr.flow",
                 "d.flow", "d.ic.flow", "exit.flow",
                 "nDx", "nDx.pos", "nDx.pos.sympt", "nDx.pos.fn",
                 "se.pp.flow", "se.pc.flow", "se.cp.flow", "se.cc.flow",
                 "meanAge", "meanClinic")
  if (at == 1) {
    for (i in 1:length(var.names)) {
      dat$epi[[var.names[i]]] <- rep(0, nsteps)
    }
  }

  # Update Outputs
  dat$epi$num[at] <- sum(active == 1)

  dat$epi$s.num[at] <- sum(active == 1 & status == "s")
  dat$epi$e.num[at] <- sum(active == 1 & status == "e")
  dat$epi$a.num[at] <- sum(active == 1 & status == "a")
  dat$epi$ip.num[at] <- sum(active == 1 & status == "ip")
  dat$epi$ic.num[at] <- sum(active == 1 & status == "ic")
  dat$epi$r.num[at] <- sum(active == 1 & status == "r")

  dat$epi$dx.pos.num[at] <- sum(active == 1 & dxStatus == 2, na.rm = TRUE)

  dat$epi$i.pass.num[at] <- sum(active == 1 & status %in% c("ip", "ic", "a") & type == "p")
  dat$epi$i.crew.num[at] <- sum(active == 1 & status %in% c("ip", "ic", "a") & type == "c")

  dat$epi$meanAge[at] <- mean(dat$attr$age, na.rm = TRUE)
  dat$epi$meanClinic[at] <- mean(dat$attr$clinical, na.rm = TRUE)

  return(dat)
}

#' @rdname moduleset-corporate
#' @export
prevalence_covid_corporate <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  vax <- get_attr(dat, "vax")

  nsteps <- get_control(dat, "nsteps")

  # Initialize Outputs
  var.names <- c("num", "s.num", "e.num", "a.num", "ip.num", "ic.num", "r.num",
                 "h.num", "v1.num", "v2.num")
  if (at == 1) {
    for (i in seq_along(var.names)) {
      dat <- add_epi(dat, var.names[i])
    }
  }

  # Update Outputs
  dat <- set_epi(dat, "num", at, sum(active == 1))

  dat <- set_epi(dat, "s.num", at, sum(active == 1 & status == "s"))
  dat <- set_epi(dat, "e.num", at, sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "a.num", at, sum(active == 1 & status == "a"))
  dat <- set_epi(dat, "ip.num", at, sum(active == 1 & status == "ip"))
  dat <- set_epi(dat, "ic.num", at, sum(active == 1 & status == "ic"))
  dat <- set_epi(dat, "r.num", at, sum(active == 1 & status == "r"))
  dat <- set_epi(dat, "h.num", at, sum(active == 1 & status == "h"))
  dat <- set_epi(dat, "v1.num", at, sum(active == 1 & status == "s" & vax == 1))
  dat <- set_epi(dat, "v2.num", at, sum(active == 1 & status == "s" & vax == 3))


  return(dat)
}


#' @rdname moduleset-contacttrace
#' @export
prevalence_covid_contacttrace <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  # vax <- get_attr(dat, "vax")

  nsteps <- get_control(dat, "nsteps")

  # Initialize Outputs
  var.names <- c("num", "s.num", "e.num", "a.num", "ip.num", "ic.num", "r.num",
                 "h.num", "icu.num")
  if (at == 1) {
    for (i in seq_along(var.names)) {
      dat <- add_epi(dat, var.names[i])
    }
  }

  # Update Outputs
  dat <- set_epi(dat, "num", at, sum(active == 1))

  dat <- set_epi(dat, "s.num", at, sum(active == 1 & status == "s"))
  dat <- set_epi(dat, "e.num", at, sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "a.num", at, sum(active == 1 & status == "a"))
  dat <- set_epi(dat, "ip.num", at, sum(active == 1 & status == "ip"))
  dat <- set_epi(dat, "ic.num", at, sum(active == 1 & status == "ic"))
  dat <- set_epi(dat, "r.num", at, sum(active == 1 & status == "r"))
  dat <- set_epi(dat, "h.num", at, sum(active == 1 & status == "h"))
  dat <- set_epi(dat, "icu.num", at, sum(active == 1 & status == "icu"))
  # dat <- set_epi(dat, "v1.num", at, sum(active == 1 & status == "s" & vax == 1))
  # dat <- set_epi(dat, "v2.num", at, sum(active == 1 & status == "s" & vax == 2))

  return(dat)
}

#' @rdname moduleset-boost
#' @export
prevalence_covid_boost <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  vax <- get_attr(dat, "vax")

  nsteps <- get_control(dat, "nsteps")

  # Initialize Outputs
  var.names <- c("num", "s.num", "e.num", "a.num", "ip.num", "ic.num", "r.num",
                 "h.num", "v1.num", "v2.num", "v3.num")
  if (at == 1) {
    for (i in seq_along(var.names)) {
      dat <- add_epi(dat, var.names[i])
    }
  }

  # Update Outputs
  dat <- set_epi(dat, "num", at, sum(active == 1))

  dat <- set_epi(dat, "s.num", at, sum(active == 1 & status == "s"))
  dat <- set_epi(dat, "e.num", at, sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "a.num", at, sum(active == 1 & status == "a"))
  dat <- set_epi(dat, "ip.num", at, sum(active == 1 & status == "ip"))
  dat <- set_epi(dat, "ic.num", at, sum(active == 1 & status == "ic"))
  dat <- set_epi(dat, "r.num", at, sum(active == 1 & status == "r"))
  dat <- set_epi(dat, "h.num", at, sum(active == 1 & status == "h"))
  dat <- set_epi(dat, "v1.num", at, sum(active == 1 & status == "s" & vax == 1))
  dat <- set_epi(dat, "v2.num", at, sum(active == 1 & status == "s" & vax == 3))
  dat <- set_epi(dat, "v3.num", at, sum(active == 1 & status == "s" & vax == 5))


  return(dat)
}
