
#' @rdname moduleset-ship
#' @export
progress_covid_ship <- function(dat, at) {

  ## Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  statusTime <- dat$attr$statusTime
  infTime <- dat$attr$infTime
  clinical <- dat$attr$clinical
  age <- dat$attr$age

  ## Parameters
  prop.clinical <- dat$param$prop.clinical
  ea.rate <- dat$param$ea.rate
  ar.rate <- dat$param$ar.rate
  eip.rate <- dat$param$eip.rate
  ipic.rate <- dat$param$ipic.rate
  icr.rate <- dat$param$icr.rate

  ## Determine Subclinical (E to A) or Clinical (E to Ip to Ic) pathway
  ids.newInf <- which(active == 1 & status == "e" & statusTime <= at & is.na(clinical))
  num.newInf <- length(ids.newInf)
  if (num.newInf > 0) {
    age.group <- pmin((round(age[ids.newInf], -1)/10) + 1, 8)
    prop.clin.vec <- prop.clinical[age.group]
    if (any(is.na(prop.clin.vec))) stop("error in prop.clin.vec")
    vec.new.clinical <- rbinom(num.newInf, 1, prop.clinical)
    clinical[ids.newInf] <- vec.new.clinical
  }
  if (any(status == "e" & is.na(clinical))) browser()

  ## Subclinical Pathway
  # E to A: latent move to asymptomatic infectious
  num.new.EtoA <- 0
  ids.Es <- which(active == 1 & status == "e" & statusTime < at & clinical == 0)
  num.Es <- length(ids.Es)
  if (num.Es > 0) {
    vec.new.A <- which(rbinom(num.Es, 1, ea.rate) == 1)
    if (length(vec.new.A) > 0) {
      ids.new.A <- ids.Es[vec.new.A]
      num.new.EtoA <- length(ids.new.A)
      status[ids.new.A] <- "a"
      statusTime[ids.new.A] <- at
    }
  }

  # A to R: asymptomatic infectious move to recovered
  num.new.AtoR <- 0
  ids.A <- which(active == 1 & status == "a" & statusTime < at & clinical == 0)
  num.A <- length(ids.A)
  if (num.A > 0) {
    vec.new.R <- which(rbinom(num.A, 1, ar.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.A[vec.new.R]
      num.new.AtoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
    }
  }

  ## Clinical Pathway
  # E to Ip: latent move to preclinical infectious
  num.new.EtoIp <- 0
  ids.Ec <- which(active == 1 & status == "e" & statusTime < at & clinical == 1)
  num.Ec <- length(ids.Ec)
  if (num.Ec > 0) {
    vec.new.Ip <- which(rbinom(num.Ec, 1, eip.rate) == 1)
    if (length(vec.new.Ip) > 0) {
      ids.new.Ip <- ids.Ec[vec.new.Ip]
      num.new.EtoIp <- length(ids.new.Ip)
      status[ids.new.Ip] <- "ip"
      statusTime[ids.new.Ip] <- at
    }
  }

  # Ip to Ic: preclinical infectious move to clinical infectious
  num.new.IptoIc <- 0
  ids.Ip <- which(active == 1 & status == "ip" & statusTime < at & clinical == 1)
  num.Ip <- length(ids.Ip)
  if (num.Ip > 0) {
    vec.new.Ic <- which(rbinom(num.Ip, 1, ipic.rate) == 1)
    if (length(vec.new.Ic) > 0) {
      ids.new.Ic <- ids.Ip[vec.new.Ic]
      num.new.IptoIc <- length(ids.new.Ic)
      status[ids.new.Ic] <- "ic"
      statusTime[ids.new.Ic] <- at
    }
  }

  # Ic to R: clinical infectious move to recovered (if not mortality first)
  num.new.IctoR <- 0
  ids.Ic <- which(active == 1 & status == "ic" & statusTime < at & clinical == 1)
  num.Ic <- length(ids.Ic)
  if (num.Ic > 0) {
    vec.new.R <- which(rbinom(num.Ic, 1, icr.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.Ic[vec.new.R]
      num.new.IctoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
    }
  }

  ## Save updated status attribute
  dat$attr$status <- status
  dat$attr$statusTime <- statusTime
  dat$attr$clinical <- clinical

  ## Save summary statistics
  dat$epi$ea.flow[at] <- num.new.EtoA
  dat$epi$ar.flow[at] <- num.new.AtoR

  dat$epi$eip.flow[at] <- num.new.EtoIp
  dat$epi$ipic.flow[at] <- num.new.IptoIc
  dat$epi$icr.flow[at] <- num.new.IctoR

  return(dat)
}
