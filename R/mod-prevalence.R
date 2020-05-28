
#' @rdname moduleset-ship
#' @export
prevalence_covid_ship <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  type <- dat$attr$type

  nsteps <- dat$control$nsteps

  # Initialize Outputs
  var.names <- c("num", "s.num", "e.num", "a.num", "ip.num", "ic.num", "r.num",
                 "i.pass.num", "i.crew.num",
                 "se.flow", "ea.flow", "ar.flow", "Rt",
                 "eip.flow", "ipic.flow", "icr.flow",
                 "d.flow", "d.ic.flow", "exit.flow", "nDx", "nDx.pos", "nDx.pos.sympt",
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

  dat$epi$i.pass.num[at] <- sum(active == 1 & status %in% c("ip", "ic", "a") & type == "p")
  dat$epi$i.crew.num[at] <- sum(active == 1 & status %in% c("ip", "ic", "a") & type == "c")

  dat$epi$meanAge[at] <- mean(dat$attr$age, na.rm = TRUE)
  dat$epi$meanClinic[at] <- mean(dat$attr$clinical, na.rm = TRUE)

  return(dat)
}
