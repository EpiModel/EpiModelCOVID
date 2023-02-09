
#' @rdname moduleset-netjail
#' @export
prevalence_covid_netjail <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  vax <- get_attr(dat, "vax")

  nsteps <- get_control(dat, "nsteps")
  
  # Initialize Outputs
  var.names <- c("num", "s.num", "e.num", "a.num", "ip.num", "ic.num", "r.num",
                 "h.num", "v.num.incomplete", "v.num.partial", "v.num.immune")
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
  dat <- set_epi(dat, "v.num.incomplete", at, sum(active == 1 & vax == 1))
  dat <- set_epi(dat, "v.num.partial", at, sum(active == 1 & vax == 2))
  dat <- set_epi(dat, "v.num.immune", at, sum(active == 1 & vax == 3))

  return(dat)
}
