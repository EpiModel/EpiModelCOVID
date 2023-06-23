
#' @rdname moduleset-corporate
#' @export
prevalence_covid_corporate <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  vax <- get_attr(dat, "vax")
  isolate <- get_attr(dat, "isolate")

  # Initialize Outputs
  var.names <- c("num", "s.num", "e.num", "a.num", "ip.num", "ic.num", "r.num",
                 "h.num", "v1.num", "v2.num", "v3.num", "v4.num",
                 "iso.mild.num","iso.sev.num", "mask.post.num","mask.exp.num")
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
  dat <- set_epi(dat, "v1.num", at, sum(active == 1 & vax == 1))
  dat <- set_epi(dat, "v2.num", at, sum(active == 1 & vax == 2))
  dat <- set_epi(dat, "v3.num", at, sum(active == 1 & vax == 3))
  dat <- set_epi(dat, "v4.num", at, sum(active == 1 & vax == 4))

  dat <- set_epi(dat, "iso.mild.num", at, sum(active == 1 & isolate == 1))
  dat <- set_epi(dat, "iso.sev.num", at, sum(active == 1 & isolate == 2))
  dat <- set_epi(dat, "mask.post.num", at, sum(active == 1 & isolate == 3))
  dat <- set_epi(dat, "mask.exp.num", at, sum(active == 1 & isolate == 4))

  return(dat)
}
