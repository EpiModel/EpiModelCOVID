
#' @rdname moduleset-corporate
#' @export
prevalence_covid_corporate <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  vax <- get_attr(dat, "vax")
  isolate <- get_attr(dat, "isolate")
  non.office <- get_attr(dat, "non.office")

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
  dat <- set_epi(dat, "num.w", at, sum(active == 1 & non.office == 0))

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
  
  dat <- set_epi(dat, "s.num.w", at, sum(active == 1 & status == "s" & non.office == 0))
  dat <- set_epi(dat, "e.num.w", at, sum(active == 1 & status == "e" & non.office == 0))
  dat <- set_epi(dat, "a.num.w", at, sum(active == 1 & status == "a" & non.office == 0))
  dat <- set_epi(dat, "ip.num.w", at, sum(active == 1 & status == "ip" & non.office == 0))
  dat <- set_epi(dat, "ic.num.w", at, sum(active == 1 & status == "ic" & non.office == 0))
  dat <- set_epi(dat, "r.num.w", at, sum(active == 1 & status == "r" & non.office == 0))
  dat <- set_epi(dat, "h.num.w", at, sum(active == 1 & status == "h" & non.office == 0))
  dat <- set_epi(dat, "v1.num.w", at, sum(active == 1 & vax == 1 & non.office == 0))
  dat <- set_epi(dat, "v2.num.w", at, sum(active == 1 & vax == 2 & non.office == 0))
  dat <- set_epi(dat, "v3.num.w", at, sum(active == 1 & vax == 3 & non.office == 0))
  dat <- set_epi(dat, "v4.num.w", at, sum(active == 1 & vax == 4 & non.office == 0))

  dat <- set_epi(dat, "iso.mild.num", at, sum(active == 1 & isolate == 1, na.rm = TRUE))
  dat <- set_epi(dat, "iso.sev.num", at, sum(active == 1 & isolate == 2, na.rm = TRUE))
  dat <- set_epi(dat, "mask.post.num", at, sum(active == 1 & isolate == 3, na.rm = TRUE))
  dat <- set_epi(dat, "mask.exp.num", at, sum(active == 1 & isolate == 4, na.rm = TRUE))
  
  dat <- set_epi(dat, "iso.mild.num.w", at, sum(active == 1 & isolate == 1 & non.office == 0, na.rm = TRUE))
  dat <- set_epi(dat, "iso.sev.num.w", at, sum(active == 1 & isolate == 2 & non.office == 0, na.rm = TRUE))
  dat <- set_epi(dat, "mask.post.num.w", at, sum(active == 1 & isolate == 3 & non.office == 0, na.rm = TRUE))
  dat <- set_epi(dat, "mask.exp.num.w", at, sum(active == 1 & isolate == 4 & non.office == 0, na.rm = TRUE))

  return(dat)
}
