#' @rdname moduleset-gmc19
#' @export
nwupdate_gmc19 <- function(dat, at) {
  # Network membership. The departures module marks this step's deaths
  # (active = 0, exitTime = at) and the arrivals module appends the newborns'
  # attributes; this module, which runs right after the arrivals module, adds
  # the newborns to every layer and then removes the dead from every layer. On
  # a clique layer arrive_nodes() also places the newborns under the layer's
  # arrival rule.
  #
  # This is EpiModel's nwupdate.net() under tergmLite without its attribute
  # fill-in. get_default_attrs() sets every attribute for arrivals, so there is
  # nothing to fill in, but nwupdate.net() builds the fill-in distribution on
  # every step with arrivals by tabling every nodal attribute (continuous age
  # and the event times included). At N = 117,810 that took about 0.35 s per
  # step, which is as long as the rest of the step.
  if (!isTRUE(get_control(dat, "tergmLite"))) {
    stop("nwupdate_gmc19 supports tergmLite = TRUE only")
  }
  active <- get_attr(dat, "active")
  entrTime <- get_attr(dat, "entrTime")
  exitTime <- get_attr(dat, "exitTime")

  if (length(unique(lengths(get_attr_list(dat)))) != 1) {
    stop("nodal attributes have unequal lengths after the arrivals module; ",
         "every attribute must be set for new nodes in get_default_attrs()")
  }

  arrivals <- which(active == 1 & entrTime == at)
  departures <- which(active == 0 & exitTime == at)
  dat <- arrive_nodes(dat, length(arrivals))
  dat <- depart_nodes(dat, departures)

  return(dat)
}
