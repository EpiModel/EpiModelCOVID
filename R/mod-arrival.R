

#' @rdname moduleset-gmc19
#' @export
arrival <- function(dat, at) {
  # if (at>200) browser()
  
  ## Input
  # Attributes
  
  ## Process
  # Stationary-population demography: replace this step's departures so the
  # population size stays constant. The age-structured death rate makes a fixed
  # per-capita a.rate drift (Phase 2 review), so arrivals track departures, which
  # the departures module records as d.flow before this module runs.
  nDep <- get_epi(dat, "d.flow", at = at)
  nNew <- if (is.null(nDep) || is.na(nDep)) 0L else as.integer(nDep)
  
  

  
  # Update Attr. The new nodes join the networks in nwupdate.FUN, which runs
  # after this module, so this module only appends their attributes. On the
  # household clique layer, nwupdate.FUN places them with the layer's arrival
  # rule (hh_join_young) and connects them to their co-residents.
  if (nNew > 0) {
   dat <- init_new_nodes_attrs(dat, nNew)
  }

  ## Output
  dat <- set_epi(dat, "nNew", at, nNew)
  
  return(dat)
}



init_new_nodes_attrs <- function(dat, n_new) {
  current_timestep <- get_current_timestep(dat)
  # Core attributes (necessary for EpiModel)
  dat <- append_core_attr(dat, current_timestep, n_new)
  def_attrs <- get_default_attrs(dat)
  for (attr_name in names(def_attrs)) {
    dat <- append_attr(dat, attr_name, def_attrs[[attr_name]], n_new)
  }
  dat <- make_computed_attrs(dat, n_new, post_init = TRUE)
  return(dat)
}

#' @param new_ids Positional ids of the nodes arriving at this time step.
#' @param network Index of the household clique layer.
#'
#' @details
#' \code{hh_join_young} is the \code{arrivals.FUN} of the household layer, built
#' with \code{netclique(nw, "hh.ids", arrivals = "join", arrivals.FUN =
#' hh_join_young)}. Each newborn joins the household of a uniformly drawn
#' existing 0-9y member, so households are drawn in proportion to their number
#' of 0-9y members. EpiModel then connects the newborn to every active member of
#' that household, including other newborns placed there in the same step.
#'
#' @rdname moduleset-gmc19
#' @export
hh_join_young <- function(dat, at, new_ids, network) {
  active <- get_attr(dat, "active")
  age.grp <- get_attr(dat, "age.grp")
  hh.ids <- get_attr(dat, "hh.ids")

  # Existing households that contain a young (0-9y) member; new arrivals (age 0)
  # are placed into one of these. The pool leaves out this step's newborns
  # (their hh.ids is not set yet) and this step's dead, who are still in the
  # attribute list until nwupdate.FUN removes them after this call. Guard the
  # degenerate case (issue #39): if no 0-9y node exists the pool is empty, so
  # return NA and leave the new nodes in no household (isolates on the layer)
  # rather than misassigning them. Index via sample.int so a length-1 pool is
  # drawn from correctly (base sample() would treat a single numeric value as
  # 1:value).
  existing <- active == 1
  existing[new_ids] <- FALSE
  hh_pool <- hh.ids[existing & age.grp == "0-9y" & !is.na(hh.ids)]
  if (length(hh_pool) == 0) {
    return(rep(NA, length(new_ids)))
  }
  hh_pool[sample.int(length(hh_pool), length(new_ids), replace = TRUE)]
}
