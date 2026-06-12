

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
  
  

  
  # Update Attr
  if (nNew > 0) {
   dat <- init_new_nodes_attrs(dat, nNew)
   # After the default above, add houshold id to new nodes and create home edgelists
   dat <- set_home_attr_el(dat, at, nNew)

  }
  
  # Update Networks
  dat <- arrive_nodes(dat, nNew)

  
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

# Assign new nodes to households and create home edgelist for these nodes, run in arrivals module
set_home_attr_el <- function(dat, at, nNew) {
  new_nodes_pid <- length(get_attr(dat, "active")) - nNew + seq_len(nNew) # data frame position of new nodes
  
  age.grp <- get_attr(dat, "age.grp") # age of all nodes, include the new nodes (0)
  hh.ids <- get_attr(dat, "hh.ids" #  hh.ids, 0 for new nodes
                        )

  # correspond to all new nodes, for each of these nodes, randomly assign hh.ids
  newHH <- 
    sample(
      hh.ids[which(age.grp == "0-9y")], # existing hh.ids with nodes in "0-9y"
      nNew , # total number of new nodes in particular age grp
      replace = TRUE
      )
  
  # add hh.ids of newly arrivals to run$attr
  dat <- set_attr(dat, "hh.ids", newHH, posit_ids = new_nodes_pid)

  # Update the household edgelist: connect each new node to the EXISTING members
  # of its assigned household. Built with pure matrix ops. The previous version
  # used merge() + rbind() on data.frames, which made R call make.unique() over
  # the row names of the whole ~N-row home edgelist every step (O(N) per step;
  # ~26 s/step / ~56% of all sim time at N=117,810). The edge SET produced is the
  # same; only the data type changes.
  sel  <- which(hh.ids %in% newHH)        # existing nodes in the chosen households
  memb <- split(sel, hh.ids[sel])         # hh.id (chr) -> existing member node ids
  new.edges <- do.call(rbind, lapply(seq_len(nNew), function(i) {
    m <- memb[[as.character(newHH[i])]]   # existing members of this new node's household
    if (length(m)) cbind(new_nodes_pid[i], m) else NULL
  }))                                     # 2-col matrix (.head = new node, .tail = existing), or NULL
  if (!is.null(new.edges)) {
    dat$run$el[[dat$num.nw]] <-
      rbind(dat$run$el[[dat$num.nw]], unname(new.edges))  # matrix rbind: no data.frame, no make.unique
  }

  return(dat)

}


