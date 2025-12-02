

#' @rdname moduleset-gmc19
#' @export
arrival <- function(dat, at) {
  # if (at>200) browser(), helper for debugging at a later time step
  
  ## Input
  # Attributes
  
  # Parameters
  a.rate   <- get_param(dat, "a.rate")
  
  ## Process
  num <- get_epi(dat, "num", at = 1)
  nNew <- rpois(1, a.rate * num) 
  
  

  
  ## Update Attr
  if (nNew > 0) {
   # dat <- init_new_nodes_attrs(dat, nNew)
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
  dat <- append_attr(dat, "hh.ids", newHH, nNew)

  # update household edgelist, arrivals module
  heads <- cbind((length(hh.ids) + 1):(length(hh.ids) + nNew), newHH) # arrival nodes + their hh.ids
  tails <- cbind(which(hh.ids %in% newHH), hh.ids[which(hh.ids %in% newHH)]) # existing nodes living with new nodes + their hh.ids
  new.edges <- merge(heads, tails, by.x = 2, by.y = 2)[, 2:3] # join by the 2nd column. connect each new node to existing nodes
  new.edgelist <- as.matrix(rbind(dat$el[[dat$num.nw]], setNames(new.edges, c(".head", ".tail")))) 
  attr(new.edgelist, 'n') <- attr(dat$el[[dat$num.nw]], 'n')
  dat$run$el[[dat$num.nw]] <- new.edgelist
  
  # update dat$run$net_attr
  
  
  return(dat)
  
}


