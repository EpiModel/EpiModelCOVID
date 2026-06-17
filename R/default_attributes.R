#' Get the list the attributes used by the model with their default value
get_default_attrs <- function(dat) { # all attributes should be listed here
  list(
    # network attributes, binary
    deg_work=0,
    deg_school=0,
    deg_nonhome=0,
    # network attributes, integer
    degree_work=0,
    degree_school=0,
    degree_nonhome=0,
    degree_hh=0,
    degree_total=0,
    degree_xlayer=0,
    degree_quartile=NA,
    # cross-layer bridge indicators
    is_bridge = FALSE,
    n_layers_active =0,
    # houshold id and age
    hh.ids =0, 
    age = get_param(dat, "arrival.age"), 
    age.grp = NA,
    # disease-related attributes
    vax.age.group = NA, 
    statusTime = NA,
    infTime = NA,
    clinical= NA,
    hospit= NA,
    status="s",
    dxStatus= 0,
    dxTime= NA,
    vax= 0,
    vax1Time= NA,
    vax2Time= NA,
    vax3Time= NA,
    vax4Time= NA,
    last.dose.time = NA,
    isolate= NA,
    isoTime= NA
   
  )
}

#' Update the attributes requiring computation for new nodes
#' This function returns the `dat` object with updated attributes for the new
#' nodes.
make_computed_attrs <- function(dat, n_new, post_init) {
  new_nodes_pid <- length(get_attr(dat, "active")) - n_new + seq_len(n_new)

  
  n_attr <- list()
  if (post_init) { # after the initialization
    # Disease status and related
    #  n_attr$status <- rep("s",n_new)
    
  } else {  # at the initialization
    # Disease status 
    e.num <- get_init(dat, "e.num")
    
    active <- get_attr(dat, "active")
    num <- sum(active)
    
    status <- get_attr(dat, "status")
    if (e.num > 0) {
      status[sample(which(active == 1), size = e.num)] <- "e"
    }
    
    n_attr$status <- status
    
    
    # Infection Time 
    idsInf <- which(status == "e")
    statusTime <- get_attr(dat, "statusTime")
    statusTime[idsInf] <- 1
    
    n_attr$statusTime <- statusTime
  }
  
  age <- get_attr(dat, "age", posit_ids = new_nodes_pid)
  age.breaks <- get_param(dat, "age.breaks")
  age.grps <- get_param(dat, "age.grps")

  # Infant flag: an age-resolution refinement (issue #23) that does NOT re-bin the
  # network age groups. Infants stay in the youngest contact band (their contacts
  # are household-dominated) but are flagged so the disease modules can apply
  # infant-specific severity (RSV) and the outputs can report the infant burden
  # separately rather than diluting it across the coarse youngest band. Absent
  # param -> infant.age.max = 1 year. Recomputed for arrivals so it stays correct.
  infant.age.max <- get_param(dat, "infant.age.max", override.null.error = TRUE)
  if (is.null(infant.age.max)) infant.age.max <- 1

  n_attr <- c(n_attr, list(
    age.grp     =   cut(age,
                        age.breaks,
                        labels = age.grps,
                        right = FALSE
                        ) |> as.character(),
    is_infant   =   age < infant.age.max
  ))
  
  for (attr_name in names(n_attr)) {
    dat <- set_attr(dat, attr_name, n_attr[[attr_name]], posit_ids = new_nodes_pid)
  }
  
  return(dat)
}

