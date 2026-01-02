#' Get the list the attributes used by the model with their default value
#'
#' @inheritParams aging_msm
#'
#' @details
#' This list must be exhaustive. All attributes should get a default value (even
#' NA) here.
#'
#' @return A named list of all the attributes with default values.
#'
#' TODO: document all the attributes
get_default_attrs <- function(dat) { # all attributes should be listed here
  list(
    # network attributes
    status="s",
    deg_work=0,
    deg_school=0,
    deg_nonhome=0,
    hh.ids =0,
    age = get_param(dat, "arrival.age"), 
    age.grp = NA,
    # disease-related attributes
    vax.age.group = NA, 
    statusTime = NA,
    infTime = NA,
    clinical= NA,
    hospit= NA,
    dxStatus= 0,
    dxTime= NA,
    vax= 0,
    vax1Time= NA,
    vax2Time= NA,
    vax3Time= NA,
    vax4Time= NA,
    isolate= NA,
    isoTime= NA
   
  )
}

#' Update the attributes requiring computation for new nodes
#'
#' @inheritParams aging_msm
#' @param n_new The number of new nodes to update
#' @param post_init logical flag, TRUE if not called by the initialization
#' module.
#'
#' @details
#' This function takes care of all the attributes requiring some computations.
#' This includes the random assignment for race or the calculation of the
#' age groups. New attributes that need special assignment should be set here
#' AS WELL AS in the `get_default_attrs` function.
#'
#' @return
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
  
  
  n_attr <- c(n_attr, list(
    age.grp     =   cut(age, 
                        age.breaks, 
                        labels = age.grps, 
                        right = FALSE
                        ) |> as.character()
  ))
  
  for (attr_name in names(n_attr)) {
    dat <- set_attr(dat, attr_name, n_attr[[attr_name]], posit_ids = new_nodes_pid)
  }
  
  return(dat)
}

# add no.contact, nonhome
# init_attr  <- get_init(dat, "attr")
# no_contact <- 1L - init_attr$contact_attribute_Nonhome
# dat <- set_attr(dat, "no.contact",         no_contact)

# Generate the `late.tester` attributes
make_late_tester <- function(dat, race) {
  rates <- get_param(dat, "hiv.test.late.prob")[race]
  runif(length(rates)) < rates
}

# Generate the `ins.quot` attributes
make_ins_quot <- function(role.class) {
  ins_quot <- numeric(length(role.class))
  ins_quot[role.class == 0]  <- 1
  ins_quot[role.class == 1]  <- 0
  ins_quot[role.class == 2]  <- runif(sum(role.class == 2))
  return(ins_quot)
}

# Generate the `tt.traj` attributes
make_tt_traj <- function(dat, race, race_lvls) {
  partial <- get_param(dat, "tt.partial.supp.prob")
  full    <- get_param(dat, "tt.full.supp.prob")
  durable <- get_param(dat, "tt.durable.supp.prob")
  race.new <- race
  tt_traj <- numeric(length(race.new))
  for (r in race_lvls) {
    ids.race <- which(race == r)
    tt_traj[ids.race] <- sample(
      seq_len(3), length(ids.race), TRUE,
      c(partial[r], full[r], durable[r])
    )
  }
  return(tt_traj)
}

# Generate the `circ` attributes
make_circ <- function(dat, race, race_lvls) {
  hiv.circ.prob <- get_param(dat, "hiv.circ.prob")
  circ <- numeric(length(race))
  for (r in race_lvls) {
    ids.race <- which(race == r)
    circ[ids.race] <- runif(length(ids.race)) < hiv.circ.prob[r]
  }
  return(circ)
}

# Generate the `role.class` attributes
make_role_class <- function(dat, race, race_lvls) {
  ns <- get_param(dat, "netstats")$attr
  role_class <- numeric(length(race))
  for (r in race_lvls) {
    ids.race <- which(race == r)
    rc.probs <- prop.table(table(ns$role.class[ns$race == r]))
    role_class[ids.race] <- sample(0:2, length(ids.race), TRUE, rc.probs)
  }
  return(role_class)
}

make_prep_class <- function(dat, n_new) {
  prep.class <- sample(3, n_new, TRUE, get_param(dat, "prep.adhr.dist"))
  return(prep.class)
}