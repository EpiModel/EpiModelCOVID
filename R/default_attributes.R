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
    deg_work=0,
    deg_school=0,
    deg_nonhome=0,
    no.contact=0, 
    age = get_param(dat, "arrival.age"), 
    age.grp = NA
   
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

  # ns <- get_param(dat, "netstats")
  # age_breaks <- ns$demog$age.breaks
  # race_dist <- prop.table(table(ns$attr$race))
  # race_lvls <- as.numeric(names(race_dist))
  
  n_attr <- list()
  if (post_init) { # after the initialization
  
   n_attr$deg_work <- rep(0,n_new)
   n_attr$deg_school <- rep(0,n_new)
   n_attr$deg_nonhome <- rep(0,n_new)
   n_attr$no.contact<- rep(1,n_new)
    
  } else {  # at the initialization
    n_attr$deg_work <- get_attr(dat, "deg_work", posit_ids = new_nodes_pid)
    n_attr$deg_school <- get_attr(dat, "deg_school", posit_ids = new_nodes_pid)
    n_attr$deg_nonhome <- get_attr(dat, "deg_nonhome", posit_ids = new_nodes_pid)
    n_attr$no.contact <- get_attr(dat, "no.contact", posit_ids = new_nodes_pid)
    
    
    # n_attr$race <- get_attr(dat, "race", posit_ids = new_nodes_pid)
    # n_attr$role.class <- get_attr(dat, "role.class", posit_ids = new_nodes_pid)
    # n_attr$risk.grp <- get_attr(dat, "risk.grp", posit_ids = new_nodes_pid)
  }
  
  age <- get_attr(dat, "age", posit_ids = new_nodes_pid)
  
  n_attr <- c(n_attr, list(
    # late.tester = make_late_tester(dat, n_attr$race),
    # circ        = make_circ(dat, n_attr$race, race_lvls),
    age.grp     =  rep("0-9y", length(age)) #cut(age, age_breaks, labels = FALSE, right = FALSE),
    # ins.quot    = make_ins_quot(n_attr$role.class),
    # tt.traj     = make_tt_traj(dat, n_attr$race, race_lvls),
    # prep.class   = make_prep_class(dat, length(n_attr$race)),
    # last.neg.test = get_attr(dat, "entrTime", posit_ids = new_nodes_pid)
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