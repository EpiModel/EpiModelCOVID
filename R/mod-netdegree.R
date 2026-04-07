netdegree <- function(dat, at) {
  #if (at>20) browser()
  if (isFALSE(dat$param$compute.degree)) return(dat)
  
  degree_work    <- get_degree(dat$run$el[["work"]])
  degree_school  <- get_degree(dat$run$el[["school"]])
  degree_nonhome <- get_degree(dat$run$el[["nonhome"]])
  
  degree_hh <-
    tabulate(c(dat$run$el[[4]][,1], dat$run$el[[4]][,2]),
             nbins=attr(dat$run$el$school, "n") # number of active nodes at this time step, pulled from school layer
             )

  degree_total <- degree_work + degree_school + degree_nonhome + degree_hh
  
  # cross-layer bridging
  n_layers_active <- (degree_work > 0) + (degree_school > 0) + (degree_nonhome > 0) + (degree_hh > 0) # bridge count
  is_bridge <- n_layers_active >= 2 # bridge indicator
  
  # mean degree by age group 
  age.grp <- get_attr(dat, item = "age.grp")
  mean_degree_age <- aggregate(degree_total ~ age.grp, FUN = mean)
  
  # proportion of nodes that are cross-layer bridges  
  ## overall
  prop_bridge <- mean(is_bridge)
  ## by age group
  prop_bridge_age <-  prop.table(table(age.grp, is_bridge))
   
  # degree variance by age group 
  degree_var_age <-  aggregate(degree_total ~ age.grp, FUN = var)
  
  
  dat <- set_attr(dat, "degree_work", degree_work)
  dat <- set_attr(dat, "degree_school", degree_school)
  dat <- set_attr(dat, "degree_nonhome", degree_nonhome)
  dat <- set_attr(dat, "degree_hh", degree_hh)
  dat <- set_attr(dat, "degree_total", degree_total)
  dat <- set_attr(dat, "is_bridge", is_bridge)
  dat <- set_attr(dat, "n_layers_active", n_layers_active)
  
  dat <- set_epi(dat, "mean_degree_work", at, mean(degree_work))
  dat <- set_epi(dat, "mean_degree_school", at, mean(degree_school))
  dat <- set_epi(dat, "mean_degree_nonhome", at, mean(degree_nonhome))
  dat <- set_epi(dat, "mean_degree_hh", at, mean(degree_hh))
  dat <- set_epi(dat, "mean_degree_total", at, mean(degree_total))
  
  dat <- set_epi(dat, "prop_bridge", at, prop_bridge) 
  dat <- set_epi(dat, "prop_bridge_age1", at, prop_bridge_age[1,2])
  dat <- set_epi(dat, "prop_bridge_age2", at, prop_bridge_age[2,2])
  dat <- set_epi(dat, "prop_bridge_age3", at, prop_bridge_age[3,2])
  dat <- set_epi(dat, "prop_bridge_age4", at, prop_bridge_age[4,2])
  dat <- set_epi(dat, "prop_bridge_age5", at, prop_bridge_age[5,2])
  dat <- set_epi(dat, "prop_bridge_age6", at, prop_bridge_age[6,2])
  
  dat <- set_epi(dat, "mean_degree_age1", at, mean_degree_age[1,2])
  dat <- set_epi(dat, "mean_degree_age2", at, mean_degree_age[2,2])
  dat <- set_epi(dat, "mean_degree_age3", at, mean_degree_age[3,2])
  dat <- set_epi(dat, "mean_degree_age4", at, mean_degree_age[4,2])
  dat <- set_epi(dat, "mean_degree_age5", at, mean_degree_age[5,2])
  dat <- set_epi(dat, "mean_degree_age6", at, mean_degree_age[6,2])
  
  dat <- set_epi(dat, "degree_var_age1", at, degree_var_age[1,2])
  dat <- set_epi(dat, "degree_var_age2", at, degree_var_age[2,2])
  dat <- set_epi(dat, "degree_var_age3", at, degree_var_age[3,2])
  dat <- set_epi(dat, "degree_var_age4", at, degree_var_age[4,2])
  dat <- set_epi(dat, "degree_var_age5", at, degree_var_age[5,2])
  dat <- set_epi(dat, "degree_var_age6", at, degree_var_age[6,2])
  
  dat
}