#' @rdname moduleset-gmc19
#' @export
netdegree <- function(dat, at) {
  if (isFALSE(dat$param$compute.degree)) return(dat)
  
  degree_work    <- get_degree(dat$run$el[["work"]])
  degree_school  <- get_degree(dat$run$el[["school"]])
  degree_nonhome <- get_degree(dat$run$el[["nonhome"]])
  
  degree_hh <-
    tabulate(c(dat$run$el[[4]][,1], dat$run$el[[4]][,2]),
             nbins=attr(dat$run$el$school, "n") # number of active nodes at this time step, pulled from school layer
             )

  degree_total <- degree_work + degree_school + degree_nonhome + degree_hh

  # cross-layer (non-household) degree: prioritization metrics that exclude the
  # dense household clique, so "degree" is not dominated by household size.
  degree_xlayer <- degree_work + degree_school + degree_nonhome

  # degree quartile (rank-based, 4 equal-size groups) for outcome stratification
  n_nodes_deg <- length(degree_total)
  degree_quartile <- pmin(pmax(
    ceiling(rank(degree_total, ties.method = "first") / n_nodes_deg * 4), 1L), 4L)

  # cross-layer bridging
  n_layers_active <- (degree_work > 0) + (degree_school > 0) + (degree_nonhome > 0) + (degree_hh > 0) # bridge count
  is_bridge <- n_layers_active >= 2 # bridge indicator
  
  # mean degree by age group 
  age.grp <- get_attr(dat, item = "age.grp")
  mean_degree_age <- tapply(degree_total, age.grp, mean)
  
  # proportion of nodes that are cross-layer bridges  
  ## overall
  prop_bridge <- mean(is_bridge)
  ## by age group
  prop_bridge_age <-  tapply(is_bridge, age.grp, mean)
   
  # degree variance by age group 
  var_degree_age <-  tapply(degree_total, age.grp, var)
  
  
  dat <- set_attr(dat, "degree_work", degree_work)
  dat <- set_attr(dat, "degree_school", degree_school)
  dat <- set_attr(dat, "degree_nonhome", degree_nonhome)
  dat <- set_attr(dat, "degree_hh", degree_hh)
  dat <- set_attr(dat, "degree_total", degree_total)
  dat <- set_attr(dat, "degree_xlayer", degree_xlayer)
  dat <- set_attr(dat, "degree_quartile", degree_quartile)
  dat <- set_attr(dat, "is_bridge", is_bridge)
  dat <- set_attr(dat, "n_layers_active", n_layers_active)
  
  dat <- set_epi(dat, "mean_degree_work", at, mean(degree_work))
  dat <- set_epi(dat, "mean_degree_school", at, mean(degree_school))
  dat <- set_epi(dat, "mean_degree_nonhome", at, mean(degree_nonhome))
  dat <- set_epi(dat, "mean_degree_hh", at, mean(degree_hh))
  dat <- set_epi(dat, "mean_degree_total", at, mean(degree_total))
  
  dat <- set_epi(dat, "prop_bridge", at, prop_bridge) 
  dat <- set_epi(dat, "prop_bridge_age1", at, prop_bridge_age[[1]])
  dat <- set_epi(dat, "prop_bridge_age2", at, prop_bridge_age[[2]])
  dat <- set_epi(dat, "prop_bridge_age3", at, prop_bridge_age[[3]])
  dat <- set_epi(dat, "prop_bridge_age4", at, prop_bridge_age[[4]])
  dat <- set_epi(dat, "prop_bridge_age5", at, prop_bridge_age[[5]])
  dat <- set_epi(dat, "prop_bridge_age6", at, prop_bridge_age[[6]])
  
  dat <- set_epi(dat, "mean_degree_age1", at, mean_degree_age[[1]])
  dat <- set_epi(dat, "mean_degree_age2", at, mean_degree_age[[2]])
  dat <- set_epi(dat, "mean_degree_age3", at, mean_degree_age[[3]])
  dat <- set_epi(dat, "mean_degree_age4", at, mean_degree_age[[4]])
  dat <- set_epi(dat, "mean_degree_age5", at, mean_degree_age[[5]])
  dat <- set_epi(dat, "mean_degree_age6", at, mean_degree_age[[6]])
  
  dat <- set_epi(dat, "var_degree_age1", at, var_degree_age[[1]])
  dat <- set_epi(dat, "var_degree_age2", at, var_degree_age[[2]])
  dat <- set_epi(dat, "var_degree_age3", at, var_degree_age[[3]])
  dat <- set_epi(dat, "var_degree_age4", at, var_degree_age[[4]])
  dat <- set_epi(dat, "var_degree_age5", at, var_degree_age[[5]])
  dat <- set_epi(dat, "var_degree_age6", at, var_degree_age[[6]])
  
  dat
}
