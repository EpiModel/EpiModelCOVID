netdegree <- function(dat, at) {
  #if (at>=1) browser()
  if (isFALSE(dat$param$compute.degree)) return(dat)
  
  degree_work    <- get_degree(dat$run$el[["work"]])
  degree_school  <- get_degree(dat$run$el[["school"]])
  degree_nonhome <- get_degree(dat$run$el[["nonhome"]])
  
  if (is.null(dat$degree_hh_cache)) {
    dat$degree_hh_cache <-
    tabulate(c(dat$run$el[[4]][,1], dat$run$el[[4]][,2]), 
             nbins=sum(dat$run$attr$active)
             )
     
  }
  degree_hh <- dat$degree_hh_cache
  
  degree_total <- degree_work + degree_school + degree_nonhome + degree_hh
  
  dat <- set_attr(dat, "degree_work", degree_work)
  dat <- set_attr(dat, "degree_school", degree_school)
  dat <- set_attr(dat, "degree_nonhome", degree_nonhome)
  dat <- set_attr(dat, "degree_hh", degree_hh)
  dat <- set_attr(dat, "degree_total", degree_total)
  
  dat
}