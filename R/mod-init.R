
source("~/Documents/GitHub/EpiModelCOVID/R/default_attributes.R")

#' @rdname moduleset-gmc19
#' @export
init_gmc19 <- function(x, param, init, control, s) {
  
  ## Master Data List Setup ##
  dat <- create_dat_object(param, init, control)

  ## network and stats initialization
  dat <- init_nets(dat, x)   
  
  ## Initialize all remaining attributes
  dat <- init_attrs(dat)
  dat <- overwrite_attrs(dat)
  
  # Add household network edgelist
  ## network index
  dat$num.nw <- dat$num.nw + 1
  ## edgelist 
  dat$run$el[[dat$num.nw]] <- as.matrix(dat$param$hh.pairs) 
  ## net_attr
  dat$run$net_attr[[dat$num.nw]] <- list() 
  dat$run$net_attr[[dat$num.nw]][["n"]] <- dat$run$num
  ## control
  dat$control[["tergmLite.track.duration"]][[dat$num.nw]] <- FALSE
  
  # simulate first time step
  dat$num.nw <- 3
  dat <- sim_nets_t1(dat)
  dat <- summary_nets(dat, at = 1L)
  dat$num.nw <- 4
  
  num <- sum(get_attr(dat, "active") == 1)
  
  # Time Unit
  time.unit <- get_param(dat, "time.unit")
  # time.unit <- param$epistats$time.unit
  dat <- set_param(dat, "time.unit",  time.unit #time.unit
                   )
  
  dat[["temp"]] <- list()
  
  # Prevalence Tracking
  dat <- set_epi(dat, "num", at = 1,  num)
  
  return(dat)

 }


init_attrs <- function(dat) { 
  
  n_nodes <- sum(get_attr(dat, "active") ==  1)
  def_attrs <- get_default_attrs(dat)
  cur_attrs <- get_attr_list(dat)
  missing <- names(def_attrs)[!names(def_attrs) %in% names(cur_attrs)]
  for (attr_name in missing) {
    dat <- append_attr(dat, attr_name, def_attrs[[attr_name]], n_nodes)
  }
  
  dat <- make_computed_attrs(dat, n_nodes, post_init = FALSE) 
  
  
  return(dat)
}


