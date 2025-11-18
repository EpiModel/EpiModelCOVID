
#' @rdname moduleset-gmc19
#' @export
resim_nets_gmc19_x_layer <- function(dat, at) { # adapted from cruiseship function 
  
  nms  <- names(dat$run$el)
  idx_school <- which(nms == "school")
  idx_work <- which(nms == "work")
  idx_nonhome <- which(nms == "nonhome")
  
  dat$num.nw <- dat$num.nw - 1 # disregard household nw for this module
  
  ## Edges correction
  dat <- edges_correct(dat, at)
  
  ## network resimulation
  dat.updates <- NVL(get_control(dat, "dat.updates"), function(dat, ...) dat)
  
  # nonhome simulation
  dat <- simulate_dat(dat = dat, at = at, network = idx_nonhome)
  dat <- dat.updates(dat = dat, at = at, network = idx_nonhome)
  
  # school simulation at t+1 (on/after simulate_date), depends on work's degree at t (before simulate_dat)
  work_layer_t <-   get_network(x =dat, network = idx_work)
  deg_work_t <-   get_degree(work_layer_t)
  
  deg_bi_work_t <- #  binarize degree, still numeric
    ifelse(deg_work_t>0, yes=1, no=0)
  
  school_layer_t <-   get_network(x =dat, network = idx_school) # get the school layer's degree
  
  school_layer_t <- # assign binarized degree at work (at t) to school layer
  set_vertex_attribute(school_layer_t, attrname = "deg.x_layer" , 
                       value =deg_bi_work_t # this is contact status at work layer 
  )
  dat <- set_network(x = dat, network = idx_school, nw = school_layer_t)
  
  dat  <- simulate_dat(dat = dat, at = at, network = idx_school)
  dat  <- dat.updates(dat = dat, at = at, network = idx_school)
  
  # work simulation at t, depends on school's degree at t+1 ----
  school_layer_t <-   get_network(x =dat, network = idx_school)
  deg_school_t <-   get_degree(school_layer_t)
  
  deg_bi_school_t <- #  binarize degree, still numeric
    ifelse(deg_school_t>0, yes=1, no=0)
  
  work_layer_t <-   get_network(x =dat, network = idx_work) # get the work layer
  
  work_layer_t <- # assign binarized degree at school (at t) to work layer
    set_vertex_attribute(work_layer_t, attrname = "deg.x_layer" , 
                         value =deg_bi_school_t # this is contact status at school layer 
    )
  dat <- set_network(x = dat, network = idx_work, nw = work_layer_t)
  
  dat  <- simulate_dat(dat = dat, at = at, network = idx_work)
  dat  <- dat.updates(dat = dat, at = at, network = idx_work)

  dat$num.nw <- dat$num.nw + 1 # disregard household nw for this module
  
  return(dat)
  
}