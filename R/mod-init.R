
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
  
  # simulate first time step
  
  dat <- sim_nets_t1(dat)
  dat <- summary_nets(dat, at = 1L)
  
  num <- sum(get_attr(dat, "active") == 1)
  
  # Add household network edgelist
  ## network index
  dat$num.nw <- dat$num.nw + 1
  ## edgelist 
  dat$run$el[[dat$num.nw]] <- NULL
  dat$run$el[[dat$num.nw]] <- as.matrix(dat$param$hh.pairs) 
  ## net_attr
  dat$run$net_attr[[dat$num.nw]] <- list() 
  dat$run$net_attr[[dat$num.nw]][["n"]] <- dat$run$num
  ## control
  dat$control[["tergmLite.track.duration"]][[dat$num.nw]] <- FALSE
  
  
  
  
  # Time Unit
  # time.unit <- param$epistats$time.unit
  dat <- set_param(dat, "time.unit", 1 #time.unit
                   )
  
  dat[["temp"]] <- list()
  
  # Prevalence Tracking
  dat <- set_epi(dat, "num", at = 1,  num)
  
  return(dat)


#   # define initial nodal attributes, added as we need this for the x-layer sim
#   node_age_grp <- as.character(init$attr$node.age.grp)
#   node_age     <- as.numeric  (init$attr$node.age)
#   deg_work   <- as.integer  (init$attr$contact_attribute_Work)
#   deg_school <- as.integer  (init$attr$contact_attribute_School)
#   deg_nonhome   <- as.integer  (init$attr$contact_attribute_Nonhome)
# 
#   no_contact <- 1L - deg_nonhome
#   
#   
#   ## set attributes for age and contact
#   dat <- set_attr(dat, "age.grp", node_age_grp)
#   dat <- set_attr(dat, "age",     node_age)
#   dat <- set_attr(dat, "deg_school", deg_school)
#   dat <- set_attr(dat, "deg_work",   deg_work)
#   dat <- set_attr(dat, "no.contact",         no_contact)
# 
#   # simulate first time step, time 1 to 2
#   dat <- sim_nets_t1(dat)
#   dat <- summary_nets(dat, at = 1L)
# 
#   # ## add home layer edgelist
#    dat$num.nw <- dat$num.nw + 1L
#    dat$el[[dat$num.nw]] <- as.matrix(dat$param$hh.pairs) # design-time containers
#    attr(dat$el[[dat$num.nw]], 'n') <-  dat$run$num # muted code doesn't work- get_epi(dat, "sim.num", at = 1)
#    dat$net_attr[[dat$num.nw]] <- list(n = dat$run$num) # muted code doesn't work - list(n = get_epi(dat, "sim.num", at = 1))
#    dat$control[["tergmLite.track.duration"]][[dat$num.nw]] <- FALSE
#    
#   
#   # add hh_id attribute to each node, load from run from data
#   #dat$attr$hh.ids <-  dat$attr <- list()
#   #dat$attr$hh.ids <-  dat$run$attr$hh.ids
#   
#   # code for tergmLite. This is needed because home layer is added after the initiation is completed
#   # this is necessary to make the  depart_nodes work 
#   dat$run$el[[dat$num.nw]]        <- as.matrix(dat$param$hh.pairs)  
#   dat$run$net_attr[[dat$num.nw]]  <- dat$net_attr[[dat$num.nw]]
# 
# 
#   ## Infection Status and Time Modules
#   dat <- init_status_covid_corporate(dat)
# 
#   ## Get initial prevalence
#   dat <- prevalence_covid_corporate(dat, at = 1)
# 
#   return(dat)
 }


init_attrs <- function(dat) { #no.contact needs to be defined here.
  

   
  n_nodes <- sum(get_attr(dat, "active") ==  1)
  def_attrs <- get_default_attrs(dat)
  cur_attrs <- get_attr_list(dat)
  missing <- names(def_attrs)[!names(def_attrs) %in% names(cur_attrs)]
  for (attr_name in missing) {
    dat <- append_attr(dat, attr_name, def_attrs[[attr_name]], n_nodes)
  }
  
  dat <- make_computed_attrs(dat, n_nodes, post_init = FALSE) 
  
  # # HIV-related attributes
  # dat <- init_status_msm(dat)
  # 
  # # STI Status
  # dat <- gono_chla_init(dat, "gono")
  # dat <- gono_chla_init(dat, "chla")
  # dat <- syph_init(dat)
  
  return(dat)
}


# init_status_gmc19 <- function(dat) {
#   
#   e.num <- get_init(dat, "e.num")
#   
#   active <- get_attr(dat, "active")
#   num <- sum(active)
#   
#   ## Disease status
#   status <- rep("s", num)
#   if (e.num > 0) {
#     status[sample(which(active == 1), size = e.num)] <- "e"
#   }
#   
#   dat <- set_attr(dat, "status", status)
#   
#   # Age group for vaccination processes
#   age <- get_attr(dat, "age")
#   
#   vax.age.group <- rep(NA, length(age))
#   vax.age.group[age < 5] <- 1
#   vax.age.group[age >= 5 & age < 18] <- 2
#   vax.age.group[age >= 18 & age < 50] <- 3
#   vax.age.group[age >= 50 & age < 65] <- 4
#   vax.age.group[age >= 65] <- 5
#   
#   dat <- set_attr(dat, "vax.age.group", vax.age.group)
#   
#   # Infection Time and related attributes
#   idsInf <- which(status == "e")
#   infTime <- rep(NA, num)
#   clinical <- rep(NA, num)
#   hospit <- rep(NA, num)
#   statusTime <- rep(NA, num)
#   statusTime[idsInf] <- 1
#   dxStatus <- rep(0, num)
#   dxTime <- rep(NA, num)
#   vax <- rep(0, num)
#   vax1Time <- rep(NA, num)
#   vax2Time <- rep(NA, num)
#   vax3Time <- rep(NA, num)
#   vax4Time <- rep(NA, num)
#   isolate <- rep(NA, num)
#   isoTime <- rep(NA, num)
#   
#   dat <- set_attr(dat, "statusTime", statusTime)
#   dat <- set_attr(dat, "infTime", infTime)
#   dat <- set_attr(dat, "clinical", clinical)
#   dat <- set_attr(dat, "hospit", hospit)
#   dat <- set_attr(dat, "dxStatus", dxStatus)
#   dat <- set_attr(dat, "dxTime", dxTime)
#   dat <- set_attr(dat, "vax", vax)
#   dat <- set_attr(dat, "vax1Time", vax1Time)
#   dat <- set_attr(dat, "vax2Time", vax2Time)
#   dat <- set_attr(dat, "vax3Time", vax3Time)
#   dat <- set_attr(dat, "vax4Time", vax4Time)
#   dat <- set_attr(dat, "isolate", isolate)
#   dat <- set_attr(dat, "isoTime", isoTime)
#   
#   return(dat)
# }
  
