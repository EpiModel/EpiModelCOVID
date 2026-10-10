
#' @rdname moduleset-gmc19
#' @export
init_gmc19 <- function(x, param, init, control, s) {
  
  ## Master Data List Setup ##
  dat <- create_dat_object(param, init, control)

  # The household layer used to be an edgelist passed as hh.pairs and appended
  # here; it is now a netclique layer in `x` (GlobalMix India #85). Stop rather
  # than run without households if a caller still passes the old form.
  if (!is.null(get_param(dat, "hh.pairs", override.null.error = TRUE))) {
    stop("hh.pairs is no longer read: pass the household layer to netsim() as ",
         "netclique(nw, \"hh.ids\", arrivals = \"join\", arrivals.FUN = hh_join_young)")
  }

  ## network and stats initialization
  dat <- init_nets(dat, x)   

  
  ## Initialize all remaining attributes in default values
  dat <- init_attrs(dat)
  dat <- overwrite_attrs(dat)
  
  # # Infection status and related
  # dat <- init_status_gmc19(dat)
  
  # simulate first time step. The household layer is a netclique layer in the
  # list passed to netsim(), which sim_nets_t1() and summary_nets() handle
  # themselves: it has no network model to simulate.
  dat <- sim_nets_t1(dat)
  dat <- summary_nets(dat, at = 1L)
  
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


# Index of the household layer: the netclique layer grouped by hh.ids. Found by
# its grouping attribute rather than its position, so the layer list passed to
# netsim() may be in any order. NA when the model has no household layer.
hh_layer <- function(dat) {
  hh <- which(vapply(seq_len(dat$num.nw), function(network) {
    identical(get_nwparam(dat, network = network)$group.attr, "hh.ids")
  }, logical(1)))
  if (length(hh) == 0) NA_integer_ else hh[1]
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

# init_status_gmc19 <- # adapted from init_status_covid_corporate
#   function(dat) {
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
#   # # Age group for vaccination processes
#   # age <- get_attr(dat, "age")
#   # 
#   # vax.age.group <- rep(NA, length(age))
#   # vax.age.group[age < 5] <- 1
#   # vax.age.group[age >= 5 & age < 18] <- 2
#   # vax.age.group[age >= 18 & age < 50] <- 3
#   # vax.age.group[age >= 50 & age < 65] <- 4
#   # vax.age.group[age >= 65] <- 5
#   # 
#   # dat <- set_attr(dat, "vax.age.group", vax.age.group)
#   # 
#   # # Infection Time and related attributes
#    idsInf <- which(status == "e")
#   # infTime <- rep(NA, num)
#   # clinical <- rep(NA, num)
#   # hospit <- rep(NA, num)
#   # statusTime <- rep(NA, num)
#    statusTime <- get_attr(dat, "statusTime")
#    statusTime[idsInf] <- 1
#    # dxStatus <- rep(0, num)
#   # dxTime <- rep(NA, num)
#    # vax <- rep(0, num)
#   # vax1Time <- rep(NA, num)
#   # vax2Time <- rep(NA, num)
#   # vax3Time <- rep(NA, num)
#   # vax4Time <- rep(NA, num)
#   # isolate <- rep(NA, num)
#   # isoTime <- rep(NA, num)
#   # 
#   dat <- set_attr(dat, "statusTime", statusTime)
#   # dat <- set_attr(dat, "infTime", infTime)
#   # dat <- set_attr(dat, "clinical", clinical)
#   # dat <- set_attr(dat, "hospit", hospit)
#    # dat <- set_attr(dat, "dxStatus", dxStatus)
#   # dat <- set_attr(dat, "dxTime", dxTime)
#    # dat <- set_attr(dat, "vax", vax)
#   # dat <- set_attr(dat, "vax1Time", vax1Time)
#   # dat <- set_attr(dat, "vax2Time", vax2Time)
#   # dat <- set_attr(dat, "vax3Time", vax3Time)
#   # dat <- set_attr(dat, "vax4Time", vax4Time)
#   # dat <- set_attr(dat, "isolate", isolate)
#   # dat <- set_attr(dat, "isoTime", isoTime)
#   
#   return(dat)
# }
# 
# 
