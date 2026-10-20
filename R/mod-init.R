
#' @rdname moduleset-gmc19
#' @export
init_gmc19 <- function(x, param, init, control, s) {
  
  ## Master Data List Setup ##
  dat <- create_dat_object(param, init, control)
  # GlobalMix India #67: a daily time step, and inf.prob from the pathogen scenario
  check_required_params(dat)
  # #64: births at the crude birth rate, and the background mortality vector
  check_demography(dat)

  # The household layer used to be an edgelist passed as hh.pairs and appended
  # here; it is now a netclique layer in `x` (GlobalMix India #85). Stop rather
  # than run without households if a caller still passes the old form.
  if (!is.null(get_param(dat, "hh.pairs", override.null.error = TRUE))) {
    stop("hh.pairs is no longer read: pass the household layer to netsim() as ",
         "netclique(nw, \"hh.ids\", arrivals = \"join\", arrivals.FUN = hh_join_birth)")
  }

  ## network and stats initialization
  dat <- init_nets(dat, x)   

  # GlobalMix India #66: no retired CorporateMix parameter in the table, and
  # the per-layer parameters line up with the layers passed to netsim(). #67:
  # one warning naming every parameter that nothing on this design reads.
  check_retired_params(dat)
  check_known_params(dat)
  dat <- check_layer_params(dat)
  # #69, #70, #71: severity per infection by severity band, post-infection
  # immunity, imports and prior immunity (mod-immunity.R)
  check_severity_params(dat)
  check_immunity_params(dat)
  check_import_params(dat)
  check_prior_immunity(dat)
  check_vax_schedule(dat)

  # GlobalMix India #61: sex and occupation class come with the population, as
  # vertex attributes of the networks, like age and hh.ids
  dat <- init_pop_class(dat)

  # GlobalMix India #63: the relationship traits come from the fitted networks'
  # vertex attributes. Without them the default attributes would be appended
  # below and the formation offsets would run on zeros.
  if (activation_on(dat)) {
    traits <- c("never", "log_aF", "log_aI")
    missing <- setdiff(traits, names(get_attr_list(dat)))
    if (length(missing) > 0) {
      stop("the relationship networks carry no vertex attribute ", toString(missing),
           "; fit them with the traits from draw_relationship_traits()", call. = FALSE)
    }
    if (any(vapply(traits, function(a) anyNA(get_attr(dat, a)), logical(1)))) {
      stop("a relationship trait (never, log_aF, log_aI) is NA for some node", call. = FALSE)
    }
    check_newborn_classes(dat)
    # #65: births placed with a mother, and a household layer that keeps it
    check_birth_placement(dat)
    # #64: the class distribution by age, checked once here rather than at the
    # first transition
    normalize_occ_age(get_param(dat, "act.occ.age", override.null.error = TRUE))
  }

  
  ## Initialize all remaining attributes in default values
  dat <- init_attrs(dat)
  dat <- overwrite_attrs(dat)
  # #69, #70: the seeds are the run's first infections, and init.immune makes
  # a share of each severity band immune for the run
  dat <- init_immunity(dat)
  # #65: the mothers of the youngest children, who cannot give birth again
  # within the birth interval
  if (births_placed(dat)) dat <- init_last_birth(dat)
  
  # # Infection status and related
  # dat <- init_status_gmc19(dat)
  
  # simulate first time step. The household layer is a netclique layer in the
  # list passed to netsim(), which sim_nets_t1() and summary_nets() handle
  # themselves: it has no network model to simulate.
  dat <- sim_nets_t1(dat)
  # The relationship networks start from an exact draw of their cross-section
  # (GlobalMix India #63): ergm's MCMC start under-placed the ties of the most
  # active people (2.7% too few among the top 1% by activity), and the exact
  # daily step would carry that error for years
  if (activation_on(dat)) dat <- start_relationships(dat)
  dat <- summary_nets(dat, at = 1L)
  
  num <- sum(get_attr(dat, "active") == 1)

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

# Per-layer parameters (GlobalMix India #66). A vector such as inf.prob is
# read by position, one value per layer, and positions carry no meaning of
# their own: these values came from the CorporateMix settings, where layer 1
# was a different layer, and the school layer inherited 5 acts per day. So the
# parameter table names its layers (layer.names), the run stops unless those
# are the names of the layer list passed to netsim(), in that order, and every
# per-layer vector is then named by layer.
#
# With the activation step (GlobalMix India #63) the per-layer parameters are
# read per tie type: the two relationship networks, the new contacts the engine
# draws each day (tie type "new", which is not a network), and the household
# layer, in the order of layer.names. The networks passed to netsim() are then
# layer.names without "new", in the same order.
#
# sympt.contact.rr (#67, #68) is the per-contact hazard multiplier while the
# infector is symptomatic (status "ic"), one value per tie type; absent, 1.
layer_params <- c("inf.prob", "exposure.rate", "sympt.contact.rr")

layer_names <- function(dat) {
  names(if (!is.null(dat$run$el)) dat$run$el else dat$run$nw)
}

tie_types <- function(dat) {
  if (!is.null(dat$run$tie.types)) dat$run$tie.types else layer_names(dat)
}

check_layer_params <- function(dat) {
  layers <- layer_names(dat)
  if (is.null(layers) || any(layers == "") || anyDuplicated(layers) > 0) {
    stop("name every layer in the list passed to netsim(), e.g. ",
         "list(freq = ..., infreq = ..., home = ...), so the ",
         "per-layer parameters can be checked against it", call. = FALSE)
  }
  declared <- get_param(dat, "layer.names", override.null.error = TRUE)
  if (is.null(declared)) {
    stop("the parameters have no layer.names; set it to the names of the ",
         "layers passed to netsim(), in order (", toString(layers), ")",
         call. = FALSE)
  }
  declared <- as.character(declared)
  if (activation_on(dat)) {
    if (!all(rel_layers %in% layers)) {
      stop("with the activation step the relationship networks must be named ",
           toString(rel_layers), "; the layers are ", toString(layers), call. = FALSE)
    }
    if (!new_type %in% declared || !identical(setdiff(declared, new_type), layers)) {
      stop("layer.names in the parameters (", toString(declared), ") must be the ",
           "layers passed to netsim() (", toString(layers), ") with the new-contact ",
           "tie type \"", new_type, "\" added", call. = FALSE)
    }
    if (is.null(dat$run$el)) {
      stop("the activation step needs tergmLite = TRUE", call. = FALSE)
    }
    activation_params(dat)  # stops here, not mid-run, on a malformed parameter
    check_trait_params(dat)
    dat$run$tie.types <- declared
  } else if (any(rel_layers %in% layers) || new_type %in% declared) {
    stop("the relationship networks and the new-contact tie type need the activation ",
         "parameters (act.o, act.kappa, ...); without them every relationship would be ",
         "a contact every day", call. = FALSE)
  } else if (!identical(declared, layers)) {
    stop("layer.names in the parameters (", toString(declared), ") are not ",
         "the layers passed to netsim() (", toString(layers), ")", call. = FALSE)
  }
  for (name in intersect(layer_params, names(dat$param))) {
    dat <- set_param(dat, name, layer_param(dat, name))
  }
  # Checked here so a malformed value stops the run before the first step
  if (!is.null(get_param(dat, "contact.weight", override.null.error = TRUE))) {
    w <- contact_weights(dat)
    if (!is.null(w) && !all(vapply(w, function(m) isSymmetric(unname(m)), logical(1)))) {
      stop("every `contact.weight` matrix must be symmetric in the age pair", call. = FALSE)
    }
  }
  score <- as.character(param_or(dat, "vax.degree.score", "total"))
  if (!score %in% degree_scores) stop("unknown vax.degree.score \"", score, "\"", call. = FALSE)
  if (!score %in% c("total", "xlayer") && !activation_on(dat)) {
    stop("vax.degree.score \"", score, "\" needs the relationship model", call. = FALSE)
  }
  dat
}

# A per-layer parameter with one value per layer (tie type), named by layer. An
# unnamed vector is read in layer order (check_layer_params() matched
# layer.names to the layers at initialization); a named one must name the
# layers in order. Called on every step, so a parameter updater cannot
# misalign a vector either.
layer_param <- function(dat, name) {
  value <- get_param(dat, name)
  layers <- tie_types(dat)
  if (length(value) != length(layers)) {
    stop("`", name, "` has ", length(value), " values for ", length(layers),
         " layers (", toString(layers), ")", call. = FALSE)
  }
  if (is.null(names(value))) {
    names(value) <- layers
  } else if (!identical(names(value), layers)) {
    stop("`", name, "` is named ", toString(names(value)), " but the layers ",
         "are ", toString(layers), call. = FALSE)
  }
  value
}

# Parameters of the CorporateMix modifiers removed from the GMC19 modules
# (GlobalMix India #66): the inf.prob windows on days 40-80 and 120-160, the
# 40% transmission cut from day 15, the hospitalization multiplier from day 211,
# and the testing and isolation pathway (the dx and isolate modules, whose
# attributes the model no longer carries). Nothing reads them now, so a table
# or updater that still sets one would run as if it had an effect.
retired_params <- c(
  "inf.add", "inf.sub", "inf.boost.start", "inf.boost.stop",
  "inf.supp.start", "inf.supp.stop",
  "inf.prob.inter.rr", "inf.prob.inter.time",
  "act.rate.inter.rr", "act.rate.inter.time",
  "hosp.boost.mult", "hosp.boost.start", "hosp.boost.stop",
  "hosp.supp.mult", "hosp.supp.start", "hosp.supp.stop",
  "dx.rate.sympt", "dx.rate.other", "pcr.sens", "allow.rescreen",
  "iso.prob.office", "iso.prob.other", "notif.prob",
  "act.rate.iso.inter.rr", "act.rate.iso.inter.time", "inf.prob.mask.rr",
  "act.rate.dx.inter.rr", "act.rate.dx.inter.time",
  "act.rate.sympt.inter.rr", "act.rate.sympt.inter.time"
)

# Parameters the redo replaced with others (GlobalMix India #69, #71), each with
# what replaced it. A table that still sets one would run as if it had an effect.
replaced_params <- c(
  rs.rate = "immune.full and immune.peak/halflife/floor.* (#69): a full protection period, then partial protection",
  prop.hospit = "ihr, admissions per infection by severity band (#71)",
  prop.hospit.infant = "ihr on the infant severity bands (#71)",
  prop.hospit.infant6 = "ihr on the infant severity bands (#71)",
  dis.death.rate = "ifr, deaths per infection by severity band, drawn once per episode (#71)",
  dis.death.rate.infant = "ifr on the infant severity bands (#71)",
  dis.death.rate.infant6 = "ifr on the infant severity bands (#71)"
)

check_retired_params <- function(dat) {
  updated <- unlist(lapply(dat$param$.param.updater.list,
                           function(u) names(u$param)))
  set <- c(names(dat$param), updated)
  found <- intersect(retired_params, set)
  if (length(found) > 0) {
    stop("the model no longer reads ", toString(found), ": the CorporateMix ",
         "modifiers they set were removed (GlobalMix India #66). Remove them ",
         "from the parameters.", call. = FALSE)
  }
  found <- intersect(names(replaced_params), set)
  if (length(found) > 0) {
    stop("the model no longer reads ",
         paste0(found, " (replaced by ", replaced_params[found], ")", collapse = "; "),
         ". Remove them from the parameters (GlobalMix India).", call. = FALSE)
  }
  invisible(TRUE)
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
