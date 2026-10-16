#' Get the list the attributes used by the model with their default value
get_default_attrs <- function(dat) { # all attributes should be listed here
  defaults <- list(
    # network attributes, binary
    deg_work=0,
    deg_school=0,
    deg_nonhome=0,
    # non-home structural-isolate flag (issue #43): drives the fixed
    # offset(nodefactor("no.contact")) in the non-home netest, so it must be
    # maintained for new arrivals or resim_nets desyncs. New arrivals are infants,
    # which have no non-home contact, hence no.contact = 1 (consistent with
    # deg_nonhome = 0 above; no.contact == 1 - deg_nonhome).
    no.contact=1,
    # network attributes, integer
    degree_work=0,
    degree_school=0,
    degree_nonhome=0,
    degree_hh=0,
    degree_total=0,
    degree_xlayer=0,
    degree_quartile=NA,
    # relationship traits and activation (GlobalMix India #63). never, log_aF and
    # log_aI come from the fitted networks at initialization and are drawn for
    # newborns in make_computed_attrs(); o_i is drawn there for everyone.
    never = 0,
    log_aF = 0,
    log_aI = 0,
    log_aN = 0,
    o_i = 0,
    act.eligible = 1,
    degree_freq = 0,
    degree_infreq = 0,
    degree_active = 0,
    degree_hazard = 0,
    degree_score = 0,
    act.out.days = 0,
    act.contact.days = 0,
    act.new.contacts = 0,
    act.recent = 0,
    act.recent.new = 0,
    # step on which the node last crossed a band top (GlobalMix India #64)
    rel.cross = NA,
    # cross-layer bridge indicators
    is_bridge = FALSE,
    n_layers_active =0,
    # household id and age. NA is "in no household"; a newborn takes its
    # mother's (place_births(), GlobalMix India #65) or, on the legacy model,
    # the household clique layer's arrival rule (hh_join_birth) draws one.
    hh.ids = NA,
    age = get_param(dat, "arrival.age"), 
    age.grp = NA,
    # disease-related attributes
    vax.age.group = NA, 
    statusTime = NA,
    infTime = NA,
    clinical= NA,
    hospit= NA,
    status="s",
    vax= 0,
    vax1Time= NA,
    vax2Time= NA,
    vax3Time= NA,
    vax4Time= NA,
    last.dose.time = NA,
    # Person-level cumulative flags (issue #46). Everything else the model emits
    # is an episode flow, so with reinfection the per-100 infection numbers
    # exceed the number of infected PEOPLE and cannot be compared to any
    # seroprevalence or cohort attack rate. These two make the person the unit:
    # ever.inf is set once on a node's first S->E transition (and for the seeds
    # at initialization) and never cleared, ever.sympt once on its first entry
    # to "ic". New arrivals start at 0: they are born into the simulation
    # susceptible.
    ever.inf = 0,
    ever.sympt = 0,
    # Episodes and immunity (GlobalMix India #69, #70, #71): whether the
    # current episode is fatal (drawn at onset with hospit), the step of the
    # last recovery (the clock of post-infection protection), the number of
    # infections in the run, and immunity from before the run. Newborns have
    # none of them.
    fatal = NA,
    recTime = NA,
    inf.count = 0,
    prior.immune = 0

  )
  # The propensity attribute prior immunity is drawn with (init.immune.attr,
  # #70), carried by the population; newborns are born after the draw
  w_attr <- immune_weight_attr(dat)
  if (!is.null(w_attr)) defaults[[w_attr]] <- NA_real_
  # Sex and occupation class (GlobalMix India #61), when the population carries
  # them: newborns are children, and make_computed_attrs() draws their sex
  if (has_pop_class(dat)) {
    defaults$sex <- NA_character_
    defaults$occ <- "child"
  }
  # Step of a woman's most recent birth, NA if none is on record (GlobalMix
  # India #65): place_births() sets it, and init_last_birth() for the mothers of
  # the starting population's youngest children
  if (births_placed(dat)) {
    defaults$last.birth <- NA_real_
  }
  defaults
}

# Sex and occupation class (GlobalMix India #61) -------------------------------
#
# The population carries each person's sex and occupation class as vertex
# attributes of the networks, which init_nets() copies to the attribute list
# with age and hh.ids. Newborns are children, male with probability
# birth.prop.male. A population without them (the four-layer model's) carries
# neither, its newborns get neither, and no random number is drawn for them.
pop_class_attrs <- c("sex", "occ")
sex_levels <- c("F", "M")
occ_levels <- c("child", "student", "worker", "home")

has_pop_class <- function(dat) {
  all(pop_class_attrs %in% names(get_attr_list(dat)))
}

# Checks the population's sex and occ at initialization and stores them as
# character vectors
init_pop_class <- function(dat) {
  present <- intersect(pop_class_attrs, names(get_attr_list(dat)))
  if (length(present) == 0) {
    if (!is.null(get_param(dat, "act.class.mult", override.null.error = TRUE))) {
      stop("act.class.mult needs each person's sex and occ, as vertex attributes of the ",
           "relationship networks", call. = FALSE)
    }
    return(dat)
  }
  if (length(present) == 1) {
    stop("the population carries ", present, " but not ", setdiff(pop_class_attrs, present),
         "; set both as vertex attributes of the networks", call. = FALSE)
  }
  sex <- as.character(get_attr(dat, "sex"))
  occ <- as.character(get_attr(dat, "occ"))
  if (!all(sex %in% sex_levels)) {
    stop("sex must be F or M; the population has ", toString(setdiff(unique(sex), sex_levels)),
         call. = FALSE)
  }
  if (!all(occ %in% occ_levels)) {
    stop("occ must be one of ", toString(occ_levels), "; the population has ",
         toString(setdiff(unique(occ), occ_levels)), call. = FALSE)
  }
  p <- get_param(dat, "birth.prop.male", override.null.error = TRUE)
  if (!isTRUE(length(p) == 1 && p >= 0 && p <= 1)) {
    stop("birth.prop.male, the share of newborns who are male, must be one value in [0, 1] ",
         "when the population carries sex", call. = FALSE)
  }
  dat <- set_attr(dat, "sex", sex)
  set_attr(dat, "occ", occ)
}

#' Update the attributes requiring computation for new nodes
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

  # Infant flag: an age-resolution refinement (issue #23) that does NOT re-bin the
  # network age groups. Infants stay in the youngest contact band (their contacts
  # are household-dominated) but are flagged so the disease modules can apply
  # infant-specific severity (RSV) and the outputs can report the infant burden
  # separately rather than diluting it across the coarse youngest band. Absent
  # param -> infant.age.max = 1 year. Recomputed for arrivals so it stays correct.
  infant.age.max <- get_param(dat, "infant.age.max", override.null.error = TRUE)
  if (is.null(infant.age.max)) infant.age.max <- 1

  n_attr <- c(n_attr, list(
    age.grp     =   cut(age,
                        age.breaks,
                        labels = age.grps,
                        right = FALSE
                        ) |> as.character(),
    is_infant   =   age < infant.age.max
  ))
  
  # Sex of newborns (GlobalMix India #61), drawn only when the population
  # carries sex and occupation class
  pop_class <- has_pop_class(dat)
  if (post_init && pop_class) {
    p.male <- get_param(dat, "birth.prop.male")
    n_attr$sex <- ifelse(stats::runif(n_new) < p.male, "M", "F")
  }

  # Relationship traits and the propensity to go out (GlobalMix India #63).
  # Newborns draw their traits from the youngest band's parameters, as the
  # initial 0-9y population's were drawn, so the two are exchangeable (#65),
  # and with their own sex and occupation class (#61). o_i is drawn here for
  # everyone, so act.o and act.kappa can be set by a scenario.
  if (activation_on(dat)) {
    pars <- activation_params(dat)
    band_new <- match(n_attr$age.grp, age.grps)
    if (anyNA(band_new)) {
      stop("age.grp outside age.grps for ", sum(is.na(band_new)), " nodes", call. = FALSE)
    }
    if (post_init) {
      tr <- draw_relationship_traits(band_new, pars, sex = if (pop_class) n_attr$sex,
                                     occ = if (pop_class) get_attr(dat, "occ", posit_ids = new_nodes_pid))
      n_attr$never <- tr$never
      n_attr$log_aF <- tr$log_aF
      n_attr$log_aI <- tr$log_aI
    }
    never_new <- if (post_init) n_attr$never else get_attr(dat, "never", posit_ids = new_nodes_pid)
    n_attr$o_i <- draw_out_propensity(band_new, pars$o, pars$kappa)
    n_attr$act.eligible <- 1 - never_new
    # A persistent propensity for new contacts, mean 1: zero for a share act.pi.N
    # of people, and gamma with shape act.theta.N among the rest (new contacts
    # are concentrated in few people in the diary)
    aN <- if (is.finite(pars$theta_N)) stats::rgamma(n_new, shape = pars$theta_N, rate = pars$theta_N) else rep(1, n_new)
    aN <- aN * (stats::runif(n_new) >= pars$pi_N) / (1 - pars$pi_N)
    n_attr$log_aN <- log(pmax(aN, 1e-8))
  }

  for (attr_name in names(n_attr)) {
    dat <- set_attr(dat, attr_name, n_attr[[attr_name]], posit_ids = new_nodes_pid)
  }
  
  return(dat)
}

