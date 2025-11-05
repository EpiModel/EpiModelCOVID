
#' @rdname moduleset-common
#' @export
aging_covid <- function(dat, at) {

  age <- get_attr(dat, "age")
  age <- age + 1 / 365
  dat <- set_attr(dat, "age", age)

  return(dat)
}


#' @rdname moduleset-gmc19
#' @export
deaths_covid_gmc19 <- function(dat, at) {
  
  ## Attributes ##
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")
  
  # these doesn't work 
  # active <- dat$attr$active
  # age <- dat$attr$age
  # status <- dat$attr$status
  
  ## Parameters ##
  mort.rates <- dat$param$mort.rates
  mort.dis.mult <- dat$param$mort.dis.mult
  
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <-
    nDeathsH <- 0
  
  if (nElig > 0) {
    
    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 86) # take upper bound of numeric age as integer age
    death_rates_of_elig <- mort.rates[whole_ages_of_elig] # age-dependent death rate
    
    idsElig.inf <- which(status[idsElig] == "h") # if hospitalized, scale up the mortality rate
    death_rates_of_elig[idsElig.inf] <- death_rates_of_elig[idsElig.inf] *mort.dis.mult
    
    vecDeaths <- which(rbinom(nElig, 1, death_rates_of_elig) == 1) # among the number of active nodes, use death rate of each individual to determine who dies
    idsDeaths <- idsElig[vecDeaths] # id of individual who dies
    nDeaths <- length(idsDeaths) # number of people who dies
    nDeathsH <- length(intersect(idsDeaths, idsElig.inf)) # number of people who dies and is hospitalized
    
    if (nDeaths > 0) { # record the ids of those who dies at this time step
      # dat$attr$active[idsDeaths] <- 0, this doesn't work
      #inactive <- which(dat$attr$active == 0) this doesn't work
      dat <- depart_nodes(dat, departures = idsDeaths) # only the things in "run" are updated
      

      
    }
  }
  
  ## Summary statistics ##
  dat <- set_epi(dat, "d.flow", at, nDeaths)
  dat <- set_epi(dat, "d.h.flow", at, nDeathsH)
  
  return(dat)
}



#' @rdname moduleset-gmc19
#' @export
arrival_covid_gmc19 <- function(dat, at) {
  # if (at>200) browser()
  
  # Parameters
  a.rate   <- get_param(dat, "a.rate")
  
  ## Process
  num <- dat$epi$num[1]
  nNew <- rpois(1, a.rate * num) # number of people entering the model, based on the arrival rate
  
  ## current node id, self added
  #old_node_ids <- dat$run$attr$unique_id
  
  ## Update Attr
  if (nNew > 0) {
    dat <- setNewAttr_covid_gmc19(dat, at, nNew) 
  }
  
  # self added
  # new_node_ids <- setdiff( dat$run$attr$unique_id, old_node_ids )
  # 
  # dat$temp$new_node_ids <- list()
  # 
  # dat$temp$new_node_ids <- new_node_ids
  
  # Update Networks
  dat <- arrive_nodes(dat, nNew)
  
  ## Output
  dat <- set_epi(dat, "nNew", at, nNew)
  
  return(dat)
}


setNewAttr_covid_gmc19 <- function(dat, at, nNew) {
  
  dat <- append_core_attr(dat, at, nNew) # add nNew nodes to dat
  
  # Age related
  arrival.age <- get_param(dat, "arrival.age") # those arrive were 0 year old
  newAges <- rep(arrival.age, nNew) # make a vector of arrival age
  dat <- append_attr(dat, "age", newAges, nNew) # add to dat
  
  # age.breaks <- seq(0, 200, 10)
  # attr_age.grp <- cut(newAges, age.breaks, labels = FALSE, right = FALSE)
  AGE_LABELS <- c("0-9y","10-19y","20-29y","30-39y","40-59y","60+y")
  age.breaks  <- c(0, 10, 20, 30, 40, 60, Inf)
  attr_age.grp <- cut(newAges,
                        breaks =  age.breaks ,
                        labels = AGE_LABELS,
                        right  = FALSE) %>% as.character() # different from corportate which use numeric grp
  dat <- append_attr(dat, "age.grp",attr_age.grp, nNew) # add categorical ages to new nodes
  
  ## age at vax
  vax.age.group <- rep(NA, length(newAges))
  vax.age.group[newAges < 10] <- 1
  vax.age.group[newAges >= 10 & newAges < 20] <- 2
  vax.age.group[newAges >= 20 & newAges < 30] <- 3
  vax.age.group[newAges >= 30 & newAges < 40] <- 4
  vax.age.group[newAges >= 40 & newAges < 60] <- 5
  vax.age.group[newAges >= 60 ] <- 6
  dat <- append_attr(dat, "vax.age.group", vax.age.group, nNew)
  
  # Assign new nodes to households
  age.grp <- get_attr(dat, "age.grp") # age of all nodes, include the new nodes
  household <- get_attr(dat, "hh.ids" # all households, not include new nodes 
                        #"household"
                        )
  
  newHH <- rep(NA, length(newAges))
  for (i in seq_along(unique(attr_age.grp))
       ) { # since all new ages ==0, the first age group selected is for all new nodes
    newHH[attr_age.grp == unique(attr_age.grp)[i]
          ] <- # correspond to all new nodes, for each of these nodes, randomly assign households with nodes at age group i
      sample(household[which(age.grp == unique(attr_age.grp)[i])], # hh.ids of old and new nodes of with nodes at age grp i, less than number of nodes below
             sum(attr_age.grp == unique(attr_age.grp)[i]), # total number of new nodes in particular age grp
             replace = TRUE)
  }
  
  # Update household edgelist
  heads <- cbind((length(household) + 1):(length(household) + nNew), newHH)
  tails <- cbind(which(household %in% newHH), household[which(household %in% newHH)])
  new.edges <- merge(heads, tails, by.x = 2, by.y = 2)[, 2:3]
  new.edgelist <- as.matrix(rbind(dat$el[[dat$num.nw]], setNames(new.edges, c(".head", ".tail"))))
  attr(new.edgelist, 'n') <- attr(dat$el[[dat$num.nw]], 'n')
  dat$el[[dat$num.nw]] <- new.edgelist
  
  
  
  # Disease status and related
  dat <- append_attr(dat, "status", "s", nNew)
  dat <- append_attr(dat, "infTime", NA, nNew)
  
  dat <- append_attr(dat, "statusTime", 0, nNew)
  dat <- append_attr(dat, "clinical", NA, nNew)
  dat <- append_attr(dat, "hospit", NA, nNew)
  dat <- append_attr(dat, "dxStatus", NA, nNew)
  dat <- append_attr(dat, "dxTime", NA, nNew)
  dat <- append_attr(dat, "vax", 0, nNew)
  dat <- append_attr(dat, "vax1Time", NA, nNew)
  dat <- append_attr(dat, "vax2Time", NA, nNew)
  dat <- append_attr(dat, "vax3Time", NA, nNew)
  dat <- append_attr(dat, "vax4Time", NA, nNew)
  dat <- append_attr(dat, "isolate", NA, nNew)
  dat <- append_attr(dat, "isoTime", NA, nNew)
  
  # specific to my project
  dat <- append_attr(dat, "deg.x_layer", rep(0L, nNew),         nNew)
  dat <- append_attr(dat, "deg_school",  rep(0L, nNew),         nNew)
  dat <- append_attr(dat, "deg_work",    rep(0L, nNew),         nNew)
  dat <- append_attr(dat, "no.contact",  rep(0L, nNew),         nNew)
  # dat <- append_attr(dat, "hh.ids",      sample(unique(dat$attr$hh.ids), nNew, replace = TRUE), # replace = TRUE allows different people go to same hh
  #                    nNew
  #                    ) # self added
  
  return(dat)
}

#' @rdname moduleset-corporate
#' @export
arrival_covid_corporate <- function(dat, at) {
  
  # Parameters
  a.rate   <- get_param(dat, "a.rate")
  
  ## Process
  num <- dat$epi$num[1]
  nNew <- rpois(1, a.rate * num)
  
  ## Update Attr
  if (nNew > 0) {
    dat <- setNewAttr_covid_corporate(dat, at, nNew)
  }
  
  # Update Networks
  dat <- arrive_nodes(dat, nNew)
  
  ## Output
  dat <- set_epi(dat, "nNew", at, nNew)
  
  return(dat)
}

setNewAttr_covid_corporate <- function(dat, at, nNew) {
  
  dat <- append_core_attr(dat, at, nNew)
  
  # Age related
  arrival.age <- get_param(dat, "arrival.age")
  newAges <- rep(arrival.age, nNew)
  dat <- append_attr(dat, "age", newAges, nNew)
  
  age.breaks <- seq(0, 200, 10)
  attr_age.grp <- cut(newAges, age.breaks, labels = FALSE, right = FALSE)
  dat <- append_attr(dat, "age.grp", attr_age.grp, nNew)
  
  vax.age.group <- rep(NA, length(newAges))
  vax.age.group[newAges < 5] <- 1
  vax.age.group[newAges >= 5 & newAges < 18] <- 2
  vax.age.group[newAges >= 18 & newAges < 50] <- 3
  vax.age.group[newAges >= 50 & newAges < 65] <- 4
  vax.age.group[newAges >= 65] <- 5
  dat <- append_attr(dat, "vax.age.group", vax.age.group, nNew)
  
  # Assign new nodes to households
  age.grp <- get_attr(dat, "age.grp")
  household <- get_attr(dat, "household")
  
  newHH <- rep(NA, length(newAges))
  for (i in seq_along(unique(attr_age.grp))) {
    newHH[attr_age.grp == unique(attr_age.grp)[i]] <-
      sample(household[which(age.grp == unique(attr_age.grp)[i])],
             sum(attr_age.grp == unique(attr_age.grp)[i]), replace = TRUE)
  }
  
  # Update household edgelist
  heads <- cbind((length(household) + 1):(length(household) + nNew), newHH)
  tails <- cbind(which(household %in% newHH), household[which(household %in% newHH)])
  new.edges <- merge(heads, tails, by.x = 2, by.y = 2)[, 2:3]
  new.edgelist <- as.matrix(rbind(dat$el[[dat$num.nw]], setNames(new.edges, c(".head", ".tail"))))
  attr(new.edgelist, 'n') <- attr(dat$el[[dat$num.nw]], 'n')
  dat$el[[dat$num.nw]] <- new.edgelist
  
  dat <- append_attr(dat, "household", newHH, nNew)
  
  # Disease status and related
  dat <- append_attr(dat, "status", "s", nNew)
  dat <- append_attr(dat, "infTime", NA, nNew)
  
  dat <- append_attr(dat, "statusTime", 0, nNew)
  dat <- append_attr(dat, "clinical", NA, nNew)
  dat <- append_attr(dat, "hospit", NA, nNew)
  dat <- append_attr(dat, "dxStatus", NA, nNew)
  dat <- append_attr(dat, "dxTime", NA, nNew)
  dat <- append_attr(dat, "vax", 0, nNew)
  dat <- append_attr(dat, "vax1Time", NA, nNew)
  dat <- append_attr(dat, "vax2Time", NA, nNew)
  dat <- append_attr(dat, "vax3Time", NA, nNew)
  dat <- append_attr(dat, "vax4Time", NA, nNew)
  dat <- append_attr(dat, "isolate", NA, nNew)
  dat <- append_attr(dat, "isoTime", NA, nNew)
  dat <- append_attr(dat, "non.office", 1, nNew)
  
  return(dat)
}

