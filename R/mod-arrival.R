

#' @rdname moduleset-gmc19
#' @export
arrival <- function(dat, at) {
  # if (at>200) browser(), helper for debugging at a later time step
  
  ## Input
  # Attributes
  
  # Parameters
  a.rate   <- get_param(dat, "a.rate")
  
  ## Process
  num <- get_epi(dat, "num", at = 1)
  nNew <- rpois(1, a.rate * num) 

  
  ## Update Attr
  if (nNew > 0) {
    dat <- init_new_nodes_attrs(dat, nNew)
  }
  
  # Update Networks
  dat <- arrive_nodes(dat, nNew)
  
  ## Output
  dat <- set_epi(dat, "nNew", at, nNew)
  
  return(dat)
}



init_new_nodes_attrs <- function(dat, n_new) {
  current_timestep <- get_current_timestep(dat)
  # Core attributes (necessary for EpiModel)
  dat <- append_core_attr(dat, current_timestep, n_new)
  def_attrs <- get_default_attrs(dat)
  for (attr_name in names(def_attrs)) {
    dat <- append_attr(dat, attr_name, def_attrs[[attr_name]], n_new)
  }
  dat <- make_computed_attrs(dat, n_new, post_init = TRUE)
  return(dat)
}


# setNewAttr_covid_gmc19 <- function(dat, at, nNew) {
#   
#   dat <- append_core_attr(dat, at, nNew) # add nNew nodes to dat
#   
#   # Age related
#   arrival.age <- get_param(dat, "arrival.age") # those arrive were 0 year old
#   newAges <- rep(arrival.age, nNew) # make a vector of arrival age
#   dat <- append_attr(dat, "age", newAges, nNew) # add to dat
#   
#   # age.breaks <- seq(0, 200, 10)
#   # attr_age.grp <- cut(newAges, age.breaks, labels = FALSE, right = FALSE)
#   AGE_LABELS <- c("0-9y","10-19y","20-29y","30-39y","40-59y","60+y")
#   age.breaks  <- c(0, 10, 20, 30, 40, 60, Inf)
#   attr_age.grp <- cut(newAges,
#                         breaks =  age.breaks ,
#                         labels = AGE_LABELS,
#                         right  = FALSE) %>% as.character() # different from corportate which use numeric grp
#   dat <- append_attr(dat, "age.grp",attr_age.grp, nNew) # add categorical ages to new nodes
#   
#   ## age at vax
#   vax.age.group <- rep(NA, length(newAges))
#   vax.age.group[newAges < 10] <- 1
#   vax.age.group[newAges >= 10 & newAges < 20] <- 2
#   vax.age.group[newAges >= 20 & newAges < 30] <- 3
#   vax.age.group[newAges >= 30 & newAges < 40] <- 4
#   vax.age.group[newAges >= 40 & newAges < 60] <- 5
#   vax.age.group[newAges >= 60 ] <- 6
#   dat <- append_attr(dat, "vax.age.group", vax.age.group, nNew)
#   
#   # Assign new nodes to households
#   age.grp <- get_attr(dat, "age.grp") # age of all nodes, include the new nodes
#   household <- get_attr(dat, "hh.ids" # all households, not include new nodes 
#                         #"household"
#                         )
#   
#   newHH <- rep(NA, length(newAges))
#   for (i in seq_along(unique(attr_age.grp))
#        ) { # since all new ages ==0, the first age group selected is for all new nodes
#     newHH[attr_age.grp == unique(attr_age.grp)[i]
#           ] <- # correspond to all new nodes, for each of these nodes, randomly assign households with nodes at age group i
#       sample(household[which(age.grp == unique(attr_age.grp)[i])], # hh.ids of old and new nodes of with nodes at age grp i, less than number of nodes below
#              sum(attr_age.grp == unique(attr_age.grp)[i]), # total number of new nodes in particular age grp
#              replace = TRUE)
#   }
#   
#   # Update household edgelist
#   heads <- cbind((length(household) + 1):(length(household) + nNew), newHH)
#   tails <- cbind(which(household %in% newHH), household[which(household %in% newHH)])
#   new.edges <- merge(heads, tails, by.x = 2, by.y = 2)[, 2:3]
#   new.edgelist <- as.matrix(rbind(dat$el[[dat$num.nw]], setNames(new.edges, c(".head", ".tail"))))
#   attr(new.edgelist, 'n') <- attr(dat$el[[dat$num.nw]], 'n')
#   dat$el[[dat$num.nw]] <- new.edgelist
#   
#   
#   
#   # Disease status and related
#   dat <- append_attr(dat, "status", "s", nNew)
#   dat <- append_attr(dat, "infTime", NA, nNew)
#   
#   dat <- append_attr(dat, "statusTime", 0, nNew)
#   dat <- append_attr(dat, "clinical", NA, nNew)
#   dat <- append_attr(dat, "hospit", NA, nNew)
#   dat <- append_attr(dat, "dxStatus", NA, nNew)
#   dat <- append_attr(dat, "dxTime", NA, nNew)
#   dat <- append_attr(dat, "vax", 0, nNew)
#   dat <- append_attr(dat, "vax1Time", NA, nNew)
#   dat <- append_attr(dat, "vax2Time", NA, nNew)
#   dat <- append_attr(dat, "vax3Time", NA, nNew)
#   dat <- append_attr(dat, "vax4Time", NA, nNew)
#   dat <- append_attr(dat, "isolate", NA, nNew)
#   dat <- append_attr(dat, "isoTime", NA, nNew)
#   
#   # specific to my project
#   dat <- append_attr(dat, "deg.x_layer", rep(0L, nNew),         nNew)
#   dat <- append_attr(dat, "deg_school",  rep(0L, nNew),         nNew)
#   dat <- append_attr(dat, "deg_work",    rep(0L, nNew),         nNew)
#   dat <- append_attr(dat, "no.contact",  rep(0L, nNew),         nNew)
#   # dat <- append_attr(dat, "hh.ids",      sample(unique(dat$attr$hh.ids), nNew, replace = TRUE), # replace = TRUE allows different people go to same hh
#   #                    nNew
#   #                    ) # self added
#   
#   return(dat)
# }
# 
# 
