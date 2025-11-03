
#' @rdname moduleset-gmc19
#' @export

home_update_gmc19 <- function(dat, at) {
  
  # which layer is home? we added it last in init
  home_ix <- dat$num.nw
  
  # current home edgelist
  el_home <- dat$el[[home_ix]]
  
  # number of current nodes
  cur_n <- dat$run$num # 11781 here, may have problem as death not counted
  
  # how many nodes the home layer has last time
  old_n <- attr(el_home, "n")
  if (is.null(old_n)) old_n <- cur_n  # first call
  
  
  # attrs we need
  active <- get_attr(dat, "active") # NA means the node has died
  hh.ids <- get_attr(dat, "hh.ids")
  
  # drop edges with inactive nodes
  keep <- # status of whether an edge contain inactive node
    active[el_home[, 1]] == 1L &  # for head side of each edge, identify the active status of the head node (1 - active, 0 - inactive)
    active[el_home[, 2]] == 1L # for tail side of each edge, identify the active status of the tail node (1 - active, 0 - inactive)
  
  el_home <- el_home[keep, ]
  
  # add edge list of new nodes
 
  ## identify new ids
  new_node_ids <- dat$temp$new_node_ids
  
  # 3a) connect each new node to *existing* members of its HH
  for (id in  new_node_ids) {
    h <- hh.ids[id]
    if (is.na(h)) next
    
    # existing members = nodes that were already there last step
    old_members <- which(hh.ids[seq_len(old_n)] == h & active[seq_len(old_n)] == 1L)
    
    if (length(old_members) > 0L) {
      k <- k + 1L
      add_list[[k]] <- cbind(id, old_members)
    } # 10/31 reach here. nned to better define cur_n, old_n, to see if available alternative can be used.
  }
  
  # 3b) ALSO connect new nodes among themselves if multiple new nodes
  #     landed in the same HH this step
  hh_new <- hh_ids[new_ids]
  for (h in unique(hh_new[!is.na(hh_new)])) {
    new_in_h <- new_ids[hh_ids[new_ids] == h]
    if (length(new_in_h) > 1L) {
      k <- k + 1L
      add_list[[k]] <- t(combn(new_in_h, 2))
    }
  }
  
  

}


# home_update_gmc19 <- function(dat, at) {
#   
#   # at t=1, return dat, as home layer would be undertaken in resim_nets.FUN 
#   if (at == 1) return(dat)
#   
#   # pull current home network
#   nw_home <- dat$nw$home
#   
#   # wipe all existing household edges
#   dat$el[[4]]
#   ec <- network::network.edgecount(nw_home)
#   if (ec > 0) {
#     nw_home <- network::delete.edges(nw_home, seq_len(ec))
#   }
#   
#   active <- which(dat$attr$active == 1)
#   hh_id  <- dat$attr$hh_id
#   
#   # for each household, make a clique among active members
#   hhs <- unique(hh_id[active])
#   hhs <- hhs[!is.na(hhs)]
#   
#   for (h in hhs) {
#     members <- active[hh_id[active] == h]
#     if (length(members) > 1) {
#       el <- t(combn(members, 2))
#       nw_home <- network::add.edges(nw_home, el[,1], el[,2])
#     }
#   }
#   
#   dat$nw$home <- nw_home
#   return(dat)
# }