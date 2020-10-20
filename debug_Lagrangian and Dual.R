library(network)

##################### EDGE LIST ###############################


setGeneric("edgelist", function(x) { stop("Not implemented.") })

setMethod("edgelist", c(x = "InfinitySparseMatrix"), function(x) {
  
  return(data.frame(i = x@rownames[x@rows],  j = x@colnames[x@cols], dist = x@.Data))
})

setMethod("edgelist", c(x = "matrix"), function(x) {
  return(edgelist(as.InfinitySparseMatrix(x)))
})

setMethod("edgelist", c(x = "data.frame"), function(x){
  stopifnot(ncol(x)==3, setequal(colnames(x), c('control', 'treated', 'distance')),
            is.numeric(x$distance))
  data.frame(i = as.character(x$treated), j=as.character(x$control), dist=x$distance)
})

###############################################################


evaluate_dual <- function(distances, solution) {
  stopifnot(is(solution, "FullmatchMCFSolutions"),
            nrow(solution@subproblems)==1 || !any(solution@subproblems[["flipped"]]) )
  flipped  <- solution@subproblems[1L, "flipped"]

  sum_supply_price <- sum(solution@nodes$supply * solution@nodes$price)
  
  suppressWarnings(
    bookkeeping_ij <- left_join(solution@arcs@bookkeeping,
                                as.data.frame(unclass(solution@nodes)),
                                by = c("start" = "name", "groups"="groups")) %>%
      left_join(y = subset(solution@nodes,
                           is.na(upstream_not_down)), #assumes bookkeeping arcs...
                by = c("end" = "name", "groups"="groups"), #... terminate only in bookkeeping nodes
                suffix = c(x = ".i", y = ".j"))
  )
  
  nonpositive_flowcosts_bookkeeping  <-
    pmin(0,
         bookkeeping_ij$price.j - bookkeeping_ij$price.i
    ) * bookkeeping_ij$capacity
  
  eld <- edgelist(distances)
  

  
  
  if (!flipped) {
    cantadd <- unique(eld$i)
    canadd <- unique(eld$j)
  } else {
    cantadd <- unique(eld$j)
    canadd <- unique(eld$i)
  }
  
  cantadd <- as.character(cantadd)
  canadd <- as.character(canadd)
  
  upstream <- split(solution@nodes, solution@nodes$upstream_not_down)
  
  ## can't impute a node price for these missing node prices
  if (any(!(cantadd %in% upstream[["TRUE"]]$name))) {
    stop("Cannot impute node price for upstream nodes (usually treatment) that were not included in original matching problem.")
  }

  impute_price <- min(solution@nodes$price[is.na(solution@nodes$upstream_not_down)])
  
  newnames <- canadd[!(canadd %in% solution@nodes$name)]
  k <- length(newnames)
  upstream[['FALSE']] <- rbind(upstream[['FALSE']],
                               new("NodeInfo", data.frame(stringsAsFactors = FALSE,
                                                          name = newnames,
                                                          price = rep(impute_price, k),
                                                          upstream_not_down = rep(FALSE, k),
                                                          supply = rep(0L, k),
                                                          groups = as.factor(rep(NA, k))))) # TODO: get any group labels from the distance?
  
  suffices  <-
    if (!flipped) c(x =".i", y =".j") else c(x =".j", y =".i")
  suppressWarnings(
    matchable_ij <- left_join(eld,
                              upstream[['TRUE']],
                              by = c("i" = "name")) %>%
      left_join(y = upstream[["FALSE"]],
                by = c("j" = "name"),
                suffix = suffices) # if nec., flip right at end.
  )
  
  
  ### DEBUG ###
#  print(bookkeeping_ij,max=10000)
  
#  print(matchable_ij)
  
  #############
  
  
  
  
  
  nonpositive_flowcosts_matchables <-
    pmin(0,
         matchable_ij$dist - (matchable_ij$price.i - matchable_ij$price.j)
    )
  
#   nonpositive_flowcosts_matchables = c()
#   index=1
#   for (index in 1:99) {
# 
#     # NOTE: Should this be min(dist), or min(dist - (p_i - p_j))?
#     matchables_j = filter(matchable_ij, j==index)
# #    print(paste0("index: ", index, "ncontrol: ", nrow(matchables_j)))
#     flowcost =  min(matchables_j$dist - (matchables_j$price.i - matchables_j$price.j))
#     flowcost = min(0, flowcost)
#     nonpositive_flowcosts_matchables =
#       c(nonpositive_flowcosts_matchables, flowcost)
# #    print(index)
#   }

    
#   print("*************Dual Value Debug***************")
# #  print(nonpositive_flowcosts_matchables)
#   print(paste0("supply price: ",sum_supply_price))
#   print(paste0("flowcost_matchables: ",sum(nonpositive_flowcosts_matchables)))
#   print(paste0("flowcost_bookkeeping: ",sum(nonpositive_flowcosts_bookkeeping)))
#   print("********************************************")
  
  
  return(sum_supply_price +
           sum(nonpositive_flowcosts_bookkeeping) +
           sum(nonpositive_flowcosts_matchables))
}

#########################################
############# Lagrangian ################
#########################################

evaluate_lagrangian <- function(distances, solution) {
  stopifnot(is(solution, "MCFSolutions"),
            nrow(solution@subproblems)==1 || !any(solution@subproblems[["flipped"]]))
  flipped  <- solution@subproblems[1L, "flipped"]
  
  suppressWarnings(# re factor conversion
    main_ij <- left_join(solution@arcs@matches,
                         subset(solution@nodes, upstream_not_down),
                         by = c("upstream" = "name")) %>%
      left_join(y = subset(solution@nodes, !upstream_not_down),
                by = c("downstream" = "name"),
                suffix = c(x = ".i", y = ".j"))
  )
  
  eld <- edgelist(distances)
  suppressWarnings(
    if (!flipped) {
      main_ij <- left_join(main_ij,
                           eld,
                           by = c("upstream" = "i", "downstream"= "j"),
                           suffix = c(x = "", y = ".dist"))
    } else {
      main_ij <- left_join(main_ij,
                           eld,
                           by = c("upstream" = "j", "downstream"= "i"),
                           suffix = c(x = "", y = ".dist"))
    }
  )
  
  suppressWarnings(
    bookkeeping_ij <- left_join(solution@arcs@bookkeeping,
                                as.data.frame(unclass(solution@nodes)),
                                by = c("start" = "name", "groups"="groups")) %>%
      left_join(y = subset(solution@nodes,
                           is.na(upstream_not_down)),#assumes bookkeeping arcs...
                by = c("end" = "name", "groups"="groups"),#...terminate only in bookkeeping nodes
                suffix = c(x = ".i", y = ".j"))
  )
  
  sum_supply_price <- sum(solution@nodes$supply * solution@nodes$price)
  
  sum_flow_cost_main <- sum(main_ij$dist - (main_ij$price.i - main_ij$price.j))
  sum_flow_cost_bookkeeping <- sum(bookkeeping_ij$flow * (0 - (bookkeeping_ij$price.i - bookkeeping_ij$price.j)))
  
  
#  print(main_ij)
  
#   print("**********Lagrangian Value Debug************")
# #  print(paste0("eld length: ", length(eld)))
#   print(paste0("supply price: ",sum_supply_price))
#   print(paste0("flowcost_matchables: ",sum_flow_cost_main))
#   print(paste0("flowcost_bookkeeping: ",sum_flow_cost_bookkeeping))
#   print("********************************************")
  
  return(sum_flow_cost_main + sum_flow_cost_bookkeeping + sum_supply_price)
}