IsingEntrophy <- function(
  graph,
  thresholds,
  beta = 1,
  conditional = numeric(0), # Indices of nodes to condition on
  marginalize = numeric(0),
  base = 2,
  responses = c(0L, 1L)
  ){
  stopifnot(isSymmetric(graph))
  checkResponses(responses)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }

  Lik <- IsingLikelihood(graph, thresholds, beta, responses)

  varNames <- names(Lik)[-1]

  if (any(marginalize %in% conditional)){
    stop("can not marginalize over nodes to condition on")
  }

  if (length(marginalize) > 0){
    Lik <- Lik %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(varNames[-marginalize]))) %>%
      dplyr::summarise(Probability = sum(.data$Probability), .groups = "drop")
  }

  if (length(conditional) > 0){
    Lik <- Lik %>% dplyr::group_by(dplyr::across(dplyr::all_of(varNames[conditional])))
  } else {
    Lik <- Lik %>% dplyr::ungroup()
  }

  condLik <- Lik %>%
    dplyr::summarise(
      P = sum(.data$Probability),
      Entrophy = -sum(.data$Probability / sum(.data$Probability) *
                        log(.data$Probability / sum(.data$Probability), base)),
      .groups = "drop"
      )

  Ent <- sum(condLik$P * condLik$Entrophy)

  return(Ent)
}


NodeInformation <- function(
  graph,
  thresholds,
  beta = 1,
  base = 2,
  responses = c(0L, 1L)
){
  stopifnot(isSymmetric(graph))
  checkResponses(responses)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }

  # Shannon information of node to whole graph
  
  sapply(seq_len(ncol(graph)),
         function(Node){
           graphEntrophy <- IsingEntrophy(graph, thresholds, beta, responses=responses,base=base, marginalize = Node)
           nodeEntrophy <- IsingEntrophy(graph, thresholds, beta, responses=responses,base=base, conditional = Node)
           return(graphEntrophy - nodeEntrophy)
         })
}
