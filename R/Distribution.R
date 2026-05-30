IsingSumLikelihood <- function(graph, thresholds, beta, responses = c(0L,1L), delta = 0)
{
  stopifnot(isSymmetric(graph))
  checkResponses(responses)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  N <- nrow(graph)
  if (length(thresholds) == 1L) thresholds <- rep(thresholds, N)
  if (length(delta) == 1L) delta <- rep(delta, N)
  Allstates <- do.call(expand.grid,lapply(1:N,function(x)responses))
  P <- exp(- beta * apply(Allstates,1,function(s)H(graph,s,thresholds,delta)))
  SumScores <- rowSums(Allstates)

  df <- plyr::ddply(data.frame(Sum = SumScores, P = P),"Sum",plyr::summarize,P=sum(P))
  df$P <- df$P / sum(df$P)
  return(df)
}

IsingLikelihood <- function(graph, thresholds, beta, responses = c(0L,1L), potential = FALSE, delta = 0)
{
  stopifnot(isSymmetric(graph))
  checkResponses(responses)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  N <- nrow(graph)
  if (length(thresholds) == 1L) thresholds <- rep(thresholds, N)
  if (length(delta) == 1L) delta <- rep(delta, N)
  Allstates <- do.call(expand.grid,lapply(1:N,function(x)responses))
  P <- exp(- beta * apply(Allstates,1,function(s)H(graph,s,thresholds,delta)))
  if (potential){
    df <- cbind(Potential = P, Allstates)    
  } else {
    df <- cbind(Probability = P / sum(P), Allstates)
  }

  return(df)
}

IsingStateProb <- function(s,graph,thresholds,beta,responses=c(0L,1L), delta = 0)
{
  checkResponses(responses)
  if (!is.list(s)) s <- list(s)
  N <- length(s[[1]])
  if (length(thresholds) == 1L) thresholds <- rep(thresholds, N)
  if (length(delta) == 1L) delta <- rep(delta, N)
  Allstates <- do.call(expand.grid,lapply(1:N,function(x)responses))
  Dist <- exp(- beta * apply(Allstates,1,function(s)H(graph,s,thresholds,delta)))
  Z <- sum(Dist)

  sapply(s, function(x)exp(-beta*H(graph,x,thresholds,delta))/Z)
}