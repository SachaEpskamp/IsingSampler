# Wrapper function:
IsingSampler <- function(n, graph, thresholds, beta = 1, nIter = 100, responses = c(0L, 1L), method = c("MH","CFTP","direct"))
{
  stopifnot(!missing(graph)|!missing(thresholds))
  stopifnot(isSymmetric(graph))  
  stopifnot(length(responses)==2)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  method <- method[1]
  if (! method %in% c("MH","CFTP","direct")) stop("method must be 'MH', 'CFTP', or 'direct'")
  
  if (method %in% c("MH","CFTP"))
  {
    Res <- IsingSamplerCpp(as.integer(n), graph, thresholds, beta, as.integer(nIter), as.integer(responses), as.logical(method == "CFTP"))                 
    
    if (any(is.na(Res)) & method == "CFTP")
    {
      warning("NA's detected, this means that no CFTP sample was drawn after 100 couplings to the past. Use higher nIter value or method='MH' for inexact sampling.")
    }
  } else 
  {
    Res <- IsingDir(n, graph, thresholds, beta, responses)
  }
  
  return(Res)
}

## direct sampling function:
IsingDir <- function(n, graph, thresholds, beta,responses = c(0L,1L))
{
  stopifnot(isSymmetric(graph))  
  stopifnot(length(responses)==2)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  N <- nrow(graph)
  Allstates <- do.call(expand.grid,lapply(1:N,function(x)c(responses[1],responses[2])))
  P <- exp(- beta * apply(Allstates,1,function(s)H(graph,s,thresholds)))  
  return(Allstates[sample(1:nrow(Allstates),n,TRUE,prob=P),])
}