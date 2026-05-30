# Internal helper: validate the 'responses' argument.
# Requires at least two numeric response options. Real-valued responses such as
# seq(-1, 1, by = 0.5) are allowed; states are treated as their numeric values.
checkResponses <- function(responses)
{
  if (length(responses) < 2)
  {
    stop("'responses' must contain at least two response options.")
  }
  if (anyNA(responses) || !is.numeric(responses))
  {
    stop("'responses' must be numeric.")
  }
  invisible(responses)
}

# Wrapper function:
IsingSampler <- function(n, graph, thresholds, beta = 1, nIter = 100, responses = c(0L, 1L), method = c("MH","CFTP","direct"), CFTPretry = 10, constrain)
{
  stopifnot(!missing(graph)|!missing(thresholds))
  stopifnot(isSymmetric(graph))
  checkResponses(responses)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  method <- method[1]
  if (! method %in% c("MH","CFTP","direct")) stop("method must be 'MH', 'CFTP', or 'direct'")

  # Exact sampling (CFTP) is only implemented for binary (two-level) responses:
  if (method == "CFTP" && length(responses) != 2)
  {
    stop("method = 'CFTP' (exact sampling) is currently only implemented for two response options. Use method = 'MH' or method = 'direct' for more than two response options.")
  }

  # Check constrains:
  if (missing(constrain))
  {
    constrain <- matrix(NA_real_,n,ncol(graph))
  }
  # responses and constrain are passed to C++ as doubles so that non-integer
  # response options (e.g. seq(-1, 1, by = 0.5)) are preserved:
  storage.mode(constrain) <- "double"

  if (method %in% c("MH","CFTP"))
  {
    try <- 1
    
    repeat{
      Res <- IsingSamplerCpp(as.integer(n), graph, thresholds, beta, as.integer(nIter), as.numeric(responses), as.logical(method == "CFTP"), constrain)

      if (any(is.na(Res)) & method == "CFTP")
      {
        if (try < CFTPretry)
        {
          cat("\n Restarting CFTP chain, attempt: ",try,"\n")
          try <- try + 1
        } else
        {
          warning(paste("NA's detected, this means that no CFTP sample was drawn after 100 couplings to the past and",CFTPretry,"resets of the chain. Use higher nIter value or method='MH' for inexact sampling."))
          break
        }
      } else break
    }

    # Preserve integer storage when all responses are integer-valued (keeps the
    # binary / integer behaviour identical to previous versions):
    if (all(responses == round(responses))) storage.mode(Res) <- "integer"
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
  checkResponses(responses)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  N <- nrow(graph)
  Allstates <- do.call(expand.grid,lapply(1:N,function(x)responses))
  P <- exp(- beta * apply(Allstates,1,function(s)H(graph,s,thresholds)))  
  return(Allstates[sample(1:nrow(Allstates),n,TRUE,prob=P),])
}