# Univariate:
EstimateIsingUni <- function(data, responses, beta = 1, ...){
  if (missing(responses)){
    responses <- sort(unique(c(data)))
  }
  
  if (length(responses) != 2){
    stop("Binary data required")
  }
  
  # Rescale data to binary:
  binarize <- function(x, responses){
    x2 <- x
    x2[x==responses[1]] <- 0
    x2[x==responses[2]] <- 1
    x2
  }
  data <- binarize(data, responses)
  
  
}




# Bivariate:
EstimateIsingBi <- function(...){stop("Not implemented")}