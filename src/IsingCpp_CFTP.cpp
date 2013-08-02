#include <Rcpp.h>
#include <climits>
using namespace Rcpp;

// FUNCTIONS FOR EXACT SAMPLING //

// Inner function to resize list:
List resize( const List& x, int n ){
    int oldsize = x.size() ;
    List y(n) ;
    for( int i=0; i<oldsize; i++) y[i] = x[i] ;
    return y ;
}

// Inner function to simulate random uniforms in a matrix:
NumericMatrix RandMat(int nrow, int ncol)
 {
  int N = nrow * ncol;
  NumericMatrix Res(nrow,ncol);
  NumericVector Rands  = runif(N);
   for (int i = 0; i < N; i++) 
  {
    Res[i] = Rands[i];
  }
  return(Res);
 }

// Computes maximal and minimal probability of node flipping:
NumericVector PplusMinMax(int i, NumericMatrix J, IntegerVector s, NumericVector h, double beta, IntegerVector responses)
{
  // The function computes the probability that node i is in Response 1 instead of 0, given all other nodes, which might be missing.
  // Output: minimal and maximal probablity
  
  NumericVector H0(2, h[i] * responses[0]); // relevant part of the Hamiltonian for state = 0
  NumericVector H1(2, h[i] * responses[1]); // relevant part of the Hamiltonian for state = 1
  
  NumericVector Res(2);
  
  int N = J.nrow();
  NumericVector TwoOpts(2);
  
  for (int j=0; j<N; j++)
  {
    if (i != j)
    {
      if (s[j] != INT_MIN)
      {
       H0[0] += J(i,j) * responses[0] * s[j];
       H0[1] += J(i,j) * responses[0] * s[j];
       H1[0] += J(i,j) * responses[1] * s[j];
       H1[1] += J(i,j) * responses[1] * s[j]; 
      } else 
      {
               
        TwoOpts[0] = J(i,j) * responses[1] * responses[0];
        TwoOpts[1] = J(i,j) * responses[1] * responses[1];

        if (TwoOpts[1] > TwoOpts[0])
        {
          H1[0] += TwoOpts[0];
          H1[1] += TwoOpts[1];
          
          H0[0] += J(i,j) * responses[0] * responses[0];
          H0[1] += J(i,j) * responses[0] * responses[1];
        } else 
        {
          H1[0] += TwoOpts[1];
          H1[1] += TwoOpts[0];          
          
          H0[0] += J(i,j) * responses[0] * responses[1];
          H0[1] += J(i,j) * responses[0] * responses[0];
        }
      }
    }
  }

  Res[0] = exp(beta * H1[0]) / ( exp(beta * H0[0]) + exp(beta * H1[0]) );
  Res[1] = exp(beta * H1[1]) / ( exp(beta * H0[1]) + exp(beta * H1[1]) );
  
  
  return(Res);
}
       
// Inner function:
IntegerVector IsingEx(NumericMatrix graph, NumericVector thresholds, double beta, int nIter, IntegerVector responses, bool exact)
{
  // Parameters and results vector:
  int N = graph.nrow();
  IntegerVector state(N, INT_MIN);
  double u;
  NumericVector P(2);
  int maxChain = 100;
  List U(1);
  int minT = 0;
  bool anyNA = true;
    
  do
  { 
    // Resize U if needed:
    if (minT > 0)
    {
      U = resize(U, minT+1);
    }
    
    // Generate new random numbers:
    U[minT] = RandMat(nIter, N);
    
    // Initialize states:
    for (int i=0; i<N; i++)
    {
      if (exact)
      {
        state[i] = INT_MIN;
      } else 
      {
        state[i] = ifelse(runif(1) < 0.5, responses[1], responses[0])[0];
      }
    }    

    // START ALGORITHM
    for (int t=minT; t > -1;  t--)
    {
      for (int it=0;it<nIter;it++)
      {
        NumericMatrix Ucur = U[t];
        for (int node=0;node<N;node++)
        {
          u = Ucur(it, node);
          P = PplusMinMax(node, graph, state, thresholds, beta, responses);
          if (u < P[0])
          {
            state[node] = responses[1];
          } else if (u >= P[1])
          {
            state[node] = responses[0];
          } else 
          {
            state[node] = INT_MIN;
          }
        }
      }
    }
    
    anyNA = false;
    if (exact)
    {
      if (minT < maxChain)
      {
       for (int i=0; i<N; i++)
       {
        if (state[i] == INT_MIN)
        {
          anyNA = true;
        }
        } 
      } 
    }    
    minT++;
    
  } while (anyNA);

  // Rf_PrintValue(wrap(minT));
  return(state);
}


// FUNCTIONS FOR METROPOLIS SAMPLER //
double Pplus(int i, NumericMatrix J, IntegerVector s, NumericVector h, double beta, IntegerVector responses)
{
  // The function computes the probability that node i is in Response 1 instead of 0, given all other nodes, which might be missing.
  // Output: minimal and maximal probablity
  
  double H0 = h[i] * responses[0]; // relevant part of the Hamiltonian for state = 0
  double H1 = h[i] * responses[1]; // relevant part of the Hamiltonian for state = 1
  double Res;
  
  int N = J.nrow();
  
  for (int j=0; j<N; j++)
  {
    if (i != j)
    {
       H0 += J(i,j) * responses[0] * s[j];
       H1 += J(i,j) * responses[1] * s[j];
    }
  }
  
  return(exp(beta * H1) / ( exp(beta * H0) + exp(beta * H1) ));
}


IntegerVector IsingMet(NumericMatrix graph, NumericVector thresholds, double beta, int nIter, IntegerVector responses)
{
  // Parameters and results vector:
  int N = graph.nrow();
  IntegerVector state =  ifelse(runif(N) < 0.5, responses[1], responses[0]);
  double u;
  double P;
    
    // START ALGORITHM
    for (int it=0;it<nIter;it++)
    {
      for (int node=0;node<N;node++)
      {
        u = runif(1)[0];
        P = Pplus(node, graph, state, thresholds, beta, responses);
        if (u < P)
        {
          state[node] = responses[1];
        } else 
        {
          state[node] = responses[0];
        }
      }
    }
   
  return(state);
}

// OVERAL FUNCTION //
// [[Rcpp::export]]
IntegerMatrix IsingSamplerCpp(int n, NumericMatrix graph, NumericVector thresholds, double beta, int nIter, IntegerVector responses, bool exact)
{
  int Ni = graph.nrow();
  IntegerMatrix Res(n,Ni);
  IntegerVector state(Ni);
  if (exact)
  {
    for (int s=0;s<n;s++)
    {
      state = IsingEx(graph, thresholds, beta, nIter, responses, exact);
      for (int i=0;i<Ni;i++) Res(s,i) = state[i];
    }
  } else 
  {
    for (int s=0;s<n;s++)
    {
      state = IsingMet(graph, thresholds, beta, nIter, responses);
      for (int i=0;i<Ni;i++) Res(s,i) = state[i];
    }
  }
  
  return(Res);
}



// Hamiltonian:
// [[Rcpp::export]]
double H(NumericMatrix J, IntegerVector s, NumericVector h)
{
  double Res = 0;
  int N = J.nrow();
  for (int i=0;i<N;i++)
  {
    Res -= h[i] * s[i];
    for (int j=i;j<N;j++)
    {
      if (j!=i) Res -= J(i,j) * s[i] * s[j];
    }
  }
  return(Res);
}

// R wrappers:
/*** R
IsingSampler <- function(n, graph, thresholds, beta = 1, nIter = 100, responses = c(0L, 1L), method = c("metropolis","exact","direct"))
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
  if (! method %in% c("metropolis","exact","direct")) stop("method must be 'metropolis', 'exact', or 'direct'")
  
  if (method %in% c("metropolis","exact"))
  {
    Res <- IsingSamplerCpp(as.integer(n), graph, thresholds, beta, as.integer(nIter), as.integer(responses), as.logical(method == "exact"))                 
    
    if (any(is.na(Res)) & method == "exact")
    {
      warning("NA's detected, this means that no exact sample was drawn after 100 couplings to the past. Use higher nIter value or method='metropolis' for inexact sampling.")
    }
  } else 
  {
      Res <- IsingDir(n, graph, thresholds, beta, responses)
  }
  
  return(Res)
}

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

IsingSumLikelihood <- function(graph, thresholds, beta, responses = c(0L,1L))
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
  SumScores <- rowSums(1*(Allstates==1))
  stopifnot(require("plyr"))
  df <- ddply(data.frame(Sum = SumScores, P = P),.(Sum),summarize,P=sum(P))
  df$P <- df$P / sum(df$P)
  return(df)
}

IsingLikelihood <- function(graph, thresholds, beta, responses = c(0L,1L))
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
  df <- cbind(Probability = P / sum(P), Allstates)
  return(df)
}
  
IsingStateProb <- function(s,graph,thresh,beta,responses=c(0L,1L))
{
  if (!is.list(s)) s <- list(s)
  
  Allstates <- do.call(expand.grid,lapply(1:N,function(x)responses))
  Dist <- exp(- beta * apply(Allstates,1,function(s)H(graph,s,thresh)))  
  Z <- sum(Dist)  
  
  sapply(s, function(x)exp(-beta*H(graph,x,thresh))/Z)
}
*/

