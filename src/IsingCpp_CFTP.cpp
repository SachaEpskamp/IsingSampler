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
NumericVector PplusMinMax(int i, NumericMatrix J, NumericVector s, NumericVector h, double beta, NumericVector responses)
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
      if (!ISNA(s[j]))
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
NumericVector IsingEx(NumericMatrix graph, NumericVector thresholds, double beta, int nIter, NumericVector responses, bool exact,
NumericVector constrain)
{
  // Parameters and results vector:
  int N = graph.nrow();
  NumericVector state(N, NA_REAL);
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
        state[i] = NA_REAL;
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
            state[node] = NA_REAL;
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
        if (ISNA(state[i]))
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
// Computes the full-conditional categorical distribution of node i over all
// response options, given the (complete) states of all other nodes.
// Returns a probability vector of length responses.size().
// For two response options this reduces exactly to the binary Ising conditional.
NumericVector Pcat(int i, NumericMatrix J, NumericVector s, NumericVector h, double beta, NumericVector responses)
{
  int N = J.nrow();
  int K = responses.size();

  // Local field acting on node i: F = h_i + sum_{j != i} J_ij s_j.
  // The conditional energy of response r is r * F, so the conditional only
  // depends on this single scalar regardless of the number of responses.
  double F = h[i];
  for (int j=0; j<N; j++)
  {
    if (i != j)
    {
      F += J(i,j) * s[j];
    }
  }

  // Softmax over beta * responses[k] * F, shifted by the maximum for numerical
  // stability (prevents overflow when beta or the field is large).
  NumericVector P(K);
  double maxH = beta * responses[0] * F;
  for (int k=1; k<K; k++)
  {
    double Hk = beta * responses[k] * F;
    if (Hk > maxH) maxH = Hk;
  }
  double sum = 0;
  for (int k=0; k<K; k++)
  {
    P[k] = exp(beta * responses[k] * F - maxH);
    sum += P[k];
  }
  for (int k=0; k<K; k++) P[k] /= sum;

  return(P);
}

// Draws a response index from categorical probabilities P given uniform u.
// Iterates from the last response downward so that, for two response options,
// the mapping is identical to the original binary sampler (u < P[1] -> response 1).
int drawResponse(NumericVector P, double u)
{
  int K = P.size();
  double cum = 0;
  for (int k=K-1; k>0; k--)
  {
    cum += P[k];
    if (u < cum) return(k);
  }
  return(0);
}

// Initializes a random state with each node drawn uniformly from all response
// options. For two response options this reproduces the original binary
// initialization (uniform over the two responses); for more options it avoids
// biasing the starting state toward the first two responses, which otherwise
// skews the realized mode at high beta when the chain mixes slowly.
NumericVector randomState(int N, NumericVector responses)
{
  int K = responses.size();
  NumericVector state(N);
  if (K == 2)
  {
    state = ifelse(runif(N) < 0.5, responses[1], responses[0]);
  } else
  {
    for (int i=0; i<N; i++) state[i] = responses[(int) floor(R::runif(0, K))];
  }
  return(state);
}


NumericVector IsingMet(NumericMatrix graph, NumericVector thresholds, double beta, int nIter, NumericVector responses,
NumericVector constrain)
{
  // Parameters and results vector:
  int N = graph.nrow();
  NumericVector state = randomState(N, responses);
  for (int i=0; i<N; i++)
  {
    if (!ISNA(constrain[i]))
    {
      state[i] = constrain[i];
    }
  }
  double u;
  NumericVector P;

    // START ALGORITHM
    for (int it=0;it<nIter;it++)
    {
      for (int node=0;node<N;node++)
      {
        if (ISNA(constrain[node]))
        {
         u = runif(1)[0];
         P = Pcat(node, graph, state, thresholds, beta, responses);
         state[node] = responses[drawResponse(P, u)];
        }
      }
    }

  return(state);
}


///ISING PROCESS SAMPLER:
// [[Rcpp::export]]
NumericMatrix IsingProcess(int nSample, NumericMatrix graph, NumericVector thresholds, double beta, NumericVector responses)
{
  // Parameters and results vector:
  int N = graph.nrow();
  NumericVector state = randomState(N, responses);
  double u;
  NumericVector P;
  NumericMatrix Res(nSample,N);
  int node;

    // START ALGORITHM
    for (int it=0;it<nSample;it++)
    {
      node = floor(R::runif(0,N));
        u = runif(1)[0];
        P = Pcat(node, graph, state, thresholds, beta, responses);
        state[node] = responses[drawResponse(P, u)];
        for (int k=0; k<N; k++) Res(it,k) = state[k];
    }
   
  return(Res);
}

// OVERAL FUNCTION //
// [[Rcpp::export]]
NumericMatrix IsingSamplerCpp(int n, NumericMatrix graph, NumericVector thresholds, double beta, int nIter, NumericVector responses, bool exact,
NumericMatrix constrain)
{
  int Ni = graph.nrow();
  NumericMatrix Res(n,Ni);
  NumericVector state(Ni);
  NumericVector constrainVec(Ni);
  if (exact)
  {
    for (int s=0;s<n;s++)
    {
      for (int i=0;i<Ni;i++) constrainVec[i] = constrain(s,i);
      state = IsingEx(graph, thresholds, beta, nIter, responses, exact, constrainVec);
      for (int i=0;i<Ni;i++) Res(s,i) = state[i];
    }
  } else 
  {
    for (int s=0;s<n;s++)
    {
      for (int i=0;i<Ni;i++) constrainVec[i] = constrain(s,i);
      state = IsingMet(graph, thresholds, beta, nIter, responses, constrainVec);
      for (int i=0;i<Ni;i++) Res(s,i) = state[i];
    }
  }
  
  return(Res);
}


// HELPER FUNCTIONS //
// Hamiltonian:
// [[Rcpp::export]]
double H(NumericMatrix J, NumericVector s, NumericVector h)
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


// Likelihood without Z
// [[Rcpp::export]]
double f(NumericMatrix Y, NumericMatrix J, NumericVector h)
{
  double Res = 1;
  int Np = Y.nrow();
  int Ni = J.ncol();
  NumericVector s(Ni);
  for (int p=0;p<Np;p++)
  {
    for (int i=0;i<Ni;i++)
    {
      s[i] = Y(p,i);
    }
    Res *= exp(-1.0 * H(J, s, h));
  }
  return(Res);
}


// VECTOR VERSIONS //


// Hamiltonian:
// [[Rcpp::export]]
double Hvec(IntegerVector s, NumericVector Theta, int N)
{
  double Res = 0;
  int c=0;
  for (int i=0;i<N;i++)
  {
    Res -= Theta[c] * s[i];
    c++;
  }
  for (int i=0;i<N;i++)
  {
    for (int j=i; j<N;j++)
    {
     if (j!=i) 
     {
       Res -= Theta[c] * s[i] * s[j]; 
       c++;
     }
    }
  }
  return(Res);
}

// Likelihood without Z
double fvec(IntegerMatrix Y, NumericVector Theta)
{
  double Res = 1;
  int Np = Y.nrow();
  int Ni = Y.ncol();
  IntegerVector s(Ni);
  for (int p=0;p<Np;p++)
  {
    for (int i=0;i<Ni;i++)
    {
      s[i] = Y(p,i);
    }
    Res *= exp(-1.0 * Hvec(s, Theta, Ni));
  }
  return(Res);
}

// Log-Likelihood without Z
double fveclog(IntegerMatrix Y, NumericVector Theta)
{
  double Res = 1;
  int Np = Y.nrow();
  int Ni = Y.ncol();
  IntegerVector s(Ni);
  for (int p=0;p<Np;p++)
  {
    for (int i=0;i<Ni;i++)
    {
      s[i] = Y(p,i);
    }
    Res -= Hvec(s, Theta, Ni);
  }
  return(Res);
}
//
//IntegerMatrix vecSampler(int n, int N, NumericVector Theta, int nIter, IntegerVector responses)
//{
//   NumericVector thresh(N);
//   for (int i=0; i<N; i++)
//   {
//     thresh[i] = Theta[i];
//   }
//   
//   NumericMatrix graph(N,N);
//   int c=N+1;
//   for (int i=0;i<N;i++)
//   {
//    for (int j=i; j<N;j++)
//    {
//     if (j!=i) 
//     {
//       graph(i,j) = Theta[c];
//      graph(j,i) = Theta[c];
//       c++;
//     }
//    }
//  }
//   
//   return(IsingSamplerCpp(n, graph, thresh, 1.0, nIter, responses, true));
//}
//
//// Uniform distribution (prior):
//double FakeUnif(NumericVector x, double lower, double upper)
//{
//  double Res = 1.0;
//  
//  for (int i=0; i < x.length(); i++)
//  {
//    if (x[i] < lower || x[i] > upper)
//    {
//      Res = 0.0;
//      break;
//    }
//  }
//  
//  return(Res);
//}


// Progress bar function:
/*
int progress_bar(double x, double N)
{
    // how wide you want the progress meter to be
    int totaldotz=40;
    double fraction = x / N;
    // part of the progressmeter that's already "full"
    int dotz = round(fraction * totaldotz);

    // create the "meter"
    int ii=0;
    printf("%3.0f%% [",fraction*100);
    // part  that's full already
    for ( ; ii < dotz;ii++) {
        printf("=");
    }
    // remaining part (spaces)
    for ( ; ii < totaldotz;ii++) {
        printf(" ");
    }
    // and back to line begin - do not forget the fflush to avoid
    // output buffering problems!
    printf("]\r");
    fflush(stdout);
}
*/

//
//// EXCHANGE ALGORTIHM //
//// [[Rcpp::export]]
//NumericMatrix ExchangeAlgo(IntegerMatrix Y, double lowerBound, double upperBound, double stepSize, int nIter, IntegerVector responses,
//    bool simAn, double tempStart, double tempEnd, NumericVector StartValues)
//{
//  int Np = Y.nrow();
//  int Ni = Y.ncol();
//  
//  // Number of parameters:
//  int Npar = Ni + (Ni*(Ni-1))/2;
//  
//  // Fantasy matrix:
//  IntegerMatrix X(Np,Ni);
//  
//  // Results matrix:
//  NumericMatrix Samples(nIter, Npar);
//  
//  // Current parameter values:
//  // NumericVector curPars = runif(Npar, lowerBound, upperBound);
//  NumericVector curPars(Npar, 0.0);
//  for (int i=0; i<Npar; i++) curPars[i] = StartValues[i];
//  NumericVector propPars(Npar);
//  
//  double a;
//  double r;
//  
//  // START ITERATING //
//  for (int it=0; it<nIter; it++)
//  {
//   // progress_bar((double)it, (double)nIter);
//    // For each parameter:
//    for (int n=0;n<Npar;n++)
//    {
//      // Propose new state:
//      for (int i=0; i<Npar; i++)
//      {
//        if (i==n)
//        {
//          propPars[i] = curPars[i] + R::rnorm(0.0,stepSize);
//        } else 
//        {
//          propPars[i] = curPars[i];
//        }
//      }
//      
//      // Simulate data with new state:
//      X =  vecSampler(Np, Ni, propPars, nIter, responses);
//      
//      // Random number:
//      r = R::runif(0,1);
//      
//      // Acceptance probability:
//      a = FakeUnif(propPars,lowerBound,upperBound)/FakeUnif(curPars,lowerBound,upperBound) * 
//           exp(fveclog(Y, propPars) + fveclog(X, curPars) - fveclog(Y,curPars) - fvec(X,propPars));
//      
//      if (!simAn)
//      {
//        if (r < a)
//        {
//          curPars[n] = propPars[n];
//        }  
//      } else {
//        if (r < exp(log(a)/ (tempStart - it * (tempStart-tempEnd)/nIter)))
//        {
//          curPars[n] = propPars[n];
//        }  
//      }
//      
//      
//      Samples(it, n) = curPars[n];
//    }
//  }
//  
//  
//  return(Samples);
//}
 
///// Broderick et al 2013:



// Function to compute expected values:
// [[Rcpp::export]]
NumericVector expvalues(NumericMatrix x){
  // Sample size:
  int N = x.nrow();
  // Number of nodes:
  int P = x.ncol();
  int nPar = P + P*(P-1)/2;
  // Results vector:
  NumericVector Res(nPar, 0.0);
  
  int par = 0;
  
  // Fill
  while (par < nPar){
    
    // Means:
    for (int j=0; j<P;j++){
      double mean = 0;
      for (int k=0; k<N; k++){
        mean += x(k,j) ;
      }
      Res[par] = mean / N ;
      par++;
    }
    
    // Covariances:
    for (int i=0; i<P; i++){
      for (int j=i; j<P;j++){
        if (i != j){
         double squared = 0;
          for (int k=0; k<N; k++){
            squared += x(k,i) * x(k,j) ;
          }
          Res[par] = squared / N;
          par++; 
        }
      }
    }
  }
  
  return(Res);
}

// Function to obtain thresholds from vector:
// [[Rcpp::export]]
NumericVector vec2Thresh(NumericVector vec, int P){
  NumericVector Res(P);
  
  for (int i=0; i<P; i++){
    Res[i] = vec[i];
  }
  
  return(Res);
}

// [[Rcpp::export]]
NumericMatrix vec2Graph(NumericVector vec, int P){
  NumericMatrix Res(P, P);
  
  int par = P;
  
  for (int i=0; i<P; i++){
      for (int j=i; j<P; j++){
        if (i != j){
           Res(i,j) = Res(j,i) = vec[par];   
           par++;
        }
      }
  }
  
  return(Res);
}

// Main optimisation function:
// [[Rcpp::export]]
NumericVector Broderick2013(
  NumericMatrix x, // Data matrix
  int M, // Number of samples  to draw
  int T, // Number of iterations
  int nIter, // Temporary: number of sequences, replace with convergence test
  NumericVector responses
  )
{
  // Number of nodes:
  int P = x.ncol();
  // Number of parameters:
  int nPar = P + P*(P-1)/2;
  // Current estimtes and new estimates:
  NumericVector curEsts(nPar, 0.0);
  NumericVector newEsts(nPar, 0.0);
  
    // Dummy constraints mat (ugly, should be removed):
  NumericMatrix cons(M, P);
  std::fill(cons.begin(), cons.end(), NA_REAL);

  // Observed statistics:
  NumericVector obsStats = expvalues(x);
  
  // Thresholds to mimic margins:
  for (int i=0; i<P; i++){
    curEsts[i] = (responses[1] - responses[0]) *log(obsStats[i]);
  }

  double step = 1;

  // Start iterating:
  for (int s=0; s<nIter; s++){
    // Set new ests:
    for (int i=0;i<nPar;i++){
      curEsts[i] = newEsts[i];
    }
    
    // Generate monte carlo samples:
    NumericMatrix Samples =  IsingSamplerCpp(M, vec2Graph(curEsts, P), vec2Thresh(curEsts, P), 1.0, 1000, responses, false,cons);
    
    // Statistics:
    NumericVector sampStats = expvalues(Samples);
    
    for (int t=0; t<T; t++){
      // For each statistic, move up if expected value too low, down if expected value too high:
      for (int par=0;par<nPar;par++){
        // Estimate sampStat:
        double stat = (sampStats[par] * exp(-(newEsts[par] - curEsts[par]) * sampStats[par]) ) / (exp(-(newEsts[par] - curEsts[par]) * sampStats[par]) );
        
        
        // too high:
        if (stat > obsStats[par]){
          newEsts[par] -= step;
        } else {
          // Too low:
          newEsts[par] += step;
        }
      }
    }

    step *= 0.5;
  }



  return(newEsts);
}
