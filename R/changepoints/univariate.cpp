#include <Rcpp.h>
#include <vector>
#include <cmath>     // for std::pow()
#include <algorithm> // for std::nth_element()
#include <string>
//#include <iostream>


/////////////////// theta-functions ///////////////

/* takes the mean over an interval [start, stop],
 The bounds start and stop assume zero-indexing */
double mean_ab_cpp (std::vector<double> x, int start, int stop)
{
  double total = 0.0;
  
  for(int i = start; i <= stop; i++)
  {
    total += x[i];
  }
  
  int n;
  if(start <= stop)
  {
    n = stop - start + 1;
  } else {
    
    // in case start > stop, return zero
    n = 1;
  }
  
  return total / n;
  
}


// calculate median of a vector x
// n > 0 is the length of x
double median(std::vector<double> x, int n)
{
  // declare final value
  double med;
  
  if(n % 2 == 0)
  {
    
    // case 1: n is even
    // in this case, the median is the mean of the n/2-th and (n/2+1)-th 
    // largest element in x
    int m1 = n / 2; 
    int m0 = m1 - 1; // account for zero-indexing
    
    // partially sort x in a way s.t. the m0-th index is the m0-th largest value
    std::nth_element(x.begin(), x.begin() + m0, x.end());
    double med0 = x[m0];
    
    // partially sort x in a way s.t. the m1-th index is the m1-th largest value
    std::nth_element(x.begin(), x.begin() + m1, x.end());
    double med1 = x[m1];
    
    // averaging yields the median
    med = 0.5 * (med0 + med1);
    
  } else
  {
    // case 2: n is odd
    // in this case, the median is the (n+1)/2-th largest element in x
    int m = (n+1)/2 - 1; // account for zero-indexing
    
    // partially sort x in a way s.t. the m-th index is the m-th largest value
    std::nth_element(x.begin(), x.begin() + m, x.end());
    
    // the median is the m-th largest value
    med = x[m];
  }
  
  return med;
}


/* takes the median over an interval [start, stop],
 The bounds start and stop assume zero-indexing */
double median_ab_cpp (std::vector<double> x, int start, int stop)
{
  // declare returned variable
  double out;
  
  if(start <= stop)
  {
    /// case 1: start index does not exceed stop index
    // length of subsetted vector
    int n = stop - start + 1;
    
    // declare and fill up subsetted vector y
    std::vector<double> y(n);
    int ct = 0;
    for(int i = start; i <= stop; i++)
    {
      y[ct] = x[i];
      ct += 1;
    }
    
    // calculate the median of y
    out = median(y, n);
    
  } else {
    /// case 2: start index exceeds stop index. Return 0 here.
    out = 0.0;
  }
  
  return out;
}



// calculate the statistic along which we expect break points
// 'theta_fun' is a pointer to either mean_ab_cpp() or median_ab_cpp()
// The bounds 'start' and 'stop' assume one-indexing
double theta_ab(std::vector<double> x, int start, int stop,
                double (*theta_fun)(std::vector<double>, int, int))
{
  // call the corresponding function. It expects 0-indexing in 
  // start and stop, so we account for this in the call below
  return theta_fun(x, start-1, stop-1);
}






//////////////// SNsingle functions ////////////



// calculate Dn(k) 
// NOTE: k assumes one-indexing and theta_fun is a pointer
// that points to the function along which we want to test for breaks
double Dn_k (std::vector<double> x, int n, int k,
             double (*theta_fun)(std::vector<double>, int, int))
{
    double theta0 = theta_ab(x, 1    , k, theta_fun); 
    double theta1 = theta_ab(x, k + 1, n, theta_fun);
    
    return (double)k * (n - k) * (theta0 - theta1) / std::pow(n, 1.5);
}

// calculate the left sum in V_n(k)
// NOTE: k assumes one-indexing and theta_fun is a pointer
// that points to the function along which we want to test for breaks
double Vn_k_left(std::vector<double> x, int n, int k,
                 double (*theta_fun)(std::vector<double>, int, int))
{
    double out = 0.0;

    for(int i = 1; i <= k; i++)
    {
        double theta0 = theta_ab(x, 1  , i, theta_fun);
        double theta1 = theta_ab(x, i+1, k, theta_fun);
        double dtheta = theta0 - theta1;
        double scale  = (double)i * (k-i);
        
        out += scale * scale * dtheta * dtheta;
    }
    
    double nk2 = (double)n * n * k * k; 
    return out / nk2;
}


// calculate the right sum in V_n(k)
// NOTE: k assumes one-indexing and theta_fun is a pointer
// that points to the function along which we want to test for breaks
double Vn_k_right(std::vector<double> x, int n, int k,
                  double (*theta_fun)(std::vector<double>, int, int))
{
    double out = 0.0;

    for(int i = k+1; i <= n; i++)
    {
        double theta0 = theta_ab(x, i  , n  , theta_fun);
        double theta1 = theta_ab(x, k+1, i-1, theta_fun);
        double dtheta = theta0 - theta1;
        double scale  = (double)(n-i+1) * (i-k-1);
        
        out += scale * scale * dtheta * dtheta;
    }
    
    double nk2 = (double)n * (n-k) * n * (n-k);
    return out / nk2;
}


// calculate V_n(k)
// NOTE: k assumes one-indexing and theta_fun is a pointer
// that points to the function along which we want to test for breaks
double Vn_k(std::vector<double> x, int n, int k,
            double (*theta_fun)(std::vector<double>, int, int))
{
    return Vn_k_left(x,n,k,theta_fun) + Vn_k_right(x,n,k,theta_fun);
}


// calculate T_n(k)
// NOTE: k assumes one-indexing and theta_fun is a pointer
// that points to the function along which we want to test for breaks
double Tn_k(std::vector<double> x, int n, int k,
            double (*theta_fun)(std::vector<double>, int, int))
{
    double D = Dn_k(x, n, k, theta_fun);
    double V = Vn_k(x, n, k, theta_fun);
    return D * D / V;
}


// calculate a vector holding T_n(k), for k = 1,...,n-1
// theta_fun is a pointer that points to the function 
// along which we want to test for breaks
std::vector<double> Tn_cpp (std::vector<double> x, int n,
                            double (*theta_fun)(std::vector<double>, int, int))
{
    int nm1 = n - 1;
    std::vector<double> out(nm1);

    for(int k = 1; k <= nm1; k++)
    {
        out[k-1] = Tn_k(x, n, k, theta_fun); // account for one-indexing in x
    }

    return out;
}


// Returns a vector holding T_n(k), for k = 1,...,n-1
// [[Rcpp::export]]
Rcpp::NumericVector Tn_univariate (Rcpp::NumericVector x, Rcpp::StringVector theta)
{
  // convert Rcpp functions STL functions, because those 
  // are expected by the declared functions 
  std::vector<double> y = Rcpp::as<std::vector<double> >(x);
  std::string theta_cpp = Rcpp::as<std::string >(theta);
  int n = y.size();
  
  
  // declare and define a function pointer that is intended 
  // to point to either 'mean_ab_cpp()' or 'median_ab_cpp()'
  double (*theta_ptr)(std::vector<double>, int, int);
  
  if(theta_cpp == "mean")
  {
    theta_ptr = &mean_ab_cpp;
  }
  else if(theta_cpp == "median")
  {
    theta_ptr = &median_ab_cpp;
    
  } else{
    Rcpp::stop("'theta' must be 'mean' or 'median'");
  }
  
  // pass to the function which expects STL input and returns STL vector
  std::vector<double> out = Tn_cpp(y, n, theta_ptr);
  
  // return R vector
  return Rcpp::wrap(out);
  
}


//////////////////////////////// unused functions ////////////////////
/*
// get index of maximum of a vector
int which_max(std::vector<double> x)
{      
  int idx = 0;
  int n = x.size();
  
  for (int i = 1; i < n; i++)
  {
    if (x[i] > x[idx])
    {
      idx = i;
    }
  }
  
  return idx;
}


// calculate elementwise cumulative sum of a vector. NOT USED
std::vector<double> cumsum (std::vector<double> x)
{
  int n = x.size();
  std::vector<double> out(n);
  double cum = 0.0;
  
  for(int i = 0; i < n; i++)
  {
    cum += x[i];
    out[i] = cum;
  }
  
  return out;
} */

/* Run SNsingle. 
// Returns index of the changepoint (one-indexed) that equals zero if no 
// index has been flagged
int SNsingle_cpp (std::vector<double> x, 
                  double Kn = 40.1)
{

    int n = x.size();
    int out;

    // calculate a vector holding T_n(k), for k = 1,...,n-1
    std::vector<double> T = Tn(x, n);

    // get the index k that maximizes Tn(k) (zero-indexed)
    int khat = which_max(T);

    if(T[khat] > Kn)
    {
        out = khat + 1; // make it one-indexed
    } else {
        out = 0;
    }

    return out;
}


// Returns index of one flagged changepoint (one-indexed and possibly empty)
// the output is a vector of length one (not int as it was easier to return
// an empty R integer objetc like this)
// [[Rcpp::export]]
Rcpp::IntegerVector SNsingle (Rcpp::NumericVector x, 
                              double Kn = 40.1)
{
  // convert Rcpp vector to STL vector, because SNsingle_cpp() expects 
  //   an object of the STL's vector class as input 
  std::vector<double> y = Rcpp::as<std::vector<double> >(x);
  
  // pass to the function which expects STL input
  int khat = SNsingle_cpp(y, Kn);
  
  Rcpp::IntegerVector out;
  
  if(khat > 0)
  {
    out = Rcpp::wrap(khat);
  }
  
  return out;
  
} */


/*
int main()
{

  std::vector<double> x = { 0, 0.1, -0.2, 0, 0.3, -0.2, 0.1, 1, 1.2, 1.1, 0.9, 1, 1.2, 1, 0.9 };
  int k = SNsingle_cpp(x);
  std::cout << k << ' ';

  return 0;
}
 */

