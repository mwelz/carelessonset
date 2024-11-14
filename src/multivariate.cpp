#include <vector>
#include <cmath>     // for std::pow()
#include <algorithm> // for std::nth_element()
#include <string>
//#include <eigen3/Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;


////////////////// theta-functions ////////////////////////

/* takes the rowwise-mean of a matrix x over an interval [start, stop],
In our context, the matrix x holds n time series observations (each of dimension d) 
in its columns, so x is d x n.
NOTE: The bounds start and stop assume zero-indexing */

VectorXd mean_ab_cpp (MatrixXd x, int start, int stop, int d)
{
  // declare
  VectorXd out;

  if(start <= stop)
  {
    // length of the observation period
    int len = stop - start + 1;

    // subet x to be [d x len] 
    MatrixXd reduced = x.block(0,start,d,len);

    // return d-dimensional mean vector
    out = reduced.rowwise().mean();

  } else
  {
    out = VectorXd::Zero(d);
  }

  return out;
  
}



/*
// calculate median of a univariate vector x (NB: this is copy-pasted from SNsingle.cpp)
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



/* takes the rowwise-median (coordinatewise median) of a matrix x over an interval [start, stop],
In our context, the matrix x holds n time series observations (each of dimension d) 
in its columns, so x is d x n. 
// NOTE: The bounds start and stop assume zero-indexing 
VectorXd median_ab_cpp (MatrixXd x, int start, int stop, int d)
{
  // initialize
  VectorXd out = VectorXd::Zero(d);

  if(start <= stop)
  {
    // length of the observation period
    int len = stop - start + 1;

    // subet x to be [d x len] 
    MatrixXd reduced = x.block(0,start,d,len);

    // calculate the median along every dimension
    for(int j = 0; j < d; j++)
    {
      VectorXd xj = reduced.row(j);

      // copy the data to a new STL vector called y
      std::vector<double> y(len);
      VectorXd::Map(&y[0], len) = xj;

      // get median along dimension j
      out[j] = median(y, len);
    }

  } else
  {
    // do nothing s.t. a zero vector will be returned
  }

  return out;
  
}
*/


// calculate the statistic along which we expect break points
// 'theta_fun' is a pointer to either mean_ab_cpp() or median_ab_cpp()
// The bounds 'start' and 'stop' assume one-indexing
VectorXd theta_ab(MatrixXd x, int start, int stop, int d,
                  VectorXd (*theta_fun)(MatrixXd, int, int, int))
{
  // call the corresponding function. It expects 0-indexing in 
  // start and stop, so we account for this in the call below
  return theta_fun(x, start-1, stop-1, d);
}



//////////////// SNsingle functions ////////////


// calculate Dn(k) 
// NOTE: k assumes one-indexing and theta_fun is a pointer
// that points to the function along which we want to test for breaks
VectorXd Dn_k (MatrixXd x, int k, int d, int n,
               VectorXd (*theta_fun)(MatrixXd, int, int, int))
{
  VectorXd theta0 = theta_ab(x, 1    , k, d, theta_fun); 
  VectorXd theta1 = theta_ab(x, k + 1, n, d, theta_fun);
    
  return (double)k * (n-k) * (theta0 - theta1) / std::pow(n, 1.5);
}


// calculate the left sum in V_n(k)
// NOTE: k assumes one-indexing and theta_fun is a pointer
// that points to the function along which we want to test for breaks
MatrixXd Vn_k_left(MatrixXd x, int k, int d, int n,
                   VectorXd (*theta_fun)(MatrixXd, int, int, int))
{
  MatrixXd out = MatrixXd::Zero(d, d);

  for(int i = 1; i <= k; i++)
  {
    VectorXd theta0 = theta_ab(x, 1  , i, d, theta_fun);
    VectorXd theta1 = theta_ab(x, i+1, k, d, theta_fun);
    VectorXd dtheta = theta0 - theta1;
    double scale    = (double)i * (k-i);

    out +=  scale * scale * dtheta * dtheta.transpose();
  }

  double nk2 = (double)n * n * k * k; // TODO: can be taken out of function (same with right Vn)
  return out / nk2;
}



// calculate the right sum in V_n(k)
// NOTE: k assumes one-indexing and theta_fun is a pointer
// that points to the function along which we want to test for breaks
MatrixXd Vn_k_right(MatrixXd x, int k, int d, int n,
                    VectorXd (*theta_fun)(MatrixXd, int, int, int))
{
  MatrixXd out = MatrixXd::Zero(d, d);

  for(int i = k+1; i <= n; i++)
  {
    VectorXd theta0 = theta_ab(x, i  , n  , d, theta_fun);
    VectorXd theta1 = theta_ab(x, k+1, i-1, d, theta_fun);
    VectorXd dtheta = theta0 - theta1;
    double scale    = (double)(n-i+1) * (i-k-1);

    out +=  scale * scale * dtheta * dtheta.transpose();
    }

    double nk2 = (double)n * (n-k) * n * (n-k);
    return out / nk2;
}

// calculate Vn_(k)
// NOTE: k assumes one-indexing and theta_fun is a pointer
// that points to the function along which we want to test for breaks
MatrixXd Vn_k(MatrixXd x, int k, int d, int n,
              VectorXd (*theta_fun)(MatrixXd, int, int, int))
{
  return Vn_k_left(x,k,d,n,theta_fun) + Vn_k_right(x,k,d,n,theta_fun);
}


// calculate T_n(k)
// NOTE: k assumes one-indexing and theta_fun is a pointer
// that points to the function along which we want to test for breaks
double Tn_k(MatrixXd x, int k, int d, int n,
            VectorXd (*theta_fun)(MatrixXd, int, int, int))
{
  VectorXd D   = Dn_k(x,k,d,n, theta_fun);
  MatrixXd V   = Vn_k(x,k,d,n, theta_fun);
  MatrixXd out = D.transpose() * V.inverse() * D; // 1x1 matrix
  return out(0,0); 
}


// calculate a vector holding T_n(k), for k = 1,...,n-1
// theta_fun is a pointer that points to the
// function along which we want to test for breaks
std::vector<double> Tn_cpp (MatrixXd x, int d, int n,
                            VectorXd (*theta_fun)(MatrixXd, int, int, int))
{
    int nm1 = n - 1;
    std::vector<double> out(nm1);

    for(int k = 1; k <= nm1; k++)
    {
        out[k-1] = Tn_k(x, k, d, n, theta_fun); // account for one-indexing in x
    }

    return out;
}


// Returns a vector holding T_n(k), for k = 1,...,n-1.
// In our context, the matrix x holds n time series observations (each of dimension d) 
// in its columns, so x is (d x n).
// [[Rcpp::export]]
Rcpp::NumericVector Tn_multivariate (Rcpp::NumericMatrix x, Rcpp::StringVector theta)
{
  // convert x to a Eigen::MatrixXd object called y
  const Map<MatrixXd> y(Rcpp::as<Map<MatrixXd> >(x));
  
  // convert theta to an STL string called theta_cpp
  std::string theta_cpp = Rcpp::as<std::string >(theta);
  
  
  // declare and define a function pointer that is intended 
  // to point to either 'mean_ab_cpp()' or 'median_ab_cpp()'
  VectorXd (*theta_ptr)(MatrixXd, int, int, int);
  
  if(theta_cpp == "mean")
  {
    theta_ptr = &mean_ab_cpp;
  }
  else if(theta_cpp == "median")
  {
    Rcpp::stop("median not yet implemented");
    
  } else{
    Rcpp::stop("'theta' must be 'mean' or 'median'");
  }
  
  
  // call the function that expects Eigen input
  int d = y.rows();
  int n = y.cols();
  std::vector<double> out = Tn_cpp(y, d, n, theta_ptr);
    
  // return R vector
  return Rcpp::wrap(out);
}

/*
int main()
{
  int d = 3; int n = 10; int k = 4; int i=2;
  //MatrixXd x = MatrixXd::Random(d,d);
  //VectorXd v(3);
  //v << 1, 2, 3;
  //double vv = v.transpose() * x.inverse() * v;
  //std::cout << vv << std::endl;

  MatrixXd x = MatrixXd::Random(d,n);
  std::vector<double> out0 = Tn_cpp(x,d,n);
  for(int i = 0; i < 2; i++) {
  std::cout << out0[i] << ' ';
}
  //std::cout << out0 << std::endl;
  //std::cout << mean_ab(x, i, k-1, d) << std::endl;

VectorXd v = VectorXd::Ones(2);
MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);

MatrixXd y = v.transpose() * m * v; // works: y is a 1x1 matrix

std::vector<double> out(2, 0.0);
out[0] = y(0,0);
for(int i = 0; i < 2; i++) {
  std::cout << out[i] << ' ';
}

//std::cout << y(0,0) << std::endl;


}*/
