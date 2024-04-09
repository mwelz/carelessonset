#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
Rcpp::IntegerVector LSP(Rcpp::IntegerVector x, int J = 1)
{
  // initialize
  int n = x.size();
  int streak = 1;
  int start_idx = 0;
  Rcpp::IntegerVector out(n);
  
  // calculate streak of consecutive items 
  for(int k = 0; k < n - J; ++k)
  {
    if(x[k] == x[k+J])
    {
      streak += 1; 
    } else
    {
      Rcpp::IntegerVector fill(k - start_idx + 1, streak);
      out[Rcpp::Range(start_idx, k)] = fill;
      start_idx = k + 1;
      streak = 1;
    }
  }
  
  // now, fill in the last J indices (there are 3 edge cases)
  
  if(streak >= J)
  {
    // case 1: current streak is ongoing and we have seen it completely at least once, so
    // assume it will last until the end and complete the current streak (which has not
    // yet been broken)
    Rcpp::IntegerVector fill(n - start_idx, streak + J - 1); 
    out[Rcpp::Range(start_idx, n-1)] = fill;
  } else if (1 < streak && streak < J)
  {
    // case 2: there is an ongoing streak, but we have no chance of seeing it completely
    // (since streak < J), so we do not know how the streak looks like. In this case,
    // fill in the current value of "streak" until the (n-J)-th item and fill in 1
    // for the last J items from (n-J+1)-th item onward (because we do not know how the 
    // complete streak looks like; we assume one-indexing throughout this comment)
    Rcpp::IntegerVector fill_streak(n-J-start_idx, streak-1); // decrement streak to ensure correct streak length 
    Rcpp::IntegerVector fill_ones(J, 1);
    out[Rcpp::Range(start_idx, n-J-1)] = fill_streak;
    out[Rcpp::Range(n-J, n-1)]         = fill_ones;
  } else
  {
    // case 3: there is no ongoing streak (streak = 1), so fill in 1 for the last J indices
    Rcpp::IntegerVector fill(J, 1); 
    out[Rcpp::Range(n-J, n-1)] = fill;
  }
  
  return out;
}
