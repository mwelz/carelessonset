#' Detects at most one changepoint in a multidimensional time series 
#' 
#' The matrix x holds n time series observations (each of dimension d) in its columns, so x is (d x n)
#' 
#' @param x A matrix of size (d x n)
#' @param Kn a vector of critical values
#' @param theta statistic along which to look for changepoints. Either "mean" or "median"
#' 
#' @export
changepoint_multivariate <- function(x, Kn, theta = "mean") 
{
  # get vector of Tn(k), k=1,...,n-1 and SN test statistic
  t <- Tn_multivariate(x, theta)
  khat <- which.max(t)
  sn <- t[khat]
  flagged_ls <- rep(list(NA_integer_), length(Kn))

  # iterate over all critical values
  for(i in seq_along(Kn))
  {
    if(sn > Kn[i])
    {
      flagged <- khat
    } else{
      flagged <- integer()
    } # IF
    
    flagged_ls[[i]] <- flagged
    
  } # FOR
  
  # prepare output
  names(flagged_ls) <- Kn
  return(list(Tn = t, flagged = flagged_ls))
}
