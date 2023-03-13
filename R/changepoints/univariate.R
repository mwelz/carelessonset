#' Detects at most one changepoint in a univariate time series 
#' 
#' @param x A vector of time series observations
#' @param Kn a vector of critical values
#' @param theta statistic along which to look for changepoints. Either "mean" or "median"
#' 
#' @export
changepoint_univariate <- function(x, Kn = 40.1, theta = "mean")
{

  # get vector of Tn(k), k=1,...,n-1 and SN test statistic
  t <- Tn_univariate(x, theta)
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
  
} # FUN

