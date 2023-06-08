
get_onset <- function(x, alpha = 0.001)
{
  ## input checks
  stopifnot(inherits(x, "carelessonset"))
  stopifnot(length(alpha) == 1)
  alphas <- x$alpha
  stopifnot(alpha %in% alphas)
  
  # get flagged changepoints
  cp <- x$changepoints[[as.character(alpha)]]
  
  ## index of respondents with a flagged changepoint at given level
  flagged <- which(!is.na(cp)) 
  
  ## onset of carelessness for flagged samples
  location <- cp[flagged]
  
  # prepare returns
  out <- cbind(flagged, location)
  colnames(out) <- c("idx flagged", "onset")
  return(out)
} # FUN


plot.carelessonset <- function(x, idx = 1, alpha = 0.001, ...)
{
  
  cp <- get_onset(x = x, alpha = alpha)
  idx_cp <- which(cp[,1] == idx)
  
  if(length(idx_cp) == 0)
  {
    onset <- NULL
  } else{
    onset <- cp[idx_cp, "onset"]
  } # IF
  
  ## make the plot
  pl <- plot_dimension(RE = x$series$RE[idx,],
                       ALSP = x$series$ALSP[idx,], 
                       time = x$series$time[idx,], 
                       onset = onset, 
                       color = "red")
  
  return(pl)
  
} # FUN
