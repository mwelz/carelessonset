## load all required functions
library("Rcpp")
Rcpp::sourceCpp("R/changepoints/univariate.cpp")
Rcpp::sourceCpp("R/changepoints/multivariate.cpp") # tons of compiler noise...
source("R/changepoints/univariate.R")
source("R/changepoints/multivariate.R")


#' Detect changepoints across multiple univariate series
#' 
#' @param data data matrix, each row of which represents a time series along which a changepoint is searched
#' @param alpha vector of significance levels
#' @param theta statistic along which to look for changepoints. Either "mean" or "median"
#' @param mc_cores number of cores to run in parallel (=1 for no parallelization). Note: Parallelization is currently only supported for Unix-like systems
#' @param matrix logical. If TRUE, then a binary matrix of the dimensionality of the input data is returned in which each 1 represents a changepoint
#' @param teststat Logical. Shall the test statistics for each period be returned?
#' 
#' @export
detect_changepoints_univariate <- function(data,
                                           alpha = c(0.005, 0.001),
                                           theta = "mean",
                                           mc_cores = parallel::detectCores(),
                                           matrix = FALSE,
                                           teststat = FALSE)
{
  ## input check 
  stopifnot(is.matrix(data))
  
  # critical values of limit distribution (Table 1 in Shao & Zhao, 2010, JASA)
  critvals <- cbind(
    level = c(0.100, 0.050, 0.025, 0.010, 0.005, 0.001),
    Kn = c(29.6, 40.1, 52.2, 68.6, 84.6, 121.9)
  )
  
  # are critical values for all passed levels available?
  stopifnot(all(alpha %in% critvals[,"level"]))
  
  ## order passed levels
  alpha_ord <- sort(alpha, decreasing = TRUE)
  
  # get critical values associated with passed levels, in increasing magnitude so that the least conservative critical value comes first
  # higher critical values are associated with more conservative levels
  Kn <- unname(critvals[critvals[,"level"] %in% alpha_ord, "Kn"])
  
  # get number of series, grid size, and dimensionality of data series
  n         <- nrow(data)
  d         <- ncol(data)
  
  
  # get flagged checkpoints: note that Kn is an increasing sequence to respect the order of the sorted alpha vector
  if(.Platform$OS.type == "unix" && mc_cores > 1L)
  {
    cp_ls <- 
      parallel::mclapply(mc.cores = mc_cores, 
                         X = seq_len(n),  FUN = function(i){
                           changepoint_univariate(data[i,], Kn, theta)
                         })
  } else{
    cp_ls <- 
      lapply(X = seq_len(n), function(i){
        changepoint_univariate(data[i,], Kn, theta)
      })
  } # IF
  
  
  if(matrix)
  {
    # get it in right shape
    CP <- lapply(seq_along(Kn), function(...) matrix(0L, nrow = n, ncol = d))
    names(CP) <- alpha_ord
    
    for(i in seq_len(n))
    {
      for(k in seq_along(Kn))
      {
        flagged <- cp_ls[[i]]$flagged[[k]]
        if(length(flagged) > 0){
          CP[[k]][i, flagged] <- 1L
        } # IF flagged
      }
    } # FOR
  } else
  {
    # get it in right shape
    CP <- lapply(seq_along(Kn), function(...) rep(NA_integer_, n))
    names(CP) <- alpha_ord
    
    for(i in seq_len(n))
    {
      for(k in seq_along(Kn))
      {
        flagged <- cp_ls[[i]]$flagged[[k]]
        if(length(flagged) > 0){
          CP[[k]][i] <- flagged
        } # IF flagged
      }
    } # FOR
  } # IF matrix
  
  if(teststat)
  {
    T. <- t(sapply(seq_len(n), function(i) cp_ls[[i]]$Tn))
    
    out <- list(changepoints = CP, 
                SNCP = list(Tn = T., Tmax = apply(T., 1L, max)))
  } else
  {
    out <- CP
  }
 
  return(out)
  
} # FUN


#' #' Detect changepoints across multiple multivariate series
#' 
#' @param data A list that contains matrices, each of which holds n time series observations (each of dimension d) in its columns, so each matrix is (d x n). We search for changepoints along these data
#' @param alpha vector of significance levels
#' @param theta statistic along which to look for changepoints. Either "mean" or "median"
#' @param mc_cores number of cores to run in parallel (=1 for no parallelization). Note: Parallelization is currently only supported for Unix-like systems
#' @param matrix logical. If TRUE, then a binary matrix of the dimensionality of the input data is returned in which each 1 represents a changepoint
#' @param teststat Logical. Shall the test statistics for each period be returned?
#' 
#' @export
detect_changepoints_multivariate <- function(data,
                                             alpha =  c(0.005, 0.001),
                                             theta = "mean",
                                             mc_cores = parallel::detectCores(),
                                             matrix = FALSE,
                                             teststat = FALSE)
{
  stopifnot(is.list(data))
  
  # get number of series, grid size, and dimensionality of series
  num_respondents        <- length(data)
  num_items              <- ncol(data[[1]])
  dimension              <- nrow(data[[1]])
  
  # critical values of limit distribution (Table 1 in Shao & Zhao, 2010, JASA)
  critvals <- 
    cbind(
      level = c(0.100, 0.050, 0.025, 0.010, 0.005, 0.001),
      Kn2 = c(56.5, 73.7, 92.2, 117.7, 135.3, 192.5),
      Kn3 = c(81.5, 103.6, 128.9, 160.0, 182.9, 246.8),
      Kn4 = c(114.7, 141.5, 171.9, 209.7, 246.6, 319.2)
    )
  
  # are critical values for all passed levels available?
  stopifnot(all(alpha %in% critvals[,"level"]))
  stopifnot(dimension > 1 & dimension < 5)
  
  # order passed levels
  alpha_ord <- sort(alpha, decreasing = TRUE)
  
  # get critical values associated with passed levels, in increasing magnitude so that the least conservative critical value comes first
  # higher critical values are associated with more conservative levels
  Kn <- unname(critvals[critvals[,"level"] %in% alpha_ord, dimension])
  
  # get flagged checkpoints: note that Kn is an increasing sequence to respect the order of the sorted alpha vector
  if(.Platform$OS.type == "unix" && mc_cores > 1L)
  {
    cp_ls <- 
      parallel::mclapply(mc.cores = mc_cores, 
                         X = seq_len(num_respondents),
                         FUN = function(i){
                           changepoint_multivariate(data[[i]], Kn, theta)
                         })
  } else{
    cp_ls <- 
      lapply(X = seq_len(num_respondents), 
             FUN = function(i){
               changepoint_multivariate(data[[i]], Kn, theta)
             })
  } # IF
  
  
  # get it in right shape
  CP <- lapply(seq_along(Kn), function(...) matrix(0L, 
                                                   nrow = num_respondents, 
                                                   ncol = num_items))
  names(CP) <- alpha_ord
  
  if(matrix)
  {
    # get it in right shape
    CP <- lapply(seq_along(Kn), function(...) matrix(0L, nrow = num_respondents, ncol = num_items))
    names(CP) <- alpha_ord
    
    for(i in seq_len(num_respondents))
    {
      for(k in seq_along(Kn))
      {
        flagged <- cp_ls[[i]]$flagged[[k]]
        if(length(flagged) > 0){
          CP[[k]][i, flagged] <- 1L
        } # IF flagged
      }
    } # FOR
  } else
  {
    # get it in right shape
    CP <- lapply(seq_along(Kn), function(...) rep(NA_integer_, num_respondents))
    names(CP) <- alpha_ord
    
    for(i in seq_len(num_respondents))
    {
      for(k in seq_along(Kn))
      {
        flagged <- cp_ls[[i]]$flagged[[k]]
        if(length(flagged) > 0){
          CP[[k]][i] <- flagged
        } # IF flagged
      }
    } # FOR
  } # IF matrix
  
  if(teststat)
  {
    T. <- t(sapply(seq_len(num_respondents), function(i) cp_ls[[i]]$Tn))

    out <- list(changepoints = CP, 
                SNCP = list(Tn = T., Tmax = apply(T., 1L, max)))
  } else
  {
    out <- CP
  }
  
  return(out)
  
} # FUN
