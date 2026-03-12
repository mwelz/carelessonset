## load all required functions
# library("Rcpp")
# Rcpp::sourceCpp("R/changepoints/univariate.cpp")
# Rcpp::sourceCpp("R/changepoints/multivariate.cpp") # tons of compiler noise...
# source("R/changepoints/univariate.R")
# source("R/changepoints/multivariate.R")


#' Detect changepoints across multiple univariate series
#' 
#' @param data data matrix, each row of which represents a time series along which a changepoint is searched
#' @param alpha vector of significance levels
#' @param theta statistic along which to look for changepoints. Either "mean" or "median"
#' @param mc_cores number of cores to run in parallel (=1 for no parallelization). Note: Parallelization is currently only supported for Unix-like systems
#' @param matrix logical. If TRUE, then a binary matrix of the dimensionality of the input data is returned in which each 1 represents a changepoint
#' @param teststat Logical. Shall the test statistics for each period be returned?
#' @param CP_at_segment_end Logical. Should the changepoint be the last period of an ending segment (TRUE, default, as in Shao et al.) or the first period of a beginning segment (FALSE)?
#' 
#' @export
changepoints_univariate <- function(data,
                                    alpha = c(0.005, 0.001),
                                    theta = "mean",
                                    mc_cores = parallel::detectCores(),
                                    matrix = FALSE,
                                    teststat = FALSE,
                                    CP_at_segment_end = TRUE)
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
                           changepoint_univariate_singleseries(data[i,], Kn, theta)
                         })
  } else{
    cp_ls <- 
      lapply(X = seq_len(n), function(i){
        changepoint_univariate_singleseries(data[i,], Kn, theta)
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
  
  ## if requested, convert the changepoint location to start of segment
  CP <- changepoint_conversion(
    CP, matrix = matrix, CP_at_segment_end = CP_at_segment_end)
  
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


#' Detect changepoints across multiple multivariate series
#' 
#' @param data A list that contains matrices, each of which holds n time series observations (each of dimension d) in its columns, so each matrix is (d x n). We search for changepoints along these data
#' @param alpha vector of significance levels
#' @param theta statistic along which to look for changepoints. Either "mean" or "median"
#' @param mc_cores number of cores to run in parallel (=1 for no parallelization). Note: Parallelization is currently only supported for Unix-like systems
#' @param matrix logical. If TRUE, then a binary matrix of the dimensionality of the input data is returned in which each 1 represents a changepoint
#' @param teststat Logical. Shall the test statistics for each period be returned?
#' @param CP_at_segment_end Logical. Should the changepoint be the last period of an ending segment (TRUE, default, as in Shao et al.) or the first period of a beginning segment (FALSE)?
#' 
#' @export
changepoints_multivariate <- function(data,
                                      alpha =  c(0.005, 0.001),
                                      theta = "mean",
                                      mc_cores = parallel::detectCores(),
                                      matrix = FALSE,
                                      teststat = FALSE,
                                      CP_at_segment_end = TRUE)
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
                           changepoint_multivariate_singleseries(data[[i]], Kn, theta)
                         })
  } else{
    cp_ls <- 
      lapply(X = seq_len(num_respondents), 
             FUN = function(i){
               changepoint_multivariate_singleseries(data[[i]], Kn, theta)
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
  
  
  ## if requested, convert the changepoint location to start of segment
  CP <- changepoint_conversion(
    CP, matrix = matrix, CP_at_segment_end = CP_at_segment_end)
  
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


## in the Shao et al. paper whose CP detection method we implement, a CP is defined
# as the last point in a segment, i.e., the last attentive response in our context.
# Here we add 1 to the location of the detected CP so that a CP corresponds
# to the first careless response (as described in the paper)
# The inout `CP` is a list where each elements is a vector or matrix of CP 
# locations for each significance level. Logical input `matrix` says if CP holds matrices
# or vectors
changepoint_conversion <- function(CP, matrix, CP_at_segment_end = TRUE)
{
  if(CP_at_segment_end)
  {
    CP_adj <- CP ## nothinng needs to be done here
  } else
  {
    if(matrix)
    {
      ## here CP is a binary matrix that is 1 if the corresponding column 
      # is a CP. To convert, we can simply shift the 1s one column to the right 
      # and replace their original position with 0.
      CP_adj <- lapply(CP, function(x) {
        ## here x is a matrix
        shifted <- cbind(0L, x[, -ncol(x)])
        return(shifted)
      })
    } else
    {
      ## here CP is vector of NAs, and non-NA values are CP locations. 
      # So adding 1 does the job
      CP_adj <- lapply(CP, function(x) x + 1L)
    }
  }
  return(CP_adj)
}
