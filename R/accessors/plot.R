library("ggplot2")

#' get critical values of SNCP test (Table 1 in Shao & Zhao, 2010, JASA)
#' @param d dimension of the initial series
#' @param alpha vector of significance levels; possibly \code{NULL}
get_Kn <- function(d, alpha)
{
  
  # critical values as in Table 1 in Shao & Zhao (2010, JASA, https://doi.org/10.1198/jasa.2010.tm10103)
  critvals <- 
    cbind(
      alpha = 1 - c(0.9, 0.95, 0.975, 0.99, 0.995, 0.999),
      Kn1 = c(29.6, 40.1, 52.2, 68.6, 84.6, 121.9),
      Kn2 = c(56.5, 73.7, 92.2, 117.7, 135.3, 192.5),
      Kn3 = c(81.5, 103.6, 128.9, 160.0, 182.9, 246.8),
      Kn4 = c(114.7, 141.5, 171.9, 209.7, 246.6, 319.2)
    )
  
  # check if critical value for choices of d and alpha is available
  stopifnot(0 < d & d < 5)
  stopifnot(all(round(alpha, 6) %in% round(critvals[,1], 6)))
  
  # return corresponding critical values
  Kn        <- unname(critvals[round(critvals[,1], 4) %in% round(alpha, 4), d+1L])
  names(Kn) <- alpha
  return(Kn)
} # FUN



#' visualize SNCP test statistic over all periods
#' @param x a vector of test statistics of length \code{p-1}
#' @param d dimension of the initial series
#' @param alpha vector of significance levels; possibly \code{NULL}
plot_teststat <- function(x, d, alpha = c(0.05, 0.025, 0.01)) 
{
  # get length of series and maximum value of x
  stopifnot(is.vector(x))
  p    <- length(x) + 1L
  khat <- which.max(x)
  
  # prepare data frame and plot
  df <- data.frame(x = seq_len(p), y = c(x, 0.0))
  gg <- ggplot(df, mapping = aes(x = x, y = y)) +
    theme_bw() +
    geom_line() +
    geom_vline(xintercept = khat, linetype = "dashed", col = "grey50") +
    scale_colour_brewer(palette = "RdBu") +
    ylab("Test statistic") +
    xlab("Item Index")
    # geom_text(aes(x=khat, y=0, label = "Tmax", vjust = 1, hjust = 1))
  
  if(!is.null(alpha))
  {
    # get critical values
    Kn <- get_Kn(d = d, alpha = alpha)
    
    gg <- gg + 
      geom_hline(data = data.frame(Kn = Kn, alpha = alpha),
                    aes(yintercept = Kn, 
                        color = factor(alpha, levels = sort(alpha, decreasing = FALSE)))) +
      theme(legend.position = "bottom")
    gg$labels$colour <- "Level"
    
  } # IF
  
  return(gg)
  
} # FUN


# gg <- gg +
#    geom_hline(yintercept = Kn, color = cols) +
#    geom_text(x=0, y=Kn[i], label = paste0("α = ", alpha[i]), vjust = -0.2, hjust = 0), color = cols[i])
#geom_text(aes(x=0, y=Kn[i], label = paste0("α = ", alpha[i]), vjust = -0.2, hjust = 0), color = cols[i])


#' visualize dimensions
#' @param RE vector of reconstruction errors
#' @param ALSP vector of ALSP values
#' @param time vector of response times
#' @param onset index of carelessness onset
#' @param color color of vertical line denoting onset
plot_dimension <- function(RE = NULL, ALSP = NULL, time = NULL, onset = NULL, color = "blue")
{
  # get information about supplied input
  dims      <- c("RE", "ALSP", "time")
  dims_used <- dims[c(!is.null(RE), !is.null(ALSP), !is.null(time))]
  if(length(dims_used) < 1L) stop("At least one dimension needs to be supplied")
  
  # get length and values to plot
  lens <-  sapply(seq_along(dims_used), function(i) length(get(dims_used[i])))
  stopifnot(all(lens == lens[1L]))
  p <- lens[1L]
  y <- c(RE, ALSP, time)
  x <- rep(seq_len(p), length(dims_used))
  
  statistic <- NULL
  dims_used[dims_used %in% "time"] <- "Time" # capitalize "time"
  
  for(i in seq_along(dims_used))
  {
    statistic <- c(statistic, rep(dims_used[i], p))
  }
  statistic <- factor(statistic, levels = dims_used)
  
  # prepare data frame and get ggplot object
  df <- data.frame(y = y, x = x, statistic = statistic)
  gg <- ggplot(df, mapping = aes(x = x, y = y)) +
    theme_bw() +
    geom_line() +
    facet_grid(rows = vars(statistic), scales = "free_y") +
    xlab("Item Index") +
    ylab("Value of Dimension")
  
  # plot onset if desired
  if(!is.null(onset))
  {
    gg <- 
      gg + geom_vline(xintercept = onset, col = color, size = 0.7) 
  } # IF
  
  return(gg)
} # FUN

