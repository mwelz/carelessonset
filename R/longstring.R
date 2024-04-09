# library("Rcpp")

# load LSP()
# sourceCpp("R/autoencoder/longstring.cpp")

#' calculate adaptive longstrongpattern
#' 
#' @param x a vector of responses of a given respondent
#' @param maxlen maximum length that a pattern can encompass
#'
#' @export
ALSP <- function(x, maxlen = 5)
{
  longstring_J <- sapply(seq_len(maxlen), function(J) LSP(x = x, J = J))
  return(apply(longstring_J, MARGIN = 1, base::max))
} # FUN
