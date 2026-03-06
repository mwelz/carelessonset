#' calculate longstrongpattern (LSP)
#' 
#' @param x a vector of responses of a given respondent
#' @param maxlen maximum length that a pattern can encompass
#'
#' @export
LSP <- function(x, maxlen = 5)
{
  longstring_J <- sapply(seq_len(maxlen), function(J) LSP_J(x = x, J = J))
  return(apply(longstring_J, MARGIN = 1, base::max))
} # FUN
