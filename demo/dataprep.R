get_prepped_big5 <- function()
{
  library("qgraph")
  data(big5) # 500 respondents and 240 items (latent dimensionality: 30). No time and page design info
  detach("package:qgraph", unload = TRUE)
  
  # restrict analysis to integer responses
  big5_red <- big5
  noninteger <- !big5 %in% c(1,2,3,4,5)
  big5_red[noninteger] <- NA
  big5_red <- na.omit(big5_red)
  
  return(big5_red)
}

big5 <- get_prepped_big5()

save(big5, file = "demo/big5.Rdata")