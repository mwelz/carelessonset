rm(list = ls()) ; cat("\014")
source("R/singleCP_plot.R")


load("application/big5/big5.Rdata")
alpha <- "0.005"
dir_save <- "application/big5/plots"


sum(!is.na(CP_LS_SC$changepoints[[alpha]])) # 4
sum(!is.na(CP_SC$changepoints[[alpha]])) # 3
sum(!is.na(CP_LS$changepoints[[alpha]])) # 5

# get estimated changepoints
cparr <- CP_LS_SC$changepoints[[alpha]]
n <- length(cparr)
p <- ncol(big5_red) # number of items
crlss_idx <- which(!is.na(cparr))
crfl_idx <- setdiff(1:n, crlss_idx)
idx_red <- rownames(big5_red)
  

# SNCP value for careless
Tn_crlss <- CP_LS_SC$SNCP$Tmax[crlss_idx]
teststats <- CP_LS_SC$SNCP$Tn

# sd of careful respondents (and careless)
sd_crfl <- mean(sapply(crfl_idx, function(j) var(big5_red[j,])))
sd_crlss <- mean(sapply(crlss_idx, function(j) var(big5_red[j,])))

# init
khat <- sd_pre_crfl <- sd_post_crfl <- sd_pre_crlss <- sd_post_crlss <- 
  rep(NA_real_, length(crlss_idx))

for(i in seq_along(crlss_idx))
{
  idx <- crlss_idx[i]
  k <-  cparr[idx]
  khat[i] <- k
  idx_red_i <- idx_red[idx]
  
  # standard deviation within periods, averaged across careful respondents...
  sd_pre_crfl[i]  <- mean(sapply(crfl_idx, function(j) var(big5_red[j, 1:(k-1)])))
  sd_post_crfl[i] <- mean(sapply(crfl_idx, function(j) var(big5_red[j, k:p])))
  
  # ... and the careless respondent
  sd_pre_crlss[i] <- var(big5_red[idx, 1:(k-1)])
  sd_post_crlss[i] <- var(big5_red[idx, k:p])
  
  # plot
  gg <- plot_dimension(RE = scores[idx,], ALSP = lngstrng[idx,], onset = k)
  
  gg_Tn <- plot_teststat(x = teststats[idx,], d = 2, alpha = 1 - c(0.975, 0.99, 0.995, 0.999)) 
  
  ggsave(filename = paste0("/appl_flagged_CP_LS_SC_idx", idx_red_i, ".pdf"), 
         plot = gg, path = dir_save, 
         device = "pdf", width = 297*0.3,  height = 210*0.3, units = "mm") 
  
  ggsave(filename = paste0("/appl_Tn_CP_LS_SC_idx", idx_red_i, ".pdf"), 
         plot = gg_Tn, path = dir_save, 
         device = "pdf", width = 297*0.35,  height = 210*0.35, units = "mm") 
  
} # FOR


out <- cbind(idx = crlss_idx, khat = khat, Tn = Tn_crlss, 
      sd_pre_crlss = sd_pre_crlss, sd_post_crlss = sd_post_crlss,
      sd_pre_crfl = sd_pre_crfl, sd_post_crfl = sd_post_crfl)
write.csv(out, file = "application/big5/results.csv")

sd_crfl
sd_crlss

