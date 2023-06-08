# library("reticulate")
# In Rstudio: global options -> Python -> select python 3.10 as interpreter [we use this as virtualenv instead of conda]
# sudo python3.10 -m pip install virtualenv tensorflow
# py_available(initialize = T) # T
tensorflow::tf_config()

rm(list = ls()) ; cat("\014")
library("qgraph")
data(big5) # 500 respondents and 240 items (latent dimensionality: 30). No time and page design info
detach("package:qgraph", unload = TRUE)
dim(big5) # 500 240
any(is.na(big5)) # FALSE

rownames(big5) <- seq_len(nrow(big5))

# for some reason, there are non-integer responses -> imputation?
table(as.numeric(big5))
sum(!big5 %in% c(1,2,3,4,5))

# restrict analysis to integer responses
big5_red <- big5
noninteger <- !big5 %in% c(1,2,3,4,5)
sum(noninteger)
big5_red[noninteger] <- NA
big5_red <- na.omit(big5_red)
dim(big5_red) # 400


# prepare autoencoder architecture
source("R/changepoints_single.R")
source("R/tf.R")

# initialize
num_constructs <- 30L
num_likert     <- 5
batch_size     <- 10L
alpha          <- 1 - c(0.9, 0.95, 0.975, 0.99, 0.995, 0.999)
mc_cores       <- 6L

# number of items
n         <- nrow(big5_red)
num_items <- ncol(big5_red)

# specify hidden layers 
hidden_layers <- c(floor(1.5 * num_items), num_constructs, floor(1.5 * num_items))

# run autoencoder on contaminated dataset
set.seed(12345)
autoencoder_ls <- autoencoder(data = big5_red, 
                              hidden_layers = hidden_layers, 
                              activation = c("tanh", "linear", "tanh"), 
                              batch_size = batch_size,
                              epochs = 100, 
                              loss = pseudo_huber_loss,
                              kernel_regularizer_HL1 = NULL,
                              optimizer = optimizer_sgd(learning_rate = 1e-04),
                              verbose = 0, 
                              seed = 12345)

# reconstruct data and get outlyingness scores
reconstructed <- autoencoder_ls$reconstructed
scores        <- ((big5_red - reconstructed) / num_likert)^2

# longstring
lngstrng <- t(sapply(1:n, function(i) longstring_varyJ(big5_red[i,], maxlen = num_likert) ))

# induce some tiny random noise in LS (otherwise CP detection crashes due to division by 0)
set.seed(1)
lngstrng0 <- lngstrng + 
  matrix(rnorm(n * num_items, mean = 0, sd = 0.01),
         n, num_items)

# combine scores and time and longstring 
LS_SC    <- lapply(1:n, function(i) rbind(scores[i,], lngstrng0[i,]))

# multidimensional CP detection
CP_LS_SC <- detect_changepoints_SNsingle_MTS(LS_SC, alpha = alpha, mc_cores = mc_cores, matrix = FALSE, teststat = TRUE)

# univariate CP detection
CP_LS <- detect_changepoints_SNsingle(lngstrng0, alpha = alpha, mc_cores = mc_cores, matrix = FALSE, teststat = TRUE)
CP_SC <- detect_changepoints_SNsingle(scores, alpha = alpha, mc_cores = mc_cores, matrix = FALSE, teststat = TRUE)


save(big5, big5_red, LS_SC, CP_LS, CP_LS_SC, CP_SC, lngstrng, scores, file = "application/big5/big5.Rdata")