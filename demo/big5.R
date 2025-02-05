## load package
# devtools::install_github("mwelz/carelessonset")

## check for correct keras specification
# library("reticulate")
# In Rstudio: global options -> Python -> select python 3.10 as interpreter [we use this as virtualenv instead of conda]
# sudo python3.10 -m pip install virtualenv tensorflow
# py_available(initialize = TRUE) # TRUE
# tensorflow::tf_config()

# load functions (may cause compiler noise)
library("carelessonset")

## get preprocessed data
# 400 participants replied to 240 NEO-PI-R items
# measures 30 facets of personality traits
# unfortunately no information on time or page membership
load("demo/big5.Rdata")

# view data
dim(big5)
View(big5)

## detect changepoints
onset <-  carelessonset:::carelessonset(responses = big5, 
                        num_scales = 30L, 
                        num_likert = 5L, 
                        seed = 12345)

# save(onset, plot_dimension, get_onset, plot.carelessonset, file = "demo/carelessonset.Rdata")
# load(file = "demo/carelessonset.Rdata")


# one changepoint flagged, for respondent no. 94
alpha <- 0.001
carelessonset:::get_onset(onset, alpha = alpha)

# plot the two dimensions for this respondent
(p <- plot.carelessonset(onset, idx = 94, alpha = alpha))


# # save
# ggsave(filename = "demo/plot_idx94.pdf", 
#        plot = p, 
#        device = "pdf", width = 297*0.3,  height = 210*0.3, units = "mm") 
