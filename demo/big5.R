## check for correct keras specification
# library("reticulate")
# In Rstudio: global options -> Python -> select python 3.10 as interpreter [we use this as virtualenv instead of conda]
# sudo python3.10 -m pip install virtualenv tensorflow
# py_available(initialize = TRUE) # TRUE
# tensorflow::tf_config()

source("demo/dataprep.R")
source("R/load.R")

## get preprocessed data
# 400 participants replied to 240 NEO-PI-R items
# measures 30 facets of personality traits
# unfortunately no information on time or page membership
load("demo/big5.Rdata")

## detect changepoints
x <- carelessonset(responses = big5$responses, 
                   num_scales = big5$num_scales,
                   num_likert = big5$num_likert, 
                   seed = 12345, 
                   mc_cores = 6) # only deterministic operations are parallelized

# save(x, plot_dimension, get_changepoints, plot.carelessonset, file = "demo/carelessonset.Rdata")
load(file = "demo/carelessonset.Rdata")


# one changepoint flagged, for respondent no. 94
get_changepoints(x)

# plot the two dimensions for this respondent
(p <- plot.carelessonset(x, idx = 94))

# save
ggsave(filename = "demo/plot_idx94.pdf", 
       plot = p, 
       device = "pdf", width = 297*0.3,  height = 210*0.3, units = "mm") 
