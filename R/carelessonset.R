# TODO: the order at which alpha is passed seems to matter, otherwise get_changpoint() output has the changepoints of the wrong dimension. This error is in R/changepoints/multivariate.r (and possibly also univariate)
carelessonset <- function(responses, 
                          num_scales,
                          num_likert,
                          time = NULL,
                          alpha = c(0.005, 0.001),
                          mc_cores = parallel::detectCores(),
                          encoder_width = floor(1.5 * ncol(responses)),
                          encoder_activation = "tanh",
                          bottleneck_activation = "linear",
                          loss = pseudo_huber_loss,
                          optimizer = optimizer_sgd(learning_rate = 1e-04),
                          kernel_regularizer_HL1 = NULL, 
                          bias_regularizer_HL1 = NULL,
                          epochs = 100L,
                          batch_size = 10L,
                          verbose = 0L,
                          seed = NULL)
{
 
  ## specify hidden layers and activation
  hidden_layers <- c(encoder_width, num_scales, encoder_width)
  activation <- c(encoder_activation, bottleneck_activation, encoder_activation)
  
  ## run the autoencoder
  ann <- autoencoder(data = responses, 
                     hidden_layers = hidden_layers,
                     verbose = verbose,
                     seed = seed, 
                     activation = activation, 
                     loss = loss,
                     optimizer = optimizer,
                     kernel_regularizer_HL1 = kernel_regularizer_HL1, 
                     bias_regularizer_HL1 = bias_regularizer_HL1, 
                     epochs = epochs, 
                     batch_size = batch_size) 
  
  ## reconstruct data and get RE
  reconstructed <- ann$reconstructed
  RE            <- ((responses - reconstructed) / num_likert)^2
  
  ## run ALSP
  n <- nrow(responses)
  num_items <- ncol(responses)
  lngstrng <- t(sapply(seq_len(n), function(i) ALSP(responses[i,], maxlen = num_likert) ))
  
  ## induce some tiny random noise in LS and time (otherwise CP detection crashes due to division by 0)
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  
  lngstrng0 <- lngstrng + 
    matrix(rnorm(n * num_items, mean = 0, sd = 0.01),
           n, num_items)
  
  if(!is.null(time))
  {
    time0 <- time + 
      matrix(rnorm(n * num_items, mean = 0, sd = 0.01),
             n, num_items)
  } else{
    time0 <- NULL
  }
  
  # reinitialize seed. TODO: make this more elegant
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  
  # combine RE and time and longstring 
  if(is.null(time))
  {
    series <- lapply(1:n, function(i) rbind(RE[i,], lngstrng0[i,]))
  } else{
    series <- lapply(1:n, function(i) rbind(RE[i,], time0[i,], lngstrng0[i,]))
  } # IF
  
  # detect changepoints
  cp <- detect_changepoints_multivariate(data = series, 
                                         alpha = alpha, 
                                         theta = "mean", 
                                         mc_cores = mc_cores, 
                                         matrix = FALSE,
                                         teststat = TRUE)
  
  # prepare output
  out <- 
    list(changepoints = cp$changepoints,
         teststatistics = cp$SNCP,
         autoencoder = ann, 
         alpha = alpha,
         series = list(RE = RE, ALSP = lngstrng, time = time0))
  class(out) <- "carelessonset"
  return(out)
} # FUN

