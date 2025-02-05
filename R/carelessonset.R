#' detect the onset of careless responding
#' TODO: the order at which alpha is passed seems to matter, otherwise get_changpoint() output has the changepoints of the wrong dimension. This error is in R/changepoints/multivariate.r (and possibly also univariate) EDIT: fixed it, but double-check! Also, autoencoder() prompts training process (verbose); change that. Also, currently a miniconda installation  of python is implicitly required. Make that explicit and add an installer function. Also, some functions are not exported properly (so cannot be called with ordinary namespace)
#' 
#' @param responses data matrix that holds the responses of a given respondent in its rows. Must be in the order as presented to each participant
#' @param num_scales number of psychometric scales in the data
#' @param num_likert number of likert-type respnse options (TODO: allow for vector-valued input)
#' @param time data matrix of per-item response time (TODO: allow for per-page time passing)
#' @param longstring shall longstring indices be computed and used?
#' @param item_order A matrix holding the item indices on the participant level. If responses ae reshuffled according to this order, then each column in the response matrix are responses to the same item
#' @param alpha significance levels
#' @param mc_cores number of cores for parallelization
#' @param encoder_width TODO
#' @param encoder_activation TODO
#' @param bottleneck_activation TODO
#' @param loss a function object for the loss function
#' @param optimizer keras object for the optimizer
#' @param kernel_regularizer_HL1 regularization imposed in hidden layer on weights
#' @param bias_regularizer_HL1 regularization imposed in hidden layer on biases
#' @param maxlen Maximum length of an LSP pattern
#' @param epochs number of epochs
#' @param batch_size batch size
#' @param verbose manage prints
#' @param seed random seed
#' 
#' @import keras
#' 
#' @export
carelessonset <- function(responses, 
                          num_scales,
                          num_likert,
                          time = NULL,
                          longstring = TRUE,
                          item_order = NULL,
                          alpha = c(0.005, 0.001),
                          mc_cores = parallel::detectCores(),
                          encoder_width = floor(1.5 * ncol(responses)),
                          encoder_activation = "tanh",
                          bottleneck_activation = "linear",
                          loss = get_pseudo_huber(),
                          optimizer = keras::optimizer_sgd(learning_rate = 1e-04),
                          kernel_regularizer_HL1 = NULL, 
                          bias_regularizer_HL1 = NULL,
                          maxlen = max(num_likert),
                          epochs = 100L,
                          batch_size = 10L,
                          verbose = 0L,
                          seed = NULL)
{
  
  ## dimensions of data
  n <- nrow(responses)
  num_items <- ncol(responses)
  responses_ordered <- responses
 
  if(!is.null(item_order))
  {
    # if item order differs for each participant, rearrange them so that every column measures responses to same item
    for(i in seq_len(n))
    {
      responses_i <- responses[i,]
      order_i <- item_order[i,]
      responses_i <- responses_i[order_i] # reorder to item order
      responses_ordered[i,] <- responses_i
    }
  } # IF
  
  ## specify hidden layers and activation
  hidden_layers <- c(encoder_width, num_scales, encoder_width)
  activation <- c(encoder_activation, bottleneck_activation, encoder_activation)
  
  ## run the autoencoder on the ordered responses
  ann <- autoencoder(data = responses_ordered, 
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
  RE_ordered    <- ((responses_ordered - reconstructed) / (num_likert - 1))^2
  
  ## put the REs back in the order in which items were presented to respondents
  # the columns in 'RE' are in the original order (in which each respondent responded) 
  if(!is.null(item_order))
  {
    RE <- matrix(NA_real_, n, num_items)
    
    for(i in seq_len(n))
    {
      RE_i <- RE_ordered[i,]
      order_i <- item_order[i,]
      RE_i <- RE_i[order(order_i)] # reshuffle to original order by ordering the ordered indices
      RE[i,] <- RE_i
    } # FOR
  } else
  {
    RE <- RE_ordered
  } # IF
  
  
  if(!is.null(time))
  {
    time0 <- time + 
      matrix(stats::rnorm(n * num_items, mean = 0, sd = 0.01),
             n, num_items)
  } else{
    time0 <- NULL
  }
  
  # reinitialize seed. TODO: make this more elegant
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  
  
  ## combine RE and time and longstring 
  # NB: these scores are expected to be in the order as presented to each participant
  # the code above ensures that this is done. Important because we need original order
  # for changepoint detection
  if(longstring)
  {
    lngstrng <- t(sapply(seq_len(n), function(i) ALSP(responses[i,], maxlen = maxlen) ))
    
    ## induce some tiny random noise in LS and time (otherwise CP detection crashes due to division by 0)
    if(!is.null(seed))
    {
      set.seed(seed)
    }
    
    lngstrng0 <- lngstrng + 
      matrix(stats::rnorm(n * num_items, mean = 0, sd = 0.01),
             n, num_items)
    
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
  } else{
    
    lngstrng <- lngstrng0 <- NULL
    
    if(is.null(time))
    {
      series <- RE
      cp <- detect_changepoints_univariate(data = series, 
                                           alpha = alpha, 
                                           theta = "mean", 
                                           mc_cores = mc_cores, 
                                           matrix = FALSE,
                                           teststat = TRUE)
      
    } else{
      series <- lapply(1:n, function(i) rbind(RE[i,], time0[i,]))
      cp <- detect_changepoints_multivariate(data = series, 
                                             alpha = alpha, 
                                             theta = "mean", 
                                             mc_cores = mc_cores, 
                                             matrix = FALSE,
                                             teststat = TRUE)
    } 
  } # IF
  
  
  
  # prepare output
  out <- 
    list(changepoints = cp$changepoints,
         teststatistics = cp$SNCP,
         autoencoder = ann, 
         alpha = alpha,
         series = list(RE = RE, ALSP = lngstrng, ALSP_noise = lngstrng0, time = time))
  class(out) <- "carelessonset"
  return(out)
} # FUN

