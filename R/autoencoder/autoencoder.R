library("keras")

#' calculates the group lasso penalty where each group consists of items that were on the same page
#' 
#' for simplicity, we assume that there is the same number of items on each page
#' 
#' @param x a weight matrix of class df.Tensor of dimenson num_items x num_nodes_L1
#' @param num_pages Number of pages
#' @param num_items_per_page number of items on each page 
#' @param num_nodes_L1 number of nodes in the first hidden layer
#' @param lambda tuning parameter 
#'
#' @export
grouplasso_page <- function(x, num_pages, num_items_per_page, num_nodes_L1, lambda = 0.01){
  
  # initialize
  out <- 0
  
  # initialize base set
  base <- 1:num_items_per_page
  
  for(j in 0:(num_pages - 1)){
    
    # calculate L2 norm
    idx <- base + j * num_items_per_page
    l2 <- k_sqrt(k_sum(k_pow(x[idx,], 2)))
    
    # add to total sum
    out <- l2 + out
  }
  
  # multiply with cardinality per group (constant here so can be outside the for-loop) 
  lambda * sqrt(num_items_per_page * num_nodes_L1) * out
} # FUN

# specify loss function
pseudo_huber_loss <- function(y_true, y_predict)
{
  delta <- k_constant(1)
  resid <- y_true - y_predict
  one   <- k_constant(1)
  k_sum(k_pow(delta, 2) * (k_sqrt(1 + k_pow(resid / delta, 2)) - 1))
} # FUN


#' fit an autoencoder
#' @param data data matrix that holds the responses of a given respondent in its rows
#' @param hidden_layers specify number of hidden layers and the number of nodes therein
#' @param activation specify activation function in each hidden layer
#' @param loss a function object for the loss function
#' @param optimizer keras object for the optimizer
#' @param kernel_regularizer_HL1 regularization imposed in hidden layer on weights
#' @param bias_regularizer_HL1 regularization imposed in hidden layer on biases
#' @param epochs number of epochs
#' @param batch_size batch size
#' @param verbose manage prints
#' @param seed random seed
#' 
#' @export
autoencoder <- function(data,
                        hidden_layers = c(10, 2, 10),
                        activation = c("tanh", "linear", "tanh"),
                        loss = pseudo_huber_loss,
                        optimizer = optimizer_sgd(learning_rate = 1e-04),
                        kernel_regularizer_HL1 = NULL, 
                        bias_regularizer_HL1 = NULL,
                        epochs = 100L,
                        batch_size = 10L,
                        verbose = 2L,
                        seed = NULL)
{
  stopifnot(!any(is.na(data)) && is.numeric(data))
  
  ## standardize data
  data        <- as.matrix(data)
  input_size  <- ncol(data)
  m           <- colMeans(data)
  s           <- sapply(1:input_size, function(j) sqrt(var(data[,j])) )
  X           <- sapply(1:input_size, function(j) (data[,j] - m[j]) / s[j] )
  colnames(X) <- colnames(data)
  
  ## prepare the layers
  num_hl <- length(hidden_layers) # number of hidden layers
  stopifnot(length(activation) == num_hl)
  layers <- lapply(1:(num_hl + 1L), function(j){
    
    if(j == 1){
      
      # input layer
      layer_dense(units = hidden_layers[j], 
                  activation = activation[j],
                  input_shape = input_size, 
                  kernel_regularizer = kernel_regularizer_HL1,
                  bias_regularizer = bias_regularizer_HL1)
      
    } else if(j == num_hl + 1L){
      
      # output layer
      layer_dense(units = input_size,
                  activation = "linear")
      
    } else {
      
      # hidden layers
      layer_dense(units = hidden_layers[j], 
                  activation = activation[j])
    } # IF
    
  }) # SAPPLY
  
  
  ## define the autoencoder
  model <- keras_model_sequential(layers)
  
  
  ## details of the optimizer
  model %>% compile(
    loss      = loss,
    optimizer = optimizer
  )
  
  
  ## specify random seed
  if(!is.null(seed)){
    tensorflow::set_random_seed(seed)
  }
  
  
  ## fit the model
  model %>% fit(X, X,
                epochs = epochs,
                batch_size = batch_size, 
                verbose = verbose)
  
  ## make in-sample prediction
  Xhat <- model %>% predict(X)
  
  ## de-standardize the predictions
  Xhat_ds <- sapply(1:input_size, function(j) Xhat[,j] * s[j] + m[j] )
  colnames(Xhat_ds) <- colnames(data)
  
  ## reconstruction error
  rec_err <-  sapply(1:nrow(X), function(i) mean((Xhat_ds[i,] - data[i,])^2) )
  
  return(list(reconstructed = Xhat_ds,
              reconstruction_error = rec_err,
              model = model))
  
} # FUN
