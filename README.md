# carelessonset 

Implementation of the working paper "I Don't Care Anymore: Identifying the Onset of Careless Responding" by Max Welz and Andreas Alfons. This paper proposes a method to identify the onset of careless responding (or an absence thereof) in lengthy questionnaires on the participant-level. 

This repository is organized as follows. The folder `R` contains two subfolders:

- `autoencoder` contains code for the autoencoder as well as the longstring sequences,
- `changepoints` contains code for the chanegpoint identification (main file: `changepoints.R`).

**We emphasize that this implementation has not yet been thoroughly tested, so we cannot yet guarantee correctness or stability!**

## Dependencies
You need to have the following `C++` installations:
- `C++17` and a `C++` compiler,
- The `eigen3` library (http://eigen.tuxfamily.org/index.php?title=Main_Page#Download),
- CRAN packages `Rcpp` and `Rcpp`; installable via the `R` command `install.packages(c("Rcpp", "RcppEigen"))`.

In addition, the code requires an installation of the Python libraries [TensorFlow](https://www.tensorflow.org/install) and [Keras](https://keras.io/), as well as an installation of [Python](https://www.python.org/downloads/) itself and required dependencies. All components can easily be installed by uncommenting and running the `R` script below:

```R
## download R interface for keras
install.packages("keras") # install keras

## Note: you may be asked to install Miniconda if there is no Python installation on your machine. Agree to this.

## install keras
keras::install_keras(method = "auto", conda = "auto", version = "2.9")  
```

In case of questions, please get in touch with Max Welz (`welz <at> ese <dot> eur <dot> nl`).

