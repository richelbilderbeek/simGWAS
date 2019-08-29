#' @useDynLib simGWAS
#' @importFrom Rcpp sourceCpp
#' @importFrom corpcor make.positive.definite
#' @importFrom combinat hcube
#' @importFrom dplyr group_by_ summarise
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rgamma rbinom runif
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

NULL
