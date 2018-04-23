#' @useDynLib simGWAS
#' @importFrom Rcpp sourceCpp
#' @importFrom corpcor make.positive.definite
#' @importFrom methods new
#' @importFrom combinat hcube
#' @importFrom dplyr group_by_ summarise
#' @importFrom mvtnorm rmvnorm

Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

NULL
