#' Compute Robust weights
#'
#' \code{robust_weights} computes the robust weights given a data matrix \code{X}, 
#' a scale parameter \code{zeta} and a parameter that controls the weights \code{delta}.
#' Namely, the lth weight \code{w[l]} is given by 
#' \deqn{
#' w[l] = exp(-\zeta{\sum_{j\in D_{1}}(X_{i'j}-X_{ij})^2+\sum_{j\in D_{2}}{\delta}^2})
#' , where the lth pair of nodes is (\code{i},\code{i'})
#' and \code{D1={j:|X_{ij}-X_{i'j}|<delta}}, \code{D2={j:|X_{ij}-X_{i'j}|>delta}}.
#' }
#' @param X The data matrix to be clustered. The rows are observations, and the columns 
#' are features.
#' @param delta The nonnegative parameter that controls the scale of robust weights 
#' when there is outlier(s) in the data.
#' @param zeta The nonnegative parameter that controls the scale of robust weights.
#' @author Chenyu Liu, Qiang Sun, Kean Ming Tan
#' @useDynLib Rcvxclustr
#' @import gdata
#' @export
#' @return A vector \cite{wt} of weights for robust convex clustering.
  robust_weights <- function(X, delta, zeta){
      sol <- robustweights(X=X,delta=delta,zeta=zeta)
      weights <- lowerTriangle(sol)
      return(weights/max(weights))
}
