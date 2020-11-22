\name{Huber_ADMM}
\alias{Huber_ADMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Alternating Direction Method of Multipleirs Algorithm for Solving Convex Clustering Objective}
\description{
ADMM algorithm for solving convex clustering objective, with Huber loss and fusion penalty. Weights in fusion penalty can be either specified or not (by default, all equal to 1). See more in reference. 
}
\usage{
Huber_ADMM(X, initV = NULL, initY = NULL, initW = NULL, initZ = NULL, maxit = 1000, rho = 1, tau = 3, lam = 0.5, eps = 0.001, wt = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The data matrix with n observations and p features.}
  \item{initV}{The initial matrix V.}
  \item{initY}{The initial matrix Y.}
  \item{initW}{The initial matrix W.}
  \item{initZ}{The initail matrix Z.}
  \item{maxit}{The maximum iterate time.}
  \item{rho}{The nonnegative tuning parameter for the ADMM algorithm.}
  \item{tau}{The robustification parameter.}
  \item{lam}{The tunning parameter.}
  \item{eps}{The tolerance level.}
  \item{wt}{The weight.}
}
\details{}
\value{
  \item{W}{W determines the desired centriod matrix U together with V.}
  \item{V}{V determines the desired centriod matrix U together with W.}
  \item{iteration}{The final iteration time.}
  \item{convergence}{The p2 norm of (W - Wold), which describes how the iteration converges. }
}
\references{
Liu, Chenyu & Sun, Qiang & Tan, Kean. (2019). Robust convex clustering: How does fusion penalty enhance robustness?. https://arxiv.org/abs/1906.09581v2
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
### TODO
}
