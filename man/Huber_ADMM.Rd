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

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, initV = NULL, initY = NULL, initW = NULL, initZ = NULL, 
    maxit = 1000, rho = 1, tau = 3, lam = 0.5, eps = 0.001, wt = NULL) 
{
    n = nrow(X)
    p = ncol(X)
    I = diag(rep(1, n))
    for (i in seq(from = n - 1, to = 1, by = -1)) {
        if (i > 1) {
            temp = diag(rep(-1, i))
        }
        else {
            temp = -1
        }
        temp1 = cbind(rep(1, i), temp)
        if (n - i - 1 == 0) {
            E = temp1
        }
        else {
            E = rbind(E, cbind(matrix(0, ncol = n - i - 1, nrow = i), 
                temp1))
        }
    }
    k = dim(E)[1]
    if (is.null(initV) == TRUE) {
        V <- matrix(rnorm(k * p), nrow = k, ncol = p)
    }
    else {
        V <- initV
    }
    if (is.null(initY) == TRUE) {
        Y <- matrix(rnorm(k * p), nrow = k, ncol = p)
    }
    else {
        Y <- initY
    }
    if (is.null(initW) == TRUE) {
        W <- matrix(rnorm(n * p), nrow = n, ncol = p)
    }
    else {
        W <- initW
    }
    if (is.null(initZ) == TRUE) {
        Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
    }
    else {
        Z <- initZ
    }
    if (is.null(wt) == TRUE) {
        wt = rep(1, k)
    }
    sigma <- lam * wt/rho
    A <- matrix(tau, nrow = n, ncol = p)
    for (t in 2:maxit) {
        Vold <- V
        Yold <- Y
        Wold <- W
        Zold <- Z
        U <- solve(t(E) \%*\% E + I, t(E) \%*\% (Vold + Yold) + Wold + 
            Zold)
        B <- abs(rho * (X - U + Zold))/(1 + rho)
        index <- which(B <= A, arr.ind = T)
        W <- X - Soft(X - U + Zold, tau/rho)
        W[index] <- (X[index] + rho * (U[index] - Zold[index]))/(1 + 
            rho)
        v <- E \%*\% U - Yold
        v.norm <- apply(v, 1, function(x) {
            norm(x, type = "2")
        })
        V <- v * positive_part(1 - sigma/v.norm)
        Y <- Yold - rho * (E \%*\% U - V)
        Z <- Zold - rho * (U - W)
        if (norm(W - Wold, type = "2") < eps) {
            break
        }
        Uold <- U
    }
    return(list(W = W, V = V, iteration = t, convergence = norm(W - 
        Wold, type = "2")))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
