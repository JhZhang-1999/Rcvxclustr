\name{weight_delta}
\alias{weight_delta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The proposed new weight function that can be applied in fusion penalty
}
\description{
Calculate a weight matrix according to design matrix X, and parameters zeta and delta. 
}
\usage{
weight_delta(X, delta, zeta)
}
\arguments{
  \item{X}{The data matrix of the size n*n.}
  \item{delta}{To judge whether two rows of X are close enough. A distance less equal than delta is considered small. }
  \item{zeta}{A nonnegative number used in the weight function.}
}
\details{
}
\value{
  \item{wt}{The weight matrix.}
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
}
