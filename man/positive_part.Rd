\name{positive_part}
\alias{positive_part}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Positive Part Function
}
\description{
Extract the positive part of a vector. All negative part will be set to zero.
}
\usage{
positive_part(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The variable which is a vector. }
}
\details{}
\value{
\item{a}{The vector returned. }
}
\references{}
\author{}
\note{}
\seealso{}
\examples{
x <- c(0,-1,3,5,-10,7)
y <- positive_part(x)
}
