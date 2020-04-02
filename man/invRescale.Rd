\name{invRescale}
\alias{invRescale}
\title{Inverse rescaling}
\description{Inverse rescaling for a set of values of a medical variable.}
\arguments{
  \item{v}{A numerical value or a numerical vector containing the values to be rescaled.}
  \item{vscale}{A vector of 4 cutpoints if the variable is continuous, or \code{NA} if the variable is binary.}
}
\value{A numerical value or a numerical vector.}
\examples{
invRescale(v=seq(-1.5,1.5,length=21), vscale=c(0,15,25,40))
}