\name{rescale}
\alias{rescale}
\title{Rescaling procedure}
\description{Rescaling for a set values of a medical variable.}
\usage{rescale(v, vscale)}
\arguments{
  \item{v}{A numerical value or a numerical vector containing the values to be rescaled.}
  \item{vscale}{A vector of 4 cutpoints if the variable is continuous, or \code{NA} if the variable is binary.}
}
\value{A numerical value or a numerical vector.}
\examples{
rescale(v=seq(0,40,length=21), vscale=c(0,15,25,40))
}