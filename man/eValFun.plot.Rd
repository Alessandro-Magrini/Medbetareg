\name{eValFun.plot}
\alias{eValFun.plot}
\title{Marginal conditional expected value functions}
\description{Plot of the conditional expected value as a function of
  a continuous explanatory variable (given that all other explanatory variables
  take their respective neutral values).}
\usage{eValFun.plot(B=1, yscale=NULL, xscale=NULL)}
\arguments{
  \item{B}{The regression coefficient associated to the explanatory variable. Default is 1.}
  \item{yscale}{A vector including the 4 cut points of the response. If omitted, rescaled cut points are used.}
  \item{xscale}{A vector including the 4 cut points of the explanatory variable. If omitted, rescaled cut points are used.}
  }
\value{An object of class \code{plot}.}
\examples{
eValFun.plot(B=1)
eValFun.plot(B=1, yscale=c(0,15,25,40), xscale=c(0,0,30,100))
}