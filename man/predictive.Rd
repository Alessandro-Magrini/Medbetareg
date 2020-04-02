\name{predictive}
\alias{predictive}
\title{Predictive distributions}
\description{Plot of a predictive distribution chosen by the user.}
\usage{predictive(x, cnfg, nrep = 500, title = NULL)}
\arguments{
  \item{x}{An object of class \code{mbr}.}
  \item{cnfg}{A numerical vector representing a configuration of explanatory variables (original scale).}
  \item{nrep}{An integer positive number representing the number of Monte Carlo replications. Default is 500.}
  \item{title}{The title of the plot. If \code{NULL} (the default), the configuration of explanatory variables is shown as a title.}
}
%\value{}
\seealso{\code{\link{newPrior}}}
\examples{
assess.test <- 'RESP RespRate (0,15,25,40);
  CEV intraShunt (0,2,5,100) 0.5 hp 0.5 hp 5;
  CEV deadSpace (0,0,30,100) 0.5 hp 0.5 hp 5;
  CEV extraShunt (0,0,5,100) 0.5 hp 0.5 hp 5;
  CEV redAlvSpace (0,0,5,100) 0.5 hp 0.5 hp 5;
  BEV Panic 0.25 hp 25;
  BEV Neuromusc 0.6 lp 100;
  INTER intraShunt deadSpace 0.9 hp 5;
  TAU 0.3 n'
###  NOT RUN: replicate the results in Magrini et al. (2018)
# set.seed(10)
# prior.test <- newPrior(assess.test, nrep=5000)
# set.seed(10)
# predictive(prior.test, xcfg1, nrep=50000, title="Configuration 1")
# set.seed(10)
# predictive(prior.test, xcfg2, nrep=50000, title="Configuration 2")
# set.seed(10)
# predictive(prior.test, xcfg3, nrep=50000, title="Configuration 3")
# set.seed(10)
# predictive(prior.test, xcfg4, nrep=50000, title="Configuration 4")
#############
prior.test <- newPrior(assess.test, nrep=100)
predictive(prior.test, c(5,30,5,5,0,0), nrep=500)
}