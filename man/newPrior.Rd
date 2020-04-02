\name{newPrior}
\alias{newPrior}
\title{Elicitation of the joint prior distribution}
\description{Computation of the joint prior distribution on parameters from expert assessments.
  Note: the original scales of variables are not required.}
\usage{newPrior(assess.code, nrep = 100, digits = 6)}
\arguments{
  \item{assess.code}{The code for expert assessments. See \code{details}.}
  \item{nrep}{The number of bootstrap replications to approximate the prior covariance matrix. Default is 100.}
  \item{digits}{The number of decimals to be shown. Default is 6.}
}
\details{The code for expert assessments must be a character string containing several instructions.
Each instruction is delimited from the other by a semicolumn or a line break, and must be of one of the types detailed below.

- Features of the response variable. Syntax:

\code{  RESP nameOfResponse scaleOfResponse}

where \code{RESP} is the keyword of the instruction, \code{nameOfResponse} is the name of the response variable,
\code{scaleOfResponse} is the set of its cutpoints indicated within round brackets and separated by commas.
There must be one and only one of this instruction;

- Features of a continuous explanatory variable. Syntax:

\code{  CEV nameOfEV scaleOfEV SelectedProp SelectedRange expectedProp expectedRange nCase}

where \code{CEV} is the keyword of the instruction, \code{nameOfEV} is the name of the continuous explanatory variable,
\code{scaleOfEV} is the set of its cutpoints indicated within round brackets and separated by commas,
\code{SelectedProp} and \code{SelectedRange} are the relative position and the range selected by the expert,
\code{expectedProp} and \code{expectedRange} are the relative position and the range of the expected value of the response assessed by the expert,
\code{nCase} is the number of patient cases on which the assessment is based.
There must be any instances of this instruction, even zero;

- Features of a binary explanatory variable. Syntax:

\code{  BEV nameOfEV expectedProp expectedRange nCase}

where \code{BEV} is the keyword of the instruction, \code{nameOfEV} is the name of the binary explanatory variable,
\code{expectedProp} and \code{expectedRange} are the relative position and the range of the expected value of the response assessed by the expert,
\code{nCase} is the number of patient cases on which the assessment is based.
There must be any instances of this instruction, even zero;

- Features of an interaction among a set of explanatory variables. Syntax:

\code{  INTER nameOfEV_first ... nameOfEV_last expectedProp expectedRange nCase}

where \code{INTER} is the keyword of the instruction, \code{nameOfEV_first ... nameOfEV_last} are the names of the interacting explanatory variables separated by spaces,
\code{expectedProp} and \code{expectedRange} are the relative position and the range of the expected value of the response assessed by the expert,
\code{nCase} is the number of patient cases on which the assessment is based.
There must be any number of instances of this instruction, even zero;
There must be any instances of this instruction, even zero;

- Assessments to determine the precision parameter. Syntax:

\code{  TAU expectedProp expectedRange}

where \code{TAU} is the keyword of the instruction,
\code{expectedProp} and \code{expectedRange} are the relative position and the range of the expected value of the response
when all explanatory variables take their respectove neutral values.
There must be one and only one of this instruction.
}
\value{An object of class \code{mbr}, that is a list with the following components:
  \itemize{\item \code{scales}: a named list containing, for each variable in the model, either a vector
  of 4 cutpoints if the variable is continuous, or \code{NA} if the variable is binary.
  The first component of the list is referred to the response variable.
  \item \code{scales}: a list containing the cutpoints for each variable.  
  \item \code{mean}: a vector containing the prior mean of parameters.            
  \item \code{vcov}: a matrix containing the prior covariance matrix of parameters.}
  S3 methods available for the \code{mbr} class:
  \itemize{\item \code{print};
  \item \code{summary}.}
  }
\seealso{\code{\link{predictive}}}
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
#############
prior.test <- newPrior(assess.test, nrep=100)
}

