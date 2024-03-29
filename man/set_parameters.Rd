% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_parameters.R
\name{set_parameters}
\alias{set_parameters}
\title{set_parameters}
\usage{
set_parameters(
  K = 200,
  lambda = 1,
  mu = 0.5,
  betaA = seq(1e-04, 0.001, 1e-04),
  betaB = seq(1e-04, 0.001, 1e-04),
  sigmaA = seq(0, 1, 0.1),
  sigmaB = seq(0, 1, 0.1),
  sigmaAB = seq(0, 1, 0.1),
  sigmaBA = seq(0, 1, 0.1),
  nuA = 0.01,
  nuB = 0.01
)
}
\arguments{
\item{K}{: carrying capacity of number of host species within the system.}

\item{lambda}{: host speciation rate}

\item{mu}{: loss rate  of host spcies}

\item{betaA}{:  transmission rate of species A}

\item{betaB}{: transmission rate of species B}

\item{sigmaA}{: priority effects of species A}

\item{sigmaB}{: priority effects of species B}

\item{sigmaAB}{: priority effects of co-infections}

\item{sigmaBA}{:priority effects of co-infections}

\item{nuA}{: loss rate of species A}

\item{nuB}{: loss rate of species B}

\item{one_species}{: If TRUE,  parameters for two sets of species will be created.}

\item{dems}{:}
}
\value{
data frame of all parameter combinations
}
\description{
set_parameters
}
\examples{
set_parameters(K = 200,lambda = 1,mu = 0.5,betaA = seq(0.0001,0.001, 0.0001),betaB = 0.001,sigmaA = 1,sigmaB = 1,sigmaAB =  0.1,sigmaBA  = 0.1,nuA = 0.01,nuB = 0.01)
}
