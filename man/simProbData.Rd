% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simConnectivity.R
\name{simProbData}
\alias{simProbData}
\title{Simulate Dirichlet-based probability table data}
\usage{
simProbData(
  psi,
  originRelAbund,
  sampleSize,
  shapes,
  captured = "target",
  requireEveryOrigin = FALSE
)
}
\arguments{
\item{psi}{Transition probabilities between B origin sites and W target
sites. B by W matrix}

\item{originRelAbund}{Vector of relative abundances at B origin sites}

\item{sampleSize}{Either the total number of data points to simulate or a
vector with the number at each target or origin site. If only the total is
provided, sampling will be done in proportion to abundance}

\item{shapes}{If captured == "target", a B by B matrix, each row of which is
the shape parameters for the Dirichlet distribution of an animal whose true
origin assignment is that row's. If captured == "origin", a W by W matrix,
each row of which is the shape parameters for the Dirichlet distribution of
an animal whose true target assignment is that row's.}

\item{captured}{Either "target" (the default) or "origin", indicating which
side animal data were collected on}

\item{requireEveryOrigin}{If TRUE, the function will throw an error if it
looks like at least one origin site has no animals released in or
migrating to it, or if it can, keep simulating until representation is
met. This helps estTransition or estMC not throw an error. Default FALSE}
}
\value{
\code{simProbData} returns a list with the elements:
\describe{
  \item{\code{originAssignment}}{Vector with true origin site of each animal}
  \item{\code{targetAssignment}}{Vector with true target site of each animal}
  \item{\code{genProbs}}{Table of assignment site probabilities for each
   animal}
  \item{\code{input}}{List containing the inputs to function}
}
}
\description{
Simulate Dirichlet-based probability table data
}
