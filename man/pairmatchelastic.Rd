\name{pairmatchelastic}
\alias{pairmatchelastic}
\alias{elastic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Optimal Subset Matching without Balance Constraints
}
\description{
	Conducts optimal subset matching as described in the reference.  
}
\usage{
pairmatchelastic(mdist, n = 0, val = 0)
elastic(mdist, n = 0, val = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdist}{
	distance matrix with rows corresponding to treated units and columns corresponding to controls.
}
  \item{n}{
	maximum number of treated units that can be excluded.
}
  \item{val}{
	cost of excluding a treated unit (i.e. we prefer to exclude a treated unit if it increases the total matched distance by more than \code{val}).
}
}
\details{
\code{pairmatchelastic} is the main function, which conducts an entire match.  \code{elastic} is a helper function which augments the original distance matrix as described in the reference.
The original versions of these functions were written by Paul Rosenbaum and distributed in the supplemental material to the reference.
}
\value{
\code{elastic} returns an augmented version of the input matrix \code{mdist}.  
\code{pairmatchelastic} returns a matrix of 1 column whose values are the column numbers of matched controls and whose rownames are the row numbers of matched treated units.
}
\references{
Rosenbaum, Paul R. (2012) "Optimal Matching of an Optimally Chosen Subset in Observational Studies."  Journal of Computational and Graphical Statistics, 21.1, 57-71.
}
\author{
Paul R. Rosenbaum (original forms),
  modifications by Luke Keele and Sam Pimentel
}

