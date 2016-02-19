\name{smahal}
\alias{smahal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Mahalanobis Distance
}
\description{
Computes robust Mahalanobis distance between treated and control units.
}
\usage{
smahal(z, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
	vector of treatment indicators (1 for treated, 0 for controls).
}
  \item{X}{
	matrix of numeric variables to be used for computing the Mahalanobis distance.  Row count must match length of \code{z}.
}
}
\details{
For an explanation of the robust Mahalanobis distance, see section 8.3 of the first reference.  This function was written by Paul Rosenbaum and distributed in the supplemental material to the second reference.
}
\value{
a matrix of robust Mahalanobis distances, with a row for each treated unit and a column for each control.
}
\references{	
Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
  Springer-Verlag.	
	
Rosenbaum, Paul R. (2012) "Optimal Matching of an Optimally Chosen Subset in Observational Studies."  Journal of Computational and Graphical Statistics, 21.1, 57-71.
}
\author{
Paul R. Rosenbaum.
}
\keyword{internal}