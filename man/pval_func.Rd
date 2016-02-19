\name{pval_func}
\alias{pval_func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Outcome analysis.
}
\description{
Calcualtes p-values for test of sharp null for treatment effect.
}
\usage{
pval_func(obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, wt = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A multiMatch object
}
  \item{out.name}{
Name of outcome covariate
}
  \item{schl_id_name}{
Name of school (group) identifier
}
  \item{treat.name}{
Name of treatment indicator
}
  \item{wt}{
Logical flag for whether p-value should weight strata by size.
}
}

\value{
A p-value for constant-additive treatment effect.
}
\references{
  Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
  Springer-Verlag.

  Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
  Springer-Verlag.
}
\author{
  Luke Keele, Penn State University, \email{ljk20@psu.edu}

  Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}

\keyword{ internal }