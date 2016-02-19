\name{ci_func}
\alias{ci_func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Outcome analysis.
}
\description{
Calculates confidence interval via grid search.
}
\usage{
ci_func(beta, obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, alpha,
alternative = "less")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
Confidence interval value
}
  \item{obj}{
a multiMatch object
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
  \item{alpha}{
Level of test for confidence interval, default is .05 for 95\% CI.
}
  \item{alternative}{
Direction of test.
}
}
\value{
The endpoint of an estimated confidence interval.
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
