\name{balanceTable}
\alias{balanceTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Create Balance Table
}
\description{
	Given an unmatched sample of treated and control units and (optionally) a matched sample from the same data, produces a table with pre- and post-match measures of covariate balance.
}
\usage{
balanceTable(df.orig, df.match = NULL, treatment, treat.wts = NULL, ctrl.wts = NULL,
mt.wts = NULL, mc.wts = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df.orig}{
	a data frame containing the data before matching
}
  \item{df.match}{
	an optional data frame containing the matched sample
}
  \item{treatment}{
	name of the binary indicator for treatment status
}
  \item{treat.wts}{
	optional weights for treated units in the original sample
}
  \item{ctrl.wts}{
	optional weights for control units in the original sample
	}
  \item{mt.wts}{

	optional weights for treated units in the matched sample
}
  \item{mc.wts}{

	optional weights for treated units in the matched sample
}
  \item{verbose}{
	a logical value indicating whether detailed output should be printed.
}
}

\value{
A matrix of balance measures, with one row for each covariate in \code{df.orig} except \code{treatment}, and columns for treated and control means, standardized differences in means, p-values from a 2-sample t-test, and p-values from either Fisher's exact test (if the covariate is binary) or a Wilcoxon signed rank test otherwise.  If \code{df.match} is specified there are twice as many columns, one set for the pre-match samples and one set for the post-match samples.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\keyword{internal}
