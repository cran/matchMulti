\name{sdiff}
\alias{sdiff}
\alias{ttest.balance}
\alias{fisher.balance}
\alias{wilc.balance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Balance Measures
}
\description{
Balance assessment for individual variables, before and after matching
}

\usage{
sdiff(varname, treatment, orig.data, match.data = NULL, 
treat.wts = NULL, ctrl.wts = NULL, mt.wts = NULL, mc.wts = NULL)	
ttest.balance(varname, treatment, orig.data, match.data = NULL,
treat.wts = NULL, ctrl.wts = NULL, mt.wts = NULL, mc.wts = NULL)	
fisher.balance(varname, treatment, orig.data, match.data = NULL,
treat.wts = NULL, ctrl.wts = NULL, mt.wts = NULL, mc.wts = NULL)
wilc.balance(varname, treatment, orig.data, match.data = NULL,
treat.wts = NULL, ctrl.wts = NULL, mt.wts = NULL, mc.wts = NULL)
}
\arguments{
  \item{varname}{
	name of the variable on which to test balance
}
  \item{treatment}{
	name of the binary indicator for treatment status
}
  \item{orig.data}{
	a data frame containing the data before matching
}
  \item{match.data}{
	an optional data frame containing the matched sample
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
}
\details{
	The \code{sdiff} function computes the standardized difference in means.  The other functions perform different kinds of balance tests: \code{t.balance} does the 2-sample t-test, \code{fisher.balance} does Fisher's exact test for binary variable, and \code{wilc.balance} does Wilcoxon's signed rank test.
}
\value{
	a labeled vector.  For \code{sdiff}, the vector has six elements if \code{match.data} is provided: treated and control means and standardized differences before and after matching.  If \code{match.data} is not provided, the vector has only the three elements corresponding to the pre-match case.
	
	For the other functions, if \code{match.data} is provided, the vector contains p-values for the test before and after matching. Otherwise a single p-value is given for the pre-match data.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
