\name{handleNA}
\alias{handleNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Handle Missing Values 
}
\description{
Preprocesses a dataframe of matching covariates so the Mahalanobis distance can be calculated.}
\usage{
handleNA(X, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a matrix or dataframe of covariates to be used for matching
}
  \item{verbose}{
logical value indicating whether detailed output should be provided.
}
}
\details{
Preprocessing involves three main steps: (1) converting factors to matrices of dummy variables (2) for any variable with NAs, adding an additional binary variable indicating whether it is missing (3) imputing all NAs with the column mean.  This follows the recommendations of Rosenbaum in section 9.4 of the referenced text.
}
\value{
a matrix containing the preprocessed data.
}
\references{
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
