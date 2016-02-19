\name{agg}
\alias{agg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract School-Level Covariates 
}
\description{
	Given a vector of variables of interest for students in a single school, extracts a single value for the school
	}
\usage{
agg(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	a vector containing student-level observations for a school. If it is a factor it must contain only a single level.
}
}
\details{
If the input is numeric, \code{agg} returns the mean; if the input is not numeric, an error will be thrown unless all values are the same, in which case the single unique value will be returned.
}
\value{
	A single value of the same type as the input vector.
}
\author{
  Luke Keele, Penn State University, \email{ljk20@psu.edu}

  Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
