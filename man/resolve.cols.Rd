\name{resolve.cols}
\alias{resolve.cols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Ensure Dataframes Share Same Set Columns
}
\description{
Takes in two dataframes.  For each column name that is in the second frame but not in the first frame, a new column of zeroes is added to the first frame.
}
\usage{
resolve.cols(df1, df2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df1}{
	a dataframe.
}
  \item{df2}{
	a dataframe.
}
}
\value{
	a dataframe
}
\author{
  Luke Keele, Penn State University, \email{ljk20@psu.edu}

  Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}

\keyword{internal }
