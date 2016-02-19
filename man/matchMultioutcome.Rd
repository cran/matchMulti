\name{matchMultioutcome}
\alias{matchMultioutcome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs an outcome analysis after multilevel matching.
}
\description{
This function returns a point estimate, 95\% confidence interval, and p-values for the matched multilevel data. All results are based on randomization inference.
}
\usage{
matchMultioutcome(obj, out.name = NULL, schl_id_name = NULL, 
treat.name = NULL, end.1 = -1000, end.2 = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A multilevel match object.
}
  \item{out.name}{
Outcome variable name
}
  \item{schl_id_name}{
Level 2 ID variabel name. This variable identifies the clusters in the data that you want to match.
}
  \item{treat.name}{
Treatment variable name, must be zero or one.
}
  \item{end.1}{
Lower bound for point estimate search, default is -1000.
}
  \item{end.2}{
Upper bound for point estimate search, default is 1000.
}
}
\details{
It may be necessary to adjust the lower and upper bounds if one expects the treatment effect confidence interval to be outside
the range of -1000 or 1000.
}
\value{
\item{pval.c }{One-sided approximate p-value for test of the sharp null.}
\item{pval.p }{One-sided approximate p-value for test of the sharp null assuming treatment effects vary with cluster size}
\item{ci1}{Lower bound for 95\% confidence interval.}
\item{ci2}{Upper bound for 95\% confidence interval.}
\item{p.est}{Point estimate for the group level treatment effect.}
}
\references{
Rosenbaum, Paul R. (2002) Observational Studies. Springer-Verlag.
}
\author{
Luke Keele, Penn State University, ljk20@psu.edu

  Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}

\seealso{
See Also as \code{\link{matchMulti}}, \code{\link{matchMultisens}}
}
\examples{
# Load Catholic school data
data(catholic_schools)

student.cov <- c('minority','female','ses','mathach')

# Check balance student balance before matching
balanceTable(catholic_schools[c(student.cov,'sector')],  treatment = 'sector')

#Match schools but not students within schools
match.simple <- matchMulti(catholic_schools, treatment = 'sector',
school.id = 'school', match.students = FALSE)

#Check balance after matching - this checks both student and school balance
balanceMulti(match.simple, student.cov = student.cov)

#Estimate treatment effect
output <- matchMultioutcome(match.simple, out.name = "mathach",
schl_id_name = "school",  treat.name = "sector")

  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ outcomes }
