\name{matchMultisens}
\alias{matchMultisens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rosenbaum Bounds after Multilevel Matching
}
\description{
Function to calculate Rosenbaum bounds for continuous outcomes after multilevel matching.
}
\usage{
matchMultisens(obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, Gamma = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A multilevel match object
}
  \item{out.name}{
Outcome variable name
}
  \item{schl_id_name}{
Level 2 ID variable name, that is this variable identifies clusters matched in the data.
}
  \item{treat.name}{
Treatment indicator name
}
  \item{Gamma}{
Sensitivity analysis parameter value. Default is one.
}
}
\details{
This function returns a single p-value, but actually conducts two tests.  The first assumes that the treatment effect does not vary with cluster size.  The second allows the treatment effect to vary with cluster size.  The function returns a single p-value that is corrected for multiple testing.  This p-value is the upper bound for a single Gamma value
}
\value{

\item{pval}{Upper bound on one-sided approximate p-value for test of the sharp null.}

}
\references{
Rosenbaum, Paul R. (2002) Observational Studies. Springer-Verlag.
}
\author{
Luke Keele, Penn State University, ljk20@psu.edu

Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{matchMulti}}, \code{\link{matchMultioutcome}}
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

# Perform sensitivity analysis using Rosenbaum bound -- increase Gamma to increase effect of
# possible hidden confounder 
         
matchMultisens(match.simple, out.name = "mathach",
          schl_id_name = "school", 
          treat.name = "sector", Gamma=1.3)
          
          }