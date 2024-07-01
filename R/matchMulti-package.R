

#' 1980 and 1982 High School and Beyond Data
#'
#' These data are a subset of the data used in Raudenbush and Bryk (1999) for
#' multilevel modeling.
#'
#'
#' @name catholic_schools
#' @docType data
#' @format A \code{data.frame} with 1595 observations on the following
#'   variables.
#'
#'   school: unique school level identifier
#'
#'   ses: student level socio-economic status scale ranges from approx. -3.578
#'   to 2.692
#'
#'   mathach: senior year mathematics test score, outcome measure
#'
#'   female: student level indicator for sex
#'
#'   minority: student level indicator for minority
#'
#'   minority_mean: school level measure of percentage of student body that is
#'   minority
#'
#'   female_mean: school level measure of percentage of student body that is
#'   female
#'
#'   ses_mean: school level measure of average level of student socio-economic
#'   status
#'
#'
#'   sector: treatment indicator 1 if catholic 0 if public
#'
#'   size: school level measure of total number of enrolled students
#'
#'   acad: school level measure of the percentage of students on the academic
#'   track
#'
#'   discrm: school level measure of disciplinary climate ranges from approx.
#'   -2.4 to 2.7
#'
#'   size_large: school level indicator for schools with more than 1000 students
#'
#'   minority_mean_large: school level indicator for schools with more than ten
#'   percent minority
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom coin wilcox_test pvalue
#' @importFrom MASS ginv
#' @importFrom stats  as.formula cov fisher.test mahalanobis model.matrix pnorm
#'   quantile uniroot var binomial glm predict sd
#'
#' @references United States Department of Education. National Center for
#'   Education Statistics.  High School and Beyond, 1980: Sophomore and Senior
#'   Cohort First Follow-Up (1982).
#' @source Raudenbush, S. W. and Bryk, A. (2002).  \emph{Hierarchical Linear
#'   Models: Applications and Data Analysis Methods}.  Thousand Oaks, CA: Sage.
#' @keywords datasets
NULL



#' @title Mini-data set for illustration
#' @name minischool 
#' @docType data
#' @description The Catholic schools dataset subset to a smaller number of
#'   schools (with only 6 Catholic schools).  See full dataset documentation for
#'   more information.
#' @format A data frame with 1500 rows and 12 variables, as described in the
#'   `catholic_schools` dataset.
#' @seealso catholic_schools
#' @source See documentation page for `catholic_schools` dataset.
NULL




#' matchMulti Package
#' 
#' \code{matchMulti} provides and easy to use set of functions to do matching
#' with multilevel data.  It is designed for use with grouped data such as
#' students in schools, where the user wishes to match a set of treated groups
#' to control groups to make the two groups more comparable.
#' 
#' This package will match treated groups to control groups, but allows for
#' trimming of both units and groups to increase balance.  There are also
#' functions for assessing balance after matching, estimating treatment effects
#' and performing sensitivity analysis for hidden confounders.
#' 
#' @name matchMulti-package
#' @docType package
#' 
#' @seealso See also \code{\link{matchMulti}}, \code{\link{matchMultisens}},
#' \code{\link{balanceMulti}}, \code{\link{matchMultioutcome}},
#' \code{\link{rematchSchools}}
#' @keywords matchMulti
#' @examples
#' 
#' \dontrun{
#' # Load Catholic school data
#' data(catholic_schools)
#' 
#' student.cov <- c('minority','female','ses','mathach')
#' 
#' # Check balance student balance before matching
#' balanceTable(catholic_schools[c(student.cov,'sector')],  treatment = 'sector')
#' 
#' #Match schools but not students within schools
#' match.simple <- matchMulti(catholic_schools, treatment = 'sector',
#' school.id = 'school', match.students = FALSE)
#' 
#' #Check balance after matching - this checks both student and school balance
#' balanceMulti(match.simple, student.cov = student.cov)
#' 
#' #Estimate treatment effect
#' output <- matchMultioutcome(match.simple, out.name = "mathach",
#' schl_id_name = "school",  treat.name = "sector")
#' 
#' # Perform sensitivity analysis using Rosenbaum bound -- increase Gamma to increase effect of
#' # possible hidden confounder 
#'          
#' matchMultisens(match.simple, out.name = "mathach",
#'           schl_id_name = "school", 
#'           treat.name = "sector", Gamma=1.3)
#'           
#' # Now match both schools and students within schools          
#' match.out <- matchMulti(catholic_schools, treatment = 'sector',
#' school.id = 'school', match.students = TRUE, student.vars = student.cov)
#' 
#' # Check balance again
#' bal.tab <- balanceMulti(match.out, student.cov = student.cov)
#' 
#' # Now match with fine balance constraints on whether the school is large 
#' # or has a high percentage of minority students
#' match.fb <- matchMulti(catholic_schools, treatment = 'sector', school.id = 'school', 
#' match.students = TRUE, student.vars = student.cov, 
#' school.fb = list(c('size_large'),c('size_large','minority_mean_large')))
#' 
#' # Estimate treatment effects
#' matchMultioutcome(match.fb, out.name = "mathach", schl_id_name = "school",  treat.name = "sector")
#' 
#' #Check Balance
#' balanceMulti(match.fb, student.cov = student.cov)
#' }
#' 
NULL



