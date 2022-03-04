# This is the balanceMulti function.




#' Performs balance checking after multilevel matching.
#'
#' This function checks balance after multilevel balance. It checks balance on
#' both level-one (student) and level-two (school) covariates.
#'
#' This function returns a list which include balance checks for before and
#' after matching for both level-one and level-two covariates. Balance
#' statistics include treated and control means, standardized differences, which
#' is the difference in means divided by the pooled standard deviation before
#' matching, and p-values for mean differences. It extracts the matched data and
#' calls `balanceTable` for student and school level covariates.
#'
#' @param match.obj A multilevel match object
#' @param student.cov Names of student level covariates that you want to check
#'   balance
#' @param school.cov Names of school level covariates for which you want to
#'   check balance, if any.
#' @param include.tests If TRUE include tests for balance.  FALSE just report
#'   the means and differences.
#' @param single.table If FALSE include a list of student and school covariates
#'   separately.  TRUE means single balance table.
#'
#' @return \item{students}{Balance table for student level covariates, as a
#'   dataframe.} \item{schools}{Balance table for school level covariates, as a
#'   dataframe.}
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu} Sam
#'   Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @seealso See also \code{\link{matchMulti}}, \code{\link{matchMultisens}},
#'   \code{\link{matchMultioutcome}}, \code{\link{rematchSchools}}
#' @keywords balance
#' @examples
#'
#' 	\dontrun{
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
#' }
#'
#' @export balanceMulti
balanceMulti <- function(match.obj, student.cov = NULL, school.cov = NULL,
                         include.tests = TRUE,
                         single.table = FALSE ) {
  
	if(is.null(match.obj)){
		 warning('match.obj has value NULL.  Did you make sure optmatch was available before generating matched samples?')
		return(NULL)
	}
	treatment <- match.obj$treatment
	school.id <- match.obj$school.id
	
	# if no student covariates are provided, compute balance on all variables in
	# dataset
	if (is.null(student.cov)) {
	  student.cov = setdiff(colnames(match.obj$raw), c(treatment, match.obj$school.id))
	}
	
	validate_inputs( match.obj$raw, treatment = treatment, 
	                 school.id = match.obj$school.id,
	                 student.vars = student.cov,
	                 school.vars = school.cov )
	
	#student balance
	student.bal <- balanceTable( df.orig = match.obj$raw[c(student.cov, school.id, treatment)], 
	                            df.match = match.obj$matched[c(student.cov, school.id, treatment)], 
	                            treatment, school.id,
	                            include.tests = include.tests )
	
	#school balance
	school.bal <- NULL
	if(!is.null(school.cov)){
		schools.raw <- students2schools(match.obj$raw, c(school.cov, treatment), match.obj$school.id)
		schools.matched <- schools.raw[which(schools.raw[[match.obj$school.id]] %in% as.vector(match.obj$school.match)),]
		
		id.idx <- which(colnames(schools.raw) == match.obj$school.id)
	
		school.bal <- balanceTable(schools.raw, 
		                           schools.matched, 
		                           treatment, school.id,
		                           include.tests = include.tests )
	}
	
	if ( single.table ) {
	  bal = rbind( student.bal, school.bal )
	  return(bal)
	} else {
  	return(list('students' = student.bal, 'schools' = school.bal))
	}
}
