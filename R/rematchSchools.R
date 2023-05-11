#' Repeat School Match Only
#' 
#' After \code{matchMulti} has been called, repeats the school match (with
#' possibly different parameters) without repeating the more computationally
#' intensive student match.
#' 
#' The \code{school.fb} argument encodes a refined covariate balance
#' constraint: the matching algorithm optimally balances the interaction of the
#' variables in the first list element, then attempts to further balance the
#' interaction in the second element, and so on.  As such variables should be
#' added in order of priority for balance.
#' 
#' The \code{keep.target} and \code{school.penalty} parameters allow optimal
#' subset matching within the school match. When the \code{keep.target}
#' argument is specified, the school match is repeated for different values of
#' the \code{school.penalty} parameter in a form of binary search until an
#' optimal match is obtained with the desired number of treated schools or a
#' stopping rule is reached.  The \code{tol} parameter controls the stopping
#' rule; smaller values provide a stronger guarantee of obtaining the exact
#' number of treated schools desired but may lead to greater computational
#' costs.
#' 
#' It is not recommended that users specify the \code{school.penalty} parameter
#' directly in most cases.  Instead the \code{keep.target} parameter provides
#' an easier way to consider excluding schools.
#' 
#' @param match.out an object returned by a call to \code{matchMulti}.
#' @param students a dataframe containing student and school covariates, with a
#' different row for each student.
#' @param school.fb an optional list of character vectors, each containing a
#' subset of the column names of \code{students}.  Each element of the list
#' should contain all the names in previous elements (producing a nested
#' structure).
#' @param verbose a logical value indicating whether detailed output should be
#' printed.
#' @param keep.target an optional numeric value specifying the number of
#' treated schools desired in the final match.
#' @param school.penalty an optional numeric value, treated as the cost (to the
#' objective function in the underlying optimization problem) of excluding a
#' treated school.  If it is set lower, more schools will be excluded.
#' @param tol a numeric tolerance value for comparing distances.  It may need
#' to be raised above the default when matching with many levels of refined
#' balance.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of California, Berkeley, \email{spi@@berkeley.edu}
#' @seealso \code{\link{matchMulti}}.
#' @references Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
#' Springer-Verlag.
#' 
#' Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
#' Springer-Verlag.
#' 
#' Rosenbaum, Paul R. (2012) "Optimal Matching of an Optimally Chosen Subset in
#' Observational Studies."  Journal of Computational and Graphical Statistics,
#' 21.1, 57-71.
#' @examples
#' 
#' \dontrun{
#' # Load Catholic school data
#' data(catholic_schools)
#' 
#' student.cov <- c('minority','female','ses')
#' school.cov <- c('minority_mean','female_mean', 'ses_mean', 'size', 'acad')
#' 
#' #Match schools but not students within schools
#' match.simple <- matchMulti(catholic_schools, treatment = 'sector',
#' school.id = 'school', match.students = FALSE)
#' 
#' #Check balance after matching - this checks both student and school balance
#' balanceMulti(match.simple, student.cov = student.cov, school.cov = school.cov)
#' 
#' #now rematch excluding 2 schools
#' match.trimmed <- rematchSchools(match.simple, catholic_schools, keep.target = 13)
#' match.trimmed$dropped$schools.t
#' }
#' 
#' @export rematchSchools
rematchSchools <-
function(match.out, students, school.fb = NULL, verbose = FALSE, keep.target = NULL, school.penalty = NULL, tol = 1e-3){

	#Are school fine balance constraints nested appropriately?
	if(!is.null(school.fb) && length(school.fb) > 1){
		for(i in c(1:(length(school.fb)-1))) {
			if(!all(school.fb[[i]] %in% school.fb[[i+1]])){
					stop('Each element of school.fb must contain all variables listed in previous elements')
				}	
		}	
	}
	
	if(!is.null(school.penalty) && school.penalty <= 0) stop('School penalty must be positive') 	

	if(is.null(match.out$student.matches)) stop('match.out has no first-stage information: rerun matchMulti with save.first.stage = TRUE')

	treatment <- match.out$treatment
	school.id <- match.out$school.id

	########## REMATCH ###########	

	school.match <- matchSchools(match.out$student.matches$schools.matrix, students, treatment, school.id, school.fb, school.penalty, verbose, tol = tol) 

	#if keep.target is provided, iterate until we get a match that keeps (close to) the desired number of schools
	if(!is.null(keep.target)) {
		treat.schools <- unique(students[[school.id]][students[[treatment]] == 1])
		if (keep.target <= 0 || keep.target > length(treat.schools) || keep.target %% 1 != 0) stop("keep.target must be a positive integer no greater than the total number of treated schools") 
		STARTVAL <- 1000
		MAXITER <- 1000
		SCALE_FACTOR <- 10
		STOPRULE <- 1
		ubound <- Inf
		lbound <- 0
		cur <- school.penalty
		if (is.null(cur)) {
			cur <- STARTVAL	
			school.match <- matchSchools(match.out$student.matches$schools.matrix, students, treatment, school.id, school.fb, penalty = cur, verbose, tol = tol) 
		} 
		next.match <- school.match		
		for (i in 1:MAXITER) {
			nkeep <- length(intersect(treat.schools, next.match[,1]))
			if (nkeep == keep.target) {
				if (verbose) print(paste('Iteration',i,': reached target number of schools with penalty of',cur))
				break
			}
			if (nkeep > keep.target) {
				if (verbose) print(paste('Iteration',i,': too many schools retained with penalty of',cur))			
			ubound <- cur
				cur <- lbound + (ubound - lbound)/2
			}
			if (nkeep < keep.target) {
				if (verbose) print(paste('Iteration',i,': not enough schools retained with penalty of',cur))			
				lbound <- cur
				if (is.finite(ubound)) {
					cur <- lbound + (ubound - lbound)/2
				} else {
					cur <- SCALE_FACTOR*cur
				}
			}
			if(ubound - lbound < STOPRULE) break
			next.match <- matchSchools(match.out$student.matches$schools.matrix, students, treatment, school.id, school.fb, cur, verbose,tol = tol) 
		}
		school.match <- next.match
	}
	
	########### OUTPUT ###########

	out.match <- assembleMatch(match.out$student.matches$student.matches, school.match, school.id, treatment)
	
	#record dropped schools and students
	drop.obj <- list()
	dropped.schools <- setdiff(unique(students[[school.id]]), school.match)
	treat.table <- table(students[[school.id]], students[[treatment]])
	treat.table <- treat.table[match(as.character(dropped.schools),rownames(treat.table)),]
	drop.z <- apply(treat.table, 1, function(x) which(x >0 )-1)
	drop.obj$schools.t <- dropped.schools[drop.z == 1]
	drop.obj$schools.c <- dropped.schools[drop.z == 0]
	
 	student.count.df <- data.frame('school.id' = c(students[[school.id]],out.match[[school.id]]), 'before.after' = c(rep(0,nrow(students)), rep(1, nrow(out.match))))
 	schooltab <- table(student.count.df$school.id, student.count.df$before.after)
 	schooltab[,2] <- schooltab[,1] - schooltab[,2]
 	colnames(schooltab) <- c('Student Count', 'Dropped')
	drop.obj$students.by.school <- schooltab
	
	out.obj <- match.out
	out.obj$matched <- out.match
	out.obj$school.match <- school.match
	out.obj$dropped <- drop.obj
	out.obj
}
