#' A function that performs multilevel matching.
#'
#' This is the workhorse function in the package which matches groups and units
#' within groups. For example, it will match both schools and students in
#' schools, where the goal is to make units more comparable to estimate
#' treatment effects.
#'
#' \code{matchMulti} first matches students (or other individual units) within
#' each pairwise combination of schools (or other groups); based on these
#' matches a distance matrix is generated for the schools.  Then schools are
#' matched on this distance matrix and the student matches for the selected
#' school pairs are combined into a single matched sample.
#'
#' School covariates are not used to compute the distance matrix for schools
#' (since it is generated from the student match).  Instead imbalances in school
#' covariates should be addressed through the\code{school.fb} argument, which
#' encodes a refined covariate balance constraint. School covariates in
#' \code{school.fb} should be given in order of priority for balance, since the
#' matching algorithm optimally balances the variables in the first list
#' element, then attempts to further balance the those in the second element,
#' and so on.
#'
#' @param data A data frame for use in matching.
#' @param treatment Name of covariate that defines treated and control groups.
#' @param school.id Identifier for groups (for example schools)
#' @param match.students Logical flag for whether units within groups should
#'   also be matched.  If set to \code{FALSE}, all units will be retained in
#'   both groups.
#' @param student.vars Names of student level covariates on which to measure
#'   balance.  School-level distances will be penalized when student mathces are
#'   imbalanced on these variables. In addition, when \code{match.students} is
#'   \code{TRUE}, students are matched on a distance computed from these
#'   covariates.
#' @param school.caliper matrix with one row for each treated school and one
#'   column for each control school, containing zeroes for pairings allowed by
#'   the caliper and \code{Inf} values for forbidden pairings.  When \code{NULL}
#'   no caliper is imposed.
#' @param school.fb A list of discrete group-level covariates on which to
#'   enforce fine balance, i.e., ensure marginal distributions are balanced.
#'   First group is most important, second is second most, etc.  If a simple
#'   list of variable names, one group is assumed.  A list of list will give
#'   this hierarchy.
#' @param verbose Logical flag for whether to give detailed output.
#' @param keep.target an optional numeric value specifying the number of treated
#'   schools desired in the final match.
#' @param student.penalty.qtile This helps exclude students if they are
#'   difficult to match. Default is 0.05, which implies that in the match we
#'   would prefer to exclude students rather than match them at distances larger
#'   than this quantile of the overall student-student robust Mahalanobis
#'   distance distribution
#' @param min.keep.pctg Minimum percentage of students (from smaller school) to
#'   keep when matching students in each school pair.
#' @param school.penalty A penalty to remove groups (schools) in the group
#'   (school) match
#' @param save.first.stage Should first stage matches be saved.
#' @param tol a numeric tolerance value for comparing distances, used in the
#'   school match.  It may need to be raised above the default when matching 
#'   with many levels of refined balance or in very large problems (when these
#'   distances will often be  at least on the order of the tens of 
#'   thousands).
#' @param solver Name of package used to solve underlying network flow problem
#'   for the school match, one of 'rlemon' and 'rrelaxiv'.  rrelaxiv carries an 
#'   academic license and is not hosted on CRAN so it must be installed 
#'   separately.
#' @return \item{raw}{The unmatched data before matching.} \item{matched}{The
#'   matched dataset of both units and groups. Outcome analysis and balance
#'   checks are peformed on this item.} \item{school.match}{Object with two
#'   parts. The first lists which treated groups (schools) are matched to which
#'   control groups. The second lists the population of groups used in the
#'   match.} \item{school.id}{Name of school identifier} \item{treatment}{Name
#'   of treatment variable}
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#'
#'   Sam Pimentel, University of California, Berkeley, \email{spi@@berkeley.edu}
#' @seealso See also \code{\link{matchMulti}}, \code{\link{matchMultisens}},
#'   \code{\link{balanceMulti}}, \code{\link{matchMultioutcome}},
#'   \code{\link{rematchSchools}}
#' @examples
#'
#'
#' #toy example with short runtime
#' library(matchMulti)
#'
#' #Load Catholic school data
#' data(catholic_schools)
#'
#' # Trim data to speed up example
#' catholic_schools <- catholic_schools[catholic_schools$female_mean >.45 &
#'  catholic_schools$female_mean < .60,]
#'
#' #match on a single covariate
#' student.cov <- c('minority')
#'
#' \dontshow{if (requireNamespace("optmatch", quietly = TRUE))} match.simple <- 
#' matchMulti(catholic_schools, treatment = 'sector',
#'                              school.id = 'school', match.students = FALSE,
#'                              student.vars = student.cov, verbose=TRUE, tol=.01)
#'
#' #Check balance after matching - this checks both student and school balance
#' \dontshow{if (requireNamespace("optmatch", quietly = TRUE))}  balanceMulti(match.simple, student.cov = student.cov)
#'
#'
#' \dontrun{
#' #larger example
#' data(catholic_schools)
#'
#' student.cov <- c('minority','female','ses')
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
#' matchMultisens(match.simple, out.name = "mathach",
#'           schl_id_name = "school",
#'           treat.name = "sector", Gamma = 1.3)
#'
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
#' school.fb = list( c('size_large'), c('minority_mean_large') )
#'
#' # Estimate treatment effects
#' matchMultioutcome(match.fb, out.name = "mathach", schl_id_name = "school",  treat.name = "sector")
#'
#' #Check Balance
#' balanceMulti(match.fb, student.cov = student.cov)
#'
#' }
#'
#'
#' @export matchMulti
matchMulti <- function(data, treatment, school.id, match.students = TRUE, 
                       student.vars = NULL, school.caliper = NULL,  school.fb = NULL, 
                       verbose = FALSE, keep.target = NULL, student.penalty.qtile = 0.05,
                       min.keep.pctg = 0.8, school.penalty = NULL, save.first.stage = TRUE, 
                       tol = 1e1, solver = 'rlemon'){
	
	students <- data
	
	##### Validate input #####
	
  validate_inputs( students, treatment, school.id, 
                   student.vars = student.vars, school.vars = school.fb )

	#do schools nest within treatment categories?
	treat.tab <- table(students[[school.id]], students[[treatment]]) 
	if(any(apply(treat.tab, 1, min) > 0)) {
		stop('Some schools contain both treated and control students')
	}
	
	# Nest the school fine balance constraints, if there are multiple lists
	if(!is.null(school.fb)) {
	  if ( is.character( school.fb ) ) {
	    school.fb = list( school.fb )
	  } else if ( length(school.fb) > 1){
	    for(i in c(1:(length(school.fb)-1))) {
	      school.fb[[i+1]] = union( school.fb[[i]], school.fb[[i+1]] )
	      #			if(!all(school.fb[[i]] %in% school.fb[[i+1]])){
	      #					stop('Each element of school.fb must contain all variables listed in previous elements')
	      #				}	
	    }
	  }
	}
	
	if(student.penalty.qtile < 0 || student.penalty.qtile > 1)
	  stop('Student penalty quantile must be in [0,1]') 

	if(!is.null(school.penalty) && school.penalty <= 0) 
	  stop('School penalty must be positive') 	
	
	if(match.students && is.null(student.vars)) 
	  stop('Cannot match students unless student variables are specified')
	
	
	########### MATCHING ###########
		
	#student matches in all school pairings	
	student.matches <- matchStudents(students, treatment, school.id, match.students, student.vars,
	                                 school.caliper, verbose, student.penalty.qtile, min.keep.pctg)
	if(is.null(student.matches)) return(NULL)

	school.match <- matchSchools(student.matches$schools.matrix, students, treatment, 
	                             school.id, school.fb, school.penalty, verbose, 
	                             tol = tol, solver = solver) 

	#if keep.target is provided, iterate until we get a match that keeps (close to) the desired number of schools
	if( nrow(school.match) > 0 && !is.null(keep.target)) {
		treat.schools <- unique(students[[school.id]][students[[treatment]] == 1])
		if (keep.target <= 0 || keep.target > length(treat.schools) || keep.target %% 1 != 0) {
		  stop("keep.target must be a positive integer no greater than the total number of treated schools") 
		}
		STARTVAL <- 1000
		MAXITER <- 1000
		SCALE_FACTOR <- 10
		ubound <- Inf
		lbound <- 0
		cur <- school.penalty
		if (is.null(cur)) {
			cur <- STARTVAL	
			school.match <- tryCatch({
			  matchSchools(student.matches$schools.matrix, students, treatment, 
			                             school.id, school.fb, penalty = cur, verbose, 
			                             tol = tol, solver = solver) 
			}, warning = function(cond){
			  if(cond$message == 'No matched pairs formed. Try using a higher exclusion penalty?'){
			    return.val <- suppressWarnings(matchSchools(student.matches$schools.matrix, students, treatment, 
			                                                school.id, school.fb, cur, verbose,tol = tol,
			                                                solver = solver) )
			  }else{
			    return.val <- matchSchools(student.matches$schools.matrix, students, treatment, 
			                               school.id, school.fb, cur, verbose,tol = tol,
			                               solver = solver) 
			  }
			  return(return.val)
			})
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
					if(is.na(as.integer((SCALE_FACTOR^2)*cur))){
						print('Cannot increase penalty further, terminating search')
						break
					}
					cur <- SCALE_FACTOR*cur
				}
			}
			if(ubound - lbound < tol/2){
			  print('Penalty cannot be varied at a finer resolution, terminating search.
			        Consider using lower tolerance.')
			  break 
			}
			next.match <- tryCatch({
			  matchSchools(student.matches$schools.matrix, students, treatment, 
			                           school.id, school.fb, cur, verbose,tol = tol,
			                           solver = solver) 
			  }, warning = function(cond){
			    if(cond$message == 'No matched pairs formed. Try using a higher exclusion penalty?'){
			      return.val <- suppressWarnings(matchSchools(student.matches$schools.matrix, students, treatment, 
			                                                  school.id, school.fb, cur, verbose,tol = tol,
			                                                  solver = solver) )
			    }else{
			     return.val <- matchSchools(student.matches$schools.matrix, students, treatment, 
			                                      school.id, school.fb, cur, verbose,tol = tol,
			                                      solver = solver) 
			    }
			    return(return.val)
			  })
		}
		school.match <- next.match
	}


	########### OUTPUT ###########
	
	if(nrow(school.match) == 0) stop('No schools matched. Try lower school exclusion penalty?')

	out.match <- assembleMatch(student.matches$student.matches, school.match, school.id, treatment)
	
	#record dropped schools and students
	drop.obj <- list()
	dropped.schools <- setdiff(unique(students[[school.id]]), school.match)
	treat.table <- table(students[[school.id]], students[[treatment]])
	treat.table <- treat.table[match(as.character(dropped.schools),rownames(treat.table)),]
	drop.z <- apply(treat.table, 1, function(x) which(x >0 )-1)
	drop.obj$schools.t <- dropped.schools[drop.z == 1]
	drop.obj$schools.c <- dropped.schools[drop.z == 0]
	
 	student.count.df <- data.frame('school.id' = c(students[[school.id]],out.match[[school.id]]), 
 	                               'before.after' = c(rep(0,nrow(students)), rep(1, nrow(out.match))))
 	schooltab <- table(student.count.df$school.id, student.count.df$before.after)
 	schooltab[,2] <- schooltab[,1] - schooltab[,2]
 	colnames(schooltab) <- c('Student Count', 'Dropped')
	drop.obj$students.by.school <- schooltab
	
	out.obj <- list('raw' = students, 
	                'matched' = out.match, 
	                'school.match' = school.match, 
	                'dropped' = drop.obj, 
	                'school.id' = school.id, 
	                'treatment' = treatment)
	
	if(save.first.stage){ 
		out.obj$student.matches <- student.matches
	}
	
	class( out.obj ) <- "matchMultiResult"
	
	out.obj$student.vars = student.vars
	out.obj$school.fb = school.fb
	
	out.obj
}

