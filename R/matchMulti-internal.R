#' Extract School-Level Covariates
#' 
#' Given a vector of variables of interest for students in a single school,
#' extracts a single value for the school
#' 
#' If the input is numeric, \code{agg} returns the mean; if the input is not
#' numeric, an error will be thrown unless all values are the same, in which
#' case the single unique value will be returned.
#' 
#' @param x a vector containing student-level observations for a school. If it
#' is a factor it must contain only a single level.
#' @return A single value of the same type as the input vector.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export agg
agg <-
function(x){
	if(length(unique(x)) == 1) return(x[1])
	return(mean(x, na.rm = TRUE))	
}



#' Collect Matched Samples
#' 
#' After students and schools have both been matched separately, assembles the
#' matched student samples corresponding to the school match into a single
#' dataframe of student-level data.
#' 
#' 
#' @param student.matches a list of lists object produced by
#' \code{matchStudents}, with each element of the second list containing a
#' dataframe composed of a matched sample for a different treated-control
#' school pairing.
#' @param school.match a dataframe, produced by \code{matchSchools}, with two
#' columns, one containing treated school IDs and the other containing matched
#' control school IDs.
#' @param school.id the name of the column storing the unique school identifier
#' (in the dataframes stored in \code{student.matches})
#' @param treatment the name of the column storing the binary treatment status
#' indicator (in the dataframes stored in \code{student.matches})
#' @return a dataframe containing the full set of matched samples for the
#' multilevel match.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export assembleMatch
assembleMatch <-
function(student.matches, school.match, school.id, treatment){
	#school.match <- school.match.list$matches[,1,drop = FALSE]
	#final.match <- student.matches[[1]][[1]]
	#final.match$pair.id <- rep(NA, nrow(final.match))
	#final.match <- final.match[-c(1:nrow(final.match)),]
	final.match <- NULL
	for(i in 1:nrow(school.match)){
		bind.obj <- student.matches[[school.match[i,1]]][[school.match[i,2]]]
		bind.obj$pair.id <- i
		if (is.null(final.match)) {
			final.match <- bind.obj
		} else {
			final.match <- rbind(final.match, bind.obj)			
		}
	}	
	return(final.match)
}

#' @export
#'
#' @rdname pairmatchelastic
elastic <- function (mdist, n = 0, val = 0) {
    st <- max(as.numeric(c(rownames(mdist), colnames(mdist))))
    h <- matrix(val, dim(mdist)[1], n)
    if(n > 0){
	    colnames(h) <- (st + 1):(st + n)
	}
    cbind(mdist, h)
}




fisher.balance <- function(varname, treatment, orig.data, match.data = NULL){
		
	orig.v <- orig.data[,which(colnames(orig.data) == varname)]
	
	tab.orig <- table(orig.v, orig.data[[treatment]])
	if(nrow(tab.orig) == 1) {
		t.orig <- list('p.value' = 1)
	} else {
		t.orig <- fisher.test(tab.orig)	
	}
	
	if(is.null(match.data)) return(c('Fisher Test Pvalue' = t.orig$p.value))

	match.v  <- match.data[,which(colnames(match.data) == varname)]	
	
	tab.match <- table(match.v, match.data[[treatment]])
	if(nrow(tab.match) == 1) {
		t.match <- list('p.value' = 1)
	} else {
		t.match <- fisher.test(tab.match)	
	}
	return(c('Fisher Test Pvalue Before' = t.orig$p.value, 'Fisher Test Pvalue After' = t.match$p.value))
}




#' Handle Missing Values
#' 
#' Preprocesses a dataframe of matching covariates so the Mahalanobis distance
#' can be calculated.
#' 
#' Preprocessing involves three main steps: (1) converting factors to matrices
#' of dummy variables (2) for any variable with NAs, adding an additional
#' binary variable indicating whether it is missing (3) imputing all NAs with
#' the column mean.  This follows the recommendations of Rosenbaum in section
#' 9.4 of the referenced text.
#' 
#' @param X a matrix or dataframe of covariates to be used for matching
#' @param verbose logical value indicating whether detailed output should be
#' provided.
#' @return a matrix containing the preprocessed data.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @references Rosenbaum, Paul R. (2010). \emph{Design of Observational
#' Studies}.  Springer-Verlag.
#' @keywords internal
#' @importFrom plyr laply
#' @export handleNA
handleNA <-
function(X, verbose = FALSE){
	if (is.data.frame(X)) {
        X.chars <- which(plyr::laply(X, class) == "character")
        if (verbose && length(X.chars) > 0) {
            print("character variables found in X, converting to factors")
            for (i in X.chars) {
                X[, i] <- factor(X[, i])
            }
        }
        X.factors <- which(plyr::laply(X, class) == "factor")
        for (i in which(plyr::laply(X, function(x) any(is.na(x))))) {
			if (verbose) print(paste('Missing values found in variable ', colnames(X)[i] ,'; imputing and adding missingness indicator'))        	
            if (i %in% X.factors) {
                X[, i] <- addNA(X[, i])
            }
            else {
                X[[paste(colnames(X)[i], "NA", sep = "")]] <- is.na(X[, 
                  i])
                X[which(is.na(X[, i])), i] <- mean(X[, i], na.rm = TRUE)
            }
        }
        for (i in rev(X.factors)) {
            X <- cbind(X[, -i], model.matrix(as.formula(paste("~", 
                colnames(X)[i], "-1")), data = X))
        }
    } else {
		if(is.null(colnames(X))) colnames(X) <- 1:ncol(X)
        for (i in c(1:ncol(X))) {
            if (any(is.na(X[, i]))) {
                X <- cbind(X, is.na(X[, i]))
                colnames(X)[ncol(X)] <- paste(colnames(X)[i], 
                  "NA", sep = "")
                X[which(is.na(X[, i])), i] <- mean(X[, i], na.rm = TRUE)
            }
        }
    }
	X
}




#' Check if a variable is binary
#' 
#' Examines a vector that is not coded as a logical to see if it contains only
#' 0s and 1s.
#' 
#' 
#' @param x A vector.
#' @return a logical value, \code{TRUE} if the vector contains only 0s and 1s
#' and \code{FALSE} otherwise.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export is.binary
is.binary <-
function(x){
	if (max(x) != 1 || min(x) != 0) return(FALSE)
	if (length(unique(x)) == 2) return(TRUE)
	#for now leave things with imputed means continuous
	#if (length(unique(x)) == 3){
	#	other.val <- setdiff(unique(x),c(0,1))
		#check if this other value is an imputed mean (due to NAs) - if so round and treat as binary	
	#}
	return(FALSE)
}




#' Compute School Distance from a Student Match
#' 
#' Defines a distance between two schools whose students have been matched
#' based on the size of the resulting matched sample and on the student-level
#' covariate balance.
#' 
#' The distance is computed by (1) subtracting the harmonic mean of the treated
#' and control counts in the matched sample from \code{largeval} (2) adding
#' \code{largeval} for each covariate among \code{studentvars} that has an
#' absolute standardized difference exceeding 0.2.  This encourages the school
#' match to choose larger schools with better balance.
#' 
#' @param matchFrame dataframe containing all matched students.
#' @param treatFrame dataframe containing all students from the treated school.
#' @param ctrlFrame dataframe containing all students from the control school.
#' @param student.vars names of variables on which to evaluate balance in the
#' matched sample.  Must be present in the column names of each of
#' \code{matchFrame}, \code{treatFrame} and \code{ctrlFrame}.
#' @param treatment name of the treatment variable. Must be present in the
#' column names of each of \code{matchFrame}, \code{treatFrame} and
#' \code{ctrlFrame}.
#' @param largeval a large penalty value to be added to the distance for each
#' student-level imbalance.
#' @return a numeric distance.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export match2distance
match2distance <-
function(matchFrame, treatFrame, ctrlFrame, student.vars, treatment, largeval){
	treat.v <- matchFrame[[treatment]]
	#compute harmonic mean of sample sizes and invert it
	harm.mean <- 2/(1/sum(treat.v) + 1/sum(1-treat.v) )
	#check balance and count imbalances
	if (!is.null(student.vars)) {
		sdiff.student <- sapply(student.vars, function(x) sdiff(x, treatment = treatment, orig.data= rbind(treatFrame, ctrlFrame), match.data = matchFrame)[6])
		return(largeval - harm.mean + largeval*sum(abs(sdiff.student) > 0.2))
	} else {
		return(largeval - harm.mean)		
	}
}



#' Match Schools on Student-based Distance
#' 
#' Takes in a school distance matrix created using information from the
#' first-stage student match and matches schools optimally, potentially
#' 
#' The \code{school.fb} argument encodes a refined covariate balance
#' constraint: the matching algorithm optimally balances the interaction of the
#' variables in the first list element, then attempts to further balance the
#' interaction in the second element, and so on.  As such variables should be
#' added in order of priority for balance.
#' 
#' @param dmat a distance matrix for schools, with a row for each treated
#' school and a column for each control school.
#' @param students a dataframe containing student and school covariates, with a
#' different row for each student.
#' @param treatment the column name of the binary treatment status indicator in
#' the \code{students} dataframe.
#' @param school.id the column name of the unique school ID in the
#' \code{students} dataframe.
#' @param school.fb an optional list of character vectors, each containing a
#' subset of the column names of \code{students}.  Each element of the list
#' should contain all the names in previous elements (producing a nested
#' structure).
#' @param penalty a numeric value, treated as the cost to the objective
#' function of excluding a treated school.  If it is set lower, more schools
#' will be excluded.
#' @param verbose a logical value indicating whether detailed output should be
#' printed.
#' @param tol a numeric tolerance value for comparing distances.  It may need
#' to be raised above the default when matching with many levels of refined
#' balance.
#' @return a dataframe with two columns, one containing treated school IDs and
#' the other containing matched control school IDs.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @importFrom plyr ddply
#' @importFrom rcbsubset rcbsubset
#' @keywords internal
#' @export matchSchools
matchSchools <-
function(dmat, students, treatment, school.id,
 school.fb, penalty, verbose, tol){
 	#DIST_TOLERANCE <- 1e-3
	
	balance.covs <- unique(c(unlist(school.fb)))
	school.df <- students2schools(students, c(balance.covs, treatment), school.id)
	
	#reorder school.df to ensure schools are in same order as implied by dimnames on dmat
	reord <- match(unlist(dimnames(dmat)), as.character(school.df[[school.id]]))
	if (any(is.na(reord))) stop('matchSchools found schools not handled by matchStudents')
	school.df <- school.df[reord,,drop = FALSE]
	if(is.null(school.fb)){
		match.out <- rcbsubset::rcbsubset(dmat,
		                       exclude.penalty = penalty, tol = tol) #	DIST_TOLERANCE)
	} else {
		match.out <- rcbsubset::rcbsubset(dmat, fb.list = school.fb, 
		                       treated.info = school.df[school.df[[treatment]] == 1,
		                                                ,drop = FALSE], 
		                       control.info = school.df[school.df[[treatment]] == 0,
		                                                ,drop = FALSE], 
		                       exclude.penalty = penalty, tol = tol)# tol =  	DIST_TOLERANCE)
	}
	out.frame <- cbind('TreatID' = school.df[[school.id]][which(school.df[[treatment]] == 1)][as.numeric(rownames(match.out$matches))], 'CtrlID' = school.df[[school.id]][which(school.df[[treatment]] == 0)][match.out$matches])
	
	return(out.frame)
}

#NB: when empty matches are formed between students, an NA goes into the distance
#matrix.  This works fine downstream because rcbsubset trims out NAs from
#distance matrix as forbidden matches because they are not finite.


#' Compute Student Matches for all Pairs of Schools
#' 
#' Iterates over all possible treated-control school pairs, optionally computes
#' and stores an optimal student match for each one, and generates a distance
#' matrix for schools based on the quality of each student match.
#' 
#' The \code{penalty.qtile} and \code{min.keep.pctg} control the rate at which
#' students are trimmed from the match.  If the quantile is high enough no
#' students should be excluded in any match; if the quantile is very low the
#' \code{min.keep.pctg} can still ensure a minimal sample size in each match.
#' 
#' @param students a dataframe containing student covariates, with a different
#' row for each student.
#' @param treatment the column name of the binary treatment status indicator in
#' the \code{students} dataframe.
#' @param school.id the column name of the unique school ID in the
#' \code{students} dataframe.
#' @param match.students logical value.  If \code{TRUE}, students are matched
#' within school pairs and some students will be excluded.  If \code{FALSE},
#' all students will be retained in the matched sample for each school pair.
#' @param student.vars column names of variables in \code{students} on which to
#' match students and assess balance of student matches in evaluating match
#' quality.
#' @param school.caliper matrix with one row for each treated school and one
#' column for each control school, containing zeroes for pairings allowed by
#' the caliper and \code{Inf} values for forbidden pairings.  When \code{NULL}
#' no caliper is imposed.
#' @param verbose a logical value indicating whether detailed output should be
#' printed.
#' @param penalty.qtile a numeric value between 0 and 1 specifying a quantile
#' of the distribution of all student-student matching distances.  The
#' algorithm will prefer to exclude treated students rather than form pairs
#' with distances exceeding this quantile.
#' @param min.keep.pctg a minimum percentage of students in the smaller school
#' in a pair which must be retained, even when treated students are excluded.
#' @return A list with two elements: \item{student.matches}{ a list with one
#' element for each treated school.  Each element is a list with one element
#' for each control school, and each element of these secondary lists is a
#' dataframe containing the matched sample for the corresponding
#' treated-control pairing. } \item{schools.matrix}{ a matrix with one row for
#' each treated school and one column for each control school, giving matching
#' distances based on the student match. }
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export matchStudents
matchStudents <-
function(students, treatment, school.id, match.students, student.vars, school.caliper = NULL, verbose, penalty.qtile, min.keep.pctg){	
	if(match.students){
		X <- students[student.vars]
		X <- handleNA(X, verbose = verbose)
		varying <- apply(X,2, function(x) length(unique(x)) > 1)
		if(!all(varying) && verbose){
			print('Constant-value variables found, they will not be used to calculate Mahalanobis distance.')
		} 
		X <- X[,which(varying),drop = FALSE]
		#TODO: examine performance of this step, consider optimization
		maha.mat <- smahal(students[[treatment]],X)
		#TODO: decide how to handle calipers
		student.penalty <- quantile(maha.mat, penalty.qtile)
	}
	
	#TODO: add documentation to show treatment should be 0/1
	t.schools <- unique(students[[school.id]][students[[treatment]]==1])
	c.schools <- unique(students[[school.id]][students[[treatment]]==0]) 
	nt.schools <- length(t.schools)
	nc.schools <- length(c.schools)
	schoolmatch.mat <- matrix(nrow = nt.schools, ncol = nc.schools)
	matches.list <- list()
	for (i in 1:nt.schools) {
		matches.list[[t.schools[i]]] <- list()
		schli <- students[students[[school.id]] == t.schools[i],,drop = FALSE]
 		for (j in 1:nc.schools) {
			schlj <- students[students[[school.id]] == c.schools[j],,drop = FALSE] 
			if (match.students) {
				#grab correct slice of maha distance
				treat.idx <- c(students[[school.id]] == t.schools[i])[students[[treatment]] == 1]
				ctrl.idx <- c(students[[school.id]] == c.schools[j])[students[[treatment]] == 0]
				maha.slice <- maha.mat[which(treat.idx), which(ctrl.idx),drop = FALSE]
				min.keep <- floor((1-min.keep.pctg)*min(dim(maha.slice)))
				match.out <- pairmatchelastic(maha.slice, n = min.keep, val = student.penalty)
				if(is.null(match.out)) return(NULL)
				matches.list[[t.schools[i]]][[c.schools[j]]] <- rbind(schli[as.numeric(rownames(match.out)),,drop = FALSE], schlj[match.out,,drop = FALSE])
			} else {
				matches.list[[t.schools[i]]][[c.schools[j]]] <- rbind(schli,schlj)
			}
			schoolmatch.mat[i,j] <- match2distance(matches.list[[t.schools[i]]][[c.schools[j]]], schli, schlj, student.vars, treatment, largeval = nrow(students)/2)
		}
	}
	if(!is.null(school.caliper)) schoolmatch.mat <- schoolmatch.mat + school.caliper

	rownames(schoolmatch.mat) <- t.schools
	colnames(schoolmatch.mat) <- c.schools

	return(list('student.matches' = matches.list, 'schools.matrix' = schoolmatch.mat))
}





#' Optimal Subset Matching without Balance Constraints
#'
#' Conducts optimal subset matching as described in the reference.
#'
#' \code{pairmatchelastic} is the main function, which conducts an entire match.
#' \code{elastic} is a helper function which augments the original distance
#' matrix as described in the reference. 
#' 
#' The original versions of these functions were written by Paul Rosenbaum and
#' distributed in the supplemental material to the paper: "Optimal Matching of
#' an Optimally Chosen Subset in Observational Studies," Paul R. Rosenbaum,
#' Journal of Computational and Graphical Statistics, Vol. 21, Iss. 1, 2012.
#'
#' @param mdist distance matrix with rows corresponding to treated units and
#'   columns corresponding to controls.
#' @param n maximum number of treated units that can be excluded.
#' @param val cost of excluding a treated unit (i.e. we prefer to exclude a
#'   treated unit if it increases the total matched distance by more than
#'   \code{val}).
#' @return \code{elastic} returns an augmented version of the input matrix
#'   \code{mdist}.  \code{pairmatchelastic} returns a matrix of 1 column whose
#'   values are the column numbers of matched controls and whose rownames are
#'   the row numbers of matched treated units.
#' @author Paul R. Rosenbaum (original forms), modifications by Luke Keele and
#'   Sam Pimentel
#'   
#' @references Rosenbaum, Paul R. (2012) "Optimal Matching of an Optimally
#'   Chosen Subset in Observational Studies."  Journal of Computational and
#'   Graphical Statistics, 21.1, 57-71.
#'   
#' @export pairmatchelastic
pairmatchelastic <- function (mdist, n = 0, val = 0) {
    ro <- dim(mdist)[1]
    co <- dim(mdist)[2]
    k <- ro + co
    mdist <- elastic(mdist, n = n, val = val)
    m <- rcbsubset::rcbsubset(mdist)
	if(is.null(m)) return(NULL)
	mt <- m
	mt$matches <- m$matches[which(m$matches <= co),,drop=FALSE]
    mt$matches
}




#' Ensure Dataframes Share Same Set Columns
#' 
#' Takes in two dataframes.  For each column name that is in the second frame
#' but not in the first frame, a new column of zeroes is added to the first
#' frame.
#' 
#' 
#' @param df1 a dataframe.
#' @param df2 a dataframe.
#' @return a dataframe
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export resolve.cols
resolve.cols <-
function(df1, df2){
	new.cols <- setdiff(colnames(df2), colnames(df1))
	for(newc in new.cols){
		df1[[newc]] <- rep(0, nrow(df1))
	}
	df1
}



#' Balance Measures
#'
#' Balance assessment for individual variables, before and after matching
#'
#' The \code{sdiff} function computes the standardized difference in means. The
#' other functions perform different kinds of balance tests: \code{t.balance}
#' does the 2-sample t-test, \code{fisher.balance} does Fisher's exact test for
#' binary variable, and \code{wilc.balance} does Wilcoxon's signed rank test.
#'
#' @aliases sdiff ttest.balance fisher.balance wilc.balance
#' @param varname name of the variable on which to test balance
#' @param treatment name of the binary indicator for treatment status
#' @param orig.data a data frame containing the data before matching
#' @param match.data an optional data frame containing the matched sample
#'
#' @return a labeled vector.  For \code{sdiff}, the vector has six elements if
#'   \code{match.data} is provided: treated and control means and standardized
#'   differences before and after matching.  If \code{match.data} is not
#'   provided, the vector has only the three elements corresponding to the
#'   pre-match case.
#'
#'   For the other functions, if \code{match.data} is provided, the vector
#'   contains p-values for the test before and after matching. Otherwise a
#'   single p-value is given for the pre-match data.
#' 
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#'
#'   Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @references Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
#'   Springer-Verlag.
#'
#'   Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
#'   Springer-Verlag.
#' @keywords internal
#' @export sdiff
sdiff <-
function(varname, treatment, orig.data, match.data = NULL) {
		
	orig.v <- orig.data[,which(colnames(orig.data) == varname)]
	m1 <- mean(orig.v[orig.data[[treatment]] == 1], na.rm = TRUE)
	m0 <- mean(orig.v[orig.data[[treatment]] == 0], na.rm = TRUE)
	
	s1 <- sqrt(var(orig.v[orig.data[[treatment]] == 1], na.rm= TRUE))
	if(length(orig.v[orig.data[[treatment]] == 1]) == 1) s1 <- 0
	s0 <- sqrt(var(orig.v[orig.data[[treatment]] == 0], na.rm = TRUE))
	if(length(orig.v[orig.data[[treatment]] == 0]) == 1) s0 <- 0
	sdiff.before <- (m1 -m0)/sqrt(0.5*(s1^2 + s0^2))
	if(s1 == 0 && s0 == 0) sdiff.before <- ifelse(m1 == m0, 0, Inf)
	if(is.null(match.data)) return(c('Treated Mean' = m1, 'Control Mean' = m0,'SDiff' = sdiff.before))
		
	#if match data is also provided
	match.v  <- match.data[,which(colnames(match.data) == varname)]
	m1.m <- mean(match.v[match.data[[treatment]] == 1], na.rm= TRUE)
	m0.m <- mean(match.v[match.data[[treatment]] == 0], na.rm = TRUE)	

	sdiff.after <- (m1.m - m0.m)/sqrt(0.5*(s1^2 + s0^2))
	if(s1 == 0 && s0 == 0) sdiff.after <- ifelse(m1.m == m0.m, 0, Inf)

	return(c('Treated Mean Before' = m1, 'Control Mean Before' = m0,'SDiff Before' = sdiff.before, 'Treated Mean After' = m1.m, 'Control Mean After' = m0.m, 'SDiff After' = sdiff.after))
}

## The original version of the following function was written by Paul Rosenbaum and distributed in the supplemental material to the paper: "Optimal Matching of an Optimally Chosen Subset in Observational Studies," Paul R. Rosenbaum, Journal of Computational and Graphical Statistics, Vol. 21, Iss. 1, 2012.


#' Robust Mahalanobis Distance
#' 
#' Computes robust Mahalanobis distance between treated and control units.
#' 
#' For an explanation of the robust Mahalanobis distance, see section 8.3 of
#' the first reference.  This function was written by Paul Rosenbaum and
#' distributed in the supplemental material to the second reference.
#' 
#' @param z vector of treatment indicators (1 for treated, 0 for controls).
#' @param X matrix of numeric variables to be used for computing the
#' Mahalanobis distance.  Row count must match length of \code{z}.
#' @return a matrix of robust Mahalanobis distances, with a row for each
#' treated unit and a column for each control.
#' @author Paul R. Rosenbaum.
#' @references Rosenbaum, Paul R. (2010). \emph{Design of Observational
#' Studies}.  Springer-Verlag.
#' 
#' Rosenbaum, Paul R. (2012) "Optimal Matching of an Optimally Chosen Subset in
#' Observational Studies."  Journal of Computational and Graphical Statistics,
#' 21.1, 57-71.
#' @keywords internal
#' @export smahal
smahal <-
function(z,X){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z)
    for (j in 1:k) X[,j]<-rank(X[,j])
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    out<-matrix(NA,m,n-m)
    Xc<-X[z==0,]
    Xt<-X[z==1,]
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
    icov<-ginv(cv)
    for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    out
}




#' Aggregate Student Data into School Data
#' 
#' Takes a dataframe of student-level covariates and aggregates selected
#' columns into a dataframe of school covariates.
#' 
#' Aggregation is either done by taking averages or by selecting the unique
#' factor value when a school has only one value for a factor.  As a result,
#' \code{school.covs} should only include variables that are numeric or do not
#' vary within schools.
#' 
#' @param students a dataframe of students.
#' @param school.cov a character vector of column names in \code{students} that
#' should be aggregated by school.
#' @param school.id the name of the column in \code{students} containing the
#' unique school identifier.
#' @return a dataframe of aggregated data, with one row for each school and
#' columns in \code{school.covs} and \code{school.id}.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @importFrom plyr ddply colwise
#' @export students2schools
students2schools <-
function(students, school.cov, school.id){
	school.df <- students[c(school.id, school.cov)]
	school.df <- plyr::ddply(school.df, school.id, plyr::colwise(agg))
	school.df	
}





#' Outcome analysis.
#' 
#' Calculates confidence interval via grid search.
#' 
#' 
#' @param beta Confidence interval value
#' @param obj a multiMatch object
#' @param out.name Name of outcome covariate
#' @param schl_id_name Name of school (group) identifier
#' @param treat.name Name of treatment indicator
#' @param alpha Level of test for confidence interval, default is .05 for 95\%
#' CI.
#' @param alternative Direction of test.
#' @return The endpoint of an estimated confidence interval.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @references Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
#' Springer-Verlag.
#' 
#' Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
#' Springer-Verlag.
#' @keywords internal
#' @export ci_func
ci_func <- function(beta, obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, alpha, alternative="less"){
	 match.data <- obj$matched
     #No of Strata
     n.s <- dim(unique(match.data[paste(schl_id_name)]))[1]
     Q.s <- matrix(NA, n.s, 1)

    #Rank Within Paired Strata
     z <- match.data[[treat.name]]
     match.data$adj.y <- match.data[[out.name]] - z*beta
     strata.ranks <- tapply(match.data$adj.y, match.data$pair.id, rank, simplify=TRUE)
     match.data$ranks <- matrix(unlist(strata.ranks), nrow(match.data), 1)
     
     sub.trt <- match.data[z==1,]
     sub.ctrl <- match.data[z==0,]
     n.s1 <- tapply(sub.trt$ranks, sub.trt$pair.id, length)
     n.s2 <- tapply(sub.ctrl$ranks, sub.ctrl$pair.id, length)
     d <- (sum(n.s1 + n.s2))
     q1 <- as.vector(tapply(sub.trt$ranks, sub.trt$pair.id, mean))
     q2 <- as.vector(tapply(sub.ctrl$ranks, sub.ctrl$pair.id, mean))

     #Constant Weights
     ws.1 <- 1
     Q.s.c <- (q1- q2)*ws.1
     T.c <- sum(Q.s.c)
     var.T.c <- sum(Q.s.c^2)
     #One-sided p-value for Test of the Sharp Null
     if(alternative=="less"){
     	 pval <- pnorm(T.c/sqrt(var.T.c))
     } else {
     	pval <- 1 - pnorm(T.c/sqrt(var.T.c))
     }
     pval - alpha
}



#' Outcome analysis.
#' 
#' Calculates Hodges-Lehmann point estimate for treatment effect.
#' 
#' 
#' @param beta Point estimate value
#' @param obj A multiMatch object
#' @param out.name Name of outcome covariate
#' @param schl_id_name Name of school (group) identifier
#' @param treat.name Name of treatment indicator
#' @return A point estimate for constant-additive treatment effect.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @references Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
#' Springer-Verlag.
#' 
#' Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
#' Springer-Verlag.
#' @keywords internal
#' @export pe_func
pe_func <- function(beta, obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL){
	 match.data <- obj$matched
     #No of Strata
     n.s <- dim(unique(match.data[paste(schl_id_name)]))[1]
     Q.s <- matrix(NA, n.s, 1)

    #Rank Within Paired Strata
     z <- match.data[[treat.name]]
     match.data$adj.y <- match.data[[out.name]] - z*beta
     strata.ranks <- tapply(match.data$adj.y, match.data$pair.id, rank, simplify=TRUE)
     match.data$ranks <- matrix(unlist(strata.ranks), nrow(match.data), 1)
     
     sub.trt <- match.data[z==1,]
     sub.ctrl <- match.data[z==0,]
     n.s1 <- tapply(sub.trt$ranks, sub.trt$pair.id, length)
     n.s2 <- tapply(sub.ctrl$ranks, sub.ctrl$pair.id, length)
     d <- (sum(n.s1 + n.s2))
     q1 <- as.vector(tapply(sub.trt$ranks, sub.trt$pair.id, mean))
     q2 <- as.vector(tapply(sub.ctrl$ranks, sub.ctrl$pair.id, mean))

     #Constant Weights
     ws.1 <- 1
     Q.s.c <- (q1- q2)*ws.1
     T.c <- sum(Q.s.c)
     var.T.c <- sum(Q.s.c^2)
     E.T <- T.c/var.T.c
     E.T - 0
     }



#' Outcome analysis.
#' 
#' Calcualtes p-values for test of sharp null for treatment effect.
#' 
#' 
#' @param obj A multiMatch object
#' @param out.name Name of outcome covariate
#' @param schl_id_name Name of school (group) identifier
#' @param treat.name Name of treatment indicator
#' @param wt Logical flag for whether p-value should weight strata by size.
#' @return A p-value for constant-additive treatment effect.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @references Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
#' Springer-Verlag.
#' 
#' Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
#' Springer-Verlag.
#' @keywords internal
#' @export pval_func
pval_func <- function(obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, wt=TRUE){
	 match.data <- obj$matched
     #No of Strata
     n.s <- dim(unique(match.data[paste(schl_id_name)]))[1]
     Q.s <- matrix(NA, n.s, 1)

    #Rank Within Paired Strata
     strata.ranks <- tapply(match.data[[out.name]], match.data$pair.id, rank, simplify=TRUE)
     match.data$ranks <- matrix(unlist(strata.ranks), nrow(match.data), 1)
     z <- match.data[[treat.name]]
     sub.trt <- match.data[z==1,]
     sub.ctrl <- match.data[z==0,]
     n.s1 <- tapply(sub.trt$ranks, sub.trt$pair.id, length)
     n.s2 <- tapply(sub.ctrl$ranks, sub.ctrl$pair.id, length)
     d <- (sum(n.s1 + n.s2))
     q1 <- as.vector(tapply(sub.trt$ranks, sub.trt$pair.id, mean))
     q2 <- as.vector(tapply(sub.ctrl$ranks, sub.ctrl$pair.id, mean))

     #Constant Weights
     if(wt==FALSE){
     	     ws.1 <- 1
     Q.s.c <- (q1- q2)*ws.1
     T.c <- sum(Q.s.c)
     var.T.c <- sum(Q.s.c^2)
     #One-sided p-value for Test of the Sharp Null
     pval.c1 <- pnorm(T.c/sqrt(var.T.c))
     pval.c2 <- 1 - pnorm(T.c/sqrt(var.T.c))
     pval.c <- min(pval.c1, pval.c2)
     return(pval.c)	
     } else {
     d <- (sum(n.s1 + n.s2))
     ws.2 <- (n.s1 + n.s2) /d
     Q.s.p <- ((q1- q2))*ws.2
     T.p <- sum(Q.s.p)
     var.T.p <- sum(Q.s.p^2)
     #One-sided p-value for Test of the Sharp Null
     pval.p1 <- pnorm(T.p/sqrt(var.T.p))
     pval.p2 <- 	1 -  pnorm(T.p/sqrt(var.T.p))
     pval.p <- min(pval.p1, pval.p2)
     return(pval.p)
     }
 	
}