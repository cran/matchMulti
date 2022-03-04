#' Performs an outcome analysis after multilevel matching.
#' 
#' This function returns a point estimate, 95\% confidence interval, and
#' p-values for the matched multilevel data. All results are based on
#' randomization inference.
#' 
#' It may be necessary to adjust the lower and upper bounds if one expects the
#' treatment effect confidence interval to be outside the range of -1000 or
#' 1000.
#' 
#' @param obj A multilevel match object.
#' @param out.name Outcome variable name
#' @param schl_id_name Level 2 ID variabel name. This variable identifies the
#' clusters in the data that you want to match.
#' @param treat.name Treatment variable name, must be zero or one.
#' @param end.1 Lower bound for point estimate search, default is -1000.
#' @param end.2 Upper bound for point estimate search, default is 1000.
#' @return \item{pval.c }{One-sided approximate p-value for test of the sharp
#' null.} \item{pval.p }{One-sided approximate p-value for test of the sharp
#' null assuming treatment effects vary with cluster size} \item{ci1}{Lower
#' bound for 95\% confidence interval.} \item{ci2}{Upper bound for 95\%
#' confidence interval.} \item{p.est}{Point estimate for the group level
#' treatment effect.}
#' @author Luke Keele, Penn State University, ljk20@@psu.edu
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @seealso See Also as \code{\link{matchMulti}}, \code{\link{matchMultisens}}
#' @references Rosenbaum, Paul R. (2002) Observational Studies.
#' Springer-Verlag.
#' @keywords outcomes
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
#' 
#' #Estimate treatment effect
#' output <- matchMultioutcome(match.simple, out.name = "mathach",
#' schl_id_name = "school",  treat.name = "sector")
#' 
#'   }
#' 
#' @export matchMultioutcome
matchMultioutcome <- function(obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, end.1=-1000, end.2=1000){
	
	
   pval.c <- pval_func(obj, out.name = out.name, schl_id_name = schl_id_name,  treat.name = treat.name, wt=FALSE)
   pval.p <- pval_func(obj, out.name = out.name, schl_id_name = schl_id_name,  treat.name = treat.name, wt=TRUE)
          
   ci1 <- uniroot(ci_func, c(end.1, end.2), obj=obj, out.name = out.name,
          schl_id_name = schl_id_name, treat.name = treat.name, alternative="less", alpha=.025)$root;
          
   ci2 <- uniroot(ci_func, c(end.1, end.2), obj=obj, out.name = out.name,
          schl_id_name = schl_id_name,  treat.name = treat.name, alternative="great", alpha=.025)$root; 
          
   pe <- uniroot(pe_func, c(end.1, end.2), obj=obj, out.name = out.name,
          schl_id_name = schl_id_name, 
          treat.name = treat.name)$root;   
          
   ci.lo <- min(ci1, ci2)
   ci.up <- max(ci1, ci2) 
   ci <- c(ci.lo, ci.up)          
     
	cat("Point Estimate is: ", pe, "\n")
	cat("95 confidence interval", ci,"\n")
	cat("One.sided p-value is: ", pval.c, "\n")
	res <- list(pval.c = pval.c, pval.p = pval.p, ci1=ci1, ci2=ci2, p.est=pe)
	return(res)
}



     
     
     

