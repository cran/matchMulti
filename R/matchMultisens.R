#' Rosenbaum Bounds after Multilevel Matching
#' 
#' Function to calculate Rosenbaum bounds for continuous outcomes after
#' multilevel matching.
#' 
#' This function returns a single p-value, but actually conducts two tests.
#' The first assumes that the treatment effect does not vary with cluster size.
#' The second allows the treatment effect to vary with cluster size.  The
#' function returns a single p-value that is corrected for multiple testing.
#' This p-value is the upper bound for a single Gamma value
#' 
#' @param obj A multilevel match object
#' @param out.name Outcome variable name
#' @param schl_id_name Level 2 ID variable name, that is this variable
#' identifies clusters matched in the data.
#' @param treat.name Treatment indicator name
#' @param Gamma Sensitivity analysis parameter value. Default is one.
#' @return
#' 
#' \item{pval}{Upper bound on one-sided approximate p-value for test of the
#' sharp null.}
#' @author Luke Keele, University of Pennsylvania, \email{luke.keele@gmail.com}
#' 
#' Sam Pimentel, University of California, Berkeley, \email{spi@@berkeley.edu}
#' @seealso See Also as \code{\link{matchMulti}},
#' \code{\link{matchMultioutcome}}
#' @references Rosenbaum, Paul R. (2002) Observational Studies.
#' Springer-Verlag.
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
#'           }
#' 
#' @export matchMultisens
matchMultisens <- function(obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, Gamma=1){
	 match.data <- obj$matched
     #No of Strata
     n.s <- dim(unique(match.data[paste(schl_id_name)]))[1]
     Q.s <- matrix(NA, n.s, 1)

    #Rank Within Paired Strata
     z <- match.data[[treat.name]]
     strata.ranks <- tapply(match.data[[out.name]], match.data$pair.id, rank, simplify=TRUE)
     match.data$ranks <- matrix(unlist(strata.ranks), nrow(match.data), 1)
     
     sub.trt <- match.data[z==1,]
     sub.ctrl <- match.data[z==0,]
     n.s1 <- tapply(sub.trt$ranks, sub.trt$pair.id, length)
     n.s2 <- tapply(sub.ctrl$ranks, sub.ctrl$pair.id, length)
     d <- (sum(n.s1 + n.s2))
     q1 <- as.vector(tapply(sub.trt$ranks, sub.trt$pair.id, mean))
     q2 <- as.vector(tapply(sub.ctrl$ranks, sub.ctrl$pair.id, mean))
     ws.1 <- 1
     Q.s.c <- (q1- q2)*ws.1
     T.c <- sum(Q.s.c)/n.s
     E.s.c <- ((Gamma-1)/(n.s*(Gamma+1))) *sum(abs(Q.s.c))
     V.s.c <- (4*Gamma/(n.s^2*(1+Gamma)^2))*sum(Q.s.c^2)
     dev.c <- (T.c - E.s.c)/sqrt(V.s.c)
     pval.sens <- 1 - pnorm(dev.c)

    cat("Upper bound on one-sided p-value is: ", pval.sens, "\n")
    res <- list(pval=pval.sens)
    }
