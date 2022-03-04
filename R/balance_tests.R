

agg.balance <- function(varname, treatment, school.id, 
                        data,
                        treat.wts = NULL, ctrl.wts = NULL ) {
  if ( is.null(treat.wts) != is.null( ctrl.wts ) ) {
    stop( "Cannot have weights for just treatment or just control groups" )
  }
  
  if ( is.null( data[[school.id]] ) ) {
    sstop( "Need to pass clustering/grouping variable. %s not found", school.id )  
  }
  
  txvec = data[[treatment]]
  
  if ( !is.null( treat.wts ) ) {
    warning( "Weights not fully implemented for aggregate balance check" )
    data$.weight[ txvec != 0 ] = treat.wts 
    data$.weight[ txvec == 0 ] = ctrl.wts
  } else {
    data$.weight = 1
  }
  
  school.ids = data[[school.id]]
  n_co = length( unique( school.ids[txvec==0] ) )
  n_tx = length( unique( school.ids[txvec!=0] ) )
  stopifnot( !is.na( n_co ) && !is.na( n_tx ) )
  
  
  data.agg <- data %>% dplyr::group_by( !!rlang::sym(treatment), !!rlang::sym(school.id) ) %>% 
    dplyr::summarize( mn = weighted.mean( !!rlang::sym(varname), w = .data$.weight ),
               .weight = sum( .data$.weight ) )
  
  
  form = paste0( "mn ~ 1 + ", treatment )
  M0 = lm( form, data=data.agg ) #, weights = data.agg$.weight )
  
  # est ATE (precision weighted)
  ATE_hat <- coef(M0)[[treatment]]
  
 # if ( sd( resid( M0 ) ) < 0.00000001 ) {
#    browser() 
 # }
 # print( sd( resid( M0 ) ) )
  
  # Heteroskedastic robust SEs (clustering at site level)
  vcov_sand <- sandwich::vcovHC(M0, type = "HC1")
  SE <-sqrt( vcov_sand[2,2] )
  
  2 * pt( -abs( ATE_hat / SE ), df = min( n_co, n_tx ) - 1 )  
}





CRVE.balance <- function(varname, treatment, school.id, 
                         data,
                         treat.wts = NULL, ctrl.wts = NULL ) {
  
  if ( is.null(treat.wts) != is.null( ctrl.wts ) ) {
    stop( "Cannot have weights for just treatment or just control groups" )
  }
  
  if ( is.null( data[[school.id]] ) ) {
    sstop( "Need to pass clustering/grouping variable. %s not found", school.id )  
  }
  
  txvec = data[[treatment]]
  school.ids = data[[school.id]]
  n_co = length( unique( school.ids[txvec==0] ) )
  n_tx = length( unique( school.ids[txvec!=0] ) )
  
  if ( !is.null( treat.wts ) ) {
    data$weight[ txvec != 0 ] = treat.wts 
    data$weight[ txvec == 0 ] = ctrl.wts
  }
  
  form = paste0( "`", varname, "` ~ 1 + ", treatment )
  M0 = lm( form, data=data, weights = data$weight )
  
  # est ATE (precision weighted)
  ATE_hat <- coef(M0)[[treatment]]
  
  # Cluster robust SEs (clustering at site level)
  vcov_clust <- sandwich::vcovCL(M0, data[[school.id]])
  SE <- sqrt(vcov_clust[2, 2])
  
  2 * pt( -abs( ATE_hat / SE ), df = min( n_co, n_tx ) - 1 )    
  
}




ttest.balance <-
  function(varname, treatment, orig.data, match.data = NULL, 
           treat.wts = NULL, ctrl.wts = NULL, 
           mt.wts = NULL, mc.wts = NULL) {
    if(is.null(treat.wts)) treat.wts = rep(1, sum(orig.data[[treatment]]))
    if(is.null(ctrl.wts)) ctrl.wts = rep(1, sum(orig.data[[treatment]] == 0))	
    
    treat.wts <- treat.wts/sum(treat.wts)*length(treat.wts)
    ctrl.wts <- ctrl.wts/sum(ctrl.wts)*length(ctrl.wts)
    orig.v <- orig.data[,which(colnames(orig.data) == varname)]
    
    t.orig <- weights::wtd.t.test(orig.v[orig.data[[treatment]] == 1], 
                                  orig.v[orig.data[[treatment]] == 0], 
                                  weight = treat.wts, weighty = ctrl.wts, samedata = FALSE)
    
    if(is.null(match.data)) return(c('T-test Pvalue' = t.orig$coefficients[[3]]))
    
    if(is.null(mc.wts)) mc.wts = rep(1, sum(match.data[[treatment]] == 0))
    if(is.null(mt.wts)) mt.wts = rep(1, sum(match.data[[treatment]] == 1))
    mt.wts <- mt.wts/sum(mt.wts)*length(mt.wts)
    mc.wts <- mc.wts/sum(mc.wts)*length(mc.wts)
    match.v  <- match.data[,which(colnames(match.data) == varname)]	
    
    t.match <- weights::wtd.t.test(match.v[match.data[[treatment]] == 1], 
                                   match.v[match.data[[treatment]] == 0], 
                                   weight = mt.wts, weighty = mc.wts, samedata = FALSE)
    
    return(c('T-test Before' = t.orig$coefficients[[3]], 'T-test After' = t.match$coefficients[[3]]))
  }





wilc.balance  <- function(varname, treatment, orig.data, match.data = NULL, 
                          treat.wts = NULL, ctrl.wts = NULL, 
                          mt.wts = NULL, mc.wts = NULL){
  if(is.null(treat.wts)) treat.wts = rep(1, sum(orig.data[[treatment]]))
  if(is.null(ctrl.wts)) ctrl.wts = rep(1, sum(orig.data[[treatment]] == 0))	
  
  treat.wts <- treat.wts/sum(treat.wts)*length(treat.wts)
  ctrl.wts <- ctrl.wts/sum(ctrl.wts)*length(ctrl.wts)
  orig.v <- orig.data[,which(colnames(orig.data) == varname)]
  
  t.orig <- wilcox_test(y ~ z, data = data.frame('y' = orig.v, 
                                                 'z' = factor(orig.data[[treatment]]))) #can add weights
  #	t.orig <- wilcox.test(orig.v[orig.data[[treatment]] == 1], orig.v[orig.data[[treatment]] == 0], exact = FALSE)	
  
  if(is.null(match.data)) return(c('Wilcoxon Test Pvalue' = pvalue(t.orig)))
  #return(c('Wilcoxon Test Pvalue' = t.orig$p.value))
  
  if(is.null(mc.wts)) mc.wts = rep(1, sum(match.data[[treatment]] == 0))
  if(is.null(mt.wts)) mt.wts = rep(1, sum(match.data[[treatment]] == 1))
  mt.wts <- mt.wts/sum(mt.wts)*length(mt.wts)
  mc.wts <- mc.wts/sum(mc.wts)*length(mc.wts)
  match.v  <- match.data[,which(colnames(match.data) == varname)]	
  
  t.match <- wilcox_test(y ~ z, data = data.frame('y' = match.v, 'z' = factor(match.data[[treatment]]))) #can add weights	
  #t.match <- wilcox.test(match.v[match.data[[treatment]] == 1], match.v[match.data[[treatment]] == 0], exact = FALSE)
  
  #	return(c('Wilcoxon Test Pvalue Before' = t.orig$p.value, 'Wilcoxon Test Pvalue After' = t.match$p.value))
  return(c('Wilcoxon Test Pvalue Before' = pvalue(t.orig), 'Wilcoxon Test Pvalue After' = pvalue(t.match)))
  
}

