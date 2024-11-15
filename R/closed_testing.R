#' The function for the closed testing procedure based on the heavy-tailed combination method
#'
#' @param p.vector a vector of p-values with length \code{n} of p-values. \code{n} is the
#' number of hypotheses.
#' @param alpha the pre-chosen FWER threshold
#' The default value is alpha=0.05.
#' @param type the type of transformation. Student t, Frechet, Pareto or Levy.
#' @return multiple testing decision with True representing rejecting the null 
#' and False representing accepting the null
#'
heavytail_ctp = function(p.vec, alpha=0.05, type='Frechet'){
  require(rmutil)
  
  if(type == 'Frechet'){
    transform_p = function(p){
      y = qfrechet(1-p)
      return(y)
    }
  }
  
  if(type == 'Pareto'){
    transform_p = function(p){
      y = qpareto(1-p)
      return(y)
    }
  }
  
  if(type == 'Levy'){
    transform_p = function(p){
      y = qlevy(1-p)
    }
  }
  
  if(type == 'Student t'){
    transform_p = function(p){
      y = qt(1-p,df=1)
    }
  }
  
  n = length(p.vec)
  
  if(n==1){
    return(p.vec<=alpha)
  }
  
  # rearrange the p value vector into ascending order
  order_idx = order(p.vec)
  p.ordered = p.vec[order_idx]
  
  # compute adjusted p values
  p_transformed = sapply(p.ordered, transform_p)
  threshold_Ik = sapply(alpha/(1:n),transform_p)-c(0,cumsum(p_transformed[n:1])[1:(n-1)])
  threshold_p = cummax(threshold_Ik)[n:1]
  j = min(which(p_transformed<threshold_p))
  decision = c(rep(T,j-1),rep(F,n-j+1))
  
  ori_decision = rep(NA,n)
  ori_decision[order_idx] = decision
  
  return(ori_decision)
}