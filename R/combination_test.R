#' Compute the global p-value using the heavy-tailed combination test
#'
#' @param p.vec an vector of base p-values. Missing values are allowed.
#' @param weights an vector of non-negative weights for each base hypothesis.
#' When it's missing, all hypotheses are assigned equal weights.
#' @param method one of "Cauchy", "Log Cauchy", "Levy", "Paret", "Frechet", "t", or "Inverse Gamma", with the default "Cauchy".
#' @param tail.idx is the shape parameter for Cauchy, Pareto, Frechet, t, and Inverse Gamma.
#' @param truncate logical value, default is FALSE.
#' If \code{truncate} is TRUE, then apply truncated Cauchy instead of standard Cauchy.
#' @param truncate.threshold truncated threshold (upper bound) for the p-values when \code{truncate} is TRUE.
#'
#' @return a p-value for the global test
#'
#'
#' @import invgamma
#' @import rmutil
#' @export


combination.test = function(p.vec, weights=NA,
                            method = 'Cauchy',
                            tail.idx = 1,
                            truncate = F,truncate.threshold = 0.99){
  # if the distribution is not in the list, raise error.
  if(!method %in% c('Cauchy','Log Cauchy','Levy',
                    'Pareto','Frechet','t','Inverse Gamma')){
    stop('No such distributions!')
  }

  # if one base p-value is missing, set its weight as 0
  weights[is.na(p.vec)] = 0
  p = length(p.vec)
  # if the number of p-values is 1, don't need to combine,,                                                                    p-values.
  if(p==1) return(p.vec[1])

  # if there's no given weights, set all hypotheses equal weights
  if(sum(is.na(weights))==p){
    weights = rep(1,p)
  }

  if(method == 'Cauchy'){
    kappa = sum(weights)
    normalized_weights = weights/kappa
    if(truncate){
      cauchy.threshold = qcauchy(1-truncate.threshold)
      S = sum(normalized_weights*qtcauchy(p.vec, threshold=cauchy.threshold)/p)
      p.global = min(ptcauchy(S, threshold=cauchy.threshold, lower.tail = F),1)
    }else{
      S = sum(normalized_weights*tan((0.5-p.vec)*pi)/p)
      p.global = min(pcauchy(S,lower.tail = F),1)
    }
  }

  if(method == 'Log Cauchy'){
    trans.p.vec = qcauchy(1-p.vec)
    trans.max = max(trans.p.vec)
    mid = trans.max + log(sum(weights*exp(trans.p.vec-trans.max)))
    p.global  = min(1,p*pcauchy(mid,lower.tail = F))
  }

  if(method == 'Levy'){
    tail.idx = 0.5
    kappa = sum(weights^tail.idx)
    require(rmutil)
    S = sum(weights*qlevy(1-p.vec))
    p.global = min(kappa*(1-plevy(S)),1)
  }

  if(method == 'Pareto'){
    kappa = sum(weights^tail.idx)
    S = sum(weights*qpareto(1-p.vec, alpha = tail.idx))
    p.global = min(kappa*ppareto(S,lower.tail = F, alpha = tail.idx),1)
  }

  if(method == 'Frechet'){
    kappa = sum(weights^tail.idx)
    S = sum(weights*qfrechet(1-p.vec, alpha = tail.idx))
    p.global = min(kappa*pfrechet(S, lower.tail = F, alpha = tail.idx),1)
  }

  if(method == 't'){
    kappa = sum(weights^tail.idx)
    S = sum(weights*qt(1-p.vec,df = tail.idx))
    p.global = min(kappa*pt(S,df = tail.idx,lower.tail = F),1)
  }

  if(method == 'Inverse Gamma'){
    kappa = sum(weights^tail.idx)
    require(invgamma)
    S = sum(weights*qinvgamma(1-p.vec,shape = tail.idx))
    p.global = min(kappa*pinvgamma(S,shape = tail.idx,lower.tail = F),1)
  }

  return(p.global)
}
