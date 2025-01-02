#' Compute the global p-value using the heavy-tailed combination test.
#'
#' @description Computes a global p-value using the heavy-tailed combination test.
#' The function \code{combination.test} realizes all three versions---standard, average, and weighted versions---of the combination test.
#' Specific choices also correspond to the Cauchy combination test and harmonic mean p-value method.
#'
#' @docType package
#' @name heavytailcombtest_package
NULL

#' The main function for the heavy-tailed combination test
#'
#' @description This function computes a global p-value using the heavy-tailed combination test.
#' When the weights are not specified, all hypotheses are assigned equal weights \code{1} and the standard combination test is applied.
#' When the \code{tail.idx=1} and all weights are normalized equal weights, i.e. the sum of the weights is \code{1}, the average version of combination test is applied.
#'
#' Specifically, when \code{method="Cauchy"}, the Cauchy combination test is used. When \code{method="Pareto"} and \code{tail.idx=1}, the harmonic mean p-value is returned.
#'
#' \code{truncate} and \code{truncate.thrrshold} are for \code{method="Cauchy"/"t"} to deal with occurrence of base p-value \code{1}.
#' When \code{truncate} is \code{TRUE}, truncated Cauchy/t instead of standard Cauchy/t is implemented as the heavy-tailed distribution for the combination test.
#' \code{truncate.thrshold} is the scaling upper bound chosen for the base p-values. Error will be raised if a base p-value 1 appears but \code{truncated=FALSE}.
#'
#'
#' @param p.vec an vector of base p-values. Missing values are allowed.
#' @param weights an vector of non-negative weights for each base hypothesis.
#' When it's missing, all hypotheses are assigned equal weights \code{1}.
#' @param method one of \code{"Cauchy"}, \code{"Log Cauchy"}, \code{"Levy"}, \code{"Pareto"}, \code{"Frechet"}, \code{"t"}, or \code{"Inverse Gamma"}, with the default \code{"Cauchy"}.
#' @param tail.idx is the shape parameter for Cauchy, Pareto, Frechet, t, and Inverse Gamma.
#' @param truncate logical value, default is \code{FALSE}.
#' If \code{truncate} is \code{TRUE}, then apply truncated Cauchy/t instead of standard Cauchy/t.
#' @param truncate.threshold truncated threshold (upper bound) for the p-values when \code{truncate} is \code{TRUE}.
#'
#' @return a p-value for the global test
#'
#' @examples
#' ## Simulate a one-sided p-value vector of size 5 with 5 base hypotheses under the global null
#' data <- data_Gen(n = 1, p = 5, mu = rep(0,5), rho = -0.2, copula = 'gaussian', one_sided = T)
#' p.vec <- as.vector(data$p.mat)
#'
#' ## Compute the global p-value using the Inverse Gamma combination test with a degree of freedom 2
#' p.global <- combination.test(p.vec, method = 'Inverse Gamma', tail.idx = 2)
#'
#' ## When special base p-value 1 occurs, apply truncated t
#' p.vec <- c(0,0,0,1,1)
#' p.global <- combination.test(p.vec, method = 't', tail.idx = 2, truncate = T, truncate.threshold = 0.99)
#'
#' @import invgamma
#' @import rmutil
#' @import mvtnorm
#' @export


combination.test = function(p.vec, weights=NA, method = 'Cauchy', tail.idx = 1, truncate = T, truncate.threshold = 0.9){
  # if the distribution is not in the list, raise error.
  if(!method %in% c('Cauchy','Log Cauchy','Levy',
                    'Pareto','Frechet','t','Inverse Gamma')){
    stop('No such distributions!')
  }

  # if there's no given weights, set all hypotheses equal weights
  p = length(p.vec)
  if(length(weights)==1 && is.na(weights)){
    weights = rep(1,p)
  }
  # if one base p-value is missing, set its weight as 0
  weights[is.na(p.vec)] = 0
  # if the number of p-values is 1, don't need to combine,,                                                                    p-values.
  if(p==1) return(p.vec[1])

  if(method == 'Cauchy'){
    kappa = sum(weights)
    normalized_weights = weights/kappa
    if(truncate){
      S = sum(normalized_weights*tan((0.5-truncate.threshold*p.vec)*pi))
      p.global = min(pcauchy(S,lower.tail = F)/truncate.threshold,1)
    }else{
      if(1%in%p.vec){
        stop('There exists base p-value 1. Use the truncated version!')
      }
      S = sum(normalized_weights*tan((0.5-p.vec)*pi))
      p.global = min(pcauchy(S,lower.tail = F),1)
    }
  }

  if(method == 'Log Cauchy'){
    trans.p.vec = qcauchy(1-p.vec)
    trans.max = max(trans.p.vec)
    mid = trans.max + log(sum(weights*exp(trans.p.vec-trans.max)))
    p.global  = min(1,p*pcauchy(mid,lower.tail = F))
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
    print(kappa)
    if(truncate){
      S = sum(weights*qt(1-truncate.threshold*p.vec,df = tail.idx))
      p.global = min(kappa*pt(S,df = tail.idx,lower.tail = F)/truncate.threshold,1)
    }else{
      if(1%in%p.vec){
        stop('There exists base p-value 1. Use the truncated version!')
      }
      S = sum(weights*qt(1-p.vec,df = tail.idx))
      p.global = min(kappa*pt(S,df = tail.idx,lower.tail = F),1)
    }
  }

  if(method == 'Inverse Gamma'){
    kappa = sum(weights^tail.idx)
    S = sum(weights*qinvgamma(1-p.vec,shape = tail.idx))
    p.global = min(kappa*pinvgamma(S,shape = tail.idx,lower.tail = F),1)
  }

  if(method == 'Levy'){
    tail.idx = 0.5
    kappa = sum(weights^tail.idx)
    require(rmutil)
    S = sum(weights*qlevy(1-p.vec))
    p.global = min(kappa*(1-plevy(S)),1)
  }

  return(p.global)
}
