# Combination test function
combination.test = function(p.vec, weights=NA,
                            method = list('Cauchy','Log Cauchy','Levy',
                                          'Pareto','Frechet','t','Inverse Gamma'),
                            tail.idx = 1,
                            truncate = F,truncate.threshold = 0.99){
  p = length(p.vec)
  # if the number of p-values is 1, don't need to combine,,                                                                    p-values.
  if(p==1) return(p.vec[1])

  # if there's no given weights, set all hypotheses equal weights
  if(sum(is.na(weights))==p){
    weights = rep(1,p)
  }

  if(!method %in% c('Cauchy','Log Cauchy','Levy',
                    'Pareto','Frechet','t','Inverse Gamma')){
    stop('No such distributions!')
  }

  if(method == 'Cauchy'){
    kappa = sum(weights)
    normalized_weights = weights/kappa
    if(truncate){
      S = sum(normalized_weights*qtcauchy(p.vec, pthreshold=truncate.threshold)/p)
      p.global = min(ptcauchy(S,pthreshold=truncate.threshold, lower.tail = F),1)
    }else{
      S = sum(normalized_weights*tan((0.5-p.vec)*pi)/p)
      p.global = min(pcauchy(S,lower.tail = F),1)
    }
  }

  if(method == 'Log Cauchy'){
    p.global  = lcauchy_trans(p.vec, weights=weights)
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
