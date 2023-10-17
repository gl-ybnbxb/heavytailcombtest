# cdf and quantile of the Pareto distribution
ppareto = function(x, x_m = 1, alpha = 1, lower.tail = T){
  if (lower.tail){
    y = 1-(x_m/x)^alpha
  }
  else{
    y = (x_m/x)^alpha
  }
  return(y)
}

qpareto = function(x, x_m = 1, alpha = 1){
  y = x_m/(1-x)^(1/alpha)
  return(y)
}


# cdf and quantile function of the Frechet distribution
pfrechet = function(x, m=0, s=1, alpha=1, lower.tail = T){
  if (lower.tail){
    y = exp(-((x-m)/s)^(-alpha))
  }
  else{
    y = 1-exp(-((x-m)/s)^(-alpha))
  }
  return(y)
}

qfrechet = function(x, m=0, s=1, alpha=1){
  m+s/(-log(x))^(1/alpha)
}



#  cdf and quantile function of the log Cauchy distribution
plcauchy = function(x, lower.tail = T){
  if (lower.tail){
    y = atan(log(x))/pi+0.5
  }
  else{
    y = 0.5-atan(log(x))/pi
  }
  return(y)
}


qlcauchy = function(x){
  y = exp(tan(pi*(x-0.5)))
  return(y)
}

# don't do the full transformation to avoid NA

lcauchy_mid= function(x){
  tan(pi*(x-0.5))
}
plcauchy_mid = function(x){
  0.5-atan(x)/pi
}

# the transformation of a z vector to a p value of the global test
# input is a p-value vector
# output is a global p-value
lcauchy_trans= function(x, weights){
  p = length(x)
  trans_x = sapply(x,function(s) lcauchy_mid(1-s))
  m = max(trans_x)
  mid = m+log(sum(weights*exp(trans_x-m)))
  p.lcauchy  = min(1,p*plcauchy_mid(mid))
  
  return(p.lcauchy)
}


# truncated Cauchy:
# the threshold is for p-values:
ptcauchy = function(x, pthreshold, lower.tail=T){
  (pcauchy(x)+pthreshold-1)/pthreshold
}

qtcauchy = function(x, pthreshold){
  qcauchy(pthreshold*x+1-pthreshold)
}

