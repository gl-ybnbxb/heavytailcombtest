#' The Pareto Distribution
#'
#' @description Distribution and quantile function for the Pareto distribution
#' with the location parameter equal to \code{x_m} and shape parameter equal to \code{alpha}.
#'
#' @usage ppareto(x, x_m = 1, alpha = 1, lower.tail = T)
#' @usage qpareto(p, x_m = 1, alpha = 1)
#' @param x vector of quantiles
#' @param p vectors of probabilities
#' @param x_m vectors of location parameters
#' @param alpha vectors of shape parameters
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X\le x)}; otherwise, \eqn{P(X>x)}.
#'
#' @export
ppareto = function(x, x_m = 1, alpha = 1, lower.tail = T){
  if (lower.tail){
    y = 1-(x_m/x)^alpha
  }
  else{
    y = (x_m/x)^alpha
  }
  return(y)
}

#' The Pareto Distribution
#'
#' @description Distribution and quantile function for the Pareto distribution
#' with the location parameter equal to \code{x_m} and shape parameter equal to \code{alpha}.
#'
#' @usage ppareto(x, x_m = 1, alpha = 1, lower.tail = T)
#' @usage qpareto(p, x_m = 1, alpha = 1)
#' @param x vector of quantiles
#' @param p vectors of probabilities
#' @param x_m vectors of location parameters
#' @param alpha vectors of shape parameters
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X\le x)}; otherwise, \eqn{P(X>x)}.
#'
#' @export
qpareto = function(p, x_m = 1, alpha = 1){
  y = x_m/(1-p)^(1/alpha)
  return(y)
}


#' The Frechet Distribution
#'
#' @description Distribution and quantile function for the Frechet distribution
#' with the location parameter equal to \code{m}, scale parameter equal to \code{s}, and shape parameter equal to \code{alpha}.
#'
#' @usage pfrechet(x, m=0, s=1, alpha=1, lower.tail = T)
#' @usage qfrechet(p, m=0, s=1, alpha=1)
#' @param x vector of quantiles
#' @param p vectors of probabilities
#' @param m vectors of location parameters
#' @param s vectors of scale parameters
#' @param alpha vectors of shape parameters
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X\le x)}; otherwise, \eqn{P(X>x)}.
#'
#' @export
pfrechet = function(x, m=0, s=1, alpha=1, lower.tail = T){
  if (lower.tail){
    y = exp(-((x-m)/s)^(-alpha))
  }
  else{
    y = 1-exp(-((x-m)/s)^(-alpha))
  }
  return(y)
}

#' The Frechet Distribution
#'
#' @description Distribution and quantile function for the Frechet distribution
#' with the location parameter equal to \code{m}, scale parameter equal to \code{s}, and shape parameter equal to \code{alpha}.
#'
#' @usage pfrechet(x, m=0, s=1, alpha=1, lower.tail = T)
#' @usage qfrechet(p, m=0, s=1, alpha=1)
#' @param x vector of quantiles
#' @param p vectors of probabilities
#' @param m vectors of location parameters
#' @param s vectors of scale parameters
#' @param alpha vectors of shape parameters
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X\le x)}; otherwise, \eqn{P(X>x)}.
#'
#' @export
qfrechet = function(p, m=0, s=1, alpha=1){
  m+s/(-log(p))^(1/alpha)
}



#' The Log Cauchy Distribution
#'
#' @description Distribution and quantile function for the log Cauchy distribution
#'
#' @usage plcauchy(x, lower.tail = T)
#' @usage qlcauchy(p)
#' @param x vector of quantiles
#' @param p vectors of probabilities
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X\le x)}; otherwise, \eqn{P(X>x)}.
#'
#' @export
plcauchy = function(x, lower.tail = T){
  if (lower.tail){
    y = atan(log(x))/pi+0.5
  }
  else{
    y = 0.5-atan(log(x))/pi
  }
  return(y)
}

#' The Log Cauchy Distribution
#'
#' @description Distribution and quantile function for the log Cauchy distribution
#'
#' @usage plcauchy(x, lower.tail = T)
#' @usage qlcauchy(p)
#' @param x vector of quantiles
#' @param p vectors of probabilities
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X\le x)}; otherwise, \eqn{P(X>x)}.
#'
#' @export
qlcauchy = function(p){
  y = exp(tan(pi*(p-0.5)))
  return(y)
}

