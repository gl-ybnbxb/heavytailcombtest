#' Generate p-value matrix for simulation
#'
#' @description Generate p-value and base test statistics matrix for simulations.
#' Each row corresponds to a base p-value vector from one sampling and each column represents p-values from one base hypothesis.
#' @param n number of repeated times
#' @param p number of base hypotheses
#' @param mu a vector of possible values of effect size of nonnull individual
#' hypotheses. When \code{mu=(0,0,...,0)}, it's the global null.
#' @param rho a common correlation between base hypotheses
#' @param copula base test statistics from multivariate gaussian or multivariate t
#' @param dof the degree of freedom of the multivariate t
#' @param one.sided either construct one-sided p-value or not
#'
#' @return a list of objects
#' \describe{
#' \item{p.mat}{the p-value matrix of size \code{n * p} where each column
#' corresponding to one base hypothesis}
#' \item{z.mat}{The corresponding base test statistics matrix of size \code{n * p},
#' which are either z-scores or t-scores}
#' }
#'
#' @export


data_Gen = function(n, p = 2,
                    mu, rho = 0,
                    copula = 'gaussian', dof = 5,
                    one_sided = F){

  z.mat = matrix(rnorm(n*p),nrow = n)
  Sigma = diag(p)

  if (rho != 0) {
    Sigma <- matrix(rep(rho, p^2), nrow = p)
    diag(Sigma) <- 1
    svd.Sigma <- svd(Sigma)
    Sigma_root <- t(t(svd.Sigma$u) * sqrt(svd.Sigma$d))

    z.mat <- t(Sigma_root %*% t(z.mat))
  }

  # t copula
  if(copula == 't'){
    z.mat = rmvt(n,sigma = Sigma,df=dof)
  }

  # z-value
  z.mat = t(apply(z.mat, 1, function(x) return(x+mu)))

  # p-value
  if(one_sided==F){
    if(copula == 'gaussian'){
      raw.p.mat <- pnorm(abs(z.mat), lower.tail = F)
    }else{
      raw.p.mat <- pt(abs(z.mat), df=dof, lower.tail = F)
    }
    p.mat = 2*raw.p.mat
  }else{
    if(copula == 'gaussian'){
      p.mat <- pnorm(z.mat, lower.tail = F)
    }else{
      p.mat <- pt(z.mat, df=dof, lower.tail = F)
    }
  }

  return(list(p.mat = p.mat, z.mat = z.mat))

}

