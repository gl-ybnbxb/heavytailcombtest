# generate the p-values
data_Gen = function(n, p = 2, mu,rho = 0, 
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
  
  return(list(p.mat = p.mat, z.mat = z.mat, mu = mu))
  
}

