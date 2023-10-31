# heavytailcombtest

This is an R package for global testing applying the heavy-tailed combination test. This (weighted) heavy-tailed combination test is proposed by [Fang et al., 2023](https://www3.stat.sinica.edu.tw/statistica/J33N21/J33N2101/J33N2101.html). Here we include the combination test with commonly-used heavy-tailed (regularly-varying) distributions with a changable `tail.idx`(shape parameter for the distribution): t, Pareto, Frechet, inverse Gamma, and Levy. 
This method also covers widely-applied [Cauchy combination test](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485) and [harmonic mean p-value](https://www.pnas.org/doi/abs/10.1073/pnas.1814092116).


## Installation
```{r}
library(devtools)

install_github("gl-ybnbxb/heavytailcombtest")
```

## Example

The `combination.test` takes in a base p-value vector and returns a global p-value. 

### Simulate base p-values
As an example, we simulate an one-sided p-value vector with the number of base hypotheses $= 5$ when the global null holds. The base test statistics are sampled from multivariate t with the location parameter $\mu=0$ and $\Sigma_\rho$ which has 1s on the diagonal and a common $\rho$ off the diagonal. Then, we apply the inverse gamma distribution with a shape parameter $2$ for the combination test. With a prior knowledge, we specify a weight vector for base hypotheses.
```
library(heavytailcombtest)

## Simulate a one-sided p-value vector of size 5 with 5 base hypotheses under the global null
data <- data_Gen(n = 1, p = 5, mu = rep(0,5), rho = -0.2, copula = 't', one_sided = T)
p.vec <- as.vector(data$p.mat)

## Compute the global p-value using the Inverse Gamma combination test with a shape parameter 2
p.global <- combination.test(p.vec, weights = 1:5, method = 'Inverse Gamma', tail.idx = 2)
```

### Truncated test
When there are 1's in base p-values and we apply Cauchy/t combination test, a truncated Cauchy/t should be used. In this case, we need to choose `truncate=True` and pick a scaling upper bound `truncate.threshold` for all base p-values. The default value for `truncate.threshold` is 0.99.

We provides an example using the truncated test as follows:
```
library(heavytailcombtest)

p.vec <- c(0,0,0,1,1)
p.global <- combination.test(p.vec, method = 't', tail.idx = 2, truncate = T, truncate.threshold = 0.99)
```



## Reference
- Wilson, D. J. (2019). The harmonic mean p-value for combining dependent tests. Proceedings of the National Academy of Sciences, 116(4), 1195-1200.[[Link](https://www.pnas.org/doi/abs/10.1073/pnas.1814092116)]
- Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures. Journal of the American Statistical Association, 115(529), 393-402.[[Link](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485)]
- Fang, Y., Chang, C., Park, Y., & Tseng, G. C. (2023). HEAVY-TAILED DISTRIBUTION FOR COMBINING DEPENDENT p-VALUES WITH ASYMPTOTIC ROBUSTNESS. Statistica Sinica, 33, 1115-1142. [[Link](https://www3.stat.sinica.edu.tw/statistica/J33N21/J33N2101/J33N2101.html)]
