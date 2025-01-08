# heavytailcombtest

This is an R package for global testing applying the heavy-tailed combination test. P-values are first transformed into heavy-tailed distributed variables and are then aggregated. Here we provide a function that users can perform heavy-tailed combination tests using a wide class of regularly varying tail distributions with different tail indexes. Specifically, distributions that we have implemented include t, left-truncated t, Pareto, Frechet, inverse Gamma, and Levy. We also allow users to provide pre-chosen weights and our method also include widely-applied [Cauchy combination test](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485) and [harmonic mean p-value](https://www.pnas.org/doi/abs/10.1073/pnas.1814092116). For details, please refer to [our paper](https://arxiv.org/abs/2310.20460). 

As a practical guidance, our default setting is to use left-truncated t with tail index 1 (which is a left-truncated Cauchy distribution) as the transformation for the combination test, with the truncation threshold $p_0 = 0.9$. As suggested in the paper, to guarantee type-I error control in practice, the tail index should be no larger than $1$.


## Installation
```{r}
library(devtools)

install_github("gl-ybnbxb/heavytailcombtest")
```

## Example

The `combination.test` takes in a base p-value vector and returns a global p-value. 

### Simulate base p-values
As an example, we simulate an one-sided p-value vector with the number of base hypotheses $= 5$ when the global null holds. The base test statistics are sampled from multivariate t with the location parameter $\mu=0$ and $\Sigma_\rho$ which has 1s on the diagonal and a common $\rho$ off the diagonal.
```
library(heavytailcombtest)

## Simulate a one-sided p-value vector of size 5 with 5 base hypotheses under the global null
data <- data_Gen(n = 1, p = 5, mu = rep(0,5), rho = -0.2, copula = 't', one_sided = T)
p.vec <- as.vector(data$p.mat)
```

### Apply the heavy-tailed combination test
By default, we recommend using the left-truncated t-distribution-based combination test. Denote the truncation threshold as $p_0$. Then, the truncated t distribution has support $[c, +\infty)$ where $c$ is the $1-p_0$ quantile of the corresponding regular student t distribution.
```
p.global <- combination.test(p.vec, method = 't', tail.idx = 1, truncate = T, truncate.threshold = 0.9)
```

We can also apply other distributions, for example, the inverse Gamma distribution for the combination test and provide a pre-chosen weight vector, proportional to the importance and efficiency of each individual hypothesis.
```
p.global <- combination.test(p.vec, weights = 1:5, method = 'Inverse Gamma', tail.idx = 1)
```

### Compute the Cauchy combination test or the harmonic mean p-value
To apply the Cauchy combination test using our code, select the student t distribution, tail index $1$ and weights as $1/n$ where $n$ is the number of individual hypotheses:
```
p.global <- combination.test(p.vec, weights = rep(1, 5)/5, method = 't', truncate = F, tail.idx = 1)
```

To apply the harmonic mean p-value using our code, select the Pareto distribution, tail index $1$ and weights as $1/n$ where $n$ is the number of individual hypotheses:
```
p.global <- combination.test(p.vec, weights = rep(1, 5)/5, method = 'pareto', tail.idx = 1)
```

