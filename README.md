# heavytailcombtest

This is an R package for global testing applying the heavy-tailed combination test. This (weighted) heavy-tailed combination test is proposed by [Fang et al., 2023](https://www3.stat.sinica.edu.tw/statistica/J33N21/J33N2101/J33N2101.html). Here we include the combination test with commonly-used heavy-tailed (regularly-varying) distributions with a changable `tail.idx`(shape parameter for the distribution): t, Pareto, Frechet, inverse Gamma, and Levy. 
This method also covers widely-applied [Cauchy combination test](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485) and [harmonic mean p-value](https://www.pnas.org/doi/abs/10.1073/pnas.1814092116).


## Installation
```{r}
library(devtools)

install_github("gl-ybnbxb/heavytailcombtest")
```

## Example


## Reference
- Wilson, D. J. (2019). The harmonic mean p-value for combining dependent tests. Proceedings of the National Academy of Sciences, 116(4), 1195-1200.[[Link](https://www.pnas.org/doi/abs/10.1073/pnas.1814092116)]
- Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures. Journal of the American Statistical Association, 115(529), 393-402.[[Link](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485)]
- Fang, Y., Chang, C., Park, Y., & Tseng, G. C. (2023). HEAVY-TAILED DISTRIBUTION FOR COMBINING DEPENDENT p-VALUES WITH ASYMPTOTIC ROBUSTNESS. Statistica Sinica, 33, 1115-1142. [[Link](https://www3.stat.sinica.edu.tw/statistica/J33N21/J33N2101/J33N2101.html)]
