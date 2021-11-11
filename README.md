# Particle Swarm Exchange Algorithm for Model-Discrimination Designs

This document is the user guide for using particle swarm exchange (PSE) algorithm (Chen et al., 2021+) and the exchange type algorithms (Meyer and Nachtsheim, 1995; Li and Wu, 1997) to generate the model-discrimination designs.  For the fundemental math expression of the model-discrimination design, please refer to the article and [**our catalog website**](https://pingyangchen.github.io/swarmexchange/).

To generate the model-discrimination designs, there are three steps:
1. Specify the design problem by the `rGetDesignInfo` function. 
2. Setup the parameters of the PSE algorithm by the  `rGetAlgInfo` function.
3. Run the PSE algorithm by the `rDiscreteDesignPSO` function.

To run the columnwise-pairwise (CP) exchange algorithm (Li and Wu, 1997) for searching the balanced design, use functions `rGetColPairInfo` and `rDiscreteDesignColPair` in Steps 2 and 3. To run the coordinate exchange algorithm (Meyer and Nachtsheim, 1995) for searching the unbalanced design, use functions `rGetCoorExInfo` and `rDiscreteDesignCoorEx` in Steps 2 and 3. 

Before implementing our `swarmexchange` codes, please run the following R script to make sure that the necessary R packages are installed.

```{r}
# Check whether the required packages have been installed or not.
# If not, then install them
pkgReq <- c("Rcpp", "RcppArmadillo")
pkgLoc <- installed.packages()[,1]
for (i in 1:length(pkgReq)) {
  if (!(pkgReq[i] %in% pkgLoc)) { install.packages(pkgReq[i]) }
}
```

#### Reference
- Chen, P.-Y., Chen, R.-B., Li, J.-P. and Li, W. W. (2021+). Particle Swarm Exchange Algorithms with Applications in Generating Optimal Model-Discrimination Designs. *Preprint*.
- Li, W., & Wu, C. F. J. (1997). Columnwise-pairwise algorithms with applications to the construction of supersaturated designs. Technometrics, 39 (2), 171–179.
- Meyer, R. K., & Nachtsheim, C. J. (1995). The coordinate-exchange algorithm for constructing exact optimal experimental designs. Technometrics, 37 (1), 60–69.


## Design Problem Specification

`rGetDesignInfo`: the function of specifying the model-discrimination problem. 

```{r}
rGetDesignInfo(typeCrit, n, m, mSpName = c("MEPI", "PMS"), g = 1, q = 1, labLevel = c(-1, 1), balance = 0)
```

| `rGetDesignInfo` <br> Input Option | Description | Default Value |
| ----------------------------- | ----------- | ------------- |
| `typeCrit` | **Model-discrimination Designs** <br> `1`: ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\overline{AF}) criterion <br> `2`: ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\overline{EPD}) criterion <br> `3`: Negative ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\overline{A^S}) criterion <br> `4`:  ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\overline{ENCP}) criterion <br> **Model-robust Designs** <br> `0`: Estimation capacity criterion <br> `5`:  Information capacity criterion | -- |
| `n`        | (Positive integer) Run size                     | -- |
| `m`        | (Positive integer) Number of factors            | -- |
| `mSpName`  | `'MEPI'`: MEPI space <br> `'PMS'`: PMS space    | -- |
| `g`        | (Positive integer) MEPI parameter               | -- |
| `q`        | (Positive integer) PMS parameter                | -- |
| `labLevel` | (vector) the labels of factor levels            | `c(-1, 1)` |
| `balance`  | `0`: unbalance design <br> `1`: balanced design | -- |



## Particle Swarm Exchange Algorithm

`rGetAlgInfo`: the function for configuring the particle swarm exchange algorithm.
```{r}
algInfo <- rGetAlgInfo(nSwarm, maxIter, MIX_C = 1, MIX_R = 1, JFO_RV = 1.0, JFO_R0 = 0.9, JFO_R1 = 0.3, JFO_RHO = 0.5)
```

| `rGetAlgInfo` <br> Input Option | Description | Default Value |
| ------------------------------- | ----------- | ------------- |
| `nSwarm`  | (Positive integer) the number of particles   | -- |
| `maxIter` | (Positive integer) the number of iterations  | -- |
| `MIX_C`   | (Positive integer) the number of columns to be exchanged in the COLMIX operator | 1 |
| `MIX_R`   | (Positive integer) the number of rows to be exchanged in the ROWMIX operator    | 1 |
| `JFO_R0`  | The value of ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\omega_{max}) &ast; | 0.9 |
| `JFO_R1`  | The value of ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\omega_{min}) &ast; | 0.3 |
| `JFO_RHO` | The value of ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\rho)        &ast; | 0.5 |
| `JFO_RV`  | the proportion of PSE updating iterations by updating  ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\omega^{(t)}) | 1.0 |

&ast; The parameters ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\omega_{max},\omega_{min}) and ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\rho) are used to compute the probabilities

![formula](https://render.githubusercontent.com/render/math?math=\color{white}\omega^{(t)}=\omega_{max}-\frac{t-1}{t_{max}-1}(\omega_{max}-\omega_{min}))

![formula](https://render.githubusercontent.com/render/math?math=\color{white}\omega^{(t)}_L=\rho(1-\omega^{(t)}),\omega^{(t)}_G=(1-\rho)(1-\omega^{(t)}))

to choose the target design to be mixed with the current design among random design, local best design and global best design. `JFO_RV` is between 0 and 1 and is used to denote the proportion of PSE updating iterations by updating ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\omega^{(t)}). For example, `JFO_RV=0.8` indicates that ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\omega^{(t)}) is updated by the above equations in the first 80% of the PSE iterations, and then, it is fixed as ![formula](https://render.githubusercontent.com/render/math?math=\color{white}\omega_{min}) for the remaining 20% iterations.



`rDiscreteDesignPSO`: the function for running the particle swarm exchange algorithm.

```{r}
rDiscreteDesignPSO(algInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE)
```

| `rDiscreteDesignPSO` <br> Input Option | Description | Default Value |
| ------------------------------- | ----------- | ------------- |
| `algInfo`     | R object generated by `rGetAlgInfo`  | -- |
| `designInfo`  | R object generated by `rGetDesignInfo`  | -- |
| `if_parallel` | (Boolean) to decide for using parallel computing | `TRUE` |
| `seed`        | (Positive integer) the random seed | `NULL` |
| `verbose`     | (Boolean) to decide for showing the updating message | `TRUE` |



## Exchange Type Algorithms

### Columnwise-pairwise Exchange Algorithm


`rGetColPairInfo`: the function for configuring the columnwise-pairwise exchange algorithm.

```{r}
rGetColPairInfo(maxIter = 50, nTry = 1, maximize = 1, CPk = 1, tol = 0)
```

| `rGetColPairInfo` <br> Input Option | Description | Default Value |
| ------------------------------- | ----------- | ------------- |
| `maxIter`  | (Positive integer) maximal number of iterations   | 50 |
| `nTry`     | (Positive integer) the number of replications  | 1 |
| `maximize` | `0`: minization <br> `1`: maximization | 1 |
| `CPk`      | (Positive integer) | 1 |
| `tol`      | (Positive real value) stopping criterion. Set `tol = 0` for disabling | 0 |


`rDiscreteDesignColPair`: the function for running the columnwise-pairwise exchange algorithm.

```{r}
rDiscreteDesignColPair(algInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE)
```

| `rDiscreteDesignColPair` <br> Input Option | Description | Default Value |
| ------------------------------- | ----------- | ------------- |
| `algInfo`     | R object generated by `rGetColPairInfo`  | -- |
| `designInfo`  | R object generated by `rGetDesignInfo`  | -- |
| `if_parallel` | (Boolean) to decide for using parallel computing | `TRUE` |
| `seed`        | (Positive integer) the random seed | `NULL` |
| `verbose`     | (Boolean) to decide for showing the updating message | `TRUE` |


### Coordinate Exchange Algorithm

`rGetCoorExInfo`: the function for configuring the coordinate exchange algorithm.

```{r}
rGetCoorExInfo(maxIter = 50, nTry = 1, maximize = 1, tol = 0)
```

| `rGetCoorExInfo` <br> Input Option | Description | Default Value |
| ------------------------------- | ----------- | ------------- |
| `maxIter`  | (Positive integer) maximal number of iterations   | 50 |
| `nTry`     | (Positive integer) the number of replications  | 1 |
| `maximize` | `0`: minization <br> `1`: maximization | 1 |
| `tol`      | (Positive real value) stopping criterion. Set `tol = 0` for disabling | 0 |


`rDiscreteDesignCoorEx`: the function for running the coordinate exchange algorithm.

```{r}
rDiscreteDesignCoorEx(algInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE)
```

| `rDiscreteDesignCoorEx` <br> Input Option | Description | Default Value |
| ------------------------------- | ----------- | ------------- |
| `algInfo`     | R object generated by `rGetCoorExInfo`  | -- |
| `designInfo`  | R object generated by `rGetDesignInfo`  | -- |
| `if_parallel` | (Boolean) to decide for using parallel computing | `TRUE` |
| `seed`        | (Positive integer) the random seed | `NULL` |
| `verbose`     | (Boolean) to decide for showing the updating message | `TRUE` |



## Illustrative Example

The R script [**run.R**](https://github.com/PingYangChen/swarmexchange/blob/master/run.R) is an example to use the PSE codes.  Suppose target on searching a balanced two-level ![formula](https://render.githubusercontent.com/render/math?math=\color{white}12\times4\overline{AF})-optimal design for discriminating among ![formula](https://render.githubusercontent.com/render/math?math=\color{white}MEPI_2) model space. 


The first step is to determine the model-discrimination design problem through the `rGetDesignInfo` function.

```{r}
designInfo <- rGetDesignInfo(typeCrit = 1, n = 12, m = 4, 
                             mSpName = "MEPI", g = 2, 
                             balance = 1)
```

Second, we choose for the PSE algorithm running with 100 iterations and 32 particles. In addition, the number of the columns for COLMIX operator is set as 1 and we disable the ROWMIX operator to maintain the balance design structure.

```{r}
algInfo <- rGetAlgInfo(nSwarm = 32, maxIter = 100, 
                       MIX_C = 1, MIX_R = 0)
```

The final step is to run PSE algorithm via the `rDiscreteDesignPSO` function. The necessary inputs are the R objects, `algInfo` and `designInfo`.  For multi-processor parallel computing, we set `if_parallel=TRUE`.  The output of the `rDiscreteDesignPSO` function is an R list object with two variables. The first one `result$RES` include the PSE searching result which is also an R list object. The variable `result$RES$GBest` is the global best design and `result$RES$fGBest` is its corresponding design criterion value. The computing time can be found in `result$CPUTIME` which is recorded in seconds.

```{r}
res <- rDiscreteDesignPSO(algInfo, designInfo, if_parallel = TRUE, 
                          seed = NULL, verbose = TRUE)
#names(res)
res$RES$fGBest # Global Best AF-bar-optimal criterion value
res$RES$GBest  # Global Best AF-bar-optimal design
```

The columnwise-pairwise (CP) exchange algorithm is also used for finding the balanced design. To use our code, the first step is the same that uses the `rGetDesignInfo` function to define the model-discrimination problem. Following the above example, the code for running CP algorithm is presented below.

```{r}
# Set Parameters for Columnwise-Pairwise (CP) algorithm
cpAlgInfo <- rGetColPairInfo(maxIter = 100, nTry = 32, 
                             CPk = 1)

# Run CP algorithm
cpRes <- rDiscreteDesignColPair(cpAlgInfo, designInfo, if_parallel = TRUE, 
                                seed = NULL, verbose = TRUE)
#names(cpRes)
cpRes$RES$DESIGN_VAL # Overall Best AF-bar-optimal criterion value among multiple trails
cpRes$RES$DESIGN     # Overall Best AF-bar-optimal design among multiple trails
```

For unbalanced model-discrimination designs, we can apply the PSE algorithm and coordinate exchange algorithm.  The workflow is similar and we present the example code as below. 

```{r}
### Unbalanced Design ###
# The example is the AF-bar-optimal unbalanced design 
# of n = 12, m = 4 for MEPI space with g = 2.
designInfo <- rGetDesignInfo(typeCrit = 1, n = 12, m = 4, 
                             mSpName = "MEPI", g = 2,
                             balance = 0)

## RUN PSE ALGORITHM ##
# Set Parameters for PSE algorithm
algInfo <- rGetAlgInfo(nSwarm = 32, maxIter = 100, 
                       MIX_C = 1, MIX_R = 1)

# Run PSE algorithm
res <- rDiscreteDesignPSO(algInfo, designInfo, if_parallel = TRUE, 
                          seed = NULL, verbose = TRUE)
#names(res)
res$RES$fGBest # Global Best AF-bar-optimal criterion value
res$RES$GBest  # Global Best AF-bar-optimal design


## RUN CE ALGORITHM ##
# Set Parameters for Coordinate Exchange (CE) algorithm
ceAlgInfo <- rGetCoorExInfo(maxIter = 100, nTry = 32)
# Run CE algorithm
ceRes <- rDiscreteDesignCoorEx(ceAlgInfo, designInfo, if_parallel = TRUE, 
                               seed = NULL, verbose = TRUE)
#names(ceRes)
ceRes$RES$DESIGN_VAL # Overall Best AF-bar-optimal criterion value among multiple trails
ceRes$RES$DESIGN     # Overall Best AF-bar-optimal design among multiple trails
```

