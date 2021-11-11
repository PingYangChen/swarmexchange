# Particle Swarm Exchange Algorithm for Model-Discrimination Designs

This document is the user guide for using particle swarm exchange (PSE) algorithm (Chen et al., 2021+) and the exchange type algorithms (Meyer and Nachtsheim, 1995; Li and Wu, 1997) to generate the model-discrimination designs.  For the fundemental math expression of the model-discrimination design, please refer to the article and [**our catalog website**](https://pingyangchen.github.io/swarmexchange/).

To generate the model-discrimination designs, there are three steps:
1. Specify the design problem by the `rGetDesignInfo` function. 
2. Setup the parameters of the PSE algorithm by the  `rGetAlgInfo` function.
3. Run the PSE algorithm by the `rDiscreteDesignPSO` function.

To run the columnwise-pairwise (CP) exchange algorithm (Li and Wu, 1997) for searching the balanced design, use functions `rGetColPairInfo` and `rDiscreteDesignColPair` in Steps 2 and 3. To run the coordinate exchange algorithm (Meyer and Nachtsheim, 1995) for searching the unbalanced design, use functions `rGetCoorExInfo` and `rDiscreteDesignCoorEx` in Steps 2 and 3. 


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

**run.R**



The R script **run.R** is an example to use the PSE codes. The first step is to determine the model-discrimination design problem through  `rGetDesignInfo` function.






In PSE algorithm, user first needs to decide the number of particles, \texttt{nSwarm}, and the number of iterations, \texttt{maxIter}. Larger values of them may result the better designs, but the computational cost is higher.
\texttt{MIX\_C} is the number of columns to be exchanged in the COLMIX operator, and its default value is . Please see Section \ref{sec:pseconf} for details.  \texttt{MIX\_R} is the number of rows to be exchanged in the ROWMIX operator and its default value is fixed as 1. The parameters, \texttt{JFO\_R1}, \texttt{JFO\_R0} and \texttt{JFO\_RHO}, are the probability values, $\omega_{max}$, $\omega_{min}$ and $\rho$ for the MIX operator shown in \eqref{eq:inertia1} and \eqref{eq:inertia2}. Their default values are set as follows, $\omega_{max} = 0.9$, $\omega_{min} = 0.3$ and $\rho = 0.5$. 
Suppose we want to implement the PSE algorithm with 1000 iterations and 100 particles. In addition, the number of the columns for COLMIX operator is set as 4 and we decide not to involve ROWMIX operator. The others are fixed as the defaulty values. Therefore, we should have the PSE information as follows.\\ 
%\noindent
\texttt{algInfo <- rGetAlgInfo(100, 1000, 4, 0)}\\



In the GitHub repository, we also
developed the codes for CE and CP algorithms in functions `rDiscreteDesignCoorEx` and `rDiscreteDesignColPair` respectively. To avoid local convergence, it is recommended using multiple initial designs to run CE and CP algorithms. The users can specify the number of initial designs via the option \texttt{nTry} into the algorithms' configuration functions, `rGetCoorExInfo` and `rGetColPairInfo`, and then, the CE and CP functions can run in parallel for these independent trails.  To use CE and CP algorithms, the users can find the example codes in **run.R** file. 


```{r}
ceInfo <- rGetCoorExInfo(maxIter = 50, nTry = 5) 
ceResult <- rDiscreteDesignCoorEx(ceInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE) 
ceResult$RES$DESIGN      # The best design across 'nTry' trails 
ceResult$RES$DESIGN_VAL # The corresponding design criterion value 
```

```{r}
cpInfo <- rGetColPairInfo(maxIter = 50, nTry = 5, CPk = 1)
cpResult <- rDiscreteDesignColPair(cpInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE) 
cpResult$RES$DESIGN     # The best design across 'nTry' trails 
cpResult$RES$DESIGN_VAL # The corresponding design criterion value
```