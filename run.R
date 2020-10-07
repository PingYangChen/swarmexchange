
# Check whether the required packages have been installed or not.
# If not, then install them
pkgReq <- c("Rcpp", "RcppArmadillo")
pkgLoc <- installed.packages()[,1]
for (i in 1:length(pkgReq)) {
  if (!(pkgReq[i] %in% pkgLoc)) { install.packages(pkgReq[i]) }
}

# Load packages and codes
library(Rcpp)
library(RcppArmadillo)
kernelPath <- "./kernel"
sourceCpp(file.path(kernelPath,"psoRcpp.cpp"), verbose = FALSE)
source(file.path(kernelPath, "rLaunchTools.r"))

# Begin Instruction

### Balanced Design ###
# The example runs for the balanced design of n = 12, m = 4 for MEPI space with g = 1.
designInfo <- rGetDesignInfo(typeCrit = 1, n = 20, m = 5, 
                             mSpName = "MEPI", g = 2,
                             balance = 1)

# Set SIDD algorithm
algInfo <- rGetAlgInfo(nSwarm = 32, maxIter = 100, PSO_UPDATE = 0,  
                       JFO_R0 = 0.9, JFO_R1 = 0.3, MIX_C = 1, MIX_R = 0,
                       HYBRIDEXALG = 1)
# Run SIDD algorithm
res <- rDiscreteDesignPSO(algInfo, designInfo, if_parallel = TRUE, 
                          seed = NULL, verbose = TRUE)
names(res)

res$RES$fGBest
res$RES$fPBest
res$RES$fGBestHist
res$RES$fPBestHist

res$RES$fPBestHist[,1:5]
res$RES$fPBestHist[,ncol(res$RES$fPBestHist)]

# Set Columnwise-Pairwise (CP) algorithm
cpAlgInfo <- rGetColPairInfo(maxIter = 100, nTry = 32, CPk = 1)
# Run CP algorithm
cpRes <- rDiscreteDesignColPair(cpAlgInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE)
names(cpRes)


### Non-regular Design ###
# The example runs for the non-regular design of n = 12, m = 4 for MEPI space with g = 1.
designInfo <- rGetDesignInfo(typeCrit = 1, n = 12, m = 4, 
                             mSpName = "MEPI", g = 1,
                             balance = 0)

# Set SIDD algorithm
algInfo <- rGetAlgInfo(nSwarm = 32, maxIter = 100, PSO_UPDATE = 0,  
                       JFO_R0 = 0.9, JFO_R1 = 0.3, MIX_C = 1, MIX_R = 1,
                       HYBRIDEXALG = 1)
# Run SIDD algorithm
res <- rDiscreteDesignPSO(algInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE)
names(res)

res$RES$fGBest

# Set JFO algorithm
jfoAlgInfo <- rGetAlgInfo(nSwarm = 32, maxIter = 100, PSO_UPDATE = 0,  
                          JFO_R0 = 0.9, JFO_R1 = 0.3)
# Run JFO algorithm
jfoRes <- rDiscreteDesignPSO(jfoAlgInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE)
names(jfoRes)

# Set Coordinate Exchange (CE) algorithm
ceAlgInfo <- rGetCoorExInfo(maxIter = 100, nTry = 32)
# Run CE algorithm
ceRes <- rDiscreteDesignCoorEx(ceAlgInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE)
names(ceRes)

ceRes$RES$DESIGN
ceRes$RES$DESIGN_VAL
