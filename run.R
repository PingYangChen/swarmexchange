
codePath <- "SPECIFY THE PATH OF THE 'swarmexchange' FOLDER"
codePath <- "D:\\rProject\\swarmexchange"

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
kernelPath <- file.path(codePath, "kernel")
sourceCpp(file.path(kernelPath,"psoRcpp.cpp"), verbose = FALSE)
source(file.path(kernelPath, "rLaunchTools.r"))

# Begin Instruction

### Balanced Design ###
# The example is the AF-bar-optimal balanced design 
# of n = 12, m = 4 for MEPI space with g = 2.
designInfo <- rGetDesignInfo(typeCrit = 1, n = 12, m = 4, 
                             mSpName = "MEPI", g = 2, 
                             balance = 1)

## RUN PSE ALGORITHM ##
# Set Parameters for PSE algorithm
algInfo <- rGetAlgInfo(nSwarm = 32, maxIter = 100, 
                       MIX_C = 1, MIX_R = 0)

# Run PSE algorithm
res <- rDiscreteDesignPSO(algInfo, designInfo, if_parallel = TRUE, 
                          seed = NULL, verbose = TRUE)
#names(res)
res$RES$fGBest # Global Best AF-bar-optimal criterion value
res$RES$GBest  # Global Best AF-bar-optimal design

## RUN CP ALGORITHM ##
# Set Parameters for Columnwise-Pairwise (CP) algorithm
cpAlgInfo <- rGetColPairInfo(maxIter = 100, nTry = 32, 
                             CPk = 1)

# Run CP algorithm
cpRes <- rDiscreteDesignColPair(cpAlgInfo, designInfo, if_parallel = TRUE, 
                                seed = NULL, verbose = TRUE)
#names(cpRes)
cpRes$RES$DESIGN_VAL # Overall Best AF-bar-optimal criterion value among multiple trails
cpRes$RES$DESIGN     # Overall Best AF-bar-optimal design among multiple trails


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
