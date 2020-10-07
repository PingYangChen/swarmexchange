
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
HYBRIDEXALG_SET <- c(0, 1)
N_SET <- c(12, 16, 20)
M_SET <- c(4, 5, 6)
G_SET <- c(1, 2)

CONFIG_MESH <- c()
for (k1 in 1:length(HYBRIDEXALG_SET)) {
  for (k2 in 1:length(N_SET)) {
    for (k3 in 1:length(M_SET)) {
      for (k4 in 1:length(G_SET)) {
        CONFIG_MESH <- rbind(CONFIG_MESH, c(HYBRIDEXALG_SET[k1], N_SET[k2], M_SET[k3], G_SET[k4]))
      }
    }
  }
}
colnames(CONFIG_MESH) <- c("HYBRIDEXALG", "N", "M", "G")

RESULT_MAT <- matrix(0, nrow(CONFIG_MESH), 4)
colnames(RESULT_MAT) <- c("fGBEST", "ITER_REACH_MAX", "N_PBEST_SPOT_GBEST", "CPU_TIME")
for (i in 1:nrow(CONFIG_MESH)) {
  i=28
  hybridexalg <- CONFIG_MESH[i,1]
  n <- CONFIG_MESH[i,2]
  m <- CONFIG_MESH[i,3]
  g <- CONFIG_MESH[i,4]
  designInfo <- rGetDesignInfo(typeCrit = 1, n = n, m = m, mSpName = "MEPI", g = g, balance = 1)
  
  # Set SIDD algorithm
  algInfo <- rGetAlgInfo(nSwarm = 32, maxIter = 100, PSO_UPDATE = 0,  
                         JFO_R0 = 0.9, JFO_R1 = 0.3, MIX_C = 1, MIX_R = 0,
                         HYBRIDEXALG = hybridexalg, HYBRIDEXALG_A = 0.5, HYBRIDEXALG_B = 25.0)
  # Run SIDD algorithm
  res <- rDiscreteDesignPSO(algInfo, designInfo, if_parallel = TRUE, 
                            seed = 1, verbose = TRUE)  
  
  ITER_REACH_MAX <- min(which(res$RES$fGBestHist == res$RES$fGBest))
  N_PBEST_SPOT_GBEST <- sum(res$RES$fPBest == res$RES$fGBest)
  RESULT_MAT[i,] <- c(res$RES$fGBest, ITER_REACH_MAX, N_PBEST_SPOT_GBEST, res$CPUTIME)
}

OUTPUT_TABLE <- cbind(CONFIG_MESH, RESULT_MAT)

testPath <- "./test"
testFileName <- paste0(testPath, "//test_useHybrid_", strftime(Sys.time(), "%Y%m%d%H%M%S"), ".csv")
write.csv(OUTPUT_TABLE, testFileName, row.names = FALSE)
