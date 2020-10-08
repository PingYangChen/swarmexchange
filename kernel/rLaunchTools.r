rGetDesignInfo <- function(typeCrit = 1, n = 12, m = 4, 
													 mSpName = c("MEPI", "PMS"), g = 2, q = 3, 
													 labLevel = c(-1, 1), balance = 0) {

	if (mSpName == "MEPI") {
		twoFiSet <- combn(m, 2); twoFiSp <- combn(ncol(twoFiSet), g)
		modelSp <- t(sapply(1:ncol(twoFiSp), function(i) c(1:m, as.vector(twoFiSet[,twoFiSp[,i]])) ))
		nMainEff <- m; nTwofi <- g
	} else if (mSpName == "PMS") {
		mainEffSet <- combn(m, q)
		modelSp <- t(sapply(1:ncol(mainEffSet), function(i) c(mainEffSet[,i], as.vector(combn(mainEffSet[,i], 2))) ))
		nMainEff <- q; nTwofi <- choose(q, 2)
	} else { stop("Model space should be one of 'MEPI' and 'PMS'.") }

	modelSp <- modelSp - 1 # for Cpp indexing
	
	modelVec <- numeric(ncol(modelSp)*nrow(modelSp)^2)
	modelDifVec <- numeric(ncol(modelSp)*nrow(modelSp)^2)
	modelPij <- rep(nMainEff + nTwofi, nrow(modelSp)^2)
	for (i in 1:nrow(modelSp)) {
		for (j in 1:nrow(modelSp)) {
			locPij <- (i - 1)*nrow(modelSp) + (j - 1) + 1
			locStart <- (i - 1)*(ncol(modelSp)*nrow(modelSp)) + (j - 1)*ncol(modelSp) + 1
			diffIdx <- modelSp[i,]
			
			for (im in 1:nMainEff) {
				tm <- 0; flag <- 0
				while((flag == 0) & (tm < nMainEff)) {
					tm <- tm + 1
					if (modelSp[i, im] == modelSp[j, tm]) { 
						diffIdx[im] <- -99; modelPij[locPij] <- modelPij[locPij] - 1; flag <- 1 
					}
				}
			}
			for (ig in 1:nTwofi) {
				tg <- 0; flag <- 0
				tarIdx <- (nMainEff+(ig-1)*2+1):(nMainEff+ig*2)
				while((flag == 0) & (tg < nTwofi)) {
					tg <- tg + 1
					scIdx <- (nMainEff+(tg-1)*2+1):(nMainEff+tg*2)
					if (all(modelSp[i, tarIdx] == modelSp[j, scIdx])) { 
						diffIdx[tarIdx] <- -99; modelPij[locPij] <- modelPij[locPij] - 1; flag <- 1 
					}
				}
			}
			modelVec[locStart:(locStart + ncol(modelSp) - 1)] <- modelSp[j,]
			modelDifVec[locStart:(locStart + ncol(modelSp) - 1)] <- diffIdx
		}
	}

	list(typeCrit = typeCrit, nRun = n, nFactor = m, nLevel = length(labLevel), labLevel = labLevel, balance = balance,
			 nModel = nrow(modelSp), nMainEff = nMainEff, nTwofi = nTwofi, 
			 modelSpace = modelSp)#, modelPij = modelPij, idxVecHj = modelVec, idxVecXij = modelDifVec)
}

rGetAlgInfo <- function(nSwarm = 8, maxIter = 50, PSO_UPDATE = 0, maximize = 1, tol = 0,
												MIX_C = 1, MIX_R = 1, JFO_RV = 1.0, JFO_R0 = 0.9, JFO_R1 = 0.3,
												JFO_RHO = 0.5, HYBRIDEXALG = 0) {
	list(nSwarm = nSwarm, maxIter = maxIter, PSO_UPDATE = PSO_UPDATE, maximize = maximize, tol = tol,
			 MIX_C = MIX_C, MIX_R = MIX_R, JFO_RV = JFO_RV, JFO_R0 = JFO_R0, JFO_R1 = JFO_R1,
			 JFO_RHO = JFO_RHO, HYBRIDEXALG = HYBRIDEXALG, 
			 version = '20201007')
}

rGetCoorExInfo <- function(maxIter = 50, nTry = 1, maximize = 1, tol = 0) {
	list(maxIter = maxIter, nTry = nTry, maximize = maximize, tol = tol, version = '20201007')
}

rGetColPairInfo <- function(maxIter = 50, nTry = 1, maximize = 1, CPk = 1, tol = 0) {
	list(maxIter = maxIter, nTry = nTry, maximize = maximize, CPk = CPk, tol = tol, version = '20201007')
}

rDiscreteDesignPSO <- function(algInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE) {
	set.seed(seed)
	cputime <- system.time(
		out <- DiscreteDesignPSO(algInfo, designInfo, if_parallel, verbose)  
	)[3]
	list(RES = out, CPUTIME = cputime)
}

rDiscreteDesignCoorEx <- function(algInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE) {
	set.seed(seed)
	cputime <- system.time(
		out <- DiscreteDesignCoorEx(algInfo, designInfo, if_parallel, verbose)  
	)[3]
	list(RES = out, CPUTIME = cputime)
}

rDiscreteDesignColPair <- function(algInfo, designInfo, if_parallel = TRUE, seed = NULL, verbose = TRUE) {
	set.seed(seed)
	if (designInfo$balance != 1) stop("Columnwise-Pairwise Algorithm works for Balance design.")
	if (length(designInfo$labLevel) != 2) stop("Columnwise-Pairwise Algorithm works for 2-level Balance Design.")
	cputime <- system.time(
		out <- DiscreteDesignColPair(algInfo, designInfo, if_parallel, verbose)  
	)[3]
	list(RES = out, CPUTIME = cputime)
}
