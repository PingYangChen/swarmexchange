
# Discrete Discrimination Test

isomorphicCheck <- function(d1) {

	n <- nrow(d1); m <- ncol(d1)
	allHammings <- getHammings(d1)
	allHammings
	#(13/12)^m - 2*(35/32)^m + ((5/4)^m/(n^2))*(n + 2*sum(allHammings))
}

getHammings <- function(d1) {
	out <- matrix(0, nrow(d1), nrow(d1))
	for (i in 1:nrow(d1)) {
		for (j in 1:i) {
			if (i != j) { 
				hamming <- sum(abs(d1[i,] - d1[j,]) > 0) 
				out[i,j] <- hamming#(4/5)^hamming
			}
		}
	}
	out
}

R_afbar <- function(DESIGN, D_INFO, verbose = TRUE) {

	MODELS <- lapply(1:nrow(D_INFO$modelSpace), function(i) {
		tmp <- matrix(0, 1 + D_INFO$nMainEff + D_INFO$nTwofi, D_INFO$nFactor)
		for (j in 1:(D_INFO$nMainEff + D_INFO$nTwofi)) {
			if (j <= D_INFO$nMainEff) {
				tmp[j+1, D_INFO$modelSpace[i,j] + 1] <- 1;
			} else {
				tmpIdx <- D_INFO$nMainEff + 2*(j - D_INFO$nMainEff) - 1
				tmp[j+1, D_INFO$modelSpace[i, tmpIdx] + 1] <- 1
				tmp[j+1, D_INFO$modelSpace[i, tmpIdx + 1] + 1] <- 1
			}
		}
		tmp
	})
	#MODELS
	
	nModel <- nrow(D_INFO$modelSpace)
	
	EYE <- diag(nrow(DESIGN))
	valMat <- dimMat <- matrix(0, nModel, nModel)
	
	for (j in 1:nModel) {
		Xj <- RgetModelMatrix(DESIGN, MODELS[[j]])
		XXj <- t(Xj) %*% Xj
		if (rcond(XXj) > 1e-12) {
			IminusHj <- EYE - (Xj %*% solve(XXj) %*% t(Xj))
			for (i in 1:nModel) {
				if (i != j) {
					modelDiff <- RgetDiffIdx(MODELS[[i]], MODELS[[j]])
					Xij <- RgetModelMatrix(DESIGN, modelDiff)
					Mij <- t(Xij) %*% IminusHj %*% Xij
					detM <- det(Mij)
					if (detM > 0) { valMat[i,j] = log(detM)/nrow(modelDiff); dimMat[i,j] <- nrow(modelDiff) } 
					if (verbose) {
						#cat(paste0("\ni = ", i, ", j = ", j, "\n"))
						#cat(paste0("\nXj = \n"))
						#print(Xj)
						#cat(paste0("\nXij = \n"))
						#print(Xij)
						#cat(paste0("\nMij = \n"))
						#print(Mij)
						#cat(paste0("\ndet(Mij) = ", detM, "\n"))
					}
				}
			}
		}
	}
	#if (verbose) { print(valMat) }
	val <- sum(valMat)/(nModel*(nModel - 1.0))
	list(val = val, each_val = valMat, each_dim = dimMat)
}

RgetDiffIdx <- function(MODEL_I, MODEL_J) {
	
	tmp <- MODEL_I
	for (i in 1:nrow(MODEL_I)) {
		for (j in 1:nrow(MODEL_J)) {
			if (all(MODEL_I[i,] == MODEL_J[j,])) { tmp[i,] <- rep(-99, ncol(tmp)) }
		}
	}
	ret <- which(tmp[,1] != -99)
	matrix(tmp[ret,], length(ret), ncol(tmp))
}

RgetModelMatrix <- function(DESIGN, MODEL) {
	modelMat <- matrix(0, nrow(DESIGN), nrow(MODEL))
	for (i in 1:nrow(MODEL)) {
		effVec <- rep(1, nrow(DESIGN))
		if (any(MODEL[i,] == 1)) {
			idx <- which(MODEL[i,] == 1)
			for (j in 1:length(idx)) {
				effVec <- effVec * DESIGN[,idx[j]]
			}
		}
		modelMat[,i] <- effVec
	}
	modelMat
}
