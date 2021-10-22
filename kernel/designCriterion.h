
// DECLARE FUNCTIONS
arma::mat DESIGNCRITERION(double &DESIGN_VAL, const mat &DESIGN, const DESIGN_INFO &D_INFO, int typeCrit);
double EstCapacityCalc(const mat &DESIGN, const DESIGN_INFO &D_INFO);
arma::mat getModelMatrix(const mat &DESIGN, const imat &MODEL, const DESIGN_INFO &D_INFO);
arma::imat getDiffIdx(double &pij, const imat &MODEL_I, const imat &MODEL_J);

// BODY
arma::mat DESIGNCRITERION(double &DESIGN_VAL, const mat &DESIGN, const DESIGN_INFO &D_INFO, int typeCrit)
{
	DESIGN_VAL = -1e20;
	int nModel = D_INFO.nModel;
	double nModel_double = (double)D_INFO.nModel;
	arma::mat EYE = arma::eye<arma::mat>(DESIGN.n_rows, DESIGN.n_rows);
	arma::mat valMat(nModel, nModel, fill::zeros);
	bool IS_NONSINGULAR;
	if (typeCrit == -1) typeCrit = D_INFO.typeCrit;
	switch (typeCrit) {
		case 0:
		{
			DESIGN_VAL = EstCapacityCalc(DESIGN, D_INFO);
			break;
		}
		case 1: 
		{ // AF_bar-criterion
			for (int j = 0; j < nModel; j++) {
				arma::mat Xj = getModelMatrix(DESIGN, D_INFO.modelIndices.slice(j), D_INFO);
				arma::mat XXj = Xj.t() * Xj;
				if (arma::rcond(XXj) < 1e-18) { break; }
				arma::mat XXjInv(XXj.n_rows, XXj.n_rows, fill::zeros);
				IS_NONSINGULAR = arma::inv(XXjInv, XXj);
				if (!IS_NONSINGULAR) { break; }
  			arma::mat IminusHj = EYE - (Xj * XXjInv * Xj.t());
				for (int i = 0; i < nModel; i++) {
					if (i != j) {
						double pij; 
						arma::imat modelDiff = getDiffIdx(pij, D_INFO.modelIndices.slice(i), D_INFO.modelIndices.slice(j));
						arma::mat Xij = getModelMatrix(DESIGN, modelDiff, D_INFO);
						arma::mat Mij = Xij.t() * IminusHj * Xij;					
						double detM = arma::det(Mij);
						if (detM > 0) { 
							valMat(i, j) = std::log(detM)/pij; 
						} /*else {
							valMat.zeros(); i += nModel; j += nModel;
						} */
					}  
				}
			}
			DESIGN_VAL = arma::accu(valMat)/(nModel_double*(nModel_double - 1.0));
			break;
		}
		case 2:
		{	// EPD bar
			double nRun_double = (double)DESIGN.n_rows;
			for (int i = 0; i < nModel; i++) {
				arma::mat Xi = getModelMatrix(DESIGN, D_INFO.modelIndices.slice(i), D_INFO);
				arma::mat XXi = Xi.t() * Xi;
				if (arma::rcond(XXi) < 1e-18) { break; }
				arma::mat XXiInv(XXi.n_rows, XXi.n_rows, fill::zeros);
				IS_NONSINGULAR = arma::inv(XXiInv, XXi);
				if (!IS_NONSINGULAR) { break; }
				arma::mat Hi = Xi * XXiInv * Xi.t();
				for (int j = 0; j < i; j++) {
					arma::mat Xj = getModelMatrix(DESIGN, D_INFO.modelIndices.slice(j), D_INFO);
					arma::mat XXj = Xj.t() * Xj;
					if (arma::rcond(XXj) < 1e-18) { break; }
					arma::mat XXjInv(XXj.n_rows, XXj.n_rows, fill::zeros);
					IS_NONSINGULAR = arma::inv(XXjInv, XXj);
					if (!IS_NONSINGULAR) { break; }
					arma::mat HiMinusHj = Hi - (Xj * XXjInv * Xj.t());
					arma::mat Dij = HiMinusHj * HiMinusHj;
					double trD = arma::trace(Dij);
					if (trD > 0) valMat(i,j) = trD/nRun_double;
					valMat(j,i) = valMat(i,j);
				}
			}
			DESIGN_VAL = arma::accu(valMat)/(nModel_double*(nModel_double - 1.0));
			break;
		}
		case 3:
		{	// Neg As bar
			valMat.fill(-1e20);
			valMat.diag().zeros();
			for (int j = 0; j < nModel; j++) {
				arma::mat Xj = getModelMatrix(DESIGN, D_INFO.modelIndices.slice(j), D_INFO);
				arma::mat XXj = Xj.t() * Xj;
				if (arma::rcond(XXj) < 1e-18) { break; }
				arma::mat XXjInv(XXj.n_rows, XXj.n_rows, fill::zeros);
				IS_NONSINGULAR = arma::inv(XXjInv, XXj);
				if (!IS_NONSINGULAR) { break; }
				arma::mat IminusHj = EYE - (Xj * XXjInv * Xj.t());
				for (int i = 0; i < nModel; i++) {
					if (i != j) {
						double pij; 
						arma::imat modelDiff = getDiffIdx(pij, D_INFO.modelIndices.slice(i), D_INFO.modelIndices.slice(j));
						arma::mat Xij = getModelMatrix(DESIGN, modelDiff, D_INFO);
						arma::mat Mij = Xij.t() * IminusHj * Xij;
						double detM = arma::det(Mij);
						if ((arma::rcond(Mij) > 1e-18) & (detM > 0)) { valMat(i,j) = (-1.0)*arma::trace(Mij.i())/pij; } 
					}
				}
			}
			DESIGN_VAL = arma::accu(valMat)/(nModel_double*(nModel_double - 1.0));
			break;
		}
		case 4:
		{	// ENCP bar
			for (int j = 0; j < nModel; j ++) {
				arma::mat Xj = getModelMatrix(DESIGN, D_INFO.modelIndices.slice(j), D_INFO);
				arma::mat XXj = Xj.t() * Xj;
				if (arma::rcond(XXj) < 1e-18) { break; }
				arma::mat XXjInv(XXj.n_rows, XXj.n_rows, fill::zeros);
				IS_NONSINGULAR = arma::inv(XXjInv, XXj);
				if (!IS_NONSINGULAR) {
  				break;
				}
				arma::mat IminusHj = EYE - (Xj * XXjInv * Xj.t());
				for (int i = 0; i < nModel; i ++) {
					if (i != j) {
						double pij; 
						arma::imat modelDiff = getDiffIdx(pij, D_INFO.modelIndices.slice(i), D_INFO.modelIndices.slice(j));
						arma::mat Xij = getModelMatrix(DESIGN, modelDiff, D_INFO);
						arma::mat Mij = Xij.t() * IminusHj * Xij;
						double trM = arma::trace(Mij);
						if (trM > 0) valMat(i,j) = trM/pij;
					}
				}
			}
			DESIGN_VAL = arma::accu(valMat)/(nModel_double*(nModel_double - 1.0));
			break;
		}
		case 5:
		{ // Information Criterion
			for (int j = 0; j < nModel; j++) {
				arma::mat Xj = getModelMatrix(DESIGN, D_INFO.modelIndices.slice(j), D_INFO);
				arma::mat XXj = Xj.t() * Xj;
				double localDetApprox = std::pow((double)Xj.n_rows, (double)Xj.n_cols);
				double detM = arma::det(XXj);
				valMat(j, 0) = detM/localDetApprox; 
			}
			DESIGN_VAL = arma::accu(valMat.col(0))/nModel_double;
			break;
		}
	}	
	return valMat;
}

double EstCapacityCalc(const mat &DESIGN, const DESIGN_INFO &D_INFO)
{
	double val = 0;
	int nModel = D_INFO.nModel;
	for (int i = 0; i < nModel; i++) {
		arma::mat Xi = getModelMatrix(DESIGN, D_INFO.modelIndices.slice(i), D_INFO);
		arma::mat inforMat = Xi.t() * Xi;
		if (arma::rcond(inforMat) > 1e-18) val++; 
	}
	val = val/((double)nModel);
	return val;
}

arma::imat getDiffIdx(double &pij, const imat &MODEL_I, const imat &MODEL_J)
{
	imat tmpModel = MODEL_I;
	for (uword i = 0; i < MODEL_I.n_rows; i++) {
		for (uword j = 0; j < MODEL_J.n_rows; j++) {
			if (all(MODEL_I.row(i) == MODEL_J.row(j))) {
				tmpModel.row(i).fill(-99); j = MODEL_J.n_rows;
			}
		}
	}
	imat modelDiff = tmpModel.rows(find(tmpModel.col(0) != -99));
	pij = (double)(modelDiff.n_rows);
	return modelDiff;
}

arma::mat getModelMatrix(const mat &DESIGN, const imat &MODEL, const DESIGN_INFO &D_INFO)
{
	mat modelMat(DESIGN.n_rows, MODEL.n_rows);
	vec effVec(DESIGN.n_rows, fill::ones);
	for (uword i = 0; i < MODEL.n_rows; i++) {
		effVec.fill(1);
		if (any(MODEL.row(i) == 1)) {
			uvec idx = find(MODEL.row(i) == 1);
			for (uword j = 0; j < idx.n_elem; j++) {
				effVec = effVec % DESIGN.col(idx(j));
			}
		}
		modelMat.col(i) = effVec;
	}
	return modelMat;
}