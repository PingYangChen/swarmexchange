
// DECLARE FUNCTIONS
double MIX_OPERATOR(const int &MIX_PROC, arma::mat &MIX, const arma::mat &DESIGN_A, const double &VAL_A, const arma::mat &DESIGN_B, 
										const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS, const PSO_DYN &PSO_DYN);
//double HAMMING_TWOLEVEL(arma::mat &DESIGN);

// BODY
double MIX_OPERATOR(const int &MIX_PROC, arma::mat &MIX, const arma::mat &DESIGN_A, const double &VAL_A, const arma::mat &DESIGN_B, 
										const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS, const PSO_DYN &PSO_DYN)
{
	int typeCrit = D_INFO.typeCrit;
	double MIX_VAL = VAL_A; double DESIGN_VAL;
	MIX = DESIGN_A;
	uword DELETE_WHO = 0; uword ADD_WHO = 0; arma::mat ADD_ONE_BACK = MIX;
	double DUP_VAL;
	if (maximize == 0) DUP_VAL = 1e20; else DUP_VAL = -1e20;
	double EXALG_PROB = PSO_DYN.EXALG_PROB; double DICE;

	switch (MIX_PROC) {
		case 0: 
		{	// Column Exchange 
			int MIX_C = PSO_OPTS.MIX_C;
			arma::vec valVec(MIX.n_cols, fill::zeros);
			arma::vec countVec(MIX.n_cols, fill::zeros);
			int nModel = D_INFO.nModel;
			for (int q = 0; q < MIX_C; q++) {
				// Deletion
				arma::mat DESIGN_VAL_MAT = DESIGNCRITERION(MIX_VAL, MIX, D_INFO, -1);
				valVec.fill(arma::accu(DESIGN_VAL_MAT));
				if ((typeCrit == 0) | (typeCrit == 5)) {
					countVec.fill((double)(D_INFO.nModel));
					for (int j = 0; j < nModel; j++) {
						arma::imat tmpModel = D_INFO.modelIndices.slice(j);
						for (uword k = 0; k < tmpModel.n_cols; k++) {
							if (arma::accu(tmpModel.col(k)) > 1) {
								valVec(k) -= DESIGN_VAL_MAT(j, 0); countVec(k) -= 1.0;
							}
						}
					}
				} else {
					countVec.fill((double)(D_INFO.nModel*(D_INFO.nModel - 1)));
					for (int j = 0; j < nModel; j++) {
						for (int i = 0; i < nModel; i++) {
							if (i != j) {
								double pij; 
								arma::imat modelDiff = getDiffIdx(pij, D_INFO.modelIndices.slice(i), D_INFO.modelIndices.slice(j));
								for (uword k = 0; k < modelDiff.n_cols; k++) {
									if (arma::any(modelDiff.col(k) == 1)) {
										valVec(k) -= DESIGN_VAL_MAT(i, j); countVec(k) -= 1.0;
									}
								}
							}
						}
					}
				}
				valVec = valVec/countVec;
				//vecPrintf(valVec); Rprintf("\n"); 
				if (maximize == 0) valVec.min(DELETE_WHO); else valVec.max(DELETE_WHO);
				// Addition
				valVec.fill(0); 
				for (uword i = 0; i < DESIGN_B.n_cols; i++) {
					ADD_ONE_BACK = MIX;
					if (arma::all(ADD_ONE_BACK.col(DELETE_WHO) == DESIGN_B.col(i))) {
						valVec(i) = DUP_VAL;
					} else {
						ADD_ONE_BACK.col(DELETE_WHO) = DESIGN_B.col(i);
						DESIGNCRITERION(DESIGN_VAL, ADD_ONE_BACK, D_INFO, -1);
						valVec(i) = DESIGN_VAL;
					}
				}
				//vecPrintf(valVec); Rprintf("\n");
				uword ADD_WHO;
				if (maximize == 0) MIX_VAL = valVec.min(ADD_WHO); else MIX_VAL = valVec.max(ADD_WHO);
				MIX.col(DELETE_WHO) = DESIGN_B.col(ADD_WHO);
				// Decide to use exchange algorithms
				DICE = arma::as_scalar(randu(1, 1));
				if (DICE < EXALG_PROB) {
					if (PSO_OPTS.HYBRIDEXALG == 1) {
						if (D_INFO.balance == 1) {
							// Use CP-algorithm
							MIX = ColumnPair_CORE(MIX, MIX_VAL, D_INFO, DELETE_WHO, maximize);
						} else {
							// Use Coordinate Exchange
							for (int row_i = 0; row_i < (int)MIX.n_rows; row_i++) {
								MIX = CoorExchange_CORE(MIX, MIX_VAL, D_INFO, row_i, DELETE_WHO, maximize);
							}
						}
					}
				}
			}
			break;
		}
		case 1: 
		{	// Row exchange
			int MIX_R = PSO_OPTS.MIX_R;
			arma::vec valVec(MIX.n_rows, fill::zeros);
			arma::mat DELETE_ONE_ROW(MIX.n_rows - 1, MIX.n_cols, fill::zeros);
			arma::uvec RID = arma::linspace<arma::uvec>(0, MIX.n_rows - 1, MIX.n_rows);		
			for (int q = 0; q < MIX_R; q++) {
				// Deletion
				valVec.fill(0);
				for (uword i = 0; i < MIX.n_rows; i++) {
					DELETE_ONE_ROW = MIX.rows(arma::find(RID != i));
					DESIGNCRITERION(DESIGN_VAL, DELETE_ONE_ROW, D_INFO, -1);
					valVec(i) = DESIGN_VAL;
				}
				if (maximize == 0) valVec.min(DELETE_WHO); else valVec.max(DELETE_WHO);
				// Addition
				valVec.fill(0); 
				for (uword i = 0; i < DESIGN_B.n_rows; i++) {
					ADD_ONE_BACK = MIX;
					if (arma::all(ADD_ONE_BACK.row(DELETE_WHO) == DESIGN_B.row(i))) {
						valVec(i) = DUP_VAL;
					} else {
						ADD_ONE_BACK.row(DELETE_WHO) = DESIGN_B.row(i);
						DESIGNCRITERION(DESIGN_VAL, ADD_ONE_BACK, D_INFO, -1);
						valVec(i) = DESIGN_VAL;
					}
				}
				if (maximize == 0) MIX_VAL = valVec.min(ADD_WHO); else MIX_VAL = valVec.max(ADD_WHO);
				MIX.row(DELETE_WHO) = DESIGN_B.row(ADD_WHO);
				// Decide to use exchange algorithms
				DICE = arma::as_scalar(randu(1, 1));
				if (DICE < EXALG_PROB) {			
					if (PSO_OPTS.HYBRIDEXALG == 1) {
						// Use Coordinate Exchange
						for (int col_j = 0; col_j < (int)MIX.n_cols; col_j++) {
							MIX = CoorExchange_CORE(MIX, MIX_VAL, D_INFO, DELETE_WHO, col_j, maximize);
						}
					}
				}				
			}
			break;
		}
	}
	return(MIX_VAL);
}

/*
double HAMMING_TWOLEVEL(arma::mat &DESIGN) 
{
	int n = (int)DESIGN.n_rows; //double n_double = (double)n;
	int m = (int)DESIGN.n_cols; //double m_double = (double)m;
	double HAMMING = 0;
	arma::rowvec rowDiff(m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			rowDiff = arma::abs(DESIGN.row(i) - DESIGN.row(j));
			uvec hasDiff = find(rowDiff > 0);
			HAMMING += (double)hasDiff.n_elem; //std::pow(0.8, (double)hasDiff.n_elem);
		}
	}
	//double v1 = std::pow(13.0/12.0, m_double) - 2.0*std::pow(35.0/32.0, m_double);
	//double v2 = (std::pow(1.25, m_double)/(n_double*n_double));
	//double v3 = n_double + 2.0*HAMMING;
	//return (v1 + v2*v3);
	return HAMMING;
}
*/