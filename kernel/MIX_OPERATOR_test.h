
// DECLARE FUNCTIONS
double MIX_OPERATOR(const int &MIX_PROC, mat &MIX, const mat &DESIGN_A, const double &VAL_A, const mat &DESIGN_B, const int &maximize, 
										const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS);
//double HAMMING_TWOLEVEL(arma::mat &DESIGN);

// BODY
double MIX_OPERATOR(const int &MIX_PROC, mat &MIX, const mat &DESIGN_A, const double &VAL_A, const mat &DESIGN_B, 
										const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS)
{
	double MIX_VAL = VAL_A; double DESIGN_VAL;
	MIX = DESIGN_A;
	uword DELETE_WHO = 0; arma::mat ADD_ONE_BACK = MIX;
	
	switch (MIX_PROC) {
		case 0: 
		{	// Column Exchange 
			int MIX_C = PSO_OPTS.MIX_C;
			arma::mat TMP = MIX;
			arma::vec valVec(MIX.n_cols, fill::zeros);
			arma::vec countVec(MIX.n_cols, fill::zeros);
			arma::uvec DELETE_ORDER(MIX.n_cols);
			int nModel = D_INFO.nModel;

			// Screening
			arma::mat DESIGN_VAL_MAT = DESIGNCRITERION(MIX_VAL, MIX, D_INFO, -1);
			valVec.fill(arma::accu(DESIGN_VAL_MAT));
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
			valVec = valVec/countVec;
			if (maximize == 0) DELETE_ORDER = arma::sort_index(valVec, "ascend");	else DELETE_ORDER = arma::sort_index(valVec, "descend");
			// Addition
			for (int q = 0; q < MIX_C; q++) {
				DELETE_WHO = DELETE_ORDER(q);
				for (uword i = 0; i < DESIGN_B.n_cols; i++) {
					ADD_ONE_BACK = MIX;
					if (!arma::all(ADD_ONE_BACK.col(DELETE_WHO) == DESIGN_B.col(i))) {
						ADD_ONE_BACK.col(DELETE_WHO) = DESIGN_B.col(i);
						DESIGNCRITERION(DESIGN_VAL, ADD_ONE_BACK, D_INFO, -1);
						if (maximize == 0) {
							if (DESIGN_VAL <= MIX_VAL) { MIX = ADD_ONE_BACK; MIX_VAL = DESIGN_VAL; i += DESIGN_B.n_cols; }
						} else {
							if (DESIGN_VAL >= MIX_VAL) { MIX = ADD_ONE_BACK; MIX_VAL = DESIGN_VAL; i += DESIGN_B.n_cols; }
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
			arma::uvec DELETE_ORDER(MIX.n_rows);
			arma::mat DELETE_ONE_ROW(MIX.n_rows - 1, MIX.n_cols, fill::zeros);
			arma::uvec RID = arma::linspace<uvec>(0, MIX.n_rows - 1, MIX.n_rows);	
			// Screening
			for (uword i = 0; i < MIX.n_rows; i++) {
				DELETE_ONE_ROW = MIX.rows(arma::find(RID != i));
				DESIGNCRITERION(DESIGN_VAL, DELETE_ONE_ROW, D_INFO, -1);
				valVec(i) = DESIGN_VAL;
			}
			if (maximize == 0) DELETE_ORDER = arma::sort_index(valVec, "ascend");	else DELETE_ORDER = arma::sort_index(valVec, "descend");
			// Addition
			for (int q = 0; q < MIX_R; q++) {
				DELETE_WHO = DELETE_ORDER(q);
				for (uword i = 0; i < DESIGN_B.n_rows; i++) {
					ADD_ONE_BACK = MIX;
					if (!arma::all(ADD_ONE_BACK.row(DELETE_WHO) == DESIGN_B.row(i))) {
						ADD_ONE_BACK.row(DELETE_WHO) = DESIGN_B.row(i);
						DESIGNCRITERION(DESIGN_VAL, ADD_ONE_BACK, D_INFO, -1);
						if (maximize == 0) {
							if (DESIGN_VAL <= MIX_VAL) { MIX = ADD_ONE_BACK; MIX_VAL = DESIGN_VAL; i += DESIGN_B.n_rows; }
						} else {
							if (DESIGN_VAL >= MIX_VAL) { MIX = ADD_ONE_BACK; MIX_VAL = DESIGN_VAL; i += DESIGN_B.n_rows; }
						}
					}
				}
			}
			break;
		}
	}
	return(MIX_VAL);
}

