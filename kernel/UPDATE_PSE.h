
void UPDATE_PSE(mat &DESIGN, double &DESIGN_VAL, const mat &GBest, const mat &PBest, 
								const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS, const PSO_DYN &PSO_DYN);

//
void UPDATE_PSE(mat &DESIGN, double &DESIGN_VAL, const mat &GBest, const mat &PBest, 
								const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS, const PSO_DYN &PSO_DYN)
{
	double JFO_R_CUR = PSO_DYN.JFO_R_CUR;
	double JFO_PL = PSO_OPTS.JFO_RHO;
	//double JFO_PG = 1.0 - JFO_PL;
	
	arma::rowvec JUMP(2);
	JUMP(0) = JFO_R_CUR;
	JUMP(1) = JUMP(0) + (1.0 - JFO_R_CUR)*JFO_PL;
	//JUMP(2) = JUMP(1) + (1.0 - JFO_R_CUR)*JFO_PG;
	
	double DICE = arma::as_scalar(randu(1, 1));
	
	arma::mat MIX_COL(DESIGN.n_rows, DESIGN.n_cols);
	arma::mat MIX_ROW(DESIGN.n_rows, DESIGN.n_cols);
	double MIX_COL_VAL, MIX_ROW_VAL;
	
	arma::mat MAT_TAR(DESIGN.n_rows, DESIGN.n_cols);
	
	if (DICE < JUMP(1)) {
		if (DICE < JUMP(0)) {
			// mix w Random Design
			initDesigns(MAT_TAR, D_INFO);
		} else {
			// mix w LB
			MAT_TAR = PBest;
		}
	} else {
		// mix w GB
		MAT_TAR = GBest;
	}
	//
	MIX_COL_VAL = MIX_OPERATOR(0, MIX_COL, DESIGN, DESIGN_VAL, MAT_TAR, maximize, D_INFO, PSO_OPTS);
	if (D_INFO.balance != 1) {
		MIX_ROW_VAL = MIX_OPERATOR(1, MIX_ROW, DESIGN, DESIGN_VAL, MAT_TAR, maximize, D_INFO, PSO_OPTS);
	} else {
		MIX_ROW = DESIGN; MIX_ROW_VAL = DESIGN_VAL;
	}
	//
	arma::rowvec moveCandidate(2, fill::zeros); uword MOVE;
	moveCandidate << MIX_COL_VAL << MIX_ROW_VAL << endr;
	if (maximize == 0) moveCandidate.min(MOVE); else moveCandidate.max(MOVE);
	switch (MOVE) {
		case 0: { DESIGN = MIX_COL; DESIGN_VAL = MIX_COL_VAL; break; }
		case 1: { DESIGN = MIX_ROW; DESIGN_VAL = MIX_ROW_VAL; break; }
	}

}