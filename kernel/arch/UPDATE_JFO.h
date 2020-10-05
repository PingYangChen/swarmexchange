
// DECLARE FUNCTIONS
void UPDATE_JFO(mat &DESIGN, double &DESIGN_VAL, const mat &GBest, const mat &PBest, 
								const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS, const PSO_DYN &PSO_DYN);

// BODY
void UPDATE_JFO(mat &DESIGN, double &DESIGN_VAL, const mat &GBest, const mat &PBest, 
								const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS, const PSO_DYN &PSO_DYN)
{
	arma::mat tmpDESIGN = DESIGN;
	
	double JFO_R_CUR = PSO_DYN.JFO_R_CUR;
	double JFO_PL = PSO_OPTS.JFO_RHO;
	//double JFO_PG = 1.0 - JFO_PL;
	
	arma::rowvec JUMP(2);
	JUMP(0) = JFO_R_CUR;
	JUMP(1) = JUMP(0) + (1.0 - JFO_R_CUR)*JFO_PL;
	//JUMP(2) = JUMP(1) + (1.0 - JFO_R_CUR)*JFO_PG;
	
	arma::mat DICE = randu(tmpDESIGN.n_rows, tmpDESIGN.n_cols);
	arma::imat labRand = randi(tmpDESIGN.n_rows, tmpDESIGN.n_cols, distr_param(0, D_INFO.nLevel - 1));
	
	for (uword i = 0; i < tmpDESIGN.n_rows; i++) {
		for (uword j = 0; j < tmpDESIGN.n_cols; j++) {
			if (DICE(i, j) < JUMP(1)) {
				if (DICE(i, j) < JUMP(0)) {
					// mix w Random Design
					tmpDESIGN(i, j) = D_INFO.labLevel(labRand(i, j));
				} else {
					// mix w LB
					tmpDESIGN(i, j) = PBest(i, j);
				}
			} else {
				// mix w GB
				tmpDESIGN(i, j) = GBest(i, j);
			}			
		}
	}
	
	DESIGNCRITERION(DESIGN_VAL, tmpDESIGN, D_INFO, -1);
	DESIGN = tmpDESIGN;
}

