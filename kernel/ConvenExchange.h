
void ConvenExchange(mat &DESIGN, double &DESIGN_VAL, const mat &GBest, const mat &PBest, 
										const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS);

//
void ConvenExchange(mat &DESIGN, double &DESIGN_VAL, const mat &GBest, const mat &PBest, 
										const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS)
{
	arma::mat MIX_COL_GB(DESIGN.n_rows, DESIGN.n_cols);
	arma::mat MIX_ROW_GB(DESIGN.n_rows, DESIGN.n_cols);
	
	double MIX_COL_GB_VAL = MIX_OPERATOR(0, MIX_COL_GB, DESIGN, GBest, maximize, D_INFO, PSO_OPTS);
	double MIX_ROW_GB_VAL = MIX_OPERATOR(1, MIX_ROW_GB, DESIGN, GBest, maximize, D_INFO, PSO_OPTS);

	arma::mat MIX_COL_LB(DESIGN.n_rows, DESIGN.n_cols);
	arma::mat MIX_ROW_LB(DESIGN.n_rows, DESIGN.n_cols);
	
	double MIX_COL_LB_VAL = MIX_OPERATOR(0, MIX_COL_LB, DESIGN, PBest, maximize, D_INFO, PSO_OPTS);
	double MIX_ROW_LB_VAL = MIX_OPERATOR(1, MIX_ROW_LB, DESIGN, PBest, maximize, D_INFO, PSO_OPTS);

	rowvec valSet(5);
	valSet << DESIGN_VAL << MIX_COL_GB_VAL << MIX_ROW_GB_VAL << MIX_COL_LB_VAL << MIX_ROW_LB_VAL << endr;

	// MOVE OPERATOR
	uword MOVE; double newVal;
	if (maximize == 0) newVal = valSet.min(MOVE); else newVal = valSet.max(MOVE);
	switch (MOVE) {
		case 0:
		{
			arma::mat MIX_RAND(DESIGN.n_rows, DESIGN.n_cols);
			double MIX_RAND_VAL = MIX_OPERATOR(2, MIX_RAND, DESIGN, GBest, maximize, D_INFO, PSO_OPTS);
			DESIGN = MIX_RAND; DESIGN_VAL = MIX_RAND_VAL;
		break;
		}
		case 1:
		{
			DESIGN = MIX_COL_GB; DESIGN_VAL = newVal;
		break;
		}
		case 2:
		{
			DESIGN = MIX_ROW_GB; DESIGN_VAL = newVal;
		break;
		}
		case 3:
		{
			DESIGN = MIX_COL_LB; DESIGN_VAL = newVal;
		break;
		}
		case 4:
		{
			DESIGN = MIX_ROW_LB; DESIGN_VAL = newVal;
		break;
		}
	}
}