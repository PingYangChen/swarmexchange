
// DECLARE FUNCTIONS
void UPDATE_RANDOM(mat &DESIGN, double &DESIGN_VAL,
									 const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS);

// BODY
void UPDATE_RANDOM(mat &DESIGN, double &DESIGN_VAL, 
									 const int &maximize, const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS)
{
	arma::imat labRand = randi(DESIGN.n_rows, DESIGN.n_cols, distr_param(0, D_INFO.nLevel - 1));
	
	for (uword i = 0; i < DESIGN.n_rows; i++) {
		for (uword j = 0; j < DESIGN.n_cols; j++) {
      DESIGN(i, j) = D_INFO.labLevel(labRand(i, j));
		}
	}
	DESIGNCRITERION(DESIGN_VAL, DESIGN, D_INFO, -1);
}

