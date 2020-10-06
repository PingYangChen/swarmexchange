// DECLARE FUNCTIONS
arma::mat CoorExchange_CORE(arma::mat DESIGN, double &DESIGN_VAL, const DESIGN_INFO &D_INFO, int row_i, int col_j, int maximize);
arma::mat ColumnPair_CORE(arma::mat DESIGN, double &DESIGN_VAL, const DESIGN_INFO &D_INFO, int col_j, int maximize);

// BODY
arma::mat CoorExchange_CORE(arma::mat DESIGN, double &DESIGN_VAL, const DESIGN_INFO &D_INFO, int row_i, int col_j, int maximize)
{
  int nLevel = (int)D_INFO.labLevel.n_elem;
  //
	arma::vec DESIGN_VAL_COMP(nLevel, fill::zeros);
	double DESIGN_VAL_TMP;
	arma::mat DESIGN_TMP = DESIGN;
	arma::uword EXCHANGE;
	//
	int level_k;
	for (level_k = 0; level_k < nLevel; level_k++) {
		DESIGN_TMP = DESIGN;
		if (DESIGN_TMP(row_i, col_j) == D_INFO.labLevel(level_k)) {
			DESIGN_VAL_COMP(level_k) = DESIGN_VAL;
		} else {
			DESIGN_TMP(row_i, col_j) = D_INFO.labLevel(level_k);
			DESIGNCRITERION(DESIGN_VAL_TMP, DESIGN_TMP, D_INFO, -1);
			DESIGN_VAL_COMP(level_k) = DESIGN_VAL_TMP;
		}
	}
	if (maximize == 0) DESIGN_VAL = DESIGN_VAL_COMP.min(EXCHANGE); else DESIGN_VAL = DESIGN_VAL_COMP.max(EXCHANGE);
	DESIGN(row_i, col_j) = D_INFO.labLevel(EXCHANGE);
	return DESIGN;
}

arma::mat ColumnPair_CORE(arma::mat DESIGN, double &DESIGN_VAL, const DESIGN_INFO &D_INFO, int col_j, int maximize)
{
	arma::uvec lab_one = arma::find(DESIGN.col(col_j) == D_INFO.labLevel(0));
	arma::uvec lab_two = arma::find(DESIGN.col(col_j) == D_INFO.labLevel(1));
	//
	arma::mat ADJ_DESIGN = DESIGN;
	double ADJ_VAL;
	arma::umat PAIR_CAND(lab_one.n_elem*lab_two.n_elem, 2);
	arma::vec ADJ_VAL_VEC(lab_one.n_elem*lab_two.n_elem);
	//
	double count = 0;
	for (uword i = 0; i < lab_one.n_elem; i++) {
		for (uword j = 0; j < lab_two.n_elem; j++) {
			ADJ_DESIGN = DESIGN;
			ADJ_DESIGN(lab_one(i), col_j) = D_INFO.labLevel(1);
			ADJ_DESIGN(lab_two(j), col_j) = D_INFO.labLevel(0);
			DESIGNCRITERION(ADJ_VAL, ADJ_DESIGN, D_INFO, -1);
			PAIR_CAND(count, 0) = lab_one(i);
			PAIR_CAND(count, 1) = lab_two(j);
			ADJ_VAL_VEC(count) = ADJ_VAL;
			count++;
		}
	}
	arma::uword SELECT_PAIR;
	if (maximize == 0) ADJ_VAL = ADJ_VAL_VEC.min(SELECT_PAIR); else ADJ_VAL = ADJ_VAL_VEC.max(SELECT_PAIR);
	if (((maximize == 0) & (ADJ_VAL <= DESIGN_VAL)) | ((maximize != 0) & (ADJ_VAL >= DESIGN_VAL))) {
		DESIGN(PAIR_CAND(SELECT_PAIR, 0), col_j) = D_INFO.labLevel(1);
		DESIGN(PAIR_CAND(SELECT_PAIR, 1), col_j) = D_INFO.labLevel(0);
		DESIGNCRITERION(DESIGN_VAL, DESIGN, D_INFO, -1);
	}
	return DESIGN;
}