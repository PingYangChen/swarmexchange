
void initDesigns(arma::mat &DESIGN_POOL, const DESIGN_INFO &D_INFO)
{
	int nRun = D_INFO.nRun;
	int nFactor = D_INFO.nFactor;
  // INITIALIZE RANDOM DESIGN
  int ND = (int)(DESIGN_POOL.n_rows/nRun);

  if (D_INFO.balance == 1) {
  	arma::vec balanceCol_0 = arma::repmat(D_INFO.labLevel.t(), (int)(nRun/D_INFO.labLevel.n_elem), 1);
  	for (int i = 0; i < ND; i++) {
  		arma::mat DESIGN(nRun, nFactor);
			for (int j = 0; j < nFactor; j++) {
				arma::vec randNum = arma::randu<vec>(nRun);
				arma::uvec randId = arma::sort_index(randNum);
				for (int k = 0; k < nRun; k++) {
					DESIGN(k, j) = balanceCol_0(randId(k));	
				}
			}
			DESIGN_POOL.rows(i*nRun, (i + 1)*nRun - 1) = DESIGN;
		}	
  } else {
	  arma::imat labRand = randi(ND*nRun, nFactor, distr_param(0, D_INFO.nLevel - 1));
		for (int i = 0; i < (ND*nRun); i++) {
			for (int j = 0; j < nFactor; j++) {
				DESIGN_POOL(i, j) = D_INFO.labLevel(labRand(i,j));
			}
		}	
  }
}

