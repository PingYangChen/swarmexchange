
// BODY
// CE MAIN FUNCTIONS
void CoorExchange_MAIN(Ptr_CE_Result Ptr_CE_Result, const CE_OPTIONS &CE_OPTS, const DESIGN_INFO &D_INFO, bool VERBOSE)
{
  // SET CE PARAMETERS
  int maxIter		= CE_OPTS.maxIter;
  int nTry			= CE_OPTS.nTry;
  int maximize 	= CE_OPTS.maximize;
  double tol 		= CE_OPTS.tol;
	
  int nRun = D_INFO.nRun;
  int nFactor = D_INFO.nFactor;
  int nLevel = (int)D_INFO.labLevel.n_elem;
	
  arma::mat DESIGN_POOL(nTry*nRun, nFactor);
  arma::vec DESIGN_VAL_POOL(nTry);
  
  arma::mat fvalHist_POOL(maxIter + 1, nTry); fvalHist_POOL.zeros();
	
	//arma::imat updateRec(nSwarm, maxIter);
  ////// INITIALIZE
  if (VERBOSE) Rprintf("Initializing .. "); 
  // INITIALIZE RANDOM DESIGN
  initDesigns(DESIGN_POOL, D_INFO);
  /*
	arma::imat labRand = randi(nTry*nRun, nFactor, distr_param(0, D_INFO.nLevel - 1));
	for (int i = 0; i < (nTry*nRun); i++) {
		for (int j = 0; j < nFactor; j++) {
			DESIGN_POOL(i,j) = D_INFO.labLevel(labRand(i,j));
		}
	}
	*/
	int iTry;
	if (VERBOSE) Rprintf("DONE \n"); 
	
	#pragma omp parallel private(iTry) 
  {
		#pragma omp for
		for (iTry = 0; iTry < nTry; iTry++) {
			//Rprintf("%d\n", iTry);
			arma::mat DESIGN = DESIGN_POOL.rows(iTry*nRun, (iTry + 1)*nRun - 1);
			double DESIGN_VAL;
			arma::vec fvalHist(maxIter + 1); fvalHist.zeros();
			//Rprintf("%d\n", iTry);
			DESIGNCRITERION(DESIGN_VAL, DESIGN, D_INFO, -1);
			fvalHist(0) = DESIGN_VAL;
			//Rprintf("%d\n", iTry);
			int t; // SET ITERATION COUNTER
			arma::vec DESIGN_VAL_COMP(nLevel, fill::zeros);
			double DESIGN_VAL_TMP;
			arma::mat DESIGN_TMP = DESIGN;
			arma::uword EXCHANGE;
	
			int row_i, col_j, level_k;
			
			////// COORDINATE EXCHANGING LOOP
			for (t = 0; t < maxIter; t++) {
				//Rprintf("%d, %d\n", iTry, t);
				/*if (VERBOSE) {
					if (t == 0) { Rprintf("Updating ..    "); }
					Rprintf("\b\b\b%2.0f%%", (double)((t+1)*100/maxIter)); 
					if (t == (maxIter - 1)) { Rprintf("\n"); }
				}*/
				for (row_i = 0; row_i < nRun; row_i++) {
					for (col_j = 0; col_j < nFactor; col_j++) {
						/*
						DESIGN_VAL_COMP.zeros();
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
						*/
						DESIGN = CoorExchange_CORE(DESIGN, DESIGN_VAL, D_INFO, row_i, col_j, maximize);
					}
				}
				//Rprintf("%d, %d\n", iTry, t);
				// RECORDING THE CURRENT VALUE
				fvalHist(t+1) = DESIGN_VAL; 
				//if ((t*2 >= maxIter) && (tol > 0)) {
				if (tol > 0) {	
					if (std::abs(DESIGN_VAL - fvalHist(t)) < tol) t = maxIter;
				}
				//Rprintf("%d, %d\n", iTry, t);
			}
			// Save Result
			DESIGN_POOL.rows(iTry*nRun, (iTry + 1)*nRun - 1) = DESIGN;
			DESIGN_VAL_POOL(iTry) = DESIGN_VAL;
			fvalHist_POOL.col(iTry) = fvalHist; 
		}
	} 
	
	arma::uword BestID; double BestVal;
	if (maximize == 0) {
		BestVal = DESIGN_VAL_POOL.min(BestID); 
	} else { 
		BestVal = DESIGN_VAL_POOL.max(BestID); 
	}
	int BestID_int = (int)BestID;
	arma::mat BestD = DESIGN_POOL.rows(BestID_int*nRun, (BestID_int + 1)*nRun - 1);
	arma::vec BestPath = fvalHist_POOL.col(BestID);
	
  // OUTPUT
  Ptr_CE_Result->DESIGN = BestD;
  Ptr_CE_Result->DESIGN_VAL = BestVal;
  Ptr_CE_Result->fvalHist = BestPath;

}

