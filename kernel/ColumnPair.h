
// BODY
// CP MAIN FUNCTIONS
void ColumnPair_MAIN(Ptr_CP_Result Ptr_CP_Result, const CP_OPTIONS &CP_OPTS, const DESIGN_INFO &D_INFO, bool VERBOSE)
{
  // SET CP PARAMETERS
  int maxIter		= CP_OPTS.maxIter;
  int nTry			= CP_OPTS.nTry;
  int maximize 	= CP_OPTS.maximize;
  int CPk				= CP_OPTS.CPk;
  double tol 		= CP_OPTS.tol;
	
	int nRun = D_INFO.nRun;
	int nFactor = D_INFO.nFactor;
	//int nLevel = (int)D_INFO.labLevel.n_elem;
	int nModel = D_INFO.nModel;
	
	if (CPk > nFactor) {
		CPk = nFactor;
		Rprintf("CPk cannot exceed nFactor. Use CPk = nFactor instead."); 
	}

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
			arma::mat DESIGN_VAL_MAT = DESIGNCRITERION(DESIGN_VAL, DESIGN, D_INFO, -1);
			fvalHist(0) = DESIGN_VAL;
			//Rprintf("%d\n", iTry);
			int t; // SET ITERATION COUNTER
	
			arma::vec valVec(DESIGN.n_cols);
			arma::vec countVec(DESIGN.n_cols);
			arma::uvec ADJUST_ORDER(DESIGN.n_cols);
			arma::mat ADJ_DESIGN = DESIGN; double ADJ_VAL;
			uword ADJUST_WHO = 0;

			////// COLUMNWISE-PAIRWISE ALGORITHM LOOP
			for (t = 0; t < maxIter; t++) {
				//Rprintf("%d, %d\n", iTry, t);
				/*if (VERBOSE) {
					if (t == 0) { Rprintf("Updating ..    "); }
					Rprintf("\b\b\b%2.0f%%", (double)((t+1)*100/maxIter)); 
					if (t == (maxIter - 1)) { Rprintf("\n"); }
				}*/
				
				// Deletion
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
				//if (maximize == 0) valVec.min(ADJUST_WHO); else valVec.max(ADJUST_WHO);
				if (maximize == 0) ADJUST_ORDER = arma::sort_index(valVec, "ascend");	else ADJUST_ORDER = arma::sort_index(valVec, "descend");

				// Pairwise-Exchange
				for (int k = 0; k < CPk; k++) {
					
					ADJUST_WHO = ADJUST_ORDER(k);
					DESIGN = ColumnPair_CORE(DESIGN, DESIGN_VAL, D_INFO, ADJUST_WHO, maximize);
					/*
					arma::uvec lab_one = arma::find(DESIGN.col(ADJUST_WHO) == D_INFO.labLevel(0));
					arma::uvec lab_two = arma::find(DESIGN.col(ADJUST_WHO) == D_INFO.labLevel(1));
					
					ADJ_DESIGN = DESIGN;
					arma::umat PAIR_CAND(lab_one.n_elem*lab_two.n_elem, 2);
					arma::vec ADJ_VAL_VEC(lab_one.n_elem*lab_two.n_elem);
					double count = 0;
					for (uword i = 0; i < lab_one.n_elem; i++) {
						for (uword j = 0; j < lab_two.n_elem; j++) {
							ADJ_DESIGN = DESIGN;
							ADJ_DESIGN(lab_one(i), ADJUST_WHO) = D_INFO.labLevel(1);
							ADJ_DESIGN(lab_two(j), ADJUST_WHO) = D_INFO.labLevel(0);
							DESIGNCRITERION(ADJ_VAL, ADJ_DESIGN, D_INFO, -1);
							PAIR_CAND(count, 0) = lab_one(i);
							PAIR_CAND(count, 1) = lab_two(j);
							ADJ_VAL_VEC(count) = ADJ_VAL;
							count++;
						}
					}
					uword SELECT_PAIR;
					if (maximize == 0) ADJ_VAL = ADJ_VAL_VEC.min(SELECT_PAIR); else ADJ_VAL = ADJ_VAL_VEC.max(SELECT_PAIR);
					if (((maximize == 0) & (ADJ_VAL <= DESIGN_VAL)) | ((maximize != 0) & (ADJ_VAL >= DESIGN_VAL))) {
						DESIGN(PAIR_CAND(SELECT_PAIR, 0), ADJUST_WHO) = D_INFO.labLevel(1);
						DESIGN(PAIR_CAND(SELECT_PAIR, 1), ADJUST_WHO) = D_INFO.labLevel(0);
						arma::mat DESIGN_VAL_MAT = DESIGNCRITERION(DESIGN_VAL, DESIGN, D_INFO, -1);		
					}
					*/
				}
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
  Ptr_CP_Result->DESIGN = BestD;
  Ptr_CP_Result->DESIGN_VAL = BestVal;
  Ptr_CP_Result->fvalHist = BestPath;
  Ptr_CP_Result->DESIGN_POOL = DESIGN_POOL;
  Ptr_CP_Result->DESIGN_VAL_POOL = DESIGN_VAL_POOL;
  Ptr_CP_Result->fvalHist_POOL = fvalHist_POOL;
}

