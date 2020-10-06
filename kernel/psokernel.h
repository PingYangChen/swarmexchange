
// BODY
// PSO MAIN FUNCTIONS
void PSO_MAIN(Ptr_PSO_Result Ptr_PSO_Result, const PSO_OPTIONS &PSO_OPTS, const DESIGN_INFO &D_INFO, bool VERBOSE)
{
  // SET PSO PARAMETERS
  int nSwarm = PSO_OPTS.nSwarm;
  int maxIter		= PSO_OPTS.maxIter;
  int maximize 	= PSO_OPTS.maximize;
  double tol = PSO_OPTS.tol;
	//int PSO_UPDATE = PSO_OPTS.PSO_UPDATE;
	
	int nRun = D_INFO.nRun;
	int nFactor = D_INFO.nFactor;
	
  arma::mat swarm(nSwarm*nRun, nFactor), PBest(nSwarm*nRun, nFactor);
  arma::mat GBest(nRun, nFactor);
  arma::vec fSwarm(nSwarm), fPBest(nSwarm);
  double fGBest;
  arma::uword GBestIdx;
  arma::rowvec fGBestHist(maxIter + 1); fGBestHist.zeros();
  arma::mat fPBestHist(nSwarm, maxIter + 1); fPBestHist.zeros();
	PSO_DYN PSO_DYN;
	
	//arma::imat updateRec(nSwarm, maxIter);
  ////// INITIALIZE
  if (VERBOSE) Rprintf("Initializing .. "); 
  // INITIALIZE RANDOM SWARM
  initDesigns(swarm, D_INFO);
  /*
	arma::imat labRand = randi(nSwarm*nRun, nFactor, distr_param(0, D_INFO.nLevel - 1));
	for (int i = 0; i < (nSwarm*nRun); i++) {
		for (int j = 0; j < nFactor; j++) {
			swarm(i,j) = D_INFO.labLevel(labRand(i,j));
		}
	}
	*/
	psoUpdateDynPara(PSO_OPTS, -1, PSO_DYN);
	
  int iSwarm, t; // SET PARTICLE COUNTER AND ITERATION COUNTER
  // --- Initialization of function evaluation
	#pragma omp parallel private(iSwarm) 
  {
		#pragma omp for
		for (iSwarm = 0; iSwarm < nSwarm; iSwarm++) {
			//Rprintf("Particle %d\n", iSwarm);
			arma::mat DESIGN = swarm.rows(iSwarm*nRun, (iSwarm + 1)*nRun - 1);
			double DESIGN_VAL;
			DESIGNCRITERION(DESIGN_VAL, DESIGN, D_INFO, -1);
			fSwarm(iSwarm) = DESIGN_VAL;
		}
	} 
  // Initialize LOCAL BEST
  fPBest = fSwarm;	PBest = swarm;
  // Initialize GLOBAL BEST
  if (maximize == 0) fGBest = fPBest.min(GBestIdx); else fGBest = fPBest.max(GBestIdx);	
	GBest = PBest.rows(GBestIdx*nRun, (GBestIdx + 1)*nRun - 1);	
	fGBestHist(0) = fGBest;
  fPBestHist.col(0) = fPBest;
	
	if (VERBOSE) Rprintf("DONE \n"); 
	
  ////// PSO LOOP
  for (t = 0; t < maxIter; t++) {
		if (VERBOSE) {
			if (t == 0) { Rprintf("Updating ..    "); }
			Rprintf("\b\b\b%2.0f%%", (double)((t+1)*100/maxIter)); 
			if (t == (maxIter - 1)) { Rprintf("\n"); }
		}
		// UPDATE PARTICLES
    psoUpdateParticle(D_INFO, PSO_OPTS, PSO_DYN, PBest, GBest, fPBest, fGBest, swarm, fSwarm);
    // UPDATING THE LOCAL BEST
    // UPDATING THE GLOBAL BEST
    if (maximize == 0) {
      if (any(fSwarm < fPBest)) {
        uvec udPB = find(fSwarm < fPBest);
				fPBest.elem(udPB) = fSwarm.elem(udPB);
				for (uword iUD = 0; iUD < udPB.n_elem; iUD++) {
					PBest.rows(udPB(iUD)*nRun, (udPB(iUD) + 1)*nRun - 1) = swarm.rows(udPB(iUD)*nRun, (udPB(iUD) + 1)*nRun - 1);
				}
      }
      if (min(fPBest) < fGBest) {
        fGBest = fPBest.min(GBestIdx);	
        GBest = PBest.rows(GBestIdx*nRun, (GBestIdx + 1)*nRun - 1);	
      }		
    } else {
      if (any(fSwarm > fPBest)) {
        uvec udPB = find(fSwarm > fPBest);
        fPBest.elem(udPB) = fSwarm.elem(udPB);
				for (uword iUD = 0; iUD < udPB.n_elem; iUD++) {
					PBest.rows(udPB(iUD)*nRun, (udPB(iUD) + 1)*nRun - 1) = swarm.rows(udPB(iUD)*nRun, (udPB(iUD) + 1)*nRun - 1);
				}
      }
      if (max(fPBest) > fGBest) {
        fGBest = fPBest.max(GBestIdx);	
        GBest = PBest.rows(GBestIdx*nRun, (GBestIdx + 1)*nRun - 1);	
      }
    }
		// UPDATE PSO DYNAMIC PARAMETERS
    psoUpdateDynPara(PSO_OPTS, t, PSO_DYN);
    // RECORDING THE CURRENT GLOBAL BEST VALUE
    fGBestHist(t+1) = fGBest; 
    fPBestHist.col(t+1) = fPBest;
		//if ((t*2 >= maxIter) && (tol > 0)) {
		if (tol > 0) {	
      //if (std::abs(fGBest - fGBestHist(t)) < tol) t = maxIter;
      if ((arma::stddev(fPBest)/arma::mean(fPBest)) < tol) t = maxIter;
    }
  }
  // OUTPUT
  Ptr_PSO_Result->GBest = GBest;
  Ptr_PSO_Result->fGBest = fGBest;
  Ptr_PSO_Result->fGBestHist = fGBestHist;
  Ptr_PSO_Result->PBest = PBest;
  Ptr_PSO_Result->fPBest = fPBest;
  Ptr_PSO_Result->fPBestHist = fPBestHist;
	//Ptr_PSO_Result->updateRec = updateRec;
}

