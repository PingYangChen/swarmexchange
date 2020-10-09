void psoUpdateDynPara(const double &fGBest, const arma::vec &fPBest, const arma::rowvec &fGBestHist, 
											const PSO_OPTIONS &PSO_OPTS, const int iter, PSO_DYN &PSO_DYN)
{
  if (iter < 0) { // INITIALIZE
	
  	int JFO_R_DUR = (int)(PSO_OPTS.JFO_RV*PSO_OPTS.maxIter);
  	PSO_DYN.JFO_R_DUR	= JFO_R_DUR;
	  PSO_DYN.JFO_R_CUR	= PSO_OPTS.JFO_R0;
	  PSO_DYN.JFO_R_DEC	= (PSO_OPTS.JFO_R0 - PSO_OPTS.JFO_R1)/JFO_R_DUR;      // Inertia weight change per iteration step
	  // TEMP_MAX = 1000.0
	  PSO_DYN.EXALG_PROB = 0.0;
	  arma::mat EXALG_TRIGGER(PSO_OPTS.maxIter + 1, 3, fill::zeros);
	  arma::rowvec extg_tmp(3, fill::zeros);
	  extg_tmp << 0.0 << 0.0 << arma::stddev(fGBest - fPBest) << endr; // prob, temp, tv, sv
	  EXALG_TRIGGER.row(0) = extg_tmp;
	  PSO_DYN.EXALG_TRIGGER = EXALG_TRIGGER;

  } else { // UPDATE
    
		if (iter <= PSO_DYN.JFO_R_DUR) {
			PSO_DYN.JFO_R_CUR = PSO_DYN.JFO_R_CUR - PSO_DYN.JFO_R_DEC; 
		}
		double fGB_timeVar = std::abs(fGBest - fGBestHist(iter))/std::abs(fGBestHist(iter)); // tv^p_
		double fGB_spaceVar = arma::stddev(fPBest)/std::abs(arma::mean(fPBest)); // sv^p_
		if ((fGB_timeVar + fGB_spaceVar) >= 2.) {
			PSO_DYN.EXALG_PROB = 0.;
		} else {
			PSO_DYN.EXALG_PROB = 1. - .5*(fGB_timeVar + fGB_spaceVar);//std::exp(-(fGB_spaceVar + fGB_timeVar));
		}
		//Rprintf("EXALG_PROB = %2.1f%% at ITERATION %d (tv = %2.2f, sv = %2.2f)\n", PSO_DYN.EXALG_PROB*100.0, iter, fGB_timeVar, fGB_spaceVar); 
		arma::rowvec extg_tmp(3, fill::zeros);
		extg_tmp << PSO_DYN.EXALG_PROB << fGB_timeVar << fGB_spaceVar << endr; // prob, temp, tv, sv
		PSO_DYN.EXALG_TRIGGER.row(iter+1) = extg_tmp;

  }
}
