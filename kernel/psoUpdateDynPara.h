void psoUpdateDynPara(const PSO_OPTIONS &PSO_OPTS, const int iter, PSO_DYN &PSO_DYN)
{
  if (iter < 0) { // INITIALIZE
	
  	int JFO_R_DUR = (int)(PSO_OPTS.JFO_RV*PSO_OPTS.maxIter);
  	PSO_DYN.JFO_R_DUR	= JFO_R_DUR;
	  PSO_DYN.JFO_R_CUR	= PSO_OPTS.JFO_R0;
	  PSO_DYN.JFO_R_DEC	= (PSO_OPTS.JFO_R0 - PSO_OPTS.JFO_R1)/JFO_R_DUR;      // Inertia weight change per iteration step
	  
  } else { // UPDATE
    
		if (iter <= PSO_DYN.JFO_R_DUR) PSO_DYN.JFO_R_CUR = PSO_DYN.JFO_R_CUR - PSO_DYN.JFO_R_DEC; 
		
  }
}
