
#include "psoheader.h"

// RCPP FUNCTIONS
//[[Rcpp::export]]
List DiscreteDesignPSO(Rcpp::List ALG_INFO_LIST, Rcpp::List D_INFO_LIST, bool PARALLEL, bool VERBOSE)
{
  //arma_rng::set_seed_random();
  int nc = omp_get_num_procs();
	if (nc < 3) nc = 2;
	if (PARALLEL) omp_set_num_threads(nc - 1); else omp_set_num_threads(1);
  
	DESIGN_INFO D_INFO = {};
	getInfoStruct(&D_INFO, D_INFO_LIST);

	PSO_OPTIONS PSO_OPT = {};
	getAlgStruct(&PSO_OPT, ALG_INFO_LIST);

  PSO_Result Result = {};
	if (VERBOSE) Rprintf("\nCalling Cpp PSO Kernel... ");
  PSO_MAIN(&Result, PSO_OPT, D_INFO, VERBOSE);
  
  return List::create(Named("GBest") = wrap(Result.GBest),
                      Named("fGBest") = wrap(Result.fGBest),
                      Named("fGBestHist") = wrap(Result.fGBestHist),
                      Named("PBest") = wrap(Result.PBest),
                      Named("fPBest") = wrap(Result.fPBest),
                      Named("fPBestHist") = wrap(Result.fPBestHist),
                      Named("exAlgTrig") = wrap(Result.exAlgTrig));//,
											//Named("updateRec") = wrap(Result.updateRec));
}

//[[Rcpp::export]]
List DiscreteDesignCoorEx(Rcpp::List CE_INFO_LIST, Rcpp::List D_INFO_LIST, bool PARALLEL, bool VERBOSE)
{
	int nc = omp_get_num_procs();
	if (nc < 3) nc = 2;
	if (PARALLEL) omp_set_num_threads(nc - 1); else omp_set_num_threads(1);
	
	DESIGN_INFO D_INFO = {};	getInfoStruct(&D_INFO, D_INFO_LIST);
	CE_OPTIONS CE_OPT = {};		getCoorExStruct(&CE_OPT, CE_INFO_LIST);
  CE_Result Result = {};
	if (VERBOSE) Rprintf("\nCalling Cpp Coordinate Exchange Kernel... ");
  CoorExchange_MAIN(&Result, CE_OPT, D_INFO, VERBOSE);
  
  return List::create(Named("DESIGN") = wrap(Result.DESIGN),
                      Named("DESIGN_VAL") = wrap(Result.DESIGN_VAL),
                      Named("fvalHist") = wrap(Result.fvalHist));
}

//[[Rcpp::export]]
List DiscreteDesignColPair(Rcpp::List CP_INFO_LIST, Rcpp::List D_INFO_LIST, bool PARALLEL, bool VERBOSE)
{
	int nc = omp_get_num_procs();
	if (nc < 3) nc = 2;
	if (PARALLEL) omp_set_num_threads(nc - 1); else omp_set_num_threads(1);
	
	DESIGN_INFO D_INFO = {};	getInfoStruct(&D_INFO, D_INFO_LIST);
	CP_OPTIONS CP_OPT = {};		getColumnPairStruct(&CP_OPT, CP_INFO_LIST);
  CP_Result Result = {};
	if (VERBOSE) Rprintf("\nCalling Cpp Columnwaise-Pairwise Kernel... ");
  ColumnPair_MAIN(&Result, CP_OPT, D_INFO, VERBOSE);
  
  return List::create(Named("DESIGN") = wrap(Result.DESIGN),
                      Named("DESIGN_VAL") = wrap(Result.DESIGN_VAL),
                      Named("fvalHist") = wrap(Result.fvalHist));
}


//[[Rcpp::export]]
List cppDesignCriterion(mat DESIGN, Rcpp::List D_INFO_LIST, int typeCrit = -1)
{
	DESIGN_INFO D_INFO = {};
	getInfoStruct(&D_INFO, D_INFO_LIST);
	
	double DESIGN_VAL;
	arma::mat VAL_MAT = DESIGNCRITERION(DESIGN_VAL, DESIGN, D_INFO, typeCrit);
	/* OUTPUT */	
  return List::create(Named("valMat") = wrap(VAL_MAT), 
  										Named("criVal") = wrap(DESIGN_VAL));	
}

