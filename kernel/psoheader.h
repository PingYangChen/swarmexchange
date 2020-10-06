// Assistant Tools for PSO

#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// DEFINE STUCTURES
typedef struct {
  int typeCrit, nRun, nFactor, nLevel, nModel, nMainEff, nTwofi, balance; 
	arma::rowvec labLevel;
	//arma::imat modelSpace;
	arma::icube modelIndices;
} DESIGN_INFO, *Ptr_DESIGN_INFO;

typedef struct {
	int PSO_UPDATE;
	int nSwarm;
	int maximize;
	int maxIter;
	double tol;
	// MIX Operator
	int MIX_C;
	int MIX_R;
	int HYBRIDEXALG;
	// Random Mixing Parameters
	double JFO_R0;
	double JFO_R1;
	double JFO_RV;
	// JFO Parameters
	double JFO_RHO;
	//double JFO_PL;
	//double JFO_PG;
} PSO_OPTIONS, *Ptr_PSO_OPTIONS;

// DEFINE STUCTURES OF PSO PARAMETERS WHICH WILL CHANGE ITERATIVELY
typedef struct {
  // inertia weight
  int JFO_R_DUR; // (int)(w_var*maxIter);
  double JFO_R_CUR; // w0;
  double JFO_R_DEC; // (w0 - w1)/w_varyfor;
} PSO_DYN, *Ptr_PSO_DYN;

typedef struct {
	int maximize;
	int nTry;
	int maxIter;
	double tol;
} CE_OPTIONS, *Ptr_CE_OPTIONS;

typedef struct {
	int maximize;
	int nTry;
	int maxIter;
	int CPk;
	double tol;
} CP_OPTIONS, *Ptr_CP_OPTIONS;


// DEFINE STUCTURES
typedef struct {
  arma::mat GBest;
  double fGBest;
  arma::rowvec fGBestHist;
  arma::mat PBest;
  arma::vec fPBest;	
	arma::imat updateRec;
	//arma::mat JumpProb;
} PSO_Result, *Ptr_PSO_Result;

typedef struct {
  arma::mat DESIGN;
  double DESIGN_VAL;
  arma::vec fvalHist;
} CE_Result, *Ptr_CE_Result;

typedef struct {
  arma::mat DESIGN;
  double DESIGN_VAL;
  arma::vec fvalHist;
} CP_Result, *Ptr_CP_Result;


// DECLARE FUNCTIONS
void matPrintf(const mat &m);
void imatPrintf(const imat &m);
void rvecPrintf(const rowvec &v);
void vecPrintf(const vec &v);
void getInfoStruct(Ptr_DESIGN_INFO Ptr_D_INFO, const Rcpp::List &D_INFO_LIST);
void getAlgStruct(Ptr_PSO_OPTIONS Ptr_PSO_OPT, const Rcpp::List &ALG_INFO_LIST);
void getModelIndices(arma::icube &modelIndices, const arma::imat &modelSpace, const int &nMainEff, const int &nTwofi);
void PSO_MAIN(Ptr_PSO_Result Ptr_PSO_Result, const PSO_OPTIONS &PSO_OPTS, const DESIGN_INFO &D_INFO, bool VERBOSE);
void psoUpdateParticle(const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS, const PSO_DYN &PSO_DYN, 
											 const arma::mat &PBest, const arma::rowvec &GBest, 
											 const arma::vec &fPBest, const double &fGBest, arma::mat &swarm, arma::vec &fSwarm);
void psoUpdateDynPara(const PSO_OPTIONS &PSO_OPTS, const int iter, PSO_DYN &PSO_DYN);
void initDesigns(arma::mat &DESIGN_POOL, const DESIGN_INFO &D_INFO);

arma::mat CoorExchange_CORE(arma::mat DESIGN, double DESIGN_VAL, const DESIGN_INFO &D_INFO);
void getCoorExStruct(Ptr_CE_OPTIONS Ptr_CE_OPT, const Rcpp::List &CE_INFO_LIST);
void CoorExchange_MAIN(Ptr_CE_Result Ptr_CE_Result, const CE_OPTIONS &CE_OPTS, const DESIGN_INFO &D_INFO, bool VERBOSE);

arma::mat ColumnPair_CORE(arma::mat DESIGN, double &DESIGN_VAL, const DESIGN_INFO &D_INFO, int col_j, int maximize);
void getColumnPairStruct(Ptr_CP_OPTIONS Ptr_CP_OPT, const Rcpp::List &CP_INFO_LIST);
void ColumnPair_MAIN(Ptr_CP_Result Ptr_CP_Result, const CP_OPTIONS &CP_OPTS, const DESIGN_INFO &D_INFO, bool VERBOSE);

// INCLUDE SUBFUNCTIONS
#include "designCriterion.h"
#include "exchangeAlgCore.h" 
#include "MIX_OPERATOR.h"
#include "UPDATE_JFO.h"
#include "UPDATE_PSE.h"
#include "UPDATE_RANDOM.h"
#include "psoUpdateParticle.h"
#include "psoUpdateDynPara.h"
#include "initDesigns.h"
#include "psokernel.h"
#include "CoorExchange.h"
#include "ColumnPair.h"

// BODY
void matPrintf(const mat &m) 
{
  for (uword i = 0; i < m.n_rows; i++) {
    for (uword j = 0; j < m.n_cols; j++) Rprintf("%4.4f\t", m(i,j));
    Rprintf("\n");
  }
	Rprintf("\n\n");
}

void imatPrintf(const imat &m) 
{
  for (uword i = 0; i < m.n_rows; i++) {
    for (uword j = 0; j < m.n_cols; j++) Rprintf("%d\t", m(i,j));
    Rprintf("\n");
  }
	Rprintf("\n\n");
}

void rvecPrintf(const rowvec &v) { for (uword i = 0; i < v.n_elem; i++) { Rprintf("%4.4f\t", v(i)); } Rprintf("\n"); }
void vecPrintf(const vec &v) { for (uword i = 0; i < v.n_elem; i++) { Rprintf("%4.4f\t", v(i)); } Rprintf("\n"); }
void irvecPrintf(const irowvec &v) { for (uword i = 0; i < v.n_elem; i++) { Rprintf("%d\t", v(i)); } Rprintf("\n"); }
void ivecPrintf(const ivec &v) { for (uword i = 0; i < v.n_elem; i++) { Rprintf("%d\t", v(i)); } Rprintf("\n"); }

//
void getInfoStruct(Ptr_DESIGN_INFO Ptr_D_INFO, const Rcpp::List &D_INFO_LIST)
{
	int nFactor = as<int>(D_INFO_LIST["nFactor"]);
	int nModel = as<int>(D_INFO_LIST["nModel"]);
	int nMainEff = as<int>(D_INFO_LIST["nMainEff"]);
	int nTwofi = as<int>(D_INFO_LIST["nTwofi"]);

	Ptr_D_INFO->typeCrit = as<int>(D_INFO_LIST["typeCrit"]);
	
	Ptr_D_INFO->nRun = as<int>(D_INFO_LIST["nRun"]);
	Ptr_D_INFO->nFactor = nFactor;
	Ptr_D_INFO->nLevel = as<int>(D_INFO_LIST["nLevel"]);
	Ptr_D_INFO->balance = as<int>(D_INFO_LIST["balance"]);

	Rcpp::NumericVector labLevelTmp = as<NumericVector>(D_INFO_LIST["labLevel"]);
	arma::rowvec labLevel(labLevelTmp.begin(), labLevelTmp.size(), false);
	Ptr_D_INFO->labLevel = labLevel;
	
	Ptr_D_INFO->nModel = nModel;
	Ptr_D_INFO->nMainEff = nMainEff;
	Ptr_D_INFO->nTwofi= nTwofi;

	Rcpp::IntegerMatrix modelSpaceTmp = as<IntegerMatrix>(D_INFO_LIST["modelSpace"]);
	arma::imat modelSpace(modelSpaceTmp.begin(), modelSpaceTmp.nrow(), modelSpaceTmp.ncol(), false);
	
	arma::icube modelIndices(nMainEff + nTwofi + 1, nFactor, nModel, fill::zeros);
	getModelIndices(modelIndices, modelSpace, nMainEff, nTwofi);
	Ptr_D_INFO->modelIndices = modelIndices;
}

//
void getAlgStruct(Ptr_PSO_OPTIONS Ptr_PSO_OPT, const Rcpp::List &ALG_INFO_LIST)
{
	Ptr_PSO_OPT->PSO_UPDATE = as<int>(ALG_INFO_LIST["PSO_UPDATE"]);
	Ptr_PSO_OPT->nSwarm = as<int>(ALG_INFO_LIST["nSwarm"]);
	Ptr_PSO_OPT->maxIter = as<int>(ALG_INFO_LIST["maxIter"]);
	Ptr_PSO_OPT->maximize = as<int>(ALG_INFO_LIST["maximize"]);
	Ptr_PSO_OPT->tol = as<double>(ALG_INFO_LIST["tol"]);
	// MIX Operator
	Ptr_PSO_OPT->MIX_C = as<int>(ALG_INFO_LIST["MIX_C"]);
	Ptr_PSO_OPT->MIX_R = as<int>(ALG_INFO_LIST["MIX_R"]);
	Ptr_PSO_OPT->HYBRIDEXALG = as<int>(ALG_INFO_LIST["HYBRIDEXALG"]);
	// Random Mixing Parameters
	Ptr_PSO_OPT->JFO_RV = as<double>(ALG_INFO_LIST["JFO_RV"]);
	Ptr_PSO_OPT->JFO_R0 = as<double>(ALG_INFO_LIST["JFO_R0"]);
	Ptr_PSO_OPT->JFO_R1 = as<double>(ALG_INFO_LIST["JFO_R1"]);
	// JFO Parameters
	Ptr_PSO_OPT->JFO_RHO = as<double>(ALG_INFO_LIST["JFO_RHO"]);
	//Ptr_PSO_OPT->JFO_PL = as<double>(ALG_INFO_LIST["JFO_PL"]);
	//Ptr_PSO_OPT->JFO_PG = as<double>(ALG_INFO_LIST["JFO_PG"]);
}

//
void getCoorExStruct(Ptr_CE_OPTIONS Ptr_CE_OPT, const Rcpp::List &CE_INFO_LIST)
{
	Ptr_CE_OPT->maxIter = as<int>(CE_INFO_LIST["maxIter"]);
	Ptr_CE_OPT->nTry = as<int>(CE_INFO_LIST["nTry"]);
	Ptr_CE_OPT->maximize = as<int>(CE_INFO_LIST["maximize"]);
	Ptr_CE_OPT->tol = as<double>(CE_INFO_LIST["tol"]);
}

//
void getColumnPairStruct(Ptr_CP_OPTIONS Ptr_CP_OPT, const Rcpp::List &CP_INFO_LIST)
{
	Ptr_CP_OPT->maxIter = as<int>(CP_INFO_LIST["maxIter"]);
	Ptr_CP_OPT->nTry = as<int>(CP_INFO_LIST["nTry"]);
	Ptr_CP_OPT->maximize = as<int>(CP_INFO_LIST["maximize"]);
	Ptr_CP_OPT->CPk = as<int>(CP_INFO_LIST["CPk"]);
	Ptr_CP_OPT->tol = as<double>(CP_INFO_LIST["tol"]);
}

//
void getModelIndices(arma::icube &modelIndices, const arma::imat &modelSpace, const int &nMainEff, const int &nTwofi)
{
	for (uword i = 0; i < modelSpace.n_rows; i++) {
		arma::imat tmpModel(modelIndices.n_rows, modelIndices.n_cols, fill::zeros);
		for (uword j = 0; j < (uword)(nMainEff + nTwofi); j++) {
			if (j < (uword)nMainEff) {
				tmpModel(j+1, modelSpace(i, j)) = 1;
			} else {
				int tmpIdx = nMainEff + 2*(j - nMainEff);
				tmpModel(j+1, modelSpace(i, tmpIdx)) = 1;
				tmpModel(j+1, modelSpace(i, tmpIdx + 1)) = 1;
			}
		}
		modelIndices.slice(i) = tmpModel;
	}
}