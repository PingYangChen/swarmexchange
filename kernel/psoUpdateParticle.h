
void psoUpdateParticle(const DESIGN_INFO &D_INFO, const PSO_OPTIONS &PSO_OPTS, const PSO_DYN &PSO_DYN, 
											 const arma::mat &PBest, const arma::mat &GBest, 
											 const arma::vec &fPBest, const double &fGBest, arma::mat &swarm, arma::vec &fSwarm,
											 arma::urowvec &WHICH_MIX_THIS_ITER)
{
	int nSwarm = PSO_OPTS.nSwarm;
  int maximize 	= PSO_OPTS.maximize;
	int PSO_UPDATE = PSO_OPTS.PSO_UPDATE;
	int nRun = D_INFO.nRun;
	int iSwarm;
	switch (PSO_UPDATE) {
		case 0: // PSE
		{
			#pragma omp parallel private(iSwarm) 
			{
				#pragma omp for
				for (iSwarm = 0; iSwarm < nSwarm; iSwarm++) {
					arma::mat DESIGN = swarm.rows(iSwarm*nRun, (iSwarm + 1)*nRun - 1);
					arma::mat tmpPBest = PBest.rows(iSwarm*nRun, (iSwarm + 1)*nRun - 1);
					double DESIGN_VAL = fSwarm(iSwarm);
					arma::uword USE_MIX;
					UPDATE_PSE(DESIGN, DESIGN_VAL, USE_MIX, GBest, tmpPBest, maximize, D_INFO, PSO_OPTS, PSO_DYN);
					swarm.rows(iSwarm*nRun, (iSwarm + 1)*nRun - 1) = DESIGN;
					fSwarm(iSwarm) = DESIGN_VAL;
					WHICH_MIX_THIS_ITER(iSwarm) = USE_MIX;
				}
			}
		break;
		}
		case 9: // JFO
		{
			#pragma omp parallel private(iSwarm) 
			{
				#pragma omp for
				for (iSwarm = 0; iSwarm < nSwarm; iSwarm++) {
					arma::mat DESIGN = swarm.rows(iSwarm*nRun, (iSwarm + 1)*nRun - 1);
					arma::mat tmpPBest = PBest.rows(iSwarm*nRun, (iSwarm + 1)*nRun - 1);
					double DESIGN_VAL = fSwarm(iSwarm);
					UPDATE_JFO(DESIGN, DESIGN_VAL, GBest, tmpPBest, maximize, D_INFO, PSO_OPTS, PSO_DYN);
					swarm.rows(iSwarm*nRun, (iSwarm + 1)*nRun - 1) = DESIGN;
					fSwarm(iSwarm) = DESIGN_VAL;
				}
			}
		break;
		}
		case 99: // RAND
		{
			#pragma omp parallel private(iSwarm) 
			{
				#pragma omp for
				for (iSwarm = 0; iSwarm < nSwarm; iSwarm++) {
					arma::mat DESIGN = swarm.rows(iSwarm*nRun, (iSwarm + 1)*nRun - 1);
					double DESIGN_VAL = fSwarm(iSwarm);
					UPDATE_RANDOM(DESIGN, DESIGN_VAL, maximize, D_INFO, PSO_OPTS);
					swarm.rows(iSwarm*nRun, (iSwarm + 1)*nRun - 1) = DESIGN;
					fSwarm(iSwarm) = DESIGN_VAL;
				}
			}
		break;
		}
		case 10 ... 19: // CSO
		{
			arma::mat GMean(GBest.n_rows, GBest.n_cols);
			imat rSwarm = arma::randi(GBest.n_rows, GBest.n_cols, distr_param(0, nSwarm - 1));
			for (int i = 0; i < (int)GMean.n_rows; i++) {
				for (int j = 0; j < (int)GMean.n_cols; j++) {
					GMean(i, j) = swarm(rSwarm(i, j)*nRun + i, j);
				}
			}
			uvec COUPLE_ID = arma::repmat(arma::linspace<arma::uvec>(0, nSwarm/2 - 1, nSwarm/2), 2, 1);
			vec RANDNUM = arma::randu(nSwarm, 1);
			uvec RANDIDX = sort_index(RANDNUM);
			COUPLE_ID = COUPLE_ID(RANDIDX);
			uword WINNER, LOSER;
			#pragma omp parallel private(iSwarm, WINNER, LOSER) 
			{
				#pragma omp for
				for (uword iSwarm = 0; iSwarm < ((uword)nSwarm/2); iSwarm++) {
					uvec COMPPAIR = find(COUPLE_ID == iSwarm);
					if (maximize == 0) {
						if (fSwarm(COMPPAIR(0)) < fSwarm(COMPPAIR(1))) { 
							WINNER = COMPPAIR(0); LOSER = COMPPAIR(1); 
						} else { 
							WINNER = COMPPAIR(1); LOSER = COMPPAIR(0);
						}
					} else {
						if (fSwarm(COMPPAIR(0)) > fSwarm(COMPPAIR(1))) {
							WINNER = COMPPAIR(0); LOSER = COMPPAIR(1);
						} else {
							WINNER = COMPPAIR(1); LOSER = COMPPAIR(0);
						}
					}
					arma::mat DESIGN = swarm.rows(LOSER*nRun, (LOSER + 1)*nRun - 1);
					arma::mat DESIGN_WIN = swarm.rows(WINNER*nRun, (WINNER + 1)*nRun - 1);
					double DESIGN_VAL = fSwarm(LOSER); double WINNER_VAL = fSwarm(WINNER); 
					int tmpUPDATE = PSO_UPDATE - 10;
					arma::uword USE_MIX;
					switch (tmpUPDATE) {
						case 0:
						{
							UPDATE_PSE(DESIGN, DESIGN_VAL, USE_MIX, GMean, DESIGN_WIN, maximize, D_INFO, PSO_OPTS, PSO_DYN);
							break;
						}
						case 9:
						{
							UPDATE_JFO(DESIGN, DESIGN_VAL, GMean, DESIGN_WIN, maximize, D_INFO, PSO_OPTS, PSO_DYN);
							break;
						}
					}
					swarm.rows(LOSER*nRun, (LOSER + 1)*nRun - 1) = DESIGN;
					fSwarm(LOSER) = DESIGN_VAL;
					WHICH_MIX_THIS_ITER(iSwarm) = USE_MIX;
				}
			}
			break;	
		}
	}	
}

