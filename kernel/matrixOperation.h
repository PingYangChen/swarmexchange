

arma::mat invRmOneRow(arma::mat &m_inv, arma::mat &m, arma::uword &rid)
{
	arma::rowvec rv = m.row(rid);
	arma::vec cv = rv.t();
	double g_inv = 1.0 - arma::as_scalar(rv * m_inv * cv);
	arma::mat inv_out = m_inv + g_inv * (m_inv * cv * rv * m_inv);
	return inv_out
}

arma::mat invAddOneRow(arma::mat &m_inv, arma::mat &m, arma::rowvec &rv)
{
	arma::vec cv = rv.t();
	double g_inv = 1.0 + arma::as_scalar(rv * m_inv * cv);
	arma::mat inv_out = m_inv - g_inv * (m_inv * cv * rv * m_inv);
	return inv_out
}

arma::mat invAddOneCol(arma::mat &m_inv, arma::mat &m, arma::vec &cv)
{
	int n = m_inv.n_rows;
	arma::vec mcv = m * cv;
	double k_inv = 1.0/(arma::as_scalar(cv.t() * cv));
	arma::mat i00 = m_inv + k_inv * (m_inv * mcv * mcv.t() * m_inv);
	arma::vec i01 = (-1.0) * k_inv * (m_inv * mcv);
	arma::mat inv_out(n+1, n+1);
	inv_out.submat(0, 0, n-1, n-1) = i00;
	inv_out.submat(0, n, n-1, n) = i01;
	inv_out.submat(n, 0, n, n-1) = i01.t();
	inv_out(n, n) = k_inv;
	return inv_out
}