namespace ExpFit {
  extern "C" int dgels_(char *trans, int *m, int *n, int *nrhs, double *a, int *lda, double *b,
			int *ldb, double *work, int *lwork, int *info);
  double dgels_norm(double gen, const double *x, const double *y, int n, char *trans, int *M,
		    int *N, int *nrhs, double *A, int *lda, double *B, int *ldb, double *work,
		    int *lwork, int *info, double affine);
  bool fit_decay(const double *x, const double *y, int n, double gen_min, double gen_max,
		 double *gen_ans, double *amp_exp, double *amp_aff, bool know_affine);
}
