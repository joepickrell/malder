#include <cmath>
#include <algorithm>
#include "ExpFit.hpp"

namespace ExpFit {
  double dgels_norm(double gen, const double *x, const double *y, int n, char *trans, int *M,
		    int *N, int *nrhs, double *A, int *lda, double *B, int *ldb, double *work,
		    int *lwork, int *info, double affine) {
    if (*N == 2) { // fit affine
      for (int i = 0; i < n; i++) {
	A[i] = 1;
	A[i+n] = exp(-gen*x[i]);
	B[i] = y[i];
      }
    }
    else { // know affine
      for (int i = 0; i < n; i++) {
	A[i] = exp(-gen*x[i]);
	B[i] = y[i] - affine;
      }
    }
    dgels_(trans, M, N, nrhs, A, lda, B, ldb, work, lwork, info);
    double norm = 0;
    for (int i = *N; i < *M; i++)
      norm += B[i]*B[i];
    return norm;
  }

  bool fit_decay(const double *x, const double *y, int n, double gen_min, double gen_max,
		 double *gen_ans, double *amp_exp, double *amp_aff, bool know_affine) {
    char trans = 'N'; // A is not transposed
    int M = n; // number of rows of A
    int N = know_affine ? 1 : 2; // number of cols of A: 2 (const and exponential)
    int nrhs = 1; // one column of B
    double A[2*n];
    int lda = n;
    double B[n];
    int ldb = n;
    double work[6*n];
    int lwork = 6*n;
    int info;
  
    double gen_multiplier = 1.2, gen_res = 0.005;
    double best_norm = 1e9, best_gen = gen_min;
    for (double gen = gen_min; gen < gen_max; gen *= gen_multiplier) {
      double norm = dgels_norm(gen, x, y, n, &trans, &M, &N, &nrhs, A, &lda, B, &ldb, work, &lwork,
			       &info, *amp_aff);
      if (norm < best_norm) {
	best_norm = norm;
	best_gen = gen;
      }
    }
    double lo = std::max(gen_min, best_gen / pow(gen_multiplier, 2));
    double hi = std::min(gen_max, best_gen * pow(gen_multiplier, 2));
    while (hi-lo > gen_res) {
      double gen1 = (2*lo+hi)/3;
      double norm1 = dgels_norm(gen1, x, y, n, &trans, &M, &N, &nrhs, A, &lda, B, &ldb, work,
				&lwork, &info, *amp_aff);
      double gen2 = (lo+2*hi)/3;
      double norm2 = dgels_norm(gen2, x, y, n, &trans, &M, &N, &nrhs, A, &lda, B, &ldb, work,
				&lwork, &info, *amp_aff);
      if (norm1 < norm2)
	hi = gen2;
      else
	lo = gen1;
    }
    *gen_ans = (lo+hi)/2;
    // do fitting to get solution coefficients A*e^-nd + C: A = B[0], C = B[1]
    dgels_norm(*gen_ans, x, y, n, &trans, &M, &N, &nrhs, A, &lda, B, &ldb, work, &lwork, &info,
	       *amp_aff);
    if (!know_affine) {
      *amp_aff = B[0];
      *amp_exp = B[1];
    }
    else
      *amp_exp = B[0];
    return gen_min + gen_res < *gen_ans && *gen_ans < gen_max - gen_res;
  }
}
