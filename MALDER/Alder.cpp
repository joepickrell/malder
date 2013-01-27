#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <omp.h>
#include <fftw3.h>

#include "nicklib.h"

#include "MiscUtils.hpp"
#include "CorrJack.hpp"
#include "ExpFitALD.hpp"
#include "Timer.hpp"
#include "Alder.hpp"
#include "AlderParams.hpp"

#define FFT_CONVOLUTION

namespace ALD {

  using std::cout;
  using std::endl;
  using std::flush;
  using std::string;
  using std::vector;
  using std::pair;
  using std::make_pair;
  using std::max;
  using std::min;
  using std::lower_bound;
  using std::accumulate;

  double sq(double x) { return x*x; }

  Alder::AffineData::AffineData(int n=0, int _num_refs=0) {
    count = 0;
    ws = ss = s2 = 0.0;
    wg = gg = gs = vector <double> (n);
    gigj = vector < vector <double> > (n, vector <double> (n));
    num_refs = _num_refs;
  }

  const int Alder::LIM_SIGNIFICANCE_FAILURES = 2;
  const double Alder::LD_COS_SIGNIF_THRESH = 0.05;
  const double Alder::PCA_VARIANCE_THRESH = 0.9;

  string Alder::format_mean_std(pair <double, double> mean_std) {
    if (isnan(mean_std.first)) return "too much noise";
    char buf[100]; 
    if (isnan(mean_std.second)) 
      sprintf(buf, "%6.3f +/-   ?  ", mean_std.first);
    else
      sprintf(buf, "%6.3f +/- %5.3f", mean_std.first, mean_std.second);
    return string(buf);
  }

  double Alder::compute_geno_mean(int s, const char *geno, int stride) {
    int sum_x = 0, n = 0;
    for (int i = 0; i < stride; i++) {
      int x = geno[s*stride+i];
      if (x != 9) {
	sum_x += x; n++;
      }
    }
    return n == 0 ? NAN : sum_x / (double) n;
  }

  double Alder::compute_ld(int s1, int s2, const char *geno, int stride) {
    int sum_x = 0, sum_y = 0, sum_xy = 0, n = 0;
    for (int i = 0; i < stride; i++) {
      int x = geno[s1*stride+i], y = geno[s2*stride+i];
      if (x != 9 && y != 9) {
	sum_x += x; sum_y += y; sum_xy += x*y; n++;
      }
    }
    return n <= 1 ? NAN : (sum_xy - sum_x * sum_y / (double) n) / (n-1);
  }

  double Alder::compute_ld(int s1, int s2) {
    return compute_ld(s1, s2, mixed_geno, num_mixed_indivs);
  }

  double Alder::compute_polyache_central_moment11sq(int s1, int s2, const char *geno, int stride) {
    int n = 0;
    double S10 = 0, S01 = 0, S20 = 0, S11 = 0, S02 = 0, S21 = 0, S12 = 0, S22 = 0;
    for (int i = 0; i < stride; i++) {
      int x = geno[s1*stride+i], y = geno[s2*stride+i];
      if (x != 9 && y != 9) {
	n++;
	S10 += x;
	S01 += y;
	S20 += sq(x);
	S11 += x*y;
	S02 += sq(y);
	S21 += sq(x)*y;
	S12 += x*sq(y);
	S22 += sq(x)*sq(y);
      }
    }
    double S0 = n;
    double S0p2 = S0*(S0-1);
    double S0p3 = S0p2*(S0-2);
    double S0p4 = S0p3*(S0-3);

    double ans = n <= 3 ? NAN :
      (sq(S10*S01) - S20*sq(S01) - S02*sq(S10) + S02*S20) / S0p4
      + 2*(S21*S01 - S10*S11*S01 + S12*S10) * (2*S0p3+S0p4) / (S0p3*S0p4)
      + (sq(S11) * (2*S0p2*S0p3 + S0p3*S0p4 + 2*S0p2*S0p4)
	 - S22 * (6*S0p2*S0p3 + S0p3*S0p4 + 4*S0p2*S0p4)) / (S0p2*S0p3*S0p4);

    return ans;
  }

  double Alder::compute_polyache_central_moment11sq(int s1, int s2) {
    return compute_polyache_central_moment11sq(s1, s2, mixed_geno, num_mixed_indivs);
  }

  bool Alder::x2_suff_accurate(pair <double, double> x2_mean_std) {
    return x2_mean_std.first > 5*x2_mean_std.second;
  }

  // stores polyache data in:
  // - test_data.count, test_data.sum_x2 and test_data.sum_y2 (same)
  // - ref_data.count, ref_data.sum_x2 and ref_data.sum_y2 (same)
  void Alder::compute_ld_corr_terms(int ref_ind, double bin_min, double bin_max,
				   CorrJack &corr_data, CorrJack &test_data, CorrJack &ref_data,
				   bool use_early_exit, bool compute_corr_data,
				   bool compute_polyache_data) {
    bool done_ld_prod = !compute_corr_data;
    bool done_polyache_test = !compute_polyache_data, done_polyache_ref = !compute_polyache_data;

    // iterate through s1 in layers:
    // at end of each offset layer, jackknife to decide if enough precision
    const int num_early_checks = 6;
    const int s1_stride = 1<<num_early_checks;
    int num_checks_left = num_early_checks+1;
    for (int s1_offset = 0; s1_offset < s1_stride; s1_offset++) {
#pragma omp parallel for schedule(static,1)
      for (int c = 0; c < num_chroms_used; c++) {
	int snp_start = chrom_start_inds[c], snp_end = chrom_start_inds[c+1];
	for (int s1 = snp_start+s1_offset; s1 < snp_end; s1 += s1_stride) {
	  for (int s2 = lower_bound(snp_pos.begin()+s1, snp_pos.begin()+snp_end,
				    snp_pos[s1] + bin_min) - snp_pos.begin();
	       s2 < snp_end && snp_pos[s2] < snp_pos[s1] + bin_max; s2++) {
	    if (!done_ld_prod) {
	      double LD_test = compute_ld(s1, s2);
	      double LD_ref = compute_ld(s1, s2, ref_genos[ref_ind], num_ref_indivs[ref_ind]);
	      corr_data.data[c].add_term(LD_test, LD_ref); // checks for nan
	    }
	    if (!done_polyache_test)
	      test_data.data[c].add_unbiased_sq_term(compute_polyache_central_moment11sq(s1, s2));
	    if (!done_polyache_ref)
	      ref_data.data[c].add_unbiased_sq_term(compute_polyache_central_moment11sq(s1, s2,
					             ref_genos[ref_ind], num_ref_indivs[ref_ind]));
	  }
	}
      }
      if (use_early_exit && !((s1_offset+1) & s1_offset)) { // power-of-2 stride
	if (!done_ld_prod) {
	  pair <double, double> cos_mean_std = corr_data.jackknife_cos();
	  if (erfc(cos_mean_std.first/cos_mean_std.second/sqrt(2.0)) * num_checks_left
	      < LD_COS_SIGNIF_THRESH)
	    done_ld_prod = true;
	}
	if (!done_polyache_test)
	  done_polyache_test = x2_suff_accurate(test_data.jackknife_x2_avg());
	if (!done_polyache_ref)
	  done_polyache_ref = x2_suff_accurate(ref_data.jackknife_x2_avg());
	if (done_ld_prod && done_polyache_test && done_polyache_ref) return;
	num_checks_left--;
      }
    }
  }

  double Alder::compute_polyache(int s1, int s2, double pAx, double pAy) {
    int n = 0;
    double S10 = 0, S01 = 0, S20 = 0, S11 = 0, S02 = 0, S21 = 0, S12 = 0, S22 = 0;
    for (int i = 0; i < num_mixed_indivs; i++) {
      int x = mixed_geno[s1*num_mixed_indivs+i], y = mixed_geno[s2*num_mixed_indivs+i];
      if (x != 9 && y != 9) {
	n++;
	S10 += x;
	S01 += y;
	S20 += sq(x);
	S11 += x*y;
	S02 += sq(y);
	S21 += sq(x)*y;
	S12 += x*sq(y);
	S22 += sq(x)*sq(y);
      }
    }
    double S0 = n;
    double S0p2 = S0*(S0-1);
    double S0p3 = S0p2*(S0-2);
    double S0p4 = S0p3*(S0-3);

    double ans = n <= 3 ? NAN :
      (4*pAx*pAy * (-S10 * S01 / S0p2
		    + S11 * (S0p2 + S0) / (S0 * S0p2))
       + 2*pAx * ((S10 * S01*S01 - S02 * S10) / S0p3
		  + (S12 - S11 * S01) * ((2*S0p2 + S0p3) / (S0p2 * S0p3)))
       + 2*pAy * ((S01 * S10*S10 - S20 * S01) / S0p3
		  + (S21 - S11 * S10) * ((2*S0p2 + S0p3) / (S0p2 * S0p3)))
       + (-S10*S10 * S01*S01 + S20 * S01*S01 + S02 * S10*S10 - S02 * S20) / S0p4
       + (S10 * S11 * S01 - S21 * S01 - S10 * S12) * (4*S0p3 + S0p4) / (S0p3 * S0p4)
       - S11*S11 * (2*S0p3 + S0p4) / (S0p3 * S0p4)
       + 2*S22 * (3*S0p3 + S0p4) / (S0p3 * S0p4)) / 4;
    return ans;
  }

  void Alder::convolve_accum(fftw_plan *plans, int Nby2, fftw_complex *z_accum,
			    fftw_complex *z1, fftw_complex *z2, double scale) {
    fftw_execute(plans[0]);
    fftw_execute(plans[1]);
    for (int b = 0; b <= Nby2; b++) { // z1bar * z2
      z_accum[b][0] += (z1[b][0] * z2[b][0] + z1[b][1] * z2[b][1]) * scale;
      z_accum[b][1] += (z1[b][0] * z2[b][1] - z1[b][1] * z2[b][0]) * scale;
    }
  }
  
  void Alder::self_convolve_accum(fftw_plan *plan, int Nby2, fftw_complex *z_accum,
				 fftw_complex *z1, double scale) {
    fftw_execute(*plan);
    for (int b = 0; b <= Nby2; b++) // z1bar * z1
      z_accum[b][0] += (sq(z1[b][0]) + sq(z1[b][1])) * scale;
  }

  // returns binned pairs: (weighted LD, count of pairs in bin)
  // also, affine_data contains info for computing affine term
  vector < pair <double, double> > Alder::run_chrom(int chrom, int num_refs,
						   const vector <double> &weights, double binsize,
						   int numbins, int mincount,
						   AffineData &affine_data) {
    
    int snp_start = chrom_start_inds[chrom], snp_end = chrom_start_inds[chrom+1];
    int numbins_chrom = (snp_pos[snp_end-1] - snp_pos[snp_start]) / binsize + 1;
    vector < pair <double, double> > ans(numbins);
    affine_data = AffineData(num_mixed_indivs, num_refs);
    affine_data.count = count(snp_num_missing.begin()+snp_start,
			      snp_num_missing.begin()+snp_end, 0);
    
    // find sufficient shift so that 2^shift > 2*num_bins_chrom (need space for convolution)
    int shift = 0;
    while ((1<<shift) < numbins_chrom) shift++;
    shift++;
    int N = 1<<shift, Nby2 = N>>1;
    double onebyN = 1.0/N;

    // allocate memory and make fftw plans
    double *fx = (double *) fftw_malloc(sizeof(double)<<shift);
    double *gy = (double *) fftw_malloc(sizeof(double)<<shift);    
    fftw_complex *fft_fx = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)<<shift);
    fftw_complex *fft_gy = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)<<shift);
    fftw_plan plans[2];

    fftw_complex *rev_c_arr = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)<<shift);
    double *rev_r_arr = (double *) fftw_malloc(sizeof(double)*(N+1));
    rev_r_arr[N] = 0; // useful for convenience later in summing stuff from right end
    fftw_plan rev_plan;
    
#pragma omp critical
    {
      plans[0] = fftw_plan_dft_r2c_1d(N, fx, fft_fx, FFTW_ESTIMATE);
      plans[1] = fftw_plan_dft_r2c_1d(N, gy, fft_gy, FFTW_ESTIMATE);
      rev_plan = fftw_plan_dft_c2r_1d(N, rev_c_arr, rev_r_arr, FFTW_ESTIMATE);
    }

    // count: this runs for both 2-ref and single-ref polyache
    memset(fx, 0, sizeof(double)<<shift);
    memset(gy, 0, sizeof(double)<<shift);
    for (int s = snp_start; s < snp_end; s++) {
      if (snp_ignore[s]) continue;

      int b = snp_bin[s];
      if (snp_num_missing[s] == 0) {
	fx[b] += 1.0;
	gy[b] += 0.5;
      }
      else {
	gy[b] += 1.0;
      }
    }
#ifdef FFT_CONVOLUTION
    memset(rev_c_arr, 0, sizeof(fftw_complex)<<shift);
    convolve_accum(plans, Nby2, rev_c_arr, fft_fx, fft_gy);
    fftw_execute(rev_plan);
    for (int b = 0; b < min(numbins, numbins_chrom); b++)
      ans[b].second = (rev_r_arr[b] + rev_r_arr[N-b]) * onebyN;
#else
    for (int b1 = 0; b1 < numbins_chrom; b1++)
      for (int b2 = max(0, b1-numbins+1); b2 < numbins_chrom && b2-b1 < numbins; b2++)
	ans[abs(b2-b1)].second += fx[b1] * gy[b2];
#endif

    // indivs and all
    memset(rev_c_arr, 0, sizeof(fftw_complex)<<shift); // clear; accum all terms before rev fft
    int n = num_mixed_indivs;
    if (num_refs == 2) {
      for (int i = 0; i <= num_mixed_indivs; i++) { // i == num_mixed_indivs is for the sum term
	memset(fx, 0, sizeof(double)<<shift);
	memset(gy, 0, sizeof(double)<<shift);
	for (int s = snp_start; s < snp_end; s++) {
	  if (snp_ignore[s]) continue;
	  int k = snp_num_missing[s];

	  int b = snp_bin[s];
	  if (i < num_mixed_indivs) { // indiv
	    int gtype = mixed_geno[s*num_mixed_indivs+i];
	    if (k == 0) affine_data.wg[i] += gtype * weights[s]; // for affine term
	    fx[b] += (k==0) * gtype * weights[s];
	    gy[b] += (gtype*(gtype!=9) + (gtype==9)*snp_sum[s]/(double) (n-k))
	      * weights[s] / ((1+(k==0)) * (n-k-1));
	  }
	  else { // sum term
	    if (k == 0) affine_data.ws += snp_sum[s] * weights[s]; // for affine term
	    fx[b] += (k==0) * snp_sum[s] * weights[s];
	    gy[b] -= snp_sum[s] * weights[s] / ((1+(k==0)) * (n-k) * (n-k-1));
	  }
	}
	  
#ifdef FFT_CONVOLUTION
	convolve_accum(plans, Nby2, rev_c_arr, fft_fx, fft_gy);
#else
	for (int b1 = 0; b1 < numbins_chrom; b1++)
	  for (int b2 = max(0, b1-numbins+1); b2 < numbins_chrom && b2-b1 < numbins; b2++)
	    ans[abs(b2-b1)].first += fx[b1] * gy[b2];
#endif
      }
    }
    else { // polyache
      /*
	original form:

      double ans = (4*pAx*pAy * (-S10 * S01 / S0p2
				 + S11 * (S0p2 + S0) / (S0 * S0p2))
		    + 2*pAx * ((S10 * S01*S01 - S02 * S10) / S0p3
			       + (S12 - S11 * S01) * ((2*S0p2 + S0p3) / (S0p2 * S0p3)))
		    + 2*pAy * ((S01 * S10*S10 - S20 * S01) / S0p3
			       + (S21 - S11 * S10) * ((2*S0p2 + S0p3) / (S0p2 * S0p3)))
		    + (-S10*S10 * S01*S01 + S20 * S01*S01 + S02 * S10*S10 - S02 * S20) / S0p4
		    + (S10 * S11 * S01 - S21 * S01 - S10 * S12) * (4*S0p3 + S0p4) / (S0p3 * S0p4)
		    - S11*S11 * (2*S0p3 + S0p4) / (S0p3 * S0p4)
		    + 2*S22 * (3*S0p3 + S0p4) / (S0p3 * S0p4)) / 4;

        corresponding affine terms are in comments after each block
	notes:
	- symmetric (x,y) <-> (y,x) terms from the original form are combined below as just 2(x,y)
	  - this is ok because we do (all no-missing vs. all no-missing) on the same chromosome
	  - but the symmetry breaks when calculating the inter-chromosome affine term!
	    (see comments in code for affine term calculation)
	- some additional combination of terms is possible
          ... but the quadratic term dominates the run time
      */
      double S0 = n;
      double S0p2 = S0*(S0-1);
      double S0p3 = S0p2*(S0-2);
      double S0p4 = S0p3*(S0-3);
      
      // -4*pAx*pAy * S10 * S01 / S0p2
      memset(fx, 0, sizeof(double)<<shift);
      for (int s = snp_start; s < snp_end; s++)
	if (!snp_ignore[s])
	  fx[snp_bin[s]] += weights[s] * snp_sum[s];
      affine_data.ws = accumulate(fx, fx+numbins_chrom, 0.0);
      self_convolve_accum(plans, Nby2, rev_c_arr, fft_fx, -4 / S0p2);
      //affine_data[c1].ws * affine_data[c2].ws * -4 / S0p2
      
      // 4*pAx*pAy * S11 * (S0p2 + S0) / (S0 * S0p2)
      for (int i = 0; i < n; i++) {
	memset(fx, 0, sizeof(double)<<shift);
	for (int s = snp_start; s < snp_end; s++)
	  if (!snp_ignore[s]) {
	    int gtype = mixed_geno[s*num_mixed_indivs+i];
	    fx[snp_bin[s]] += weights[s] * gtype;
	  }
	affine_data.wg[i] = accumulate(fx, fx+numbins_chrom, 0.0);
	//affine_data[c1].wg[i] * affine_data[c2].wg[i] * 4 * (S0p2 + S0) / (S0 * S0p2)
	self_convolve_accum(plans, Nby2, rev_c_arr, fft_fx, 4 * (S0p2 + S0) / (S0 * S0p2));
      }

      // 4*pAx * S10 * (S01*S01-S02) / S0p3   (combining sym term)
      memset(fx, 0, sizeof(double)<<shift); memset(gy, 0, sizeof(double)<<shift);
      for (int s = snp_start; s < snp_end; s++)
	if (!snp_ignore[s]) {
	  fx[snp_bin[s]] += weights[s] * snp_sum[s];
	  gy[snp_bin[s]] += sq(snp_sum[s]) - snp_sum2[s];
	}
      //affine_data[c1].ws * (affine_data[c2].ss - affine_data[c2].s2) * 4 / S0p3
      convolve_accum(plans, Nby2, rev_c_arr, fft_fx, fft_gy, 4 / S0p3);

      // 4*pAx * (S12 - S11*S01) * ((2*S0p2 + S0p3) / (S0p2 * S0p3))   (combining sym term)
      for (int i = 0; i < n; i++) {
	memset(fx, 0, sizeof(double)<<shift); memset(gy, 0, sizeof(double)<<shift);
	for (int s = snp_start; s < snp_end; s++)
	  if (!snp_ignore[s]) {
	    int gtype = mixed_geno[s*num_mixed_indivs+i];
	    fx[snp_bin[s]] += weights[s] * gtype;
	    gy[snp_bin[s]] += sq(gtype) - gtype * snp_sum[s];
	  }
	//affine_data[c1].wg[i] * (affine_data[c2].gg[i] - affine_data[c2].gs[i]) * 4*(2*S0p2+S0p3) / (S0p2*S0p3)
	convolve_accum(plans, Nby2, rev_c_arr, fft_fx, fft_gy, 4*(2*S0p2+S0p3) / (S0p2*S0p3));
      }

      // (2*S20 - S10*S10) * S01*S01 / S0p4   (combining sym term in first)
      memset(fx, 0, sizeof(double)<<shift); memset(gy, 0, sizeof(double)<<shift);
      for (int s = snp_start; s < snp_end; s++)
	if (!snp_ignore[s]) {
	  fx[snp_bin[s]] += 2*snp_sum2[s] - sq(snp_sum[s]);
	  gy[snp_bin[s]] += sq(snp_sum[s]);
	}
      affine_data.ss = accumulate(gy, gy+numbins_chrom, 0.0);
      //(2*affine_data[c1].s2 - affine_data[c1].ss) * affine_data[c2].ss * 1 / S0p4
      convolve_accum(plans, Nby2, rev_c_arr, fft_fx, fft_gy, 1 / S0p4);

      // -S02 * S20 / S0p4
      memset(fx, 0, sizeof(double)<<shift);
      for (int s = snp_start; s < snp_end; s++)
	if (!snp_ignore[s])
	  fx[snp_bin[s]] += snp_sum2[s];
      affine_data.s2 = accumulate(fx, fx+numbins_chrom, 0.0);
      //affine_data[c1].s2 * affine_data[c2].s2 * -1/S0p4
      self_convolve_accum(plans, Nby2, rev_c_arr, fft_fx, -1/S0p4);

      // (S10 * S11 * S01 - 2 * S21 * S01) * (4*S0p3 + S0p4) / (S0p3 * S0p4)   (combining sym term in second)
      for (int i = 0; i < n; i++) {
	memset(fx, 0, sizeof(double)<<shift); memset(gy, 0, sizeof(double)<<shift);
	for (int s = snp_start; s < snp_end; s++)
	  if (!snp_ignore[s]) {
	    int gtype = mixed_geno[s*num_mixed_indivs+i];
	    fx[snp_bin[s]] += gtype * snp_sum[s] - 2*sq(gtype);
	    gy[snp_bin[s]] += gtype * snp_sum[s];
	  }
	affine_data.gs[i] = accumulate(gy, gy+numbins_chrom, 0.0);
	//(affine_data[c1].gs[i] - 2*affine_data[c1].gg[i]) * affine_data[c2].gs[i] * (4*S0p3 + S0p4) / (S0p3 * S0p4)
	convolve_accum(plans, Nby2, rev_c_arr, fft_fx, fft_gy, (4*S0p3 + S0p4) / (S0p3 * S0p4));
      }

      // -S11*S11 * (2*S0p3 + S0p4) / (S0p3 * S0p4)
      for (int i = 0; i < n; i++)
	for (int j = i+1; j < n; j++) {
	  memset(fx, 0, sizeof(double)<<shift);
	  for (int s = snp_start; s < snp_end; s++)
	    if (!snp_ignore[s]) {
	      int gtype_i = mixed_geno[s*num_mixed_indivs+i];
	      int gtype_j = mixed_geno[s*num_mixed_indivs+j];
	      fx[snp_bin[s]] += gtype_i * gtype_j;
	    }
	  affine_data.gigj[i][j] = accumulate(fx, fx+numbins_chrom, 0.0);
	  //affine_data[c1].gigj[i][j] * affine_data[c2].gigj[i][j] * -2*(2*S0p3 + S0p4) / (S0p3 * S0p4)
	  // factor of 2 for sym (i,j) <-> (j,i)
	  self_convolve_accum(plans, Nby2, rev_c_arr, fft_fx, -2*(2*S0p3 + S0p4) / (S0p3 * S0p4));
	}

      // 2*S22 * (3*S0p3 + S0p4) / (S0p3 * S0p4)... along with square terms from the previous
      for (int i = 0; i < n; i++) {
	memset(fx, 0, sizeof(double)<<shift);
	for (int s = snp_start; s < snp_end; s++)
	  if (!snp_ignore[s]) {
	    int gtype = mixed_geno[s*num_mixed_indivs+i];
	    fx[snp_bin[s]] += sq(gtype);
	  }
	affine_data.gg[i] = accumulate(fx, fx+numbins_chrom, 0.0);
	//affine_data[c1].gg[i] * affine_data[c2].gg[i] * (4*S0p3 + S0p4) / (S0p3 * S0p4)
	self_convolve_accum(plans, Nby2, rev_c_arr, fft_fx, (4*S0p3 + S0p4) / (S0p3 * S0p4));
      }

      // divide the whole thing by 8 (4 for original polyache x 2 for double-count)
      for (int b = 0; b <= Nby2; b++) {
	rev_c_arr[b][0] /= 8;
	rev_c_arr[b][1] /= 8;
      }
    }
#ifdef FFT_CONVOLUTION
    fftw_execute(rev_plan);
    for (int b = 0; b < min(numbins, numbins_chrom); b++)
      ans[b].first = (rev_r_arr[b] + rev_r_arr[N-b]) * onebyN;
#endif
    fftw_destroy_plan(plans[0]); fftw_free(fx); fftw_free(fft_fx);
    fftw_destroy_plan(plans[1]); fftw_free(gy); fftw_free(fft_gy);
    fftw_destroy_plan(rev_plan); fftw_free(rev_c_arr); fftw_free(rev_r_arr);

    return ans;
  }

  vector < pair <double, double> > Alder::run_chrom_naive(int chrom, int num_refs,
							 const vector <double> &weights,
							 double binsize, int numbins,
							 int mincount) {

    int snp_start = chrom_start_inds[chrom], snp_end = chrom_start_inds[chrom+1];
    vector < pair <double, double> > ans(numbins);
    int maxmissing = num_mixed_indivs - mincount;
    for (int s1 = snp_start; s1 < snp_end; s1++) {
      if (snp_num_missing[s1] > maxmissing) continue;
      for (int s2 = s1+1; s2 < snp_end; s2++) {
	if (snp_num_missing[s2] > maxmissing) continue;
	//if ((snp_num_missing[s1]>0)==(snp_num_missing[s2]>0)) continue;
	int bin = SUBTRACT_THEN_BIN ?
	  (snp_pos[s2] - snp_pos[s1]) / binsize : snp_bin[s2] - snp_bin[s1];
	if (bin >= numbins) break;
	double wld = num_refs == 2 ? compute_ld(s1, s2) * weights[s1] * weights[s2]
	  : compute_polyache(s1, s2, weights[s1], weights[s2]);
	if (!isnan(wld)) {
	  ans[bin].first += wld;
	  ans[bin].second += 1;
	}
      }
    }
    return ans;
  }

  pair <double, double> Alder::compute_inter_chrom_affine(const vector <AffineData> &affdats) {
    int num_refs = affdats[0].num_refs;
    double aff = 0.0;
    double tot_pair_count = 0.0;
    double S0 = num_mixed_indivs;
    double S0p2 = S0*(S0-1);
    double S0p3 = S0p2*(S0-2);
    double S0p4 = S0p3*(S0-3);
    // note: this gets called with both all chroms and one left out
    for (int c1 = 0; c1 < (int) affdats.size(); c1++) {
      // need to do both (c1,c2) and (c2,c1) because of sym term simplification in polyache formula
      for (int c2 = 0; c2 < (int) affdats.size(); c2++) {
	if (c2 == c1) continue;
	tot_pair_count += affdats[c1].count * affdats[c2].count;
	if (num_refs == 2) {
	  for (int i = 0; i < num_mixed_indivs; i++)
	    aff += 1.0/(num_mixed_indivs-1)
	      * affdats[c1].wg[i] * affdats[c2].wg[i];
	  aff -= 1.0/(num_mixed_indivs*(num_mixed_indivs-1))
	    * affdats[c1].ws * affdats[c2].ws;
	}
	else { // polyache
	  aff += affdats[c1].ws * affdats[c2].ws * -4 / S0p2;
	  aff += affdats[c1].ws * (affdats[c2].ss - affdats[c2].s2) * 4 / S0p3;
	  aff += (2*affdats[c1].s2 - affdats[c1].ss) * affdats[c2].ss * 1 / S0p4;
	  aff += affdats[c1].s2 * affdats[c2].s2 * -1/S0p4;
	  for (int i = 0; i < num_mixed_indivs; i++) {
	    aff += affdats[c1].wg[i] * affdats[c2].wg[i] * 4 * (S0p2 + S0) / (S0 * S0p2);
	    aff += affdats[c1].wg[i] * (affdats[c2].gg[i] - affdats[c2].gs[i])
	      * 4*(2*S0p2+S0p3) / (S0p2*S0p3);
	    aff += (affdats[c1].gs[i] - 2*affdats[c1].gg[i]) * affdats[c2].gs[i]
	      * (4*S0p3 + S0p4) / (S0p3 * S0p4);
	    aff += affdats[c1].gg[i] * affdats[c2].gg[i] * (4*S0p3 + S0p4) / (S0p3 * S0p4);
	    for (int j = i+1; j < num_mixed_indivs; j++)
	      aff += affdats[c1].gigj[i][j] * affdats[c2].gigj[i][j]
		* -2*(2*S0p3 + S0p4) / (S0p3 * S0p4);
	  }
	}
      }
    }
    if (num_refs == 1) aff /= 4;
    aff /= tot_pair_count;
    return make_pair(aff, tot_pair_count);
  }
  
  void Alder::check_affine_amp(int num_refs, const vector <double> &weights) {
    // check affine amplitude with random sampling
    double wld_sum = 0, wld_sum2 = 0;
    for (int t = 0; ; t++) {
      if (t > 1000000 && (t & (t-1)) == 0)
	printf("sampled interchrom: %10.8f +/- %10.8f\n", wld_sum / t,
	       sqrt((t*wld_sum2 - sq(wld_sum)) / (t*(t-1.0))/t));
      int s1 = rand()%snp_pos.size();
      int s2 = rand()%snp_pos.size();
      if (snp_chrom_ind_squash[s1] == snp_chrom_ind_squash[s2]) { t--; continue; }
      double wld = num_refs == 2 ? compute_ld(s1, s2) * weights[s1] * weights[s2]
	: compute_polyache(s1, s2, weights[s1], weights[s2]);
      if (isnan(wld)) { t--; continue; }
      wld_sum += wld;
      wld_sum2 += sq(wld);
    }
  }

  ExpFitALD Alder::exp_fit_jackknife(const vector <AlderResults> &results_jackknife,
				    double mindis, double maxdis) {

    ExpFitALD fits(num_chroms_used+1, mindis, maxdis, use_jackknife);
    // jc == num_chroms_used is for all (no jackknife)
    // don't bother computing jc in [0..num_chroms_used) if not use_jackknife
#pragma omp parallel for
    for (int jc = use_jackknife ? 0 : num_chroms_used; jc <= num_chroms_used; jc++) {
      const AlderResults &results = results_jackknife[jc];
      int bmin = 0;
      while (bmin < (int) results.d_Morgans.size()
	     && results.d_Morgans[bmin] < mindis - 1e-9)
	bmin++;
      fits.do_fit(jc, &results.d_Morgans[bmin], &results.weighted_LD_avg[bmin],
		  results.d_Morgans.size()-bmin);
    }
    return fits;
  }

  vector <AlderResults> Alder::make_results(
      const vector < vector < pair <double, double> > > &results_allchrom,
      const vector <AffineData> &affine_data_allchrom, double binsize, bool use_naive_algo,
      double fit_start_dis) {
    
    bool use_inter_chrom_affine = !use_naive_algo;
    if (use_inter_chrom_affine) {
      if (use_jackknife && num_chroms_used <= 2) {
	cout << "WARNING: fitting exponential + unconstrained affine (A * exp(-n*d) + C)" << endl;
	cout << "need >= 3 chroms to use jackknife with inter-chrom affine terms" << endl << endl;
	use_inter_chrom_affine = false;
      }
      else if (num_chroms_used <= 1) {
	cout << "WARNING: fitting exponential + unconstrained affine (A * exp(-n*d) + C)" << endl;
	cout << "need >= 2 chroms to determine affine term from inter-chrom data" << endl << endl;
	use_inter_chrom_affine = false;
      }
    }
    int numbins = results_allchrom[0].size();
    vector <AlderResults> results_jackknife(num_chroms_used+1);
    for (int jc = 0; jc <= num_chroms_used; jc++) {
      results_jackknife[jc].fit_start_dis = fit_start_dis;
      results_jackknife[jc].jack_id = jack_ind_ids[jc];
      vector <double> &x = results_jackknife[jc].d_Morgans;
      vector <double> &y = results_jackknife[jc].weighted_LD_avg;
      vector <double> &count = results_jackknife[jc].bin_count;
      x = y = count = vector <double> (numbins);
      // set up the remove-one data
      vector <AffineData> affdats_minus_jc;
      for (int c = 0; c < num_chroms_used; c++)
	if (c != jc) {
	  for (int b = 1; b < numbins; b++) {
	    x[b] = b * binsize;
	    y[b] += results_allchrom[c][b].first;
	    count[b] += results_allchrom[c][b].second;
	  }
	  if (use_inter_chrom_affine)
	    affdats_minus_jc.push_back(affine_data_allchrom[c]);
	}
      for (int b = 1; b < numbins; b++) y[b] /= count[b];
      if (use_inter_chrom_affine) {
	x.push_back(INFINITY);
	pair <double, double> aff_tot_count = compute_inter_chrom_affine(affdats_minus_jc);
	y.push_back(aff_tot_count.first);
	count.push_back(aff_tot_count.second);
      }
    }    
    return results_jackknife;
  }

  vector <ExpFitALD> Alder::fit_results(const vector <AlderResults> &results_jackknife,
				       double fit_start_dis, double maxdis, int &fit_test_ind) {

    vector <ExpFitALD> fits_all_starts;
    if (fit_start_dis == INFINITY) {
      cout << "fit start = inf because of long-range LD correlation; not doing fitting" << endl;
      return fits_all_starts;
    }
    for (double mindis = fit_start_dis-0.002; mindis <= fit_start_dis+0.0021; mindis += 0.001) {
      if (fabs(mindis - fit_start_dis) < 1e-9) fit_test_ind = fits_all_starts.size();
      fits_all_starts.push_back(exp_fit_jackknife(results_jackknife, mindis, maxdis));
    }
    return fits_all_starts;
  }

  void Alder::count_alleles(const char *geno, int stride, int s, double &a, double &b) {
    a = b = 0;
    for (int i = 0; i < stride; i++) {
      int x = geno[s*stride+i];
      if (x != 9) {
	a += x;
	b += 2-x;
      }
    }
  }

  vector <double> Alder::compute_f2_jacks(const char *geno1, int stride1,
					  const char *geno2, int stride2) {
    // note: only use chromosomes specified at initialization!
    vector <double> f2_N_per_chrom(num_chroms_used), num_f2_per_chrom(num_chroms_used);
    for (int c = 0; c < num_chroms_used; c++)
      for (int s = chrom_start_inds[c]; s < chrom_start_inds[c+1]; s++) {
	if (snp_ignore[s]) continue;
	double a1, b1, a2, b2;
	count_alleles(geno1, stride1, s, a1, b1);
	count_alleles(geno2, stride2, s, a2, b2);
	if (a1+b1 <= 1 || a2+b2 <= 1) continue;
	double p1 = a1/(a1+b1);
	double N_bias1 = a1*b1/((a1+b1)*(a1+b1)*(a1+b1-1));
	double p2 = a2/(a2+b2);
	double N_bias2 = a2*b2/((a2+b2)*(a2+b2)*(a2+b2-1));
	f2_N_per_chrom[c] += ((p1-p2)*(p1-p2) - N_bias1 - N_bias2);
	num_f2_per_chrom[c]++;
      }
    vector <double> f2_jacks(num_chroms_used+1);
    for (int jc = 0; jc <= num_chroms_used; jc++) {
      double f2_N = 0, num_f2 = 0;
      for (int c = 0; c < num_chroms_used; c++)
	if (c != jc) {
	  f2_N += f2_N_per_chrom[c];
	  num_f2 += num_f2_per_chrom[c];
	}
      f2_jacks[jc] = f2_N / num_f2;
    }
    return f2_jacks;
  }

  // public functions

  Alder::Alder(char *_mixed_geno, int _num_mixed_indivs, const string &_mixed_pop_name,
	     const vector <char *> &_ref_genos, const vector <int> &_num_ref_indivs,
	     const vector <string> &_ref_pop_names, const vector < pair <int, double> > &snp_locs,
	     Timer &_timer) :
    mixed_geno(_mixed_geno), num_mixed_indivs(_num_mixed_indivs), mixed_pop_name(_mixed_pop_name),
    ref_genos(_ref_genos), num_ref_indivs(_num_ref_indivs), ref_pop_names(_ref_pop_names),
    timer(_timer) {
    
    int S = snp_locs.size();
    snp_chrom_ind_squash = snp_num_missing = snp_sum = snp_sum2 = snp_bin = vector <int> (S);
    snp_ignore = vector <char> (S);
    snp_pos = vector <double> (S);

    // set up snp tables
    // set up chromosome number remap: squash to 0, 1, 2, ...
    int mixed_geno_ind = 0;
    for (int s = 0; s < S; s++) {
      if (s > 0 && snp_locs[s] < snp_locs[s-1]) fatalx("snps must be sorted (error at %d)\n", s);
      if (s == 0 || snp_locs[s].first != snp_locs[s-1].first) { // new chromosome
	chrom_start_inds.push_back(s);
	jack_ind_ids.push_back(to_str(snp_locs[s].first));
      }	
      snp_chrom_ind_squash[s] = jack_ind_ids.size()-1;
      snp_pos[s] = snp_locs[s].second;
      for (int i = 0; i < num_mixed_indivs; i++) {
	int gtype = mixed_geno[mixed_geno_ind++];
	if (gtype == 9)
	  snp_num_missing[s]++;
	else {
	  snp_sum[s] += gtype;
	  snp_sum2[s] += gtype*gtype;
	}
      }
    }
    chrom_start_inds.push_back(S);
    num_chroms_used = jack_ind_ids.size();
    cout << "number of chromosomes (with data) used: " << num_chroms_used << endl;
    jack_ind_ids.push_back("none");

    if (num_chroms_used == 0) fatalx("no chromosomes with data\n");
    use_jackknife = num_chroms_used > 1;
  }
  
  int Alder::get_num_chroms_used(void) {
    return num_chroms_used;
  }

  vector <double> Alder::find_ld_corr_stops(double binsize0, bool use_early_exit, double mindis) {
    cout << "     *** Determining extent of correlated LD between test and ref pops ***" << endl;
    cout << endl;

    if (mindis != AlderParams::MINDIS_NOT_SET) {
      cout << "WARNING: skipping LD correlation computation: 'mindis' param is set" << endl;
      cout << "decay curves will be fit starting at d = " << mindis << " Morgans, as specified"
	   << endl << endl;
      return vector <double> (ref_pop_names.size(), mindis);
    }

    if (ref_genos.empty()) {
      cout << "WARNING: can't check LD corr; need reference genos (not just weights)" << endl;
      cout << "decay curves will be fit starting at the default min distance: "
	   << 100*AlderParams::DEFAULT_FIT_START << " cM" << endl << endl;
      return vector <double> (1, AlderParams::DEFAULT_FIT_START);
    }
    vector <double> ld_corr_stops(ref_pop_names.size(), INFINITY);

    double binsize = binsize0;
    const double binsize1 = 0.002;
    do {
      const int min_bin = 1;
      const double max_ld_dist = 0.02;
      const int numbins = max_ld_dist / binsize;
      for (int r = 0; r < (int) ref_pop_names.size(); r++) {
	cout << "Checking LD correlation of test pop " << mixed_pop_name << " with ref pop "
	     << ref_pop_names[r] << endl;
	cout << "  binsize: " << (100*binsize) << " cM" << endl;
	cout << "  (distances are rounded down to bins; bin starting at 0 is skipped)" << endl
	     << endl;
      
	bool compute_corr_data = true, compute_polyache_data = !use_early_exit;
	if (!compute_polyache_data)
	  printf("%6s%20s", "d (cM)", "LD corr (scaled)");
	else
	  printf("%6s%20s%12s%12s", "d (cM)", "unbiased LD corr", "RMS(D) test",
		 "RMS(D) ref");
	printf("%12s\n", "bin count");

	int num_significance_failures = 0;
	for (int b = min_bin; b < numbins; b++) {
	  CorrJack corr_data(num_chroms_used, jack_ind_ids),
	    test_data(num_chroms_used, jack_ind_ids), ref_data(num_chroms_used, jack_ind_ids);
	  compute_ld_corr_terms(r, b*binsize, (b+1)*binsize, corr_data, test_data, ref_data,
				use_early_exit, compute_corr_data, compute_polyache_data);
	  pair <double, double> corr_mean_std = compute_polyache_data ?
	    corr_data.jackknife_cos_polyache_denom(test_data, ref_data) :
	    corr_data.jackknife_corr();

	  // output bin, corr
	  printf("%6.3f%20s", 100*b*binsize, format_mean_std(corr_mean_std).c_str());
	  if (compute_polyache_data) {
	    // output unbiased LD in test, ref
	    printf("%12.5f%12.5f",
		   test_data.tot_sum_x2() / test_data.tot_count(),
		   ref_data.tot_sum_x2() / ref_data.tot_count());
	  }
	  // output count (same full count for all), newline
	  printf("%12d", (int) corr_data.tot_count());
	
	  if (erfc(corr_mean_std.first/corr_mean_std.second/sqrt(2.0)) > LD_COS_SIGNIF_THRESH ||
	      !(corr_mean_std.first > 0)) { // the latter check takes care of NAN and unknown std
	    num_significance_failures++;
	    if (use_early_exit)
	      cout << "   losing significance (" << num_significance_failures << ")";
	    if (num_significance_failures == LIM_SIGNIFICANCE_FAILURES) {
	      if (ld_corr_stops[r] == INFINITY || b*binsize > ld_corr_stops[r])
		ld_corr_stops[r] = b*binsize; // store as answer
	      if (use_early_exit) {
		printf("\n");
		cout << "lost significance; computing bias-corrected LD corr polyache" << endl;
		compute_corr_data = false; compute_polyache_data = true;
		compute_ld_corr_terms(r, b*binsize, (b+1)*binsize, corr_data, test_data, ref_data,
				      use_early_exit, compute_corr_data, compute_polyache_data);
		pair <double, double> cos_polyache_mean_std =
		  corr_data.jackknife_cos_polyache_denom(test_data, ref_data);

		printf("%6.3f%20s", 100*b*binsize, format_mean_std(cos_polyache_mean_std).c_str());
		printf("   <-- approx bias-corrected LD corr\n");
		/*
		  printf("SNP pairs used: %d (LD prod), %d (LD test), %d (LD ref)\n",
		  corr_data.tot_count(), test_data.tot_count(), ref_data.tot_count());
		*/
		break;
	      }
	    }
	  }
	  printf("\n");
	}
	cout << endl;
      }
      binsize *= 2;
    } while (binsize < binsize1);

    printhline();
    cout << "                 *** Summary of LD correlation results ***" << endl << endl;
    cout << "Decay curves will be fit starting at the folliowing min distances (cM):" << endl;
    cout << "  (to override, specify the 'mindis' parameter)" << endl << endl;
    for (int r = 0; r < (int) ref_pop_names.size(); r++) {
      printf("%20s%8.3f", ref_pop_names[r].c_str(), 100*ld_corr_stops[r]);
      if (ld_corr_stops[r] < AlderParams::DEFAULT_FIT_START) {
	printf(" --> replacing with default min %.3f", 100*AlderParams::DEFAULT_FIT_START);
	ld_corr_stops[r] = AlderParams::DEFAULT_FIT_START;
      }
      if (ld_corr_stops[r] == INFINITY)
	printf(" --> long-range LD corr with %s!", mixed_pop_name.c_str());
      printf("\n");
    }
    cout << endl << "==> Time to calculate LD correlation: " << timer.update_time()
	 << endl << endl;
    return ld_corr_stops;
  }

  // computes weighted LD on each chromosome; returns vector of results from jackknife runs
  // fit data is stored in fits_all_starts
  // ref_inds tells which ref pops are being used, just for the purpose of output
  //   (might be empty in the case of external weights)
  vector <AlderResults> Alder::run(
      int num_refs, const vector <int> &ref_inds, const vector <double> &weights, double maxdis,
      double binsize, int mincount, bool use_naive_algo, double fit_start_dis,
      vector <ExpFitALD> &fits_all_starts, int &fit_test_ind) {

    cout << "   *** Computing " << num_refs << "-ref weighted LD with weights";
    if (ref_inds.empty())
      cout << " from file";
    else
      for (int i = 0; i < (int) ref_inds.size(); i++)
	cout << " " << ref_pop_names[ref_inds[i]];
    cout << " ***" << endl << endl;

    if (use_naive_algo)
      printf("using naive pairwise algorithm (on all snps with >= %d non-missing values)\n",
	     mincount);
    if (num_refs == 1 && num_mixed_indivs < 4)
      fatalx("need at least 4 indivs in test pop to compute single-ref LD (polyache)\n");

    int numbins = maxdis / binsize;

    vector < vector < pair <double, double> > > results_allchrom(num_chroms_used);
    vector <AffineData> affine_data_allchrom(num_chroms_used);

    // set up snp_bin and snp_ignore tables
    // polyache requires no missing
    int maxmissing = num_refs == 2 ? num_mixed_indivs - mincount : 0;
    for (int s = 0; s < (int) snp_pos.size(); s++) {
      snp_bin[s] = (snp_pos[s] - snp_pos[chrom_start_inds[snp_chrom_ind_squash[s]]]) / binsize;
      snp_ignore[s] = snp_num_missing[s] > maxmissing || isnan(weights[s]);
    }

    // run computation on each chromosome
    cout << "analyzing chrom";
#pragma omp parallel for schedule(static,1)
    for (int c = 0; c < num_chroms_used; c++) {
#pragma omp critical
      cout << " " << jack_ind_ids[c] << flush;
      AffineData affine_data;
      vector < pair <double, double> > chrom_results = use_naive_algo ?
	run_chrom_naive(c, num_refs, weights, binsize, numbins, mincount) :
	run_chrom(c, num_refs, weights, binsize, numbins, mincount, affine_data);
      affine_data_allchrom[c] = affine_data;
      results_allchrom[c] = chrom_results;
    }
    cout << endl;

    vector <AlderResults> results_jackknife = make_results(results_allchrom, affine_data_allchrom,
							   binsize, use_naive_algo, fit_start_dis);

    cout << endl << "==> Time to run alder: " << timer.update_time() << endl << endl;
    
    fits_all_starts = fit_results(results_jackknife, fit_start_dis, maxdis, fit_test_ind);

    return results_jackknife;
  }
  
  double Alder::compute_mult_hyp_corr(const vector <bool> &use_ref) {
    char jobu = 'N', jobvt = 'N';
    int m = 0, n = 0;
    for (int r = 0; r < (int) use_ref.size(); r++) m += use_ref[r];
    double *A = new double[m * snp_pos.size()]; // max size needed: may throw out some snps
    double S[m];
    int ldu = 1, ldvt = 1;
    int lwork = 5*m;
    double work[lwork];
    int info;
    
    // create allele freq array
    for (int s = 0; s < (int) snp_pos.size(); s++) {
      int arr_pos = n*m;
      bool snp_good = true;
      for (int r = 0; r < (int) use_ref.size(); r++)
	if (use_ref[r]) {
	  double geno_mean = compute_geno_mean(s, ref_genos[r], num_ref_indivs[r]);
	  if (isnan(geno_mean))
	    snp_good = false;
	  else
	    A[arr_pos++] = geno_mean;
	}
      if (snp_good) {
	// mean-center the data for this snp
	double geno_mean_sum = 0.0;;
	for (int i = 0; i < m; i++) geno_mean_sum += A[n*m+i];
	for (int i = 0; i < m; i++) A[n*m+i] -= geno_mean_sum / m;
	n++;
      }
    }
    cout << "Calculating number of effective populations using data from "
	 << n << " snps..." << endl;
    if (n < m) fatalx("need at least as many snps as ref pops\n");
    
    dgesvd_(&jobu, &jobvt, &m, &n, A, &m, S, NULL, &ldu, NULL, &ldvt, work, &lwork, &info);
    delete[] A;

    cout << "Singular values of mean-centered allele frequency matrix:" << endl;
    for (int i = 0; i < m; i++)
      cout << "\t" << i+1 << ":\t" << S[i] << endl;

    double tot_variance = 0.0, cum_variance = 0.0;
    int num_eff_pops = 0;
    for (int i = 0; i < m-1; i++) tot_variance += sq(S[i]);
    for (; num_eff_pops < m-1 && cum_variance < PCA_VARIANCE_THRESH * tot_variance; num_eff_pops++)
      cum_variance += sq(S[num_eff_pops]);
    num_eff_pops++; // extra singular vector from mean-centering
    cout << "Effective number of ref pops (accounting for >90% of variance): "
	 << num_eff_pops << endl;
    double mult_hyp_corr = num_eff_pops * (num_eff_pops-1) / 2.0;
    cout << "Number of \"independent\" tests for multiple hypothesis correction: "
	 << mult_hyp_corr << endl;
    cout << endl;
    return mult_hyp_corr;
  }

  vector <double> Alder::compute_one_ref_f2_jacks(int ref_ind) {
    return compute_f2_jacks(mixed_geno, num_mixed_indivs,
			    ref_genos[ref_ind], num_ref_indivs[ref_ind]);
  }
}
