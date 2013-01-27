#ifndef ALDER_HPP
#define ALDER_HPP

#include <string>
#include <vector>
#include <utility>
#include <cmath>

#include <fftw3.h>

#include "Timer.hpp"
#include "CorrJack.hpp"
#include "ExpFitALD.hpp"

namespace ALD {

  extern "C" int dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s,
			 double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork,
			 int *info);

  using std::string;
  using std::vector;
  using std::pair;
  inline bool isnan(double x)
  {
   return (x != x);
  }
  struct AlderResults {
    vector <double> d_Morgans, weighted_LD_avg, bin_count;
    string jack_id;
    double fit_start_dis;
  };

  class Alder {
  
    struct AffineData {
      double count;
      double ws, ss, s2;
      vector <double> wg, gg, gs;
      vector < vector <double> > gigj;
      int num_refs;
      AffineData(int n, int _num_refs);
    };

    // constants for determining LD correlation extent
    static const int LIM_SIGNIFICANCE_FAILURES;
    static const double LD_COS_SIGNIF_THRESH;
    static const double PCA_VARIANCE_THRESH;

    static const bool SUBTRACT_THEN_BIN = false; // only for naive pairwise algorithm

    const char *mixed_geno;
    const int num_mixed_indivs;
    const string &mixed_pop_name;
    const vector <char *> &ref_genos;
    const vector <int> &num_ref_indivs;
    const vector <string> &ref_pop_names;
    bool use_jackknife;
    Timer &timer;

    vector <int> snp_num_missing, snp_sum, snp_sum2, snp_bin;
    vector <char> snp_ignore;
    vector <double> snp_pos;
    vector <int> snp_chrom_ind_squash;

    int num_chroms_used;
    vector <string> jack_ind_ids;
    vector <int> chrom_start_inds;

    string format_mean_std(pair <double, double> mean_std);
    double compute_geno_mean(int s, const char *geno, int stride);
    double compute_ld(int s1, int s2, const char *geno, int stride);
    double compute_ld(int s1, int s2);
    double compute_polyache_central_moment11sq(int s1, int s2, const char *geno, int stride);
    double compute_polyache_central_moment11sq(int s1, int s2);
    bool x2_suff_accurate(pair <double, double> x2_mean_std);

    // stores polyache data in:
    // - test_data.count, test_data.sum_x2 and test_data.sum_y2 (same)
    // - ref_data.count, ref_data.sum_x2 and ref_data.sum_y2 (same)
    void compute_ld_corr_terms(int ref_ind, double bin_min, double bin_max,
			       CorrJack &corr_data, CorrJack &test_data, CorrJack &ref_data,
			       bool use_early_exit, bool compute_corr_data,
			       bool compute_polyache_data);
    double compute_polyache(int s1, int s2, double pAx, double pAy);
    void convolve_accum(fftw_plan *plans, int Nby2, fftw_complex *z_accum,
			fftw_complex *z1, fftw_complex *z2, double scale=1.0);
    void self_convolve_accum(fftw_plan *plan, int Nby2, fftw_complex *z_accum,
			     fftw_complex *z1, double scale=1.0);
    // returns binned pairs: (weighted LD, count of pairs in bin)
    // also, affine_data contains info for computing affine term
    vector < pair <double, double> > run_chrom(int chrom, int num_refs,
					       const vector <double> &weights, double binsize,
					       int numbins, int mincount, AffineData &affine_data);
    vector < pair <double, double> > run_chrom_naive(int chrom, int num_refs,
						     const vector <double> &weights,
						     double binsize, int numbins, int mincount);
    pair <double, double> compute_inter_chrom_affine(const vector <AffineData> &affdats);
    void check_affine_amp(int num_refs, const vector <double> &weights);
    ExpFitALD exp_fit_jackknife(const vector <AlderResults> &results_jackknife,
				double mindis, double maxdis);
    vector <AlderResults> make_results(
        const vector < vector < pair <double, double> > > &results_allchrom,
        const vector <AffineData> &affine_data_allchrom, double binsize, bool use_naive_algo,
	double fit_start_dis);
    vector <ExpFitALD> fit_results(const vector <AlderResults> &results_jackknife,
				   double fit_start_dis, double maxdis, int &fit_test_ind);
    void count_alleles(const char *geno, int stride, int s, double &a, double &b);
    vector <double> compute_f2_jacks(const char *geno1, int stride1,
				     const char *geno2, int stride2);

  public:
    Alder(char *_mixed_geno, int _num_mixed_indivs, const string &_mixed_pop_name,
	 const vector <char *> &_ref_genos, const vector <int> &_num_ref_indivs,
	 const vector <string> &_ref_pop_names, const vector < pair <int, double> > &snp_locs,
	 Timer &_timer);  
    int get_num_chroms_used(void);
    vector <double> find_ld_corr_stops(double binsize, bool use_early_exit, double mindis);
    // computes weighted LD on each chromosome; returns vector of results from jackknife runs
    // fit data is stored in fits_all_starts
    // ref_inds tells which ref pops are being used, just for the purpose of output
    //   (might be empty in the case of external weights)
    vector <AlderResults> run(
	int num_refs, const vector <int> &ref_inds, const vector <double> &weights, double maxdis,
	double binsize, int mincount, bool use_naive_algo, double fit_start_dis,
	vector <ExpFitALD> &fits_all_starts, int &fit_test_ind);
    double compute_mult_hyp_corr(const vector <bool> &refs_to_use);
    vector <double> compute_one_ref_f2_jacks(int ref_ind);
  };  
}

#endif
