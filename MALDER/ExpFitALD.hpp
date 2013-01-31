#ifndef EXPFITALD_HPP
#define EXPFITALD_HPP

#include <string>
#include <vector>
#include <utility>
#include <map>
#include "ExpFit.hpp"
#include "Jackknife.hpp"


namespace ALD {

  using std::string;
  using std::vector;
  using std::pair;
  using std::map;

  class ExpFitALD {
    vector <double> gen, amp_tot, amp_exp, amp_aff;
    double mindis, maxdis;
    bool use_jackknife;

    bool use_inter_chrom_affine;
    int num_fit_failures;

    string mean_std_str(const char *varname) const;
    string make_fit_line(const char *descrip, int digits, pair <double, double> mean_std,
			 bool print_zscore) const;
    string make_fit_line(const char *varname, int digits, bool print_zscore,
			 bool print_jackknife_fits) const;
    pair <double, double> compute_diff_zscore_percent(const ExpFitALD &ref_fit,
						      const char *varname) const;
    static bool is_diff_significant(pair <double, double> diff_zscore_percent);

  public:
    static const int GEN_MIN = 2;
    static const int GEN_MAX = 500;
    static const int TEST_DECAY_DIFF_THRESH = 25;

    ExpFitALD(int _num_chroms_used_plus_1, double _mindis, double _maxdis, bool _use_jackknife);
    void do_fit(int jc, const double *x, const double *y, int bins);
    void print_fit_header(void) const;
    vector <double> get_var(const char *varname) const;
    void print_fit_diff(const ExpFitALD &ref_fit, const char *varname, int digits,
			const string &pop_name, const string &ref_name) const;
    bool print_test_fit_diff(const ExpFitALD &ref_fit, const char *varname,
			     const string &pop_name, const string &ref_name) const;
    void print_fit(bool print_jackknife_fits) const;
    double zscore(const char *varname) const;
    double p_value(double mult_hyp_corr) const;
    bool is_significant(double mult_hyp_corr) const;
    bool test_and_print_oneref_curve(void) const;
    static bool run_admixture_test(const ExpFitALD &tworef, const ExpFitALD &oneref_1,
				   const ExpFitALD &oneref_2, const string &mixed_pop_name,
				   const string &ref_pop_name_1, const string &ref_pop_name_2,
				   bool run_pretest, double mult_hyp_corr);
    pair <double, double> mix_frac_bound(const vector <double> &f2_jacks); // f2(C,A) jack reps
    static void print_data_header(void);


  };
}

#endif
