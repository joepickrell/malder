#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "ExpFit.hpp"
#include "Jackknife.hpp"
#include "MiscUtils.hpp"
#include "ExpFitALD.hpp"

namespace ALD {

  using std::string;
  using std::vector;
  using std::pair;
  using std::make_pair;
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::max;
  using std::min;

  string ExpFitALD::mean_std_str(const char *varname) const {
    pair <double, double> mean_std = Jackknife::mean_std(get_var(varname));
    int digits = strcmp(varname, "decay") == 0 ? 2 : 8;
    char format[100], line[100];
    sprintf(format, "%%.%df", digits);
    sprintf(line, format, mean_std.first);
    string ret(line);
    if (use_jackknife) {
      sprintf(format, " +/- %%.%df", digits);
      sprintf(line, format, mean_std.second);
      ret += string(line);
    }
    return ret;
  }

  string ExpFitALD::make_fit_line(const char *descrip, int digits, pair <double, double> mean_std,
				  bool print_zscore) const {
    char format[100], line[100];
    sprintf(format, "d>%%.2f%%12s:%%12.%df", digits);
    sprintf(line, format, 100*mindis, descrip, mean_std.first);
    string ret(line);
    if (use_jackknife) {
      sprintf(format, " +/- %%-10.%df", digits);
      sprintf(line, format, mean_std.second);
      ret += string(line);
      if (print_zscore) {
	double zscore = mean_std.first / mean_std.second;
	sprintf(line, "   z = %.2f", zscore);
	ret += string(line);
	if (zscore >= 2) ret += " *";
      }
    }
    return ret;
  }

  string ExpFitALD::make_fit_line(const char *varname, int digits, bool print_zscore,
				  bool print_jackknife_data) const {
    if (print_jackknife_data) {
      const vector <double> &var = get_var(varname);
      printf("d>%.2f%12s jackknife reps\t", 100*mindis, varname);
      for (int i = 0; i < (int) var.size(); i++)
	printf("\t%.8f", var[i]);
      cout << endl;
    }
    return make_fit_line(varname, digits, Jackknife::mean_std(get_var(varname)), print_zscore);
  }

  pair <double, double> ExpFitALD::compute_diff_zscore_percent(const ExpFitALD &ref_fit,
							       const char *varname) const {
    pair <double, double> diff_mean_std =
      Jackknife::diff_mean_std(get_var(varname), ref_fit.get_var(varname));
    double mean_of_means = (get_var(varname).back() + ref_fit.get_var(varname).back()) / 2;
    double diff_zscore = diff_mean_std.first / diff_mean_std.second; 
    double diff_percent = 100.0 * diff_mean_std.first / mean_of_means;
    return make_pair(diff_zscore, diff_percent);
  }

  bool ExpFitALD::is_diff_significant(pair <double, double> diff_zscore_percent) {
    return fabs(diff_zscore_percent.second) < TEST_DECAY_DIFF_THRESH;
  }

  // public functions

  ExpFitALD::ExpFitALD(int _num_chroms_used_plus_1, double _mindis, double _maxdis,
		       bool _use_jackknife)
    : mindis(_mindis), maxdis(_maxdis), use_jackknife(_use_jackknife) {
    gen = amp_tot = amp_exp = amp_aff = vector <double> (_num_chroms_used_plus_1);
    use_inter_chrom_affine = false;
    num_fit_failures = 0;
  }
  
  void ExpFitALD::do_fit(int jc, const double *x, const double *y, int bins) {
    if (x[bins-1] == INFINITY) {
      use_inter_chrom_affine = true;
      amp_aff[jc] = y[bins-1];
    }
    bool success =
      ExpFit::fit_decay(x, y, bins, GEN_MIN, GEN_MAX, &gen[jc], &amp_exp[jc], &amp_aff[jc],
			use_inter_chrom_affine);
    amp_tot[jc] = amp_exp[jc] + amp_aff[jc]/2;
    if (!success && jc != (int) gen.size()-1) { // failure on jackknife trial
      num_fit_failures++;
      gen[jc] = amp_exp[jc] = amp_aff[jc] = amp_tot[jc] = INFINITY;
    }
  }
  
  void ExpFitALD::print_fit_header(void) const {
    printf("---- fit on data from %.2f to %.2f cM (%s) ----\n", 100*mindis, 100*maxdis,
	   use_inter_chrom_affine ? "using inter-chrom affine term" : "unconstrained affine term");
    if (num_fit_failures)
      printf("WARNING: %d of %d jackknife reps failed to fit; unable to estimate std\n",
	     num_fit_failures, (int) gen.size()-1);
  }

  vector <double> ExpFitALD::get_var(const char *varname) const {
    if (strcmp(varname, "decay") == 0) return gen;
    else if (strcmp(varname, "amp_tot") == 0) return amp_tot;
    else if (strcmp(varname, "amp_exp") == 0) return amp_exp;
    else if (strcmp(varname, "amp_aff") == 0) return amp_aff;
    else {
      fprintf(stderr, "internal error: unknown fit varname %s\n", varname);
      exit(1);
    }
  }
  
  void ExpFitALD::print_fit_diff(const ExpFitALD &ref_fit, const char *varname, int digits,
				 const string &pop_name, const string &ref_name) const {
    cout << make_fit_line((varname + string(" diff")).c_str(), digits,
			  Jackknife::diff_mean_std(get_var(varname), ref_fit.get_var(varname)),
			  false);
    cout << pop_name << " - " << ref_name << endl;
  }

  bool ExpFitALD::print_test_fit_diff(const ExpFitALD &ref_fit, const char *varname,
				      const string &pop_name, const string &ref_name) const {
    pair <double, double> diff_zscore_percent = compute_diff_zscore_percent(ref_fit, varname);
    printf("   %-50s %5.2f   (%3.0f%%)\n", (pop_name.substr(0, 18) + " - " +
					    ref_name.substr(0, 18) + " z-score:").c_str(),
	   diff_zscore_percent.first, diff_zscore_percent.second);
    return is_diff_significant(diff_zscore_percent);
  }
  
  void ExpFitALD::print_fit(bool print_jackknife_fits) const {
    print_fit_header();
    cout << make_fit_line("decay", 2, true, print_jackknife_fits) << endl;
    cout << make_fit_line("amp_tot", 8, false, print_jackknife_fits) << endl;
    cout << make_fit_line("amp_exp", 8, true, print_jackknife_fits) << endl;
    cout << make_fit_line("amp_aff", 8, false, print_jackknife_fits) << endl;
    cout << endl;
  }

  double ExpFitALD::zscore(const char *varname) const {
    return Jackknife::zscore(get_var(varname));
  }
  
  double ExpFitALD::p_value(double mult_hyp_corr = 1.0) const {
    return erfc(min(zscore("decay"), zscore("amp_exp")) / sqrt(2.0)) * mult_hyp_corr;
  }

  bool ExpFitALD::is_significant(double mult_hyp_corr = 1.0) const {
    return p_value(mult_hyp_corr) < 0.05;
  }

  bool ExpFitALD::test_and_print_oneref_curve(void) const {
    bool test_success = is_significant();
    printf("   1-ref decay z-score:   %5.2f\n", zscore("decay"));
    printf("   1-ref amp_exp z-score: %5.2f\n", zscore("amp_exp"));
    printf("                                  %s\n",
	   (test_success ? "YES: curve is significant" : "NO: curve is not significant"));
    cout << endl;
    return test_success;
  }

  bool ExpFitALD::run_admixture_test(const ExpFitALD &tworef, const ExpFitALD &oneref_1,
				     const ExpFitALD &oneref_2, const string &mixed_pop_name,
				     const string &ref_pop_name_1, const string &ref_pop_name_2,
				     bool run_pretest, double mult_hyp_corr) {

    printhline();
    cout << "                    *** Admixture test summary ***" << endl << endl;

    cout << "Weighted LD curves are fit starting at " << 100*tworef.mindis << " cM"
	 << endl << endl;
    
    bool test_success = true;
    if (run_pretest) {
      cout << "Pre-test: Does " << mixed_pop_name << " have a 1-ref weighted LD curve with "
	   << ref_pop_name_1 << "?" << endl;
      test_success &= oneref_1.test_and_print_oneref_curve();
      cout << "Pre-test: Does " << mixed_pop_name << " have a 1-ref weighted LD curve with "
	   << ref_pop_name_2 << "?" << endl;
      test_success &= oneref_2.test_and_print_oneref_curve();
    }
    
    bool tworef_success = tworef.is_significant(mult_hyp_corr);
    test_success &= tworef_success;
    cout << "Does " << mixed_pop_name << " have a 2-ref weighted LD curve with "
	 << ref_pop_name_1 << " and " << ref_pop_name_2 << "?" << endl;
    printf("   2-ref decay z-score:   %5.2f\n", tworef.zscore("decay"));
    printf("   2-ref amp_exp z-score: %5.2f\n", tworef.zscore("amp_exp"));
    printf("                                  %s\n",
	   (tworef_success ? "YES: curve is significant" : "NO: curve is not significant"));
    cout << endl;

    cout << "Do 2-ref and 1-ref curves have consistent decay rates?" << endl;
    bool consistency_success = true;
    consistency_success &= oneref_1.print_test_fit_diff(tworef, "decay",
							"1-ref " + ref_pop_name_1, "2-ref");
    consistency_success &= oneref_2.print_test_fit_diff(tworef, "decay",
							"1-ref " + ref_pop_name_2, "2-ref");
    consistency_success &= oneref_2.print_test_fit_diff(oneref_1, "decay", 
							"1-ref " + ref_pop_name_2,
							"1-ref " + ref_pop_name_1);
    printf("                                  %s\n", consistency_success ?
	   "YES: decay rates are consistent" : "WARNING: decay rates are inconsistent");
    cout << endl;

    //test_success &= consistency_success; <-- this is now a warning

    double test_zscore = min(tworef.zscore("decay"), tworef.zscore("amp_exp"));
    double p_value = tworef.p_value(mult_hyp_corr);
    printf("Test %s (z=%.2f, p=%.2g) for %s with {%s, %s} weights\n",
	   (test_success ? "SUCCEEDS" : "FAILS"), test_zscore, p_value,
	   mixed_pop_name.c_str(), ref_pop_name_1.c_str(), ref_pop_name_2.c_str());
    if (mult_hyp_corr != 1.0)
      cout << "  note: p-value is multiplied by " << mult_hyp_corr
	   << " for multiple-hypothesis correction" << endl;
    cout << endl;

    //    succ/fail, C, A, B, z-score min of 2-ref decay and amp_exp, min z-score A-ref, min z-score B-ref, max percent discordance, decay +/- std, amp +/- std for all three curves
    printf("DATA:\t%s%s\t%.2g\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.0f%%\t%s\t%s\t%s\t%s\t%s\t%s\n",
	   (test_success ? "success" : "failure"),
	   (consistency_success ? "" : " (warning: decay rates inconsistent)"),
	   p_value, mixed_pop_name.c_str(), ref_pop_name_1.c_str(), ref_pop_name_2.c_str(),
	   test_zscore,
	   min(oneref_1.zscore("decay"), oneref_1.zscore("amp_exp")),
	   min(oneref_2.zscore("decay"), oneref_2.zscore("amp_exp")),
	   max(fabs(oneref_1.compute_diff_zscore_percent(tworef, "decay").second),
	       max(fabs(oneref_2.compute_diff_zscore_percent(tworef, "decay").second),
		   fabs(oneref_2.compute_diff_zscore_percent(oneref_1, "decay").second))),
	   tworef.mean_std_str("decay").c_str(), tworef.mean_std_str("amp_exp").c_str(),
	   oneref_1.mean_std_str("decay").c_str(), oneref_1.mean_std_str("amp_exp").c_str(),
	   oneref_2.mean_std_str("decay").c_str(), oneref_2.mean_std_str("amp_exp").c_str());
    cout << endl;
  
    return test_success;
  }
  
  /*
    assuming A''=A, amp is 2ab^3*f2(A,B)^2     (in reality, smaller)
    assuming A'=A, f2(C,A') = f2(C,A) = b^2*f2(A,B)     (in reality, bigger)
    => amp / f2(C,A')^2 = 2a/(1-a)     (in reality, smaller => estimate is lower bound)
    solve 2a/(1-a) with binary search; jackknife
  */
  pair <double, double> ExpFitALD::mix_frac_bound(const vector <double> &f2_jacks) {
    if (f2_jacks.size() != amp_tot.size()) {
      cerr << "internal error: f2_jacks and amp_tot do not have the same size" << endl;
      exit(1);
    }
    vector <double> alpha_jacks(f2_jacks.size());
    for (int jc = 0; jc < (int) f2_jacks.size(); jc++) {
      double y = amp_tot[jc] / pow(f2_jacks[jc], 2.0);
      double lo = 0, hi = 1;
      while (hi-lo > 1e-6) {
	double alpha = (lo+hi)/2;
	if (2*alpha/(1-alpha) < y)
	  lo = alpha;
	else
	  hi = alpha;
      }
      alpha_jacks[jc] = (lo+hi)/2;
    }
    return Jackknife::mean_std(alpha_jacks);
  }

  void ExpFitALD::print_data_header(void) {
    printf("DATA:\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	   "test status", "p-value", "test pop", "ref A", "ref B",
	   "2-ref z-score", "1-ref z-score A", "1-ref z-score B", "max decay diff %",
	   "2-ref decay", "2-ref amp_exp",
	   "1-ref decay A", "1-ref amp_exp A",
	   "1-ref decay B", "1-ref amp_exp B");
  }
}
