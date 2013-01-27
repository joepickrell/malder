#include <string>
#include <vector>
#include <utility>
#include <numeric>
#include <cmath>
#include "CorrJack.hpp"
#include "Jackknife.hpp"

namespace ALD {
  
  using std::string;
  using std::vector;
  using std::pair;

  void Corr::add_term(double x, double y) {
    if (std::isnan(x) || std::isnan(y)) return;
    count += 1.0; sum_x += x; sum_y += y;
    sum_xy += x*y; sum_x2 += x*x; sum_y2 += y*y;
  }
  
  // augments both x2 and y2 (only one will actually be used)
  void Corr::add_unbiased_sq_term(double sq_term) {
    if (std::isnan(sq_term)) return;
    count += 1.0; sum_x2 += sq_term; sum_y2 += sq_term;
  }


  // central = true for usual corr, false for non-zeroed
  pair <double, double> CorrJack::jackknife_corr(bool central) {
    vector <double> corr(C+1);
    for (int jc = 0; jc <= C; jc++) { // jc == C for all (no jackknife)
      double jcount = 0, jsum_x = 0, jsum_y = 0, jsum_xy = 0, jsum_x2 = 0, jsum_y2 = 0;
      for (int c = 0; c < C; c++)
	if (c != jc) {
	  jcount += data[c].count; jsum_x += data[c].sum_x; jsum_y += data[c].sum_y;
	  jsum_xy += data[c].sum_xy; jsum_x2 += data[c].sum_x2; jsum_y2 += data[c].sum_y2;
	}
      if (central)
	corr[jc] = (jcount*jsum_xy - jsum_x*jsum_y) /
	  sqrt((jcount*jsum_x2 - jsum_x*jsum_x) * (jcount*jsum_y2 - jsum_y*jsum_y));
      else
	corr[jc] = jsum_xy / sqrt(jsum_x2 * jsum_y2);
    }
    return Jackknife::mean_std(corr);
  }
  
  CorrJack::CorrJack(int _C, const vector <string> &_ids) : C(_C), jack_ind_ids(_ids) {
    data = vector <Corr> (C);
  }

  pair <double, double> CorrJack::jackknife_corr(void) {
    return jackknife_corr(true);
  }

  pair <double, double> CorrJack::jackknife_cos(void) {
    return jackknife_corr(false);
  }

  pair <double, double> CorrJack::jackknife_cos_polyache_denom(const CorrJack &test_data,
							       const CorrJack &ref_data) {
    vector <double> cos(C+1);
    for (int jc = 0; jc <= C; jc++) { // jc == C for all (no jackknife)
      double jcount_xy = 0, jsum_xy = 0, jcount_x2 = 0, jsum_x2 = 0, jcount_y2 = 0, jsum_y2 = 0;
      for (int c = 0; c < C; c++)
	if (c != jc) {
	  jcount_xy += data[c].count; jsum_xy += data[c].sum_xy;
	  jcount_x2 += test_data.data[c].count; jsum_x2 += test_data.data[c].sum_x2;
	  jcount_y2 += ref_data.data[c].count; jsum_y2 += ref_data.data[c].sum_y2;
	}
      cos[jc] = jsum_xy/jcount_xy / sqrt(jsum_x2/jcount_x2 * jsum_y2/jcount_y2);
    }
    return Jackknife::mean_std(cos);
  }

  pair <double, double> CorrJack::jackknife_x2_avg(void) {
    vector <double> x2_avg(C+1);
    for (int jc = 0; jc <= C; jc++) { // jc == C for all (no jackknife)
      double jcount = 0, jsum_x2 = 0;
      for (int c = 0; c < C; c++)
	if (c != jc) {
	  jcount += data[c].count;
	  jsum_x2 += data[c].sum_x2;
	}
      x2_avg[jc] = jsum_x2 / jcount;
    }
    return Jackknife::mean_std(x2_avg);
  }
  
  double CorrJack::tot_count(void) {
    double ans = 0.0;
    for (int c = 0; c < C; c++)
      ans += data[c].count;
    return ans;
  }

  double CorrJack::tot_sum_x2(void) {
    double ans = 0.0;
    for (int c = 0; c < C; c++)
      ans += data[c].sum_x2;
    return ans;
  }

}
