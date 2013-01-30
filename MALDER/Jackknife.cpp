#include <vector>
#include <utility>
#include <cmath>
#include "Jackknife.hpp"

namespace Jackknife {

  using std::vector;
  using std::pair;
  using std::make_pair;

  double stddev(const vector <double> &x, int n) {
    for (int i = 0; i < n; i++) if (std::isnan(x[i])) return NAN;
    for (int i = 0; i < n; i++) if (std::isinf(x[i])) return INFINITY;
    double s = 0.0, s2 = 0.0;
    for (int i = 0; i < n; i++) {
      s += x[i];
      s2 += x[i]*x[i];
    }
    return sqrt((s2 - s*s/n) * (n-1) / n);
  }
//
  // assumes last element of x is the leave-no-data-out estimator (previous are jackknife reps)
  pair <double, double> mean_std(const vector <double> &x) {
    int n = x.size()-1; // number of jackknife reps
    return make_pair(x[n], stddev(x, n));
  }

  // assumes last element of x is the leave-no-data-out estimator (previous are jackknife reps)
  double zscore(const vector <double> &x) {
    pair <double, double> mu_sigma = mean_std(x);
    return mu_sigma.first / mu_sigma.second;
  }

  // assumes last element of x is the leave-no-data-out estimator (previous are jackknife reps)
  pair <double, double> diff_mean_std(vector <double> x, const vector <double> &x_ref) {
    int n = x.size()-1; // number of jackknife reps
    for (int i = 0; i <= n; i++) x[i] -= x_ref[i];
    return make_pair(x[n], stddev(x, n));
  }
}
