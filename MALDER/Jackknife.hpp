#include <vector>
#include <utility>

namespace Jackknife {
  double stddev(const std::vector <double> &x, int n);
  std::pair <double, double> mean_std(const std::vector <double> &x);
  double zscore(const std::vector <double> &x);
  std::pair <double, double> diff_mean_std(std::vector <double> x,
					   const std::vector <double> &x_ref);
}
