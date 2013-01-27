#include <string>
#include <vector>
#include <utility>
#include "Alder.hpp"

namespace ALD {

  using std::string;
  using std::vector;
  using std::pair;

  void printhline(void);
  void plot_ascii_curve(const AlderResults &results, double mindis);
  void output_curve_data(const AlderResults &results);
  void write_raw_output(const char *filename, bool print_raw_jackknife,
			const vector <AlderResults> &results_jackknife);
  string to_str(double x);

}
