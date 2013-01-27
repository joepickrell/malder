#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cstdio>
#include <cmath>
#include "Alder.hpp"
#include "MiscUtils.hpp"

namespace ALD {

  using std::cout;
  using std::endl;
  using std::string;
  using std::vector;
  using std::pair;
  using std::map;
  using std::min;
  using std::max;
  using std::ostringstream;

  void printhline(void) {
    printf("-------------------------------------------------------------------------------\n\n");
  }

  void plot_ascii_curve(const AlderResults &results, double mindis) {
    int numcols = 79; double xres = 0.1; // cM
    vector <double> y_tot(numcols), count_tot(numcols);
    for (int b = 1; b < (int) results.d_Morgans.size(); b++) {
      if (results.d_Morgans[b] < mindis) continue;
      double x = 100 * results.d_Morgans[b];
      int xcol = (x+1e-9) / xres;
      if (xcol >= numcols) break;
      y_tot[xcol] += results.weighted_LD_avg[b] * results.bin_count[b];
      count_tot[xcol] += results.bin_count[b];
    }
    double y_min = INFINITY, y_max = -INFINITY;
    for (int xcol = 0; xcol < numcols; xcol++)
      if (count_tot[xcol]) {
	double y = y_tot[xcol] / count_tot[xcol];
	y_min = min(y_min, y);
	y_max = max(y_max, y);
      }
    if (y_min == INFINITY) return; // no data

    // set plot boundaries
    int yrow_max = 30;
    double y_scale_max = max(fabs(y_min), fabs(y_max)) * yrow_max / (yrow_max-1);
    double y_pow = pow(10.0, floor(log10(y_scale_max)));
    double y_sigdigs = y_scale_max / y_pow;
    if (y_sigdigs < 2)
      y_sigdigs = ceil(10*y_sigdigs) / 10;
    else
      y_sigdigs = ceil(y_sigdigs);
    //cout << y_scale_max << " " << y_sigdigs * y_pow << endl;
    y_scale_max = y_sigdigs * y_pow;
    double yres = y_scale_max / yrow_max;
    int yrow_min = min(-1.0, floor(y_min / yres + 0.5));
    map <int, string> plot;

    // make axes
    for (int yrow = yrow_min; yrow <= yrow_max; yrow++) {
      plot[yrow] = string(numcols, yrow == 0 ? '-' : ' ');
      plot[yrow][0] = (yrow == 0 || yrow == yrow_max) ? '+' : '|';
    }
    
    // plot data points
    for (int xcol = 0; xcol < numcols; xcol++)
      if (count_tot[xcol]) {
	double y = y_tot[xcol] / count_tot[xcol];
	double yrow = floor(y / yres + 0.5);
	plot[yrow][xcol] = 'x';
      }

    // print scale
    int prev_int_cM = 0;
    for (int xcol = 0; xcol < numcols; xcol++)
      if ((int) (xcol*xres + 1e-9) > prev_int_cM) {
	prev_int_cM = (xcol*xres + 1e-9);
	plot[0][xcol] = '+';
	plot[-1][xcol] = '0' + prev_int_cM;
      }

    // print axis labels
    char buf[100]; sprintf(buf, "%.1e", y_scale_max); plot[yrow_max] = '+' + string(buf);
    plot[-1][numcols-2] = 'c';
    plot[-1][numcols-1] = 'M';
    
    for (int yrow = yrow_max; yrow >= yrow_min; yrow--)
      cout << plot[yrow] << endl;
    cout << endl;
  }

  void output_curve_data(const AlderResults &results) {
    for (int b = 1; b < (int) results.d_Morgans.size(); b++) {
      if (!(b < 20 || b > (int) results.d_Morgans.size()-10)) {
	if (b < 23)
	  printf("%10s\t%11s\t%13s\n", ".", ".", ".");
	  continue;
      }
      printf("%10.3f\t%11.8f\t%13Ld", 100*results.d_Morgans[b], results.weighted_LD_avg[b],
	     (long long) results.bin_count[b]);
      if (results.d_Morgans[b] == INFINITY)
	printf("   (inter-chrom snp pairs)");
      printf("\n");
    }
    printf("\n");
  }

  void write_raw_output(const char *filename, bool print_raw_jackknife,
			const vector <AlderResults> &results_jackknife) {
    for (int jc = print_raw_jackknife ? 0 : results_jackknife.size()-1;
	 jc < (int) results_jackknife.size(); jc++) {
      const AlderResults &results = results_jackknife[jc];
      FILE *fptr = fopen(jc < (int) results_jackknife.size()-1 ? // append suffix if chrom left out
			 (filename+("-"+results.jack_id)).c_str() : filename, "w");
      fprintf(fptr, "#%9s\t%11s\t%13s\n", "d (cM)", "weighted LD", "bin count");
      for (int b = 1; b < (int) results.d_Morgans.size(); b++)
	fprintf(fptr, "%c%9.3f\t%11.8f\t%13Ld\n",
		results.d_Morgans[b] < results.fit_start_dis - 1e-9 ? '#' : ' ',
		100*results.d_Morgans[b], results.weighted_LD_avg[b],
		(long long) results.bin_count[b]);
      if (results.d_Morgans.back() == INFINITY)
	fprintf(fptr, "# last row: affine term computed from pairs of SNPs on different chroms\n");
      fclose(fptr);
    }
  }

  string to_str(double x) { ostringstream oss; oss << x; return oss.str(); }
}
