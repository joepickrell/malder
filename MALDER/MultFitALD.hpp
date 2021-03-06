#ifndef MULTFITALD_HPP
#define MULTFITALD_HPP

#undef max
#include <sstream>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "Alder.hpp"
#include "Jackknife.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "nnls.h"
using std::map;
using std::string;
using std::vector;
using ALD::AlderResults;
using std::cout;
using std::make_pair;
using std::stringstream;
using std::pair;
using std::ofstream;
using std::cerr;

class MultFitALD{
public:
	MultFitALD(int, map<string, vector<AlderResults> >*);
	int nmix;
	double ss_epsilon, phi, resphi;
	map<string, vector<AlderResults> > *curves;
	vector<double> times;
	map<string, vector<double> > expamps;
	map<string, double> affine_amps;
	double ss();
	double ss(int);
	pair< vector<double>, map <string, vector<double> > > fit_curves();
	pair< vector<double>, map <string, vector<double> > > fit_curves_nnls();
	void fit_curves_jack(int);
	void fit_curves_jack_nnls(int);
	bool fit_amps_nnls();
	bool fit_amps_nnls_jack(int);
	int golden_section_time(double, double, double, double, int);
	int golden_section_time_nnls(double, double, double, double, int);
	int golden_section_time_nnls(double, double, double, double, int, int);
	int golden_section_amp(double, double, double, double, string, int);
	int golden_section_time(double, double, double, double, int, int);
	int golden_section_amp(double, double, double, double, string, int, int);
	bool print_fitted(pair< vector<double>, map <string, vector<double> > >* , pair< vector<vector<double> >, vector<map <string, vector<double> > > >*);
	void print_fitted(string);
	pair< vector<double>, map <string, vector<double> > > add_mix();
	pair< vector<vector<double> >, vector<map <string, vector<double> > > > jackknife();
	void print_curves(const char *);

	// trying GSL optimization
	double nelder_term;
	pair< vector<double>, map<string, vector<double> > > GSL_optim();
	pair< vector<vector<double> >, vector<map <string, vector<double> > > >  GSL_jack();
	void GSL_optim(int);
	//double get_timese(int, vector<double>, vector< vector<double> >);
};

struct GSL_params{
        MultFitALD *d;
        int which;
};
extern double GSL_ss(const gsl_vector *, void *GSL_params);
extern double GSL_ss_jack(const gsl_vector *, void *GSL_params);

#endif
