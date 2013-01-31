#ifndef MULTFITALD_HPP
#define MULTFITALD_HPP

#undef max
#include <sstream>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include "Alder.hpp"
#include "asa047.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
using std::map;
using std::string;
using std::vector;
using ALD::AlderResults;
using std::cout;
using std::make_pair;
using std::stringstream;

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
	void fit_curves();
	int golden_section_time(double, double, double, double, int);
	int golden_section_amp(double, double, double, double, string, int);
	void print_fitted();
	void add_mix();
};


#endif
