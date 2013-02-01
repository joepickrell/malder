#include "MultFitALD.hpp"


MultFitALD::MultFitALD(int n, map<string, vector<AlderResults> >* ma){
	nmix = n;
	curves = ma;
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
		string ps = it->first;
		//cout << ps << "\n";
		vector<double> t;
		vector<double> amp;
		for (int i = 0; i < nmix; i++)	amp.push_back(0.001);
		expamps.insert(make_pair(ps, amp));
	}
	for (int i = 0; i < nmix; i++) times.push_back(10.0);
	ss_epsilon = 1e-6;
    phi = (1+sqrt(5))/2;
    resphi = 2-phi;
}

double MultFitALD::ss(){
	double toreturn = 0;
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
		string ps = it->first;
		AlderResults r = it->second.back();

		vector<double> amps = expamps[ps];
		double affine = r.weighted_LD_avg[r.bin_count.size()-1];
		for (int i =0; i < (int) r.bin_count.size()-1; i++){
			double d = r.d_Morgans[i];
			if (d < r.fit_start_dis) continue;
			double pred = affine;
			for (int j = 0; j < nmix; j++){
				pred += amps[j] * exp(-d* times[j]);
			}
			//if (ps == "French;Ju|'hoan_North") cout << ps << " "<< i << " "<< d << " "<< pred << " "<< r.weighted_LD_avg[i] <<" "<< amps[0]<<" "<< amps[1]<< " "<< times[0]<< " "<< times[1]<< " "<< affine << "\n";
			double diff = pred-r.weighted_LD_avg[i];
			double toadd = diff*diff;
			toreturn+= toadd;
		}
	}
	return toreturn*100000;
}


double MultFitALD::ss(int which){
	double toreturn = 0;
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
		string ps = it->first;
		AlderResults r = it->second[which];

		vector<double> amps = expamps[ps];
		double affine = r.weighted_LD_avg[r.bin_count.size()-1];
		for (int i =0; i < (int) r.bin_count.size()-1; i++){
			double d = r.d_Morgans[i];
			if (d < r.fit_start_dis) continue;
			double pred = affine;
			for (int j = 0; j < nmix; j++){
				pred += amps[j] * exp(-d* times[j]);
			}
			//if (ps == "French;Ju|'hoan_North") cout << ps << " "<< i << " "<< d << " "<< pred << " "<< r.weighted_LD_avg[i] <<" "<< amps[0]<<" "<< amps[1]<< " "<< times[0]<< " "<< times[1]<< " "<< affine << "\n";
			double diff = pred-r.weighted_LD_avg[i];
			double toadd = diff*diff;
			toreturn+= toadd;
		}
	}
	return toreturn*100000;
}

pair< vector<double>, map<string, vector<double> > > MultFitALD::fit_curves(){
	pair< vector<double>, map<string, vector<double> > > toreturn;
	double initss =  ss();
	//cout << initss << "\n";
	int nit= 0;
	bool done = false;
	while (!done){
		for (int i = 0; i < nmix; i++){
			//cout << "here1 "<< i << " "<< times[i]<< "\n";
			double lt = log(times[i]);
			golden_section_time(-30, lt, 20, 1e-10, i);
			//cout << times[i] << "\n";
		}
		for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
			string ps = it->first;

			vector<double> amps = expamps[ps];
			for (int i = 0; i < nmix; i++){
				golden_section_amp(-30, log(amps[i]), 10, 1e-10, ps, i);
				//cout << expamps[ps][i] << "\n";
			}
		}
		//cout << ss() << "\n";

		double currentss = ss();
		if (fabs(currentss - initss) < ss_epsilon) done = true;
		initss = currentss;
		//cout << ss() << "\n";
		//cout << "Iteration number: "<< nit << "\n";

		nit ++;
		//if (nit ==1) done = true;
	}
	for (vector<double>::iterator it = times.begin(); it != times.end(); it++) toreturn.first.push_back(*it);
	for (map<string, vector<double> >::iterator it = expamps.begin(); it != expamps.end(); it++) {
		vector<double> tmp;
		for (vector<double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) tmp.push_back(*it2);
		toreturn.second.insert(make_pair(it->first, tmp));
	}
	return toreturn;
}

void MultFitALD::fit_curves_jack(int which){
	double initss =  ss(which);
	//cout << initss << "\n";
	int nit= 0;
	bool done = false;
	while (!done){
		for (int i = 0; i < nmix; i++){
			//cout << "here1 "<< i << " "<< times[i]<< "\n";
			double lt = log(times[i]);
			golden_section_time(-30, lt, 20, 1e-10, i, which);
			//cout << times[i] << "\n";
		}
		for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
			string ps = it->first;

			vector<double> amps = expamps[ps];
			for (int i = 0; i < nmix; i++){
				golden_section_amp(-30, log(amps[i]), 10, 1e-10, ps, i, which);
				//cout << expamps[ps][i] << "\n";
			}
		}
		//cout << ss() << "\n";

		double currentss = ss(which);
		if (fabs(currentss - initss) < ss_epsilon) done = true;
		initss = currentss;
		//cout << ss() << "\n";
		//cout << "Iteration number: "<< nit << "\n";

		nit ++;
		//if (nit ==1) done = true;
	}

}

void MultFitALD::fit_curves_jack_nnls(int which){
	double initss =  ss(which);
	//cout << initss << "\n";
	int nit= 0;
	bool done = false;
	while (!done){
		for (int i = 0; i < nmix; i++){
			//cout << "here1 "<< i << " "<< times[i]<< "\n";
			double lt = log(times[i]);
			golden_section_time_nnls(-30, lt, 20, 1e-10, i, which);
			//cout << times[i] << "\n";
		}

		double currentss = ss(which);
		if (fabs(currentss - initss) < ss_epsilon) done = true;
		initss = currentss;

		nit ++;
	}

}
pair< vector<double>, map<string, vector<double> > > MultFitALD::fit_curves_nnls(){
	pair< vector<double>, map<string, vector<double> > > toreturn;
	double initss =  ss();
	//cout << initss << "\n";
	int nit= 0;
	bool done = false;
	while (!done){
		for (int i = 0; i < nmix; i++){
			double lt = log(times[i]);
			golden_section_time_nnls(-30, lt, 20, 1e-10, i);
			//cout << times[i] << "\n";
		}
		double currentss = ss();
		if (fabs(currentss - initss) < ss_epsilon) done = true;
		initss = currentss;

		nit ++;
	}
	for (vector<double>::iterator it = times.begin(); it != times.end(); it++) toreturn.first.push_back(*it);
	for (map<string, vector<double> >::iterator it = expamps.begin(); it != expamps.end(); it++) {
		vector<double> tmp;
		for (vector<double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) tmp.push_back(*it2);
		toreturn.second.insert(make_pair(it->first, tmp));
	}
	return toreturn;
}


int MultFitALD::golden_section_time(double min, double guess, double max, double tau, int which){
        double x;
        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_time = (min+max)/2;
                times[which] = exp(new_time);
                return 0;
        }

        times[which] = exp(x);
        double f_x = ss();

        times[which] = exp(guess);
        double f_guess = ss();
        //cout << which << " "<< exp(x) << " "<< exp(guess) << " "<< f_x << " "<< f_guess << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_time(guess, x, max, tau,  which);
                else return golden_section_time(min, x, guess, tau, which);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_time(min, guess, x, tau,  which);
                else return golden_section_time(x, guess, max, tau, which);
        }
}


int MultFitALD::golden_section_time_nnls(double min, double guess, double max, double tau, int which){
        double x;
        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_time = (min+max)/2;
                times[which] = exp(new_time);
                return 0;
        }

        times[which] = exp(x);
        fit_amps_nnls();
        double f_x = ss();

        times[which] = exp(guess);
        fit_amps_nnls();
        double f_guess = ss();
        //cout << which << " "<< exp(x) << " "<< exp(guess) << " "<< f_x << " "<< f_guess << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_time(guess, x, max, tau,  which);
                else return golden_section_time(min, x, guess, tau, which);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_time(min, guess, x, tau,  which);
                else return golden_section_time(x, guess, max, tau, which);
        }
}


int MultFitALD::golden_section_time_nnls(double min, double guess, double max, double tau, int which, int wjack){
        double x;
        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_time = (min+max)/2;
                times[which] = exp(new_time);
                return 0;
        }

        times[which] = exp(x);
        fit_amps_nnls_jack(wjack);
        double f_x = ss(wjack);

        times[which] = exp(guess);
        fit_amps_nnls_jack(wjack);
        double f_guess = ss(wjack);
        //cout << which << " "<< exp(x) << " "<< exp(guess) << " "<< f_x << " "<< f_guess << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_time(guess, x, max, tau,  which, wjack);
                else return golden_section_time(min, x, guess, tau, which, wjack);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_time(min, guess, x, tau,  which, wjack);
                else return golden_section_time(x, guess, max, tau, which, wjack);
        }
}
void MultFitALD::fit_amps_nnls(){
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
		string ps = it->first;
		AlderResults r = it->second.back();

		//vector<double> amps = expamps[ps];
		double affine = r.weighted_LD_avg[r.bin_count.size()-1];

		//
		// set up NNLS
		//

		//initialize the workspace
		int n = 0; // n is the number of weighted LD points to fit
		for (int i =0; i < (int) r.bin_count.size()-1; i++){
			double d = r.d_Morgans[i];
			if (d < r.fit_start_dis) continue;
			n++;
		}

		int p = nmix; // p is the number of times that are fixed

		//set up the workspace
		// solve Ax = b  with x >=0
		// A = exponentials (fixed, since times are fixed)
		// x = amplitudes (to be solved)
		// b = observed weighted LD
		NNLS_SOLVER nnls(n, p);
		//nnls._maxIter = 10;


		double * A = new double[n*p]; //A holds the exponentials exp(-t_i d)


		double * b  = new double[n];  // b contains the entries of the observed LD, subtracting off the affine term
		double * x = new double[p];   // x will be estimated, will contain the fitted amplitudes
		double rNorm;

		int index = 0;
		for (int i =0; i < (int) r.bin_count.size()-1; i++){
			double d = r.d_Morgans[i];
			double ld = r.weighted_LD_avg[i];
			if (d < r.fit_start_dis) continue;
			b[index] = ld - affine;
			for (int j = 0; j< nmix ; j++){
				double t = times[j];
				double e = exp(-t*d);
				A[j*n + index] = e;
			}
			index++;
		}

		//for (int i = 0; i < n; i++){
		//	cout << b[i] << " ";
		//	for(int j = 0; j < p ; j++){
		//		cout << A[j* n + i] << " ";
		//	}
		//	cout << "\n";
		//}
		//cout << "\n";
		//fit NNLS
		bool converged = nnls.solve(A, p, b, x, rNorm);


		// put back
		if (converged){
			for (int i =0 ; i < nmix; i++) {
				//cout << x[i]<< " "<< ss() << "\n";
				expamps[ps][i] = x[i];
			}
		}
		else{
			//cout << "Warning: fit did not converge\n";
			for (int i =0 ; i < nmix; i++) expamps[ps][i] = 0;
		}
	}
}


void MultFitALD::fit_amps_nnls_jack(int which){
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
		string ps = it->first;
		AlderResults r = it->second[which];

		//vector<double> amps = expamps[ps];
		double affine = r.weighted_LD_avg[r.bin_count.size()-1];

		//
		// set up NNLS
		//

		//initialize the workspace
		int n = 0; // n is the number of weighted LD points to fit
		for (int i =0; i < (int) r.bin_count.size()-1; i++){
			double d = r.d_Morgans[i];
			if (d < r.fit_start_dis) continue;
			n++;
		}

		int p = nmix; // p is the number of times that are fixed

		//set up the workspace
		// solve Ax = b  with x >=0
		// A = exponentials (fixed, since times are fixed)
		// x = amplitudes (to be solved)
		// b = observed weighted LD
		NNLS_SOLVER nnls(n, p);


		double * A = new double[n*p]; //A holds the exponentials exp(-t_i d)


		double * b  = new double[n];  // b contains the entries of the observed LD, subtracting off the affine term
		double * x = new double[p];   // x will be estimated, will contain the fitted amplitudes
		double rNorm;

		int index = 0;
		for (int i =0; i < (int) r.bin_count.size()-1; i++){
			double d = r.d_Morgans[i];
			double ld = r.weighted_LD_avg[i];
			if (d < r.fit_start_dis) continue;
			b[index] = ld - affine;
			for (int j = 0; j< nmix ; j++){
				double t = times[j];
				double e = exp(-t*d);
				A[j*n + index] = e;
			}
			index++;
		}

		//for (int i = 0; i < n; i++){
		//	cout << b[i] << " ";
		//	for(int j = 0; j < p ; j++){
		//		cout << A[j* n + i] << " ";
		//	}
		//	cout << "\n";
		//}
		//cout << "\n";
		//fit NNLS
		bool converged = nnls.solve(A, p, b, x, rNorm);


		// put back
		if (converged){
			for (int i =0 ; i < nmix; i++) {
				//cout << x[i]<< " "<< ss() << "\n";
				expamps[ps][i] = x[i];
			}
		}
		else{
			//cout << "Warning: fit did not converge\n";
			for (int i =0 ; i < nmix; i++) expamps[ps][i] = 0;
		}
	}
}

int MultFitALD::golden_section_amp(double min, double guess, double max, double tau, string pops, int which){
        double x;
        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_time = (min+max)/2;
                expamps[pops][which] = exp(new_time);
                return 0;
        }

        expamps[pops][which] = exp(x);
        double f_x = ss();

        expamps[pops][which] = exp(guess);
        double f_guess = ss();
       // cout << pops << " "<< which << " "<< x << " "<< guess << " "<< f_x << " "<< f_guess << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_amp(guess, x, max, tau, pops, which);
                else return golden_section_amp(min, x, guess, tau, pops, which);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_amp(min, guess, x, tau, pops, which);
                else return golden_section_amp(x, guess, max, tau, pops, which);
        }
}


int MultFitALD::golden_section_time(double min, double guess, double max, double tau, int which, int wjack){
        double x;
        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_time = (min+max)/2;
                times[which] = exp(new_time);
                return 0;
        }

        times[which] = exp(x);
        double f_x = ss(wjack);

        times[which] = exp(guess);
        double f_guess = ss(wjack);
        //cout << which << " "<< exp(x) << " "<< exp(guess) << " "<< f_x << " "<< f_guess << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_time(guess, x, max, tau,  which, wjack);
                else return golden_section_time(min, x, guess, tau, which, wjack);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_time(min, guess, x, tau,  which, wjack);
                else return golden_section_time(x, guess, max, tau, which, wjack);
        }
}

int MultFitALD::golden_section_amp(double min, double guess, double max, double tau, string pops, int which, int wjack){
        double x;
        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_time = (min+max)/2;
                expamps[pops][which] = exp(new_time);
                return 0;
        }

        expamps[pops][which] = exp(x);
        double f_x = ss(wjack);

        expamps[pops][which] = exp(guess);
        double f_guess = ss(wjack);
       // cout << pops << " "<< which << " "<< x << " "<< guess << " "<< f_x << " "<< f_guess << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_amp(guess, x, max, tau, pops, which, wjack);
                else return golden_section_amp(min, x, guess, tau, pops, which, wjack);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_amp(min, guess, x, tau, pops, which, wjack);
                else return golden_section_amp(x, guess, max, tau, pops, which, wjack);
        }
}

bool MultFitALD::print_fitted(pair< vector<double>, map <string, vector<double> > > * means,pair< vector<vector<double> >, vector<map<string, vector<double> > > >* reps ){
	stringstream ss;
	int njack = reps->first.size();
	ss << "RESULT_" << nmix <<  "\trefpops";
	for (int i = 0; i <  nmix; i++){
		ss << "\tamp"<< i <<"\ttime" << i;
	}
	cout << ss.str() << "\n";
	vector<double> max_ampz;
	vector<double> all_timez;
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
		string ps = it->first;
		cout << "RESULT_" << nmix << "\t"<< ps;
		for (int i = 0; i <  nmix; i++){
			if (it == curves->begin()) max_ampz.push_back(0.0);
			vector<double> timereps;
			vector<double> ampreps;
			for (int j = 0; j < njack; j++ ) {
				timereps.push_back(reps->first[j][i]);
				ampreps.push_back(reps->second[j][ps][i]);
			}
			timereps.push_back(means->first[i]);
			ampreps.push_back(means->second[ps][i]);
			pair<double, double> timedata = Jackknife::mean_std(timereps);
			pair<double, double> ampdata = Jackknife::mean_std(ampreps);
			double timez = timedata.first/timedata.second;
			if (it == curves->begin()) all_timez.push_back(timez);
			double ampz = ampdata.first/ampdata.second;
			if (max_ampz[i]< ampz && !std::isnan(ampz)) max_ampz[i] = ampz;
			cout << "\t" << expamps[ps][i] << " +/- "<< ampdata.second << " (Z="<< ampz <<")\t" << times[i] << " +/- " << timedata.second << " (Z="<< timez<< ")";
		}
		cout << "\n";
	}
	for (vector<double>::iterator it = max_ampz.begin(); it != max_ampz.end(); it++){
		//cout << *it << " amp\n";
		if (*it < 3) return true;
	}
	for (vector<double>::iterator it = all_timez.begin(); it != all_timez.end(); it++){
		if (*it < 3) return true;
	}
	return false;

}

void MultFitALD::print_fitted(string start){
	stringstream ss;
	ss << start << " RESULT_" << nmix <<  "\trefpops";
	for (int i = 0; i <  nmix; i++){
		ss << "\tamp"<< i <<"\ttime" << i;
	}
	cout << ss.str() << "\n";
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
		string ps = it->first;
		cout << start << " RESULT_" << nmix << "\t"<< ps;
		for (int i = 0; i <  nmix; i++){
			cout << "\t" << expamps[ps][i] << "\t" << times[i];
		}
		cout << "\n";
	}
}


pair< vector<double>, map <string, vector<double> > > MultFitALD::add_mix(){
	nmix++;
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
		string ps = it->first;
		expamps[ps].push_back(1e-7);
	}
	times.push_back(50.0);
	return(fit_curves_nnls());
}

pair< vector<vector<double> >, vector<map<string, vector<double> > > > MultFitALD::jackknife(){
	pair<vector<vector<double> >, vector<map<string, vector<double> > > > toreturn;
	vector<AlderResults> r = curves->begin()->second;
	int njack = r.size()-1;
	for (int i = 0; i < njack ; i++){
		//cout << "fitting "<< i << "\n"; cout.flush();
		fit_curves_jack_nnls(i);
		//cout << i << " " << njack <<" ";
		stringstream ss;
		ss << i;
		print_fitted(ss.str());
		vector<double> tmptimes;
		map<string,  vector<double> > tmpexp;
		for (vector<double>::iterator it = times.begin(); it != times.end(); it++) tmptimes.push_back(*it);
		for (map<string, vector<double> >::iterator it = expamps.begin(); it != expamps.end(); it++) {
			vector<double> tmp;
			for (vector<double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) tmp.push_back(*it2);
			tmpexp.insert(make_pair(it->first, tmp));
		}
		toreturn.first.push_back(tmptimes);
		toreturn.second.push_back(tmpexp);
	}
	return toreturn;
}

void MultFitALD::print_curves(const char* outname){
	ofstream outfile(outname);
	outfile << "d ";
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++) outfile << " "<< it->first;
	outfile << "\n";
	AlderResults r = curves->begin()->second.back();
	int nbin = (int) r.d_Morgans.size();
	double fixdist = r.fit_start_dis;
	for (int i = 0; i < nbin; i++){
		if (r.d_Morgans[i] < fixdist) outfile << "#";
		outfile << r.d_Morgans[i];
		for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
			outfile << " "<< it->second.back().weighted_LD_avg[i];
		}
		outfile << "\n";
	}
}
