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

void MultFitALD::fit_curves(){
	double initss =  ss();
	//cout << initss << "\n";
	int nit= 0;
	bool done = false;
	while (!done){
		for (int i = 0; i < nmix; i++){
			//cout << "here1 "<< i << " "<< times[i]<< "\n";
			double lt = log(times[i]);
			golden_section_time(-20, lt, 20, 1e-10, i);
			//cout << times[i] << "\n";
		}
		for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
			string ps = it->first;

			vector<double> amps = expamps[ps];
			for (int i = 0; i < nmix; i++){
				//cout << ps << " "<<i << " "<< expamps[ps][i] << "\n";
				double la = log(amps[i]);
				golden_section_amp(-20, log(amps[i]), 10, 1e-10, ps, i);
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

void MultFitALD::print_fitted(){
	//for (int i = 0; i < (int) times.size(); i++){
	//	cout << i << " "<< times[i] <<  " "<< times.size() << "\n";
	//}
	stringstream ss;
	ss << "RESULT_" << nmix <<  "\trefpops";
	for (int i = 0; i <  nmix; i++){
		ss << "\tamp"<< i <<"\ttime" << i;
	}
	cout << ss.str() << "\n";
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
		string ps = it->first;
		cout << "RESULT_" << nmix << "\t"<< ps;
		for (int i = 0; i <  nmix; i++){
			cout << "\t" << expamps[ps][i] <<"\t" << times[i];
		}
		cout << "\n";
	}
}
void MultFitALD::add_mix(){
	nmix++;
	for (map<string, vector<AlderResults> >::iterator it = curves->begin(); it != curves->end(); it++){
		string ps = it->first;
		expamps[ps].push_back(1e-7);
	}
	times.push_back(50.0);
	fit_curves();
}
