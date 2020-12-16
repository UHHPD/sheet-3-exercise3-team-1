#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>


double poisson(double mu, int k) {
    double e = exp(-mu);
    double mu_pow_k = pow(mu, k);
    double k_fac = tgamma(k+1);
    double P = (mu_pow_k * e)/k_fac;
    
    return P;
}

double prob(std::vector<int> daten, double mu){
    double likelihood = 1;
	
    for(unsigned int i=0; i<234; i++){
        likelihood = poisson(mu, daten[i]) * likelihood;
    }
    
    return likelihood;
}


int main() {
    using namespace std;
    vector<int> zaehler(11);
    vector<int> daten;

    ifstream fin("datensumme.txt");
    ofstream fout("likelihood.txt");
    ofstream fout2("nll.txt");
    ofstream fout3("deltanll.txt");
    
    int n_i;
    int out_i;
    double mu = 3.11538;
    double P;
    double N = 234;
    double P_Expectation;

    for(int i = 0 ; i < 234 ; ++i) {
        fin >> n_i;
        zaehler[n_i] += 1;
        daten.push_back(n_i);
    }
    
    double likelihood = prob(daten, mu);
    std::cout << likelihood << '\n' << std::scientific;
    
    
    double likelihood_scan;
    double nll;
    double deltanll;
    double deltanll_min;
    double mu_min;
    
    for(double mu_scan = 0 ; mu_scan <= 6; mu_scan += 0.1 ){
        likelihood_scan = prob(daten, mu_scan);
        nll = -2*log(mu_scan);
        deltanll = nll + 2*log(likelihood);
        fout << mu_scan << ' ' << likelihood_scan << std::endl;
        fout2 << mu_scan << ' ' << nll << std::endl;
        fout3 << mu_scan << ' ' << deltanll << std::endl;
        
        
        // find the best estimator of mu (i.e. the one which minimises the difference)
        if(mu_scan = 0){
            deltanll_min = deltanll;
        }
        
        if(deltanll < deltanll_min){
            deltanll_min = deltanll;
            mu_min = mu_scan;
        }
    }
    
    vector<double> interval;
    // determine the uncertainty on the estimated mu hat by scanning across the sample space
    for(double mu_scan = 0 ; mu_scan <= 6; mu_scan += 0.1 ){
        double ratio = -2 * (log(prob(daten, mu_scan)/log(prob(daten, mu))));
    	if(ratio < 1){
    	    interval.push_back(mu_scan);
    	}
    }
    
    double lower_interval = interval[0];
    double upper_interval = interval[static_cast<int>(interval.size())];
    
    
    cout << "Uncertainty from interval: " << lower_interval - upper_interval << std::endl;
    cout << "Error on the mean: "<< sqrt(mu/N) << std::endl;
    
    
    double max_L = 1;
    for(int i = 0; i < 234; i++)
    {
        max_L = max_L * poisson(daten[i], daten[i]);
    }
    
   cout << "Lambda: " << prob(daten, mu)/max_L << endl;
   cout << "z: " << (-2 * log(prob(daten, mu)/max_L) - 233)/sqrt(2*233) << endl;
   
   
    /*
    for(unsigned int k = 0 ; k < zaehler.size() ; ++k)
    {
	P = poisson(mu, k);
	P_Expectation = N * P; 
    	std::cout << k << ":" << zaehler[k] << std::endl;
    	fout << k << " " << zaehler[k] << std::endl;
	fout2.precision(16);
        fout2 << k << " " << zaehler[k] << " " << P_Expectation << std::endl;
	
    }
    */
    
    fin.close();
    fout.close();
    fout2.close();
    fout3.close();
}


/*
int main() {
    using namespace std;


    ifstream fin("datensumme.txt");
    int n_i;
    for(int i = 0 ; i < 234 ; ++i) {
        fin >> n_i;
    }
    fin.close();
}
*/
