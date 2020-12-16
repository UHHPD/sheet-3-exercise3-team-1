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

int main() {
    using namespace std;
    vector<int> zaehler(11);

    ifstream fin("datensumme.txt");
    ofstream fout("hist.txt");
    ofstream fout2("histpoi.txt");
    
    int n_i;
    int out_i;
    double mu = 3.11538;
    double P;
    double N = 234;
    double P_Expectation;

    for(int i = 0 ; i < 234 ; ++i) {
        fin >> n_i;
        zaehler[n_i] += 1;
    }
    
    for(unsigned int k = 0 ; k < zaehler.size() ; ++k)
    {
	P = poisson(mu, k);
	P_Expectation = N * P; 
    	std::cout << k << ":" << zaehler[k] << std::endl;
    	fout << k << " " << zaehler[k] << std::endl;
	fout2.precision(16);
        fout2 << k << " " << zaehler[k] << " " << P_Expectation << std::endl;
	
    }

    
    fin.close();
}
