#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib> 
#include <string> 
#include <algorithm>			    					// std::sort
#include <math.h>
#include <iomanip>
#include <ctime>
#include <time.h>
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/gamma_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/math/distributions/beta.hpp>
#include <boost/random.hpp>
#include <boost/math/special_functions/digamma.hpp> 				//needed for computing the modified Dirichlet bound
#include <boost/math/special_functions/gamma.hpp>
#include <boost/foreach.hpp>
#include "omp.h"
#include "get_k.h" 
#include "infiles.h" 

using namespace std;
using namespace boost::math::policies;
using namespace boost::math;

	// Define a specific policy:
typedef policy< 
	boost::math::policies::overflow_error<ignore_error> 
	> my_policy;


//****************************************************************************************************************************************
	vector< vector <double> > x; 						// double x[n][p]:	pdfvalues
	vector< int > y; 							// int y[n]:		number of mappings per read
	vector< vector <int> > c; 						// int c[n][p]:		original index of mapping trascripts per read.
	vector< vector <int> > c_new; 						// int c_new[n][p]:	relabelled index of mapping components per read.
	vector< int > components;
	vector< int > components_raw;
	vector < int > z;							//allocation variables
	vector < double > n_sum(K);						//n_sum(k) = alpha + sum_{i=1}^{n}I(z_i==k)
	double * theta;
	vector < double > dir_pars(K);
	vector < double > prior_pars(K);
	vector < double > probs(K);						// conditional membership probabilities
//****************************************************************************************************************************************
	int current_max = 20;							// initial max number of components
//****************************************************************************************************************************************
	double* vb;								//vb estimates which are read from file
	double	lmbeta_one_k; 							//lmbeta(rep(1,K))
	double* new_alpha;
	double* new_beta;
	int n_obs;
	double unifRand();
	double logl();
	double logl_theta();
	int map_index;
	vector < double > theta_map(K);
	vector < double > theta_mean(K);
							// number of components and number of observations
//****************************************************************************************************************************************
	double * dir(vector<double>&, int );				//function to generate Dirichlet random numbers
	double lmbeta(vector<double>&, int );					//multinomial Beta function
//****************************************************************************************************************************************


int main()
{


	ifstream inFile, trFile, log_l_values;
	int i, j, num, t, k, place = 0;
	double stheta, temp;
	string line, tid;
	srand(time(NULL)) ;
	clock_t start, end;

	//cout <<"main K = "<< K<<endl;
	inFile.open(prob_file.data(), ifstream::in);
	cout<< "Reading "<< prob_file.data()<<" file"<<endl;
	if (!inFile) {
    		cerr << "Unable to open file .prob file\n";
    		exit(1);   // call system to stop
	}
	
	i = -1;
	while (i<0){
		getline( inFile, line );
		cout <<line <<  endl;
		if (line[0] != '#'){
			i = 0;	
			inFile.seekg(place);	// this sets the ifstream to the previous line
		}
		place = inFile.tellg(); // this reads the current position of the ifstream
	}


	i = 0;
	inFile >> tid >> num;
	y.push_back(num);
	vector<double> row(num);
	vector<int> c_row(num);
	temp = 0.0;
	for ( j = 0; j < num; ++j){
		inFile >> c_row[j] >> row[j];	
		components.push_back(c_row[j]);
		components_raw.push_back(c_row[j]);
 
	}
	x.push_back(row);	
	c_new.push_back(c_row);	
	
	std::sort(components.begin(), components.end());
	inFile.ignore(10000000,'\n');
	i = i+1;		
	int paush;
	while( inFile.good() ){

		inFile >> tid >> num;
		y.push_back(num);
		vector<double> row(num);
		vector<int> c_row(num);
		for ( j = 0; j < num; ++j){
			inFile >> c_row[j] >> row[j];
		}
		x.push_back(row);
		c_new.push_back(c_row);
		inFile.ignore(10000000,'\n');
		i = i+1;
		if ( i % 1000000 == 0 )   {
		cout << "reading .prob file: line " << i<< endl;
		}
	}
	
	inFile.close();
	
	n_obs = i - 1;
	cout << "*******************************************"<<endl;
	cout << "Done: K = " << K << ", n = "<< n_obs <<"."<<endl;
	cout << "*******************************************"<<endl;
	z.resize(n_obs);
	//exit(1);
	// read the VB results and store them to pointer vb**************************************************************************
	//***************************************************************************************************************************



//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************

	int iter, max_iter;

	for ( k = 0; k < K; k++){
		prior_pars[k] = 1.0;
		theta_map[k] = 0.0;
		theta_mean[k] = 0.0;
	}

	ofstream outfile("gibbs_chib.log");
	int mcmc_runs = 10, iter_mcmc;
	for (iter_mcmc = 0;iter_mcmc<mcmc_runs;iter_mcmc++){
		ofstream log_l_values("log_values_l.log");
		//initialization
		theta = dir(prior_pars, K);
		double l_value_mean, l_value;
		//l_value = logl();
		map_index = 0;
		// gibbs sampling
		double denom, u, c_sum;
		int burn_in = 1000;
		max_iter = 11000;

		l_value = logl_theta();
		for (k = 0; k < K; k++){
			theta_map[k] += theta[k];
			theta_mean[k] += theta[k];
		}

		for ( iter = 0; iter < max_iter; iter++){
			for (k = 0;k<K;k++){
				n_sum[k] = prior_pars[k];
			}
		
			for (i = 0; i<n_obs; i++){
				denom = 0.0;
			
				for (k = 0; k < y[i]; k++){
					probs[k] = theta[c_new[i][k]]*x[i][k];
					denom += probs[k];
				}
				for (k = 0; k < y[i]; k++){
					probs[k] /= denom;
				}
				u = unifRand();
				k = 0;
				c_sum = probs[k];
				z[i] = c_new[i][k];
				while ( u>c_sum){
					k++;
					c_sum += probs[k];
					z[i] = c_new[i][k];
				}
				//if (z[i]> (K-1)){cout<<"da fuq: "<<z[i]<<endl;}
				n_sum[z[i]] += 1.0;	
			}

			u = logl_theta();
			log_l_values << setprecision(12)<< u<<endl;
			if (u>l_value){
				l_value = u;
				for (k = 0; k<K; k++){
					theta_map[k] = theta[k];
				}
			}
			if (iter > (burn_in-1)){
				for (k = 0; k < K; k++){
					theta_mean[k] += theta[k];
				}
			}
			delete[] theta;
			theta = dir(n_sum, K);
			if (iter % 1000 == 0) {
				cout<<iter<<endl;
			}
			
		}
		//exit(1);
		for (k = 0; k < K; k++){
			theta_mean[k] /= 1.0*(max_iter - burn_in);
		}
		l_value_mean = logl();


		cout <<"log_likehood(theta_map) = "<< setprecision(12)<<l_value<<endl;
		cout <<"log_likehood(theta_mean) = "<< setprecision(12)<<l_value_mean<<endl;
				//for (k = 0; k < K; k++){
				//	cout<<theta_map[k]<<endl;
				//}

		//extra sampling to estimate the log marginal likelihood (chib's method)
	
		log_l_values.close();
		double logm, log_pdf_theta_map = 0.0, log_pdf_theta_mean = 0.0, log_like, log_prior;
		max_iter = 10000;
		double logged_values[max_iter],logged_values_mean[max_iter];
		double max_logged_values = -999999990.0,max_logged_values_mean = -999999990.0, temp_mean;
		for ( iter = 0; iter < max_iter; iter++){
			for (k = 0;k<K;k++){
				n_sum[k] = prior_pars[k];
			}
			for (i = 0; i<n_obs; i++){
				denom = 0.0;
				for (k = 0; k < y[i]; k++){
					probs[k] = theta[c_new[i][k]]*x[i][k];
					denom += probs[k];
				}
				for (k = 0; k < y[i]; k++){
					probs[k] /= denom;
				}
				u = unifRand();
				k = 0;
				c_sum = probs[k];
				z[i] = c_new[i][k];
				while (u>c_sum){
					k++;
					c_sum += probs[k];
					z[i] = c_new[i][k];
				}
				n_sum[z[i]] += 1.0;	
			}

			temp = 0.0;
			temp_mean = 0.0;
			for (k = 0;k<K;k++){
				temp += (n_sum[k]-1.0)*log(theta_map[k]);
				temp_mean += (n_sum[k]-1.0)*log(theta_mean[k]);
			}
			temp -= lmbeta(n_sum, K);
			temp_mean -= lmbeta(n_sum, K);
			logged_values[iter] = temp;
			logged_values_mean[iter] = temp_mean;
			if (logged_values[iter]>max_logged_values){
				max_logged_values = logged_values[iter];
			}
			if (logged_values_mean[iter]>max_logged_values_mean){
				max_logged_values_mean = logged_values_mean[iter];
			}

			delete[] theta;
			theta = dir(n_sum, K);
			if (iter % 1000 == 0) {
				cout<<iter<<endl;
			}

	
		}
		cout<<"max_logged_values = "<<max_logged_values<<" (iterations = "<<iter<<")"<<endl;
		ofstream logvalues("log_values.txt");
		for (iter = 0; iter<max_iter; iter++){
			logvalues<<setprecision(12)<<logged_values[iter]<<endl;
			logged_values[iter] = logged_values[iter] - max_logged_values;
			logged_values_mean[iter] = logged_values_mean[iter] - max_logged_values_mean;
		}
		logvalues.close();
		log_pdf_theta_map =  0.0;
		log_pdf_theta_mean =  0.0;
		int und_num = 0, und_num_mean = 0;

		for (iter = 0; iter<max_iter; iter++){
			if (iter % 1000 == 0) {
				cout<< iter<<": iter, "<<"log = "<<logged_values[iter]<<", exp = "<<exp(logged_values[iter])<<endl;
			}
			u = exp(logged_values[iter]);
			if (u == 0.0){und_num++;}
			log_pdf_theta_map += u;
			u = exp(logged_values_mean[iter]);
			if (u == 0.0){und_num_mean++;}
			log_pdf_theta_mean += u;

		}
		cout << "number of underflows (map): "<< und_num <<endl;
		cout << "number of underflows (mean): "<< und_num_mean <<endl;
		log_pdf_theta_map = max_logged_values + log(log_pdf_theta_map) -log(0.0+max_iter);
		log_pdf_theta_mean = max_logged_values_mean + log(log_pdf_theta_mean) -log(0.0+max_iter);

		//cout <<"final = "<<log_pdf_theta_map<<endl;
		//cout <<"final = "<<log_pdf_theta_mean<<endl;
		log_prior = -lmbeta(prior_pars,K);
		logm = l_value + log_prior - log_pdf_theta_map;
		outfile<<setprecision(12)<<logm<<" ";

		cout << "************************************"<<endl;
		cout<<"Chib's log-marginal estimate (map) = "<<setprecision(12)<<logm<<endl;
		logm = l_value_mean + log_prior - log_pdf_theta_mean;
		outfile<<setprecision(12)<<logm<<endl;
		cout<<"Chib's log-marginal estimate (mean) = "<<setprecision(12)<<logm<<endl;
		cout << "************************************"<<endl;

	}

	outfile.close();


return(0);
}



/* 
*****************************************************************************************
	Function to compute the incomplete log-likelihood				*
											*
n_comps  : number of mixture components							*
n_reads  : number of observations							*
theta[]  : n_comps-dimensional array containing mixture weights 			*
*****************************************************************************************
*/

double logl(){
double temp, loglike = 0.0;
int i, k;
for (i = 0; i < n_obs; i++){ 
	temp = 0.0;
	//cout << "y[" << i << "] = " << y[i] << " :{";
	for (k = 0; k < y[i]; k++){
	//	j = c_new[i][k];	
	//	temp += theta[j]*x[i][k];
		temp += theta_mean[c_new[i][k]]*x[i][k];
	}
	//cout << "}" << endl;
	loglike += log(temp);
}

return loglike;
}



double logl_theta(){
double temp, loglike = 0.0;
int i, k;
for (i = 0; i < n_obs; i++){ 
	temp = 0.0;
	//cout << "y[" << i << "] = " << y[i] << " :{";
	for (k = 0; k < y[i]; k++){
	//	j = c_new[i][k];	
	//	temp += theta[j]*x[i][k];
		temp += theta[c_new[i][k]]*x[i][k];
	}
	//cout << "}" << endl;
	loglike += log(temp);
}

return loglike;
}



/* 
*****************************************************************************************
 	Function to simulate (x_1, x_2,...x_{K}) ~ Dirichlet(pars[1],...pars[K])	*
		(corresponds to "simulate3" function in R)				*			
*****************************************************************************************
*/
double * dir(vector< double >& pars, int n_comps){
	int i;
	double *theta_vec = new double[n_comps]; 
	double sum_theta = 0.0;
	int seed = rand();
	boost::mt19937 rng(seed);
	for (i = 0; i< n_comps; ++i){
		boost::gamma_distribution<double> gd(pars[i]);
		boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma( rng, gd );
		theta_vec[i] = var_gamma();
		sum_theta += theta_vec[i];		
	}
	for (i = 0; i< n_comps; ++i){
	theta_vec[i]/= sum_theta;
	}
	return theta_vec;
}



/*
*****************************************************************************************
 	lmbeta function: 	sum(lgamma(c(...)))-lgamma(sum(c(...)))			*
*****************************************************************************************
*/

double lmbeta(vector< double >& pars, int n_comps){
	int i;
	double x = 0.0, s = 0.0;
	for (i = 0; i<n_comps; ++i){
		s += pars[i];
		x += boost::math::lgamma(pars[i]);
	}
	x -=  boost::math::lgamma(s);
	return x;
}



/*
*****************************************************************************************
 	Function to simulate U(0,1) r.v's						*
*****************************************************************************************								
*/

double unifRand()
{
    return rand() / double(RAND_MAX);
}


/*
*****************************************************************************************
 	Function to simulate Bernoulli {-1, 1} r.v's					*
*****************************************************************************************								
*/

int bernRand()
{
	double u;
	int b = 1;
	u = rand() / double(RAND_MAX);
	if ( u<0.5 ) b = -1;
    return b;
}










