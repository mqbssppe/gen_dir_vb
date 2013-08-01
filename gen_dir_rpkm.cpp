#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib> 
#include <string> 
#include <algorithm>			    					// std::sort
#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdio>
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
#include "transposeFiles.h"

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
	vector< int > empty;							// empty components
//****************************************************************************************************************************************
	int current_max = 20;							// initial max number of components
//****************************************************************************************************************************************
	double* vb;								//vb estimates which are read from file
	double	lmbeta_one_k; 							//lmbeta(rep(1,K))
	double* new_alpha;
	double* new_beta;
	int n_obs;							// number of components and number of observations
//****************************************************************************************************************************************
	double * g_d_rpkm(vector<double>&, vector<double>&, int );				//function to generate Generalized Dirichlet random numbers

	double lmbeta(double [], int );						//multinomial Beta function
//****************************************************************************************************************************************


	vector<int> random_sample;
	int n_size;
	double* vb_means;							//the means of the standard VB
	double *theta_rpkm;

	vector <int> tr_length(K);


	vector <double> alpha_temp1(K);
	vector <double> beta_temp1(K); 
	vector <double> a_b_ab_boost_lgamma1(K);
	vector <double>  coefficients1(K);
	
	vector< vector <double> > g_temp2(n_threads,vector<double>(K,0.0));
	vector< vector <double> > sum_theta_k2(n_threads,vector<double>(K,0.0));

	vector< vector <double> > g_temp1(n_threads,vector<double>(K,0.0));
	vector< vector <double> > sum_theta_k1(n_threads,vector<double>(K,0.0));






int main()
{


	ifstream trFile, gdFile;
	int i, j, num, t, k, place = 0;
	double stheta, temp;
	string line, tid, gd_full_name;
	srand(time(NULL)) ;
	clock_t start, end;

	trFile.open(tr_info_file.data(), ifstream::in);
	cout<< "Reading "<< tr_info_file.data()<<" file"<<endl;	
	if (!trFile) {
    		cerr << "Unable to open .tr file\n";
    		exit(1);   // call system to stop
	}

	i = -1;
	while (i<0){
		getline( trFile, line );
		cout <<line <<  endl;
		if (line[0] != '#'){
			i = 0;	
			trFile.seekg(place);	// this sets the ifstream to the previous line
		}
		place = trFile.tellg(); // this reads the current position of the ifstream
	}

	while( trFile.good() ){

		trFile >> tid >> tid >> tr_length[i] >> tr_length[i];
		trFile.ignore(10000000,'\n');
		if ( i % 500 == 0 )   {
		cout << "reading .tr file: line " << i<< ", length = "<< tr_length[i]<< endl;
		}
		i = i+1;
	}

	
	trFile.close();

	gd_full_name = out_file_extension.data();
	gd_full_name += ".gen_dir_parameters";
	gdFile.open(gd_full_name.data(), ifstream::in);
	cout<< "Reading "<< gd_full_name.data()<<" file"<<endl;	
	if (!gdFile) {
    		cerr << "Unable to open gd parameters file\n";
    		exit(1);   // call system to stop
	}
	i = 0;
	while( gdFile.good() ){

		gdFile >> alpha_temp1[i] >> beta_temp1[i];
		trFile.ignore(10000000,'\n');
		if ( i % 500 == 0 )   {
		cout << "reading .gd file: line " << i<< ", alpha = "<<setprecision(10)<< alpha_temp1[i] << ", beta = "<<  beta_temp1[i]<< endl;
		}
		i = i+1;
	}

	
	gdFile.close();
	
	ofstream rpkm;
	string rpkmfile;
	rpkmfile = out_file_extension.data();
	rpkmfile += ".rpkm";
	cout<< "writing transcript logged rpkm values to "<< rpkmfile << " file... "<<endl;

	rpkm.open((rpkmfile+"TMP").c_str());
	rpkm <<"# M " << K - 1 <<endl;
	rpkm <<"# N " << 1000 <<endl;
	for (j = 0;j<1000;j++){
		theta_rpkm = g_d_rpkm(alpha_temp1, beta_temp1, K);
		for (i = 1; i<K; ++i){
			rpkm << setprecision(9)<< theta_rpkm[i] <<" ";}
		rpkm <<endl;		
		delete[] theta_rpkm;
	}
	rpkm.close();
	
	transposeFiles(vector<string>(1,rpkmfile+"TMP"), rpkmfile, true, "# L\n");
	//remove((rpkmfile+"TMP").c_str());

	cout<<"done!"<<endl;


	//exit(1);
	// read the VB results and store them to pointer vb**************************************************************************
	//***************************************************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************




return(0);
}



double * g_d_rpkm( vector<double>& a_pars, vector<double>& b_pars, int n_comps){
	int i, j ;
	double *theta_vec = new double[n_comps]; 
	double sum_theta, g1, g2; // g1 and g2 are the two gamma rv's in order to simulate a beta rv
	int seed = rand();
	boost::mt19937 rng(seed);
	i = 0;
	boost::gamma_distribution<double> gd1(a_pars[i]), gd2(b_pars[i]);
	boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma1( rng, gd1 ), var_gamma2( rng, gd2 );
	g1 = var_gamma1(); // g1 ~ Gamma(a_pars[], 1)
	g2 = var_gamma2(); // g2 ~ Gamma(b_pars[], 1)	
	theta_vec[i] = g1/(g1 + g2); // theta[] ~ Beta(a_pars[], b_pars[])
	sum_theta = theta_vec[i];
	for (i = 1; i< n_comps-1; ++i){
		boost::gamma_distribution<double> gd1(a_pars[i]), gd2(b_pars[i]);
		boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma1( rng, gd1 ), var_gamma2( rng, gd2 );
		g1 = var_gamma1(); // g1 ~ Gamma(a_pars[], 1)
		g2 = var_gamma2(); // g2 ~ Gamma(b_pars[], 1)	
		theta_vec[i] = g1*(1.0-sum_theta)/(g1 + g2); // theta[] ~ Beta(a_pars[], b_pars[])
		sum_theta += theta_vec[i];
	}
	theta_vec[n_comps-1] = 1.0-sum_theta;
	double log_theta_0; 
	double con;
	log_theta_0 = log(1.0 - theta_vec[0]);
	con = 9.0*log(10.0);
	for (i = 1;i< n_comps; ++i ){
		theta_vec[i] = con + log(theta_vec[i]) - log_theta_0 - log(0.0 + tr_length[i-1]);
		
	}

	return theta_vec;
}


/*
*****************************************************************************************
 	lmbeta function: 	sum(lgamma(c(...)))-lgamma(sum(c(...)))			*
*****************************************************************************************
*/

double lmbeta(double pars[], int n_comps){
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







