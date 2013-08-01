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
	double * gen_dir(vector<double>&, vector<double>&, int );				//function to generate Generalized Dirichlet random numbers
	double * gen_dir_tp(vector<double>&, vector<double>&, vector<double>&, vector<double>&, int);

	double lmbeta(double [], int );						//multinomial Beta function
//****************************************************************************************************************************************


	vector<int> random_sample;
	int n_size;
	double* vb_means;							//the means of the standard VB

	vector <int> vb_cluster(K);

	vector <double> alpha_temp2(K);
	vector <double> beta_temp2(K); 
	vector <double> a_b_ab_boost_lgamma2(K);
	vector <double>  coefficients2(K);

	vector <double> alpha_temp1(K);
	vector <double> beta_temp1(K); 
	vector <double> a_b_ab_boost_lgamma1(K);
	vector <double>  coefficients1(K);

	vector< vector <double> > g_temp2(n_threads,vector<double>(K,0.0));
	vector< vector <double> > sum_theta_k2(n_threads,vector<double>(K,0.0));

	vector< vector <double> > g_temp1(n_threads,vector<double>(K,0.0));
	vector< vector <double> > sum_theta_k1(n_threads,vector<double>(K,0.0));

	double gd_bound_cluster(vector<double>& ,int , int );
	double diff_bound(vector<double>& ,vector<double>& ,int , int );
	//double gd_sp_cluster (int , int , int , double , int, const std::string& );
	double gd_sp_cluster (int , int , int , double , int);
	double theta1[5], theta2[5];

	const int dirichlet = 0;



int main()
{


	ifstream inFile, trFile, vbclusters;
	int i, j, num, t, k, place = 0;
	double stheta, temp;
	string line, tid;
	srand(time(NULL)) ;
	clock_t start, end;

	double index[K];
	for (i = 0; i < K; ++i) index[i] = 1.0;
	lmbeta_one_k = lmbeta(index, K);

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

	//exit(1);
	// read the VB results and store them to pointer vb**************************************************************************
	//***************************************************************************************************************************
	double vb_temp[K], vb_temp_means[K], new_alpha_temp[K], new_beta_temp[K-1];	
	inFile.open(vb_file, ifstream::in);
	if (!inFile) {
    		cerr << "Unable to open file vb_result.txt\n";
    		exit(1);   // call system to stop
	}
	//cout<< "** reading vb_results.m_alphas file **"<<endl;

	i = -1;
	while (i<0){
		getline( inFile, line );
		//cout <<line <<  endl;
		if (line[0] != '#'){
			i = 0;	
			inFile.seekg(place);	// this sets the ifstream to the previous line
		}
		place = inFile.tellg(); // this reads the current position of the ifstream
	}


	double not_needed;
	i = 0;
	while (inFile.good()){
		inFile >> vb_temp_means[i] >> vb_temp[i] >> not_needed;
		inFile.ignore(10000000,'\n');
		++i;
	}
	inFile.close();
	vb = vb_temp;
	vb_means = vb_temp_means;


//	use this in case you want to run gen_dir without dirichlet results

	for (i = 0; i < K-1;++i){
		new_beta_temp[i] = 0.0;
		for (j = i + 1; j < K; ++j)new_beta_temp[i] += vb_temp[j]; 
	}
	new_alpha = vb_temp;
	new_beta = new_beta_temp;



	int n_clust;
	if (K>500){
		n_clust = 200;
	}else{
		n_clust = K-1;		
	}

	if (dirichlet==1){
		for (k = 0; k<K-1; k++){
			vb_cluster[k] = 0; 
			//cout << vb_cluster[k]<<endl;
		}
		cout<<"************** dirichlet *******************"<<endl;
		gd_sp_cluster (5500, s_size, K, 5.0*K/1.0, 1);
		exit(1);
	}

	vbclusters.open(vb_clusters_file.data());
	if (!vbclusters) {
    		cerr << "Unable to open "<< vb_clusters_file.data() <<" file\n";
    		exit(1);   // call system to stop
	}
	for (k = 0; k<K-1; k++){
		vbclusters >> vb_cluster[k]; 
		//cout << vb_cluster[k]<<endl;
	}

	vbclusters.close();
	
	cout<<"************** Gen Dirichlet *******************"<<endl;

	// run gd_cluster from the beginning
	gd_sp_cluster (5500, s_size, K, K/2.0, n_clust);



	//cout << "Dirichlet SP algorithm " <<endl;
	//cout << dir_sp2 (30, 4, 1, 25000.0) << endl;

	exit(1);
	
	


//*******************************************************************************************
// run gd_cluster given the results of a previous gd_cluster run
//*******************************************************************************************
//	2 run
	for (k = 0; k<K-1; k++){
		vbclusters >> vb_cluster[k]; 
		//cout << vb_cluster[k]<<endl;
	}
	inFile.open("gen_dir_pars_clust_200.txt", ifstream::in);
	if (!inFile) {
    		cerr << "Unable to open file gen_dir_pars_clust_200.txt\n";
    		exit(1);   // call system to stop
	}



	vbclusters.close();





	exit(1);


//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************




return(0);
}



double * gen_dir( vector<double>& a_pars, vector<double>& b_pars, int n_comps){
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
	return theta_vec;
}




double * gen_dir_tp(vector<double>& a_pars1, vector<double>& b_pars1, vector<double>& a_pars2, vector<double>& b_pars2, int n_comps){
	int i, j ;
	double *theta_vec = new double[2*n_comps];
	//double theta_vec[2*n_comps]; 
	double sum_theta1, g11, g21; // g1 and g2 are the two gamma rv's in order to simulate a beta rv
	double sum_theta2, g12, g22; // g1 and g2 are the two gamma rv's in order to simulate a beta rv
	int seed = rand();
	boost::mt19937 rng(seed);
	i = 0;

	boost::gamma_distribution<double> gd1(a_pars1[i]), gd2(b_pars1[i]), gd3(a_pars2[i]), gd4(b_pars2[i]);
	boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma1( rng, gd1 ), var_gamma2( rng, gd2 ), var_gamma3( rng, gd3 ), var_gamma4( rng, gd4 );
	g11 = var_gamma1(); // g1 ~ Gamma(a_pars[], 1)
	g21 = var_gamma2(); // g2 ~ Gamma(b_pars[], 1)	

	g12 = var_gamma3(); // g1 ~ Gamma(a_pars[], 1)
	g22 = var_gamma4(); // g2 ~ Gamma(b_pars[], 1)	
	theta_vec[i] = g11/(g11 + g21); // theta[] ~ Beta(a_pars[], b_pars[])
	theta_vec[n_comps + i] = g12/(g12 + g22); // theta[] ~ Beta(a_pars[], b_pars[])
	sum_theta1 = theta_vec[i];
	sum_theta2 = theta_vec[n_comps + i];
	for (i = 1; i< n_comps-1; ++i){
		boost::gamma_distribution<double> gd1(a_pars1[i]), gd2(b_pars1[i]), gd3(a_pars2[i]), gd4(b_pars2[i]);
		boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma1( rng, gd1 ), var_gamma2( rng, gd2 ), var_gamma3( rng, gd3 ), var_gamma4( rng, gd4 );
		g11 = var_gamma1(); // g1 ~ Gamma(a_pars[], 1)
		g21 = var_gamma2(); // g2 ~ Gamma(b_pars[], 1)	
		g12 = var_gamma3(); // g1 ~ Gamma(a_pars[], 1)
		g22 = var_gamma4(); // g2 ~ Gamma(b_pars[], 1)	
		theta_vec[i] = g11*(1.0-sum_theta1)/(g11 + g21); // theta[] ~ Beta(a_pars[], b_pars[])
		theta_vec[n_comps + i] = g12*(1.0-sum_theta2)/(g12 + g22); // theta[] ~ Beta(a_pars[], b_pars[])
		sum_theta1 += theta_vec[i];
		sum_theta2 += theta_vec[n_comps + i];
	}
	theta_vec[n_comps-1] = 1.0-sum_theta1;
	theta_vec[2*n_comps-1] = 1.0-sum_theta2;
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


double gd_bound_cluster(vector<double>& gammas,int total, int n_comps){
	double *theta;
	int th_id;
	//double theta_matrix[total][n_comps];
	double klest = 0.0, kl, lpdf, loglike = 0.0, temp; 
	//double alpha_temp[n_comps], beta_temp[n_comps-1], g_temp[n_comps - 1], sum_theta_k[n_comps - 1];
	int i, j, k, t, z;
	//double a_b_ab_boost_lgamma[K-1], coefficients[K-2];
	//ofstream malakia("malakia.txt");
	
	// z denotes the cluster
	//cout<<"this"<<endl;
	for (i = 0; i<n_comps - 1;++i){
		z = vb_cluster[i];
		alpha_temp1[i] = exp(gammas[z])*new_alpha[i];		
		beta_temp1[i] = exp(gammas[z])*new_beta[i];
		//malakia << z << " "<< alpha_temp[i] << " "<<beta_temp[i]<<endl;
		//if (alpha_temp[i] <= 0.0 || beta_temp[i] <= 0.0 ) cout<<"ooops " <<  i<< " : cluster = "<< z<<", gammas[z] = "<< gammas[z]<< endl;
		a_b_ab_boost_lgamma1[i] = boost::math::lgamma(alpha_temp1[i], my_policy()) + boost::math::lgamma(beta_temp1[i], my_policy()) - boost::math::lgamma(alpha_temp1[i] + beta_temp1[i], my_policy());
	}
	//malakia.close();
	for (i = 0; i<n_comps - 2;++i){
		coefficients1[i] = beta_temp1[i] - ( alpha_temp1[i+1] + beta_temp1[i+1] );
	}



	klest = 0.0;
	#pragma omp parallel num_threads(n_threads) reduction(+:klest)
	{
		//#pragma omp for private(i,j,k,temp, loglike, lpdf, theta, kl, sum_theta_k1, g_temp1)
		#pragma omp for private(i,j,k,temp, loglike, lpdf, theta, kl, th_id)
			for ( i = 0; i<total; i++){
				th_id = omp_get_thread_num();
				//sum_theta_k1.resize(K);
				//g_temp1.resize(K);
				loglike = 0.0;
				theta = gen_dir(alpha_temp1, beta_temp1,n_comps);
				sum_theta_k1[th_id][0] = theta[0];
				g_temp1[th_id][0] = coefficients1[0]*log(1.0 - sum_theta_k1[th_id][0]);
				lpdf = (alpha_temp1[0]-1.0)*log(theta[0])-a_b_ab_boost_lgamma1[0] + g_temp1[th_id][0];
				for (k = 1; k < n_comps - 2; ++k){
					sum_theta_k1[th_id][k] = sum_theta_k1[th_id][k-1] + theta[k];
					g_temp1[th_id][k] = coefficients1[k]*log(1.0 - sum_theta_k1[th_id][k]);
					lpdf += (alpha_temp1[k]-1.0)*log(theta[k])-a_b_ab_boost_lgamma1[k]+ g_temp1[th_id][k];
				} 
				k = n_comps - 2;
				g_temp1[th_id][k] = (beta_temp1[k] - 1.0)*log(theta[k+1]);
				lpdf += (alpha_temp1[k]-1.0)*log(theta[k])-a_b_ab_boost_lgamma1[k]+ g_temp1[th_id][k];
				for (j = 0; j < n_obs; j++){ 
					temp = 0.0;
					for (k = 0; k < y[j]; k++){
						temp += theta[c_new[j][k]]*x[j][k];
					}
					loglike += log(temp);
				}
				kl = - lpdf + loglike;
				klest +=  kl;
				delete[] theta;
				//cout<<"j = "<<i<<", "<<lpdf<<", "<< setprecision(16)<< loglike<<endl;
			}
	}
	klest /= total;





	klest -= lmbeta_one_k;
	return klest;
}




// u1 u2 mazi



double diff_bound(vector<double>& gammas1,vector<double>&  gammas2,int total, int n_comps){
	//double theta_matrix[total][n_comps];
	double *theta;
	double klest = 0.0, kl, lpdf, loglike = 0.0, temp1, temp2; 
	//double alpha_temp[n_comps], beta_temp[n_comps-1], g_temp[n_comps - 1], sum_theta_k[n_comps - 1];
	int i, j, k, t, z;
	int th_id;
	//double a_b_ab_boost_lgamma[K-1], coefficients[K-2];
	//ofstream malakia("malakia.txt");
	
	// z denotes the cluster
	//cout<<"this"<<endl;
	for (i = 0; i<n_comps - 1;++i){
		z = vb_cluster[i];
		alpha_temp1[i] = exp(gammas1[z])*new_alpha[i];		
		beta_temp1[i] = exp(gammas1[z])*new_beta[i];
		alpha_temp2[i] = exp(gammas2[z])*new_alpha[i];		
		beta_temp2[i] = exp(gammas2[z])*new_beta[i];
		//malakia << z << " "<< alpha_temp[i] << " "<<beta_temp[i]<<endl;
		//if (alpha_temp[i] <= 0.0 || beta_temp[i] <= 0.0 ) cout<<"ooops " <<  i<< " : cluster = "<< z<<", gammas[z] = "<< gammas[z]<< endl;
		a_b_ab_boost_lgamma1[i] = boost::math::lgamma(alpha_temp1[i], my_policy()) + boost::math::lgamma(beta_temp1[i], my_policy()) - boost::math::lgamma(alpha_temp1[i] + beta_temp1[i], my_policy());
		a_b_ab_boost_lgamma2[i] = boost::math::lgamma(alpha_temp2[i], my_policy()) + boost::math::lgamma(beta_temp2[i], my_policy()) - boost::math::lgamma(alpha_temp2[i] + beta_temp2[i], my_policy());

	}
	//malakia.close();
	for (i = 0; i<n_comps - 2;++i){
		coefficients1[i] = beta_temp1[i] - ( alpha_temp1[i+1] + beta_temp1[i+1] );
		coefficients2[i] = beta_temp2[i] - ( alpha_temp2[i+1] + beta_temp2[i+1] );
	}

	klest = 0.0;
	#pragma omp parallel num_threads(n_threads) reduction(+:klest)
	{
		//#pragma omp for private(i,j,k,temp1, temp2, loglike, lpdf, theta, kl ,sum_theta_k1,g_temp1,sum_theta_k2,g_temp2)
		#pragma omp for private(i,j,k,temp1, temp2, loglike, lpdf, theta, kl, th_id )
			for ( i = 0; i<total; i++){
				th_id = omp_get_thread_num();
				//sum_theta_k1.resize(K);
				//sum_theta_k2.resize(K);
				//g_temp1.resize(K);
				//g_temp2.resize(K);
				loglike = 0.0;
				theta = gen_dir_tp(alpha_temp1, beta_temp1, alpha_temp2, beta_temp2, n_comps);
				sum_theta_k1[th_id][0] = theta[0];
				sum_theta_k2[th_id][0] = theta[n_comps + 0];
				g_temp1[th_id][0] = coefficients1[0]*log(1.0 - sum_theta_k1[th_id][0]);
				g_temp2[th_id][0] = coefficients2[0]*log(1.0 - sum_theta_k2[th_id][0]);
				lpdf = (alpha_temp1[0]-1.0)*log(theta[0]) - (alpha_temp2[0]-1.0)*log(theta[n_comps + 0]) -(a_b_ab_boost_lgamma1[0] - a_b_ab_boost_lgamma2[0]) + g_temp1[th_id][0] - g_temp2[th_id][0];
				for (k = 1; k < n_comps - 2; ++k){
					sum_theta_k1[th_id][k] = sum_theta_k1[th_id][k-1] + theta[k];
					sum_theta_k2[th_id][k] = sum_theta_k2[th_id][k-1] + theta[n_comps + k];
					g_temp1[th_id][k] = coefficients1[k]*log(1.0 - sum_theta_k1[th_id][k]);
					g_temp2[th_id][k] = coefficients2[k]*log(1.0 - sum_theta_k2[th_id][k]);
					lpdf += (alpha_temp1[k]-1.0)*log(theta[k]) - (alpha_temp2[k]-1.0)*log(theta[n_comps +k]) - (a_b_ab_boost_lgamma1[k] - a_b_ab_boost_lgamma2[k] ) + g_temp1[th_id][k] - g_temp2[th_id][k];
				} 
				k = n_comps - 2;
				g_temp1[th_id][k] = (beta_temp1[k] - 1.0)*log(theta[k+1]);
				g_temp2[th_id][k] = (beta_temp2[k] - 1.0)*log(theta[n_comps + k+1]);
				lpdf += (alpha_temp1[k]-1.0)*log(theta[k]) - (alpha_temp2[k]-1.0)*log(theta[n_comps + k]) - (a_b_ab_boost_lgamma1[k] - a_b_ab_boost_lgamma2[k])+ g_temp1[th_id][k] - g_temp2[th_id][k];
				for (j = 0; j < n_obs; j++){ 
					temp1 = 0.0;
					temp2 = 0.0;
					for (k = 0; k < y[j]; k++){
						temp1 += theta[c_new[j][k]]*x[j][k];
						temp2 += theta[n_comps + c_new[j][k]]*x[j][k];
					}
					loglike += log(temp1) - log(temp2);
				}
				kl = - lpdf + loglike;
				klest +=  kl;
				delete[] theta;
			}
	}
	klest /= total;
	return klest;
}


////////////////////////////////////////////////////////


///////


//double gd_sp_cluster (int max_iter, int n_sim, int n_comps, double a_0, int n_clust, const std::string& out_file_extension){
double gd_sp_cluster (int max_iter, int n_sim, int n_comps, double a_0, int n_clust){
	double bound_gd = 0.0,  best_start, u;
	double alphas[max_iter+100], c_n[max_iter+100];
	int iter, i, j, k, z;
	//clock_t start, end;
	
	//stochastic approximation parameters
	const double c_0 = 6.0, ftol = pow(10.0,-8.0);
	const int n_pars = n_clust;	// number of variational parameters
	const int n_init = 1;		// number of initial iterations
	const int m_conv = 6;		// this determines that the number of runs until convergence should be (m_conv - 2), taking into account the last m_conv function evaluations.
	double obj_function[m_conv];			// this keeps the values of the last  m_conv function evaluations
	int signs[m_conv - 1], runs[m_conv - 2];	// the sign of successive evaluations. the runs of the successive evaluations.
	vector<double> pars(n_pars);
	double b_0;
	double init_bound;



	b_0 = 0.602*a_0;
	b_0 = pow(b_0,1.66113);
	for (iter = 0; iter<max_iter+100; ++iter){
		//alphas[iter] = 1.0/(a_0 + 1.0*(iter+1.0));
		alphas[iter] = 0.602/pow((b_0 + 1.0*(iter+1.0)),0.602);
		//c_n[iter] = pow(1.0*(iter+1.0),-1.0/c_0); 
		c_n[iter] = 1.0/pow(1.0*(iter + 1.0),0.101);
	}

	double unifRand();
	int bernRand();
	double graddiff;
	double condition = 500.0, gradplus, gradminus, mov_av, old_mov_av;
	vector<double> old_pars(n_pars);
	vector<double> u1(n_pars);
	vector<double> u2(n_pars);
	int delta[n_pars]; //this is the bernoulli rv's
	double start, end;


	start = omp_get_wtime( );


	// this is the file that the bound values will be written
	std::string b_values;// ("bound_values");
	b_values = out_file_extension;
	if (dirichlet ==0){	
		b_values += ".bound_values";
	}else{
		b_values += ".bound_values_dir";
	}
	// this is the file that the final generalized dirichlet parameters will be written
	std::string gd_final_pars;
	//gd_final_pars = "dirichlet_parameters";
	gd_final_pars = out_file_extension;
	if (dirichlet ==0){	
		gd_final_pars += ".gen_dir_parameters";
	}else{
		gd_final_pars += ".dir_parameters";
	}
	// this is the file that the while sequence of generalized dirichlet parameters will be written
	//ofstream seqfile ("seq_gen_dir_pars.txt");
	std::string par_values; //("sequence_vb_parameters");
	par_values = out_file_extension;
	if (dirichlet ==0){	
		par_values += ".sequence_vb_parameters";
	}else{
		par_values += ".sequence_vb_parameters_2";
	}
	ofstream seqfile (par_values.data());
	//ofstream boundfile ("bound_values.txt");
	ofstream boundfile ( b_values.data() );





	i = 0;
	j = i;	// index j will correspond to the best initial values
	iter = 0;
	for (k = 0; k<n_pars; ++k) {
		pars[k] = 0.0;	
		delta[k] = bernRand();
		u1[k] = pars[k] + c_n[iter]*delta[k];
		u2[k] = pars[k] - c_n[iter]*delta[k];

		}
	double test_val, prev_bound;
	prev_bound = gd_bound_cluster(pars,100, K);
	init_bound = prev_bound;


	boundfile<< setprecision(16) << prev_bound << endl;

	iter = 1;
	pars[n_pars-1] = 0.0;
	for (k=0; k<n_pars; ++k){
		old_pars[k] = pars[k];
		delta[k] = bernRand();
		u1[k] = pars[k] + c_n[iter]*delta[k];
		u2[k] = pars[k] - c_n[iter]*delta[k];
	}
	
	
	for (i = 0; i<n_pars; ++i){seqfile << pars[i] << " ";}
	seqfile << endl;

	
	
	for (i = 0; i < m_conv; i++){
		obj_function[i] = prev_bound;
	}
	
	vector<double> mov_av_pars(n_pars);
	for (i = 0; i < n_pars; i++){
		mov_av_pars[i] = 0.0;
	}

	double mean_change = 1.0, cur_change;
	while ( (iter < max_iter) && (condition > ftol) ){
		graddiff = diff_bound(u1, u2, n_sim, K);
		boundfile<< setprecision(16) << prev_bound << endl;
		condition = 0.0;
		test_val = graddiff/graddiff;
		
		if ( test_val == 1.0 ){
			cur_change = 0.0;
			for (k=0; k<n_pars; ++k){  // -1 for correcting malakia
				u = alphas[iter]*graddiff/(2.0*c_n[iter]*delta[k]);
				if (abs(u)>0.1){u = 0.0;}
				cur_change += u;
				pars[k] += u; 
				mov_av_pars[k] += pars[k];
				delta[k] = bernRand();
				u1[k] = pars[k] + c_n[iter]*delta[k];
				u2[k] = pars[k] - c_n[iter]*delta[k];
				seqfile << pars[k] << " ";
		
			}
		}else{
			for (k=0; k<n_pars; ++k){
				delta[k] = bernRand();
				u1[k] = pars[k] + c_n[iter]*delta[k];
				u2[k] = pars[k] - c_n[iter]*delta[k];
				mov_av_pars[k] += pars[k];
				seqfile << pars[k] << " ";
		
			}

		}
		
		++iter;
		condition = 5.0;
		if (iter % step == 0) {
			for (i=0; i<n_pars; i++){ mov_av_pars[i] /= (step - 1.0);}
			test_val = gd_bound_cluster(mov_av_pars,10*n_sim, K);
			if (test_val/test_val == 1.0){prev_bound = test_val;}
			for ( i = 1; i < m_conv; i++) {
				obj_function[m_conv - i] = obj_function[m_conv - i - 1];
			}
			obj_function[0] = test_val;
			for ( i = 0; i<m_conv - 1; i++){
				if( obj_function[i] - obj_function[i+1] >0.0 ){
					signs[i] = 1;
				}else{
					signs[i] = -1;
				}
			}
			j = 0;
			for (i = 0; i < m_conv - 2; i++){
				j += signs[i]*signs[i+1];
			}
			if (j ==  2 - m_conv )condition = ftol/10.0;
			for (i = 0; i < n_pars; i++){
				mov_av_pars[i] = 0.0;
			}
			end = omp_get_wtime( );
			cout << "Iteration " << iter <<": " <<setprecision(12)<<" bound: "<< test_val<<", c: "<< j <<", time: "<< setprecision(3)<< end - start<<" sec."<< "\n\n";


		}
		if ( condition != condition ) condition = 5.0;
		old_mov_av = mov_av;

		seqfile << endl;

	}
	seqfile.close();
	//boundfile.close();



	cout << "Last 50 iterations... ";;
//	another 100 iterations for getting the final estimate as an average

	j = iter + 50;
	//double mov_av_pars[n_pars];
	for (i = 0; i < n_pars; i++){
		mov_av_pars[i] = 0.0;
	}

	while ( (iter < j) ){

		//gradplus = gd_bound_cluster(u1,n_sim, K);
		//gradminus = gd_bound_cluster(u2,n_sim, K);
		graddiff = diff_bound(u1, u2, n_sim, K);	
		//boundfile<< setprecision(16) << gd_bound_cluster(pars,n_sim, K) << endl;
		//test_val = (gradplus*gradminus)/(gradplus*gradminus);
		test_val = graddiff/graddiff;
		if ( test_val == 1.0 ){

			condition = 0.0;
			for (k=0; k<n_pars; ++k){
				u = alphas[iter]*graddiff/(2.0*c_n[iter]*delta[k]);
				if (abs(u)>0.1){u = 0.0;}
				pars[k] += u;
				delta[k] = bernRand();
				u1[k] = pars[k] + c_n[iter]*delta[k];
				u2[k] = pars[k] - c_n[iter]*delta[k];
				mov_av_pars[k] += pars[k];
			}
		}else{
				delta[k] = bernRand();
				u1[k] = pars[k] + c_n[iter]*delta[k];
				u2[k] = pars[k] - c_n[iter]*delta[k];
				mov_av_pars[k] += pars[k];

		}
		++iter;
	}
	for (k=0; k<n_pars; ++k){ mov_av_pars[k] /= (50 + 0.0);}
	cout << " Done!"<<endl;





	ofstream outfile ( gd_final_pars.data() );	//this creates a new file
	for (i = 0; i<K-1; ++i){
		z = vb_cluster[i];
		outfile << setprecision(16)<<exp(mov_av_pars[z])*new_alpha[i] << " " << setprecision(16)<<exp(mov_av_pars[z])*new_beta[i] << std::endl;}
	//i = K-1;
	//outfile << exp(0.0)*new_alpha[i] << " " <<" -- " << std::endl;
	outfile.close();
	cout << "Generalized Dirichlet parameters written to file: " <<gd_final_pars <<endl;
	cout << "Variational parameters written to file: " <<par_values <<endl;

	bound_gd = gd_bound_cluster(pars,100, K);
	boundfile<< setprecision(16) << bound_gd << endl;
	boundfile.close();

	cout<<"Initial Log-marginal likelihood bound estimate: "<<setprecision(14)<< init_bound<<endl;
	cout<<"Final Log-marginal likelihood bound estimate: "<<setprecision(14)<< bound_gd<<"\n\n";;

	end = omp_get_wtime( );
	cout << "Total time: " <<  setprecision(3)<<end - start <<" sec."<< "\n\n";


	return bound_gd;
}














