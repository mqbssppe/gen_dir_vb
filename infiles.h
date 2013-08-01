#ifndef infiles
#define infiles
#include <string>

const char vb_file[] = "vb_result.m_alphas";	
const std::string prob_file = "test_data.prob";
const std::string tr_info_file = "test_data.tr";
const std::string vb_clusters_file = "vb_clusters_new.txt";

const int n_threads = 4;			// number of threads
const int s_size = 4;				// number of simulated gen_dir vectors
const int step = 50; 				//the bound is computed after averaging the last (step-1) iterations.
const std::string out_file_extension = "out";	// output file prefix
 
#endif

