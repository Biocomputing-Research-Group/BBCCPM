// bccmp.cpp : Defines the entry point for the console application.
//
#include "BlockRMS.h"

int main (int argc, char * const argv[]) {
	
	
	if( argc != 13 )
	{
		// cout << "MCMC(.exe) [data file] [number of samples] [number of variables] [number of iterations] [hyper parameter file] [block label file] [trace file] [thin] [change point prior] [init file]" << endl;
		cout << "bccpm_e(.exe)";
		cout << "[data file]";// argument 1
		cout << "[number of samples]"; // argument 2
		cout << "[number of variables]"; // argument 3
		cout << "[number of inner iterations]"; // argument 4
		cout << "[number of outer iterations]"; // argument 5
		cout << "[hyper parameter file]"; // argument 6
		cout << "[block label file prefix]"; // argument 7
		cout << "[trace file file prefix]"; // argument 8
		cout << "[thin]"; // argument 9
		cout << "[change point prior]"; // argument 10
		cout << "[group prior]"; // argument 11
		cout << "[init file]" << endl; // argument 12
		return 0;
	}
	
	string data_file;
	string hyper_par_file;
	string block_label_file;
	string trace_file;
	string init_block_file;
	int thin;
	
	
	data_file = argv[1];
	hyper_par_file = argv[6];
	block_label_file = argv[7];
	trace_file = argv[8];
	init_block_file = argv[12];
	thin = atoi(argv[9]);


	long num_iter_1 = atol(argv[4]);
	long num_iter_outer = atol(argv[5]);
	long num_sample = atol(argv[2]);
	long num_pos = atol(argv[3]);
	
	double change_point_prior = atof(argv[10]);

	double group_prior = atof(argv[11]);

    
/*
	//for debug
	data_file = "5ROI_1.txt";
	hyper_par_file = "hyper_par.txt";
	block_label_file = "block_label.txt";
	trace_file = "trace.txt";
	init_block_file= "init.txt";
	thin = 1;

	long num_iter_1 = 1000;
	
	long num_sample = 400;
	long num_pos = 5;
	
	float change_point_prior = 0;
*/

    

	RMS rms;
	rms.set_change_point_prior( change_point_prior );
	rms.set_group_prior(group_prior);
	rms.set_num_iter( num_iter_1);
	rms.set_num_iter_outer(num_iter_outer);
	rms.set_output_info(trace_file, block_label_file, thin);
	int i_status = rms.load_multiple_data( num_sample, num_pos, data_file, hyper_par_file, init_block_file);
	if(i_status != 0){
		cout << "Loading data failed. Please check command and input data."<<endl;
		return 1;
	}
	rms.metro_hast_run_2_levels();
	//rms.output_results( trace_file, block_label_file, thin );
    return 0;
}


