/*
 *  BlockRMS.h
 *  
 *  Created by Cong Li on 7/8/11.
 *  Revised by Zhichao Lian on Apr/10/2013.
 *  Copyright 2013 Yale University. All rights reserved.
 *
 */


#ifndef _RMS_H_
#define _RMS_H_

#include "Block.h"


class RMS
{
	HyperParSingle hyper_par_single;
	HyperParGroup hyper_par_group;
	
	int thin;
	long thinned_size;

	string trace_file_prefix;
	string block_label_file_prefix;
	
	long num_iter_1;
	// used for group subjects
	long num_iter_outer;
	
	gsl_matrix * data;
	// data of all subject
	vector<gsl_matrix *> v_data;
	// number of groups
	int num_groups;
	long num_sample;
	int num_pos;
	double change_point_prior;
	double group_prior;
	
	vector< int > init_block_lengths;
	vector< vector< Block > > samples;
	vector< double > traces;
	vector< vector< int > > block_labels;
	vector< Block > blocks;

	// group info
	vector<int> vi_group_assignment;
	vector<vector<int> > vvi_group_info;
	vector<vector<int> > vvi_group_assignment;
	vector<double> vd_group_traces;

	int propose_divide( double & new_post, vector< Block > & new_blocks, gsl_rng * rng );
	int propose_combine( double & new_post, vector< Block > & new_blocks, gsl_rng * rng );
	int propose_shift( double & new_post, vector< Block > & new_blocks, gsl_rng * rng );

	// group version of block operations
	int propose_divide_group( double & new_post, vector<vector< Block > > & vv_new_group_blocks, gsl_rng * rng, vector<int> & v_group_info );
	int propose_combine_group( double & new_post, vector<vector< Block > > & vv_new_group_blocks, gsl_rng * rng, vector<int> & v_group_info);
	int propose_shift_group( double & new_post, vector<vector< Block > > & vv_new_group_blocks, gsl_rng * rng, vector<int> & v_group_info);

	int convert_to_block_labels();

	// group level operation
	int group_init(gsl_rng * _rng);
	int get_group_info(vector<vector<int> > & vvi_group_info, vector<int> & vi_group_assignment);
	int group_divide(vector<int> & _vi_new_hyper_group, vector<vector<int> > & vvi_hyper_group_info, gsl_rng * _rng);
	int group_combine(vector<int> & _vi_new_hyper_group, vector<vector<int> > & vvi_hyper_group_info, gsl_rng * _rng);
	int group_swap(vector<int> & _vi_new_hyper_group, vector<vector<int> > & vvi_hyper_group_info, gsl_rng * _rng);
	int group_shift(vector<int> & _vi_new_hyper_group, vector<vector<int> > & vvi_hyper_group_info, gsl_rng * _rng);

public:
	int set_change_point_prior( double change_point_prior );
	int set_group_prior( double group_prior );
	int set_output_info(string trace_file_prefix, string block_label_file_prefix, int thin);
	int set_thin( int thin );
	int load_data( long num_sample, long num_pos, string data_file, string par_file,string init_block_file  );
	int load_multiple_data( long num_sample, long num_pos, string data_file, string par_file,string init_block_file  );
	int output_results( string trace_file, string block_label_file, int thin );
	int output_results_group( string trace_file, string block_label_file);
	int metro_hast_run();
	int metro_hast_run_2_levels();
	int metro_hast_inner_run(vector<int> & v_group_info);
	double validate_blocks( vector< int > & block_labels, long nIter );
	int get_num_pos();
	long get_num_sample();
	int set_num_iter( long num_iter_1);
	int set_num_iter_outer(long num_iter_outer);
	
	RMS();
	~RMS();
	
};


#endif
