/*
 *  Block.h
 *  
 *  Created by Cong Li on 7/8/11.
 *  Revised by Zhichao Lian on Apr/10/2013.
 *  Copyright 2013 Yale University. All rights reserved.
 *
 */

#ifndef _BLOCK_H_
#define _BLOCK_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <vector>
#include <math.h>
#include <map>
#include <numeric>
#include <algorithm>
#include <time.h>
#include <fstream>
#include <iostream>
#include <string>
#define PI 3.14156265359

using namespace std;

struct HyperParSingle
{
	double mu;
	double nu;
	double kappa;
	double lamda;
};

struct HyperParGroup
{
	double mu;
	double nu;
	double kappa;
	double lamda_diag;
	double lamda_off_diag;
};

class Block
{
	HyperParSingle hyper_par_single;
	HyperParGroup hyper_par_group;
	
	gsl_matrix_view block_data;
	long start;
	long num_sample;
	long num_pos;
	
//	vector< int > initial_grouping;
//	int initial_struct;
//	double initial_post;
	
//	vector< int > grouping_sample;
//	int str_sample;
	double current_post;
	
//	int is_valid_grouping( vector< int > & grouping );
//	int get_group_data( vector< int > & grouping, vector< int > & group_mbr, gsl_matrix * & group_data );
//	int calc_lik( vector< int > & grouping, vector< int > & group_mbr, double & likelihood );
	int calc_group_lik( gsl_matrix * group_data, double & likelihood );
	int calc_single_lik( gsl_vector * single_data, double & likelihood );
	
//	int propose_mutate( vector< int > & grouping, int str, double & new_post, gsl_rng * rng );
//	int propose_swap( vector< int > & grouping, int str, double & new_post, gsl_rng * rng );
//	int propose_str( vector< int > & grouping, int & str, double & new_post, gsl_rng * rng );
	
public:
	int set_data( long start, long num_sample, gsl_matrix * source );
	double get_current_post();
	int calc_post();
//	int get_last_grouping( vector< int > & grouping );
//	int get_last_struct( int & str );
	
//	int initialize();

//	int calc_post( vector< int > & grouping, int str, double & posterior );
//	int metro_hast_run( long n_iter, double temperature );
	
	long get_num_sample();
	long get_start();
	
	Block( int num_pos, HyperParGroup & hyper_par_group, HyperParSingle & hyper_par_single );
	~Block();
};

#endif