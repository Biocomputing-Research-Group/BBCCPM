/*
 *  Block.cpp
 *  
 *  Created by Cong Li on 7/8/11.
 *  Revised by Zhichao Lian on Apr/10/2013.
 *  Copyright 2013 Yale University. All rights reserved.
 *
 */

#include "Block.h"

Block::Block( int num_pos, HyperParGroup & hyper_par_group, HyperParSingle & hyper_par_single )
{
	this->num_pos = num_pos;
	this->hyper_par_group = hyper_par_group;
	this->hyper_par_single = hyper_par_single;
	//initial_grouping.resize( num_pos );
	
	num_sample = 0;
}

Block::~Block()
{}


int Block::set_data( long start, long num_sample, gsl_matrix * source )
{
	this->block_data = gsl_matrix_submatrix( source, start, 0, num_sample, num_pos );
	this->num_sample = num_sample;
	this->start = start;
	
	return 0;
}

int Block::calc_post()
{
	int i,j;
	gsl_matrix * group_data = NULL;
	gsl_vector * single_data = NULL;
	group_data = gsl_matrix_calloc( num_sample, num_pos);
	for( i = 0; i < num_sample; i ++ )
	{
		for( j = 0; j < num_pos; j ++ )
		{
				gsl_matrix_set( group_data, i, j, gsl_matrix_get( &(block_data.matrix), i, j ) );
		}
	}

	if( group_data == NULL )
	{
		current_post = 0;
		return 0;
	}
	if( group_data->size2 == 1 )
	{
		single_data = gsl_vector_calloc( num_sample );
		gsl_matrix_get_col( single_data, group_data, 0 );
		gsl_matrix_free( group_data );
		
		calc_single_lik( single_data, current_post );
		gsl_vector_free( single_data );
	}else
	{
		calc_group_lik( group_data, current_post);
		gsl_matrix_free( group_data );
	}	
	return 0;
    
}
/*
int Block::propose_mutate( vector< int > & grouping, int str, double & new_post, gsl_rng * rng )
{
	int old_group;
	int pos;
	
	pos = gsl_rng_uniform_int( rng, num_pos );
	old_group = grouping[pos];
	
	int new_group = gsl_rng_uniform_int( rng, 2 );
	if( new_group >= old_group )
	{
		new_group ++;
	}
	
	grouping[pos] = new_group;
	if( 0 == is_valid_grouping( grouping ) )
	{
		new_post = GSL_NEGINF;
		return 0;
	}
	
	if( calc_post( grouping, str, new_post ) )
	{
		return 1;
	}
	
	return 0;
}


int Block::propose_swap( vector< int > & grouping, int str, double & new_post, gsl_rng * rng )
{
	int i, group1, group2;
	
	group1 = gsl_rng_uniform_int( rng, 3 );
	group2 = gsl_rng_uniform_int( rng, 2 );
	
	if( group2 >= group1 )
		group2 ++;
	
	for( i = 0; i < grouping.size(); i ++ )
	{
		if( grouping[i] == group1 )
			grouping[i] = group2;
		if( grouping[i] == group2 )
			grouping[i] = group1;
	}
	
	if( 0 == is_valid_grouping( grouping ) )
	{
		new_post = GSL_NEGINF;
		return 0;
	}
	
	
	if( calc_post( grouping, str, new_post ) )
	{
		return 1;
	}
	
	return 0;
}

int Block::propose_str( vector< int > & grouping, int & str, double & new_post, gsl_rng * rng )
{
	str = 1 - str;
	
	if( calc_post( grouping, str, new_post ) )
	{
		return 1;
	}
	return 0;
}



int Block::calc_post( vector< int > & grouping, int str, double & post )
{
	string group_code;
	vector < int > group_mbr;
	double lik0, lik1, lik01, lik12, lik012;
	if( str == 0 )
	{
		group_mbr.clear();
		group_mbr.resize(1);
		group_mbr[0] = 0;
		
		calc_lik( grouping, group_mbr, lik0 );
		
		group_mbr.clear();
		group_mbr.resize(1);
		group_mbr[0] = 1;
		
		calc_lik( grouping, group_mbr, lik1 );
		
		group_mbr.clear();
		group_mbr.resize(2);
		group_mbr[0] = 0;
		group_mbr[1] = 1;
		
		calc_lik( grouping, group_mbr, lik01 );
		
		group_mbr.clear();
		group_mbr.resize(3);
		group_mbr[0] = 0;
		group_mbr[1] = 1;
		group_mbr[2] = 2;
		
		calc_lik( grouping, group_mbr, lik012 );
		
		post = lik0 + lik1 + lik012 - lik01;
	}else
	{
		group_mbr.clear();
		group_mbr.resize(2);
		group_mbr[0] = 0;
		group_mbr[1] = 1;
		
		calc_lik( grouping, group_mbr, lik01 );
		
		
		group_mbr.clear();
		group_mbr.resize(2);
		group_mbr[0] = 1;
		group_mbr[1] = 2;
		
		calc_lik( grouping, group_mbr, lik12 );
		
		group_mbr.clear();
		group_mbr.resize(1);
		group_mbr[0] = 1;
		
		calc_lik( grouping, group_mbr, lik1 );
		
		post = lik01 + lik12 - lik1;
	}
	
	return 0;
	
}


int Block::is_valid_grouping( vector< int > & grouping )
{
	int g0 = 0;
	int g1 = 0;
	int g2 = 0;
	int i;
	for( i = 0; i < grouping.size(); i ++ )
	{
		if( grouping[i] == 0 )
		{
			g0 = 1;
		}else if( grouping[i] == 1)
		{
			g1 = 1;
		}else
		{
			g2 = 1;
		}
	}
	if( g1 - g0 > 0 )
		return 0;
	if( g2 - g1 > 0 )
		return 0;
	if( g2 - g0 > 0 )
		return 0;
	
	if( g0 == 0 || g1 == 0 )
		return 0;
	else
		return 1;
}



int Block::calc_lik( vector< int > & grouping, vector< int > & group_mbr, double & likelihood )
{
	gsl_matrix * group_data = NULL;
	gsl_vector * single_data = NULL;
	
	get_group_data( grouping, group_mbr, group_data );
	if( group_data == NULL )
	{
		likelihood = 0;
		return 0;
	}
	if( group_data->size2 == 1 )
	{
		single_data = gsl_vector_calloc( num_sample );
		gsl_matrix_get_col( single_data, group_data, 0 );
		gsl_matrix_free( group_data );
		
		calc_single_lik( single_data, likelihood );
		gsl_vector_free( single_data );
	}else
	{
		calc_group_lik( group_data, likelihood );
		gsl_matrix_free( group_data );
	}	
	return 0;
}


int Block::get_group_data( vector< int > & grouping, vector< int > & group_mbr, gsl_matrix * & group_data )
{
	vector< int > pooled_group;
	int i, j;
	for( i = 0; i < grouping.size(); i ++ )
	{
		if( find( group_mbr.begin(), group_mbr.end(), grouping[i] ) != group_mbr.end() )
		{
			pooled_group.push_back( i );
		}
	}
	
	if( group_data != NULL )
	{
		gsl_matrix_free( group_data );
	}
	
	if( pooled_group.size() == 0 )
	{
		group_data = NULL;
		return 0;
	}else
	{
		group_data = gsl_matrix_calloc( num_sample, pooled_group.size() );
		for( i = 0; i < num_sample; i ++ )
		{
			for( j = 0; j < pooled_group.size(); j ++ )
			{
				gsl_matrix_set( group_data, i, j, gsl_matrix_get( &(block_data.matrix), i, pooled_group[j] ) );
			}
		}
		
		return 0;
	}
}

*/
int Block::calc_group_lik( gsl_matrix * group_data, double & likelihood )
{
	int i, j, k;
	double temp;
	
	int c = group_data -> size2;
	long m = group_data -> size1;
	
	gsl_vector * mu;
	double kappa, kappa_m;
	double nu, nu_m;
	gsl_matrix * lamda;
	gsl_matrix * lamda_m;
	gsl_matrix * s;
	gsl_matrix * s_m;
	lamda = gsl_matrix_calloc( c, c );
	lamda_m = gsl_matrix_calloc( c, c );
	s = gsl_matrix_calloc( c, c );
	s_m = gsl_matrix_calloc( c, c );
	mu = gsl_vector_calloc(c);
	
	kappa = hyper_par_group.kappa;
	////////
	nu = hyper_par_group.nu*c;
	////////
	kappa_m = kappa + m;
	nu_m = nu + m;
	
	vector< double > mean;
	mean.resize( c );
	for( i = 0; i < c; i ++ )
	{
		temp = 0;
		for( j = 0; j < m; j ++ )
		{
			temp += gsl_matrix_get( group_data, j, i );
		}
		mean[i] = temp/m;
	}
	
	gsl_vector_set_all( mu, hyper_par_group.mu );
	for( i = 0; i < c; i ++ )
	{
		for( j = 0; j < c; j ++ )
		{
			if( i == j)
			{
				gsl_matrix_set( lamda, i, j, hyper_par_group.lamda_diag );
				gsl_matrix_set( lamda_m, i, j, hyper_par_group.lamda_diag );
			}
			else
			{
				gsl_matrix_set( lamda, i, j, hyper_par_group.lamda_off_diag );
				gsl_matrix_set( lamda_m, i, j, hyper_par_group.lamda_off_diag );
			}
			temp = (mean[i] - gsl_vector_get(mu, i))*(mean[j] - gsl_vector_get(mu, j))*kappa*m/kappa_m;
			gsl_matrix_set( s_m, i, j, temp);
			for( k = 0; k < m; k ++ )
			{
				temp = gsl_matrix_get(s, i, j);
				temp += (gsl_matrix_get(group_data, k, i) - mean[i])*(gsl_matrix_get(group_data, k, j) - mean[j]);
				gsl_matrix_set( s, i, j, temp );
			}
		}
	}
	
	gsl_matrix_add( lamda_m, s );
	gsl_matrix_add( lamda_m, s_m );
	
	likelihood = -m*c*0.5*log(2*PI)+c*0.5*log(kappa/kappa_m);
	for( i = 0; i < c; i ++ )
	{
		likelihood += gsl_sf_lngamma(nu_m*0.5+(-i)*0.5) - gsl_sf_lngamma(nu*0.5+(-i)*0.5);
	}
	
	gsl_permutation * p = gsl_permutation_alloc( c );
	int signum;
	gsl_linalg_LU_decomp(lamda, p, &signum);
	gsl_linalg_LU_decomp(lamda_m, p, &signum);
	likelihood += nu*0.5*gsl_linalg_LU_lndet(lamda) - nu_m*0.5*gsl_linalg_LU_lndet(lamda_m)+m*c*0.5*log(2);
	
	
	gsl_matrix_free(lamda);
	gsl_matrix_free(lamda_m);
	gsl_matrix_free(s);
	gsl_matrix_free(s_m);
	gsl_vector_free(mu);
	gsl_permutation_free(p);
	
	return 0;
}





int Block::calc_single_lik( gsl_vector * single_data, double & likelihood )
{
	double mu = hyper_par_single.mu;
	double nu = hyper_par_single.nu;
	double kappa = hyper_par_single.kappa;
	double lamda = hyper_par_single.lamda;
	int i;
	
	long m = single_data -> size;
	double nu_m = nu + m;
	double kappa_m = kappa + m;
	
	double mean = 0;
	double var = 0;
	for( i = 0; i < m; i ++ )
	{
		mean += gsl_vector_get( single_data, i);
	}
	mean /= m;
	for( i = 0; i < m; i ++ )
	{
		var += pow( gsl_vector_get( single_data, i) - mean, 2 );
	}
	var /= m-1;
	
	double lamda_m = lamda + (m-1)*var + (kappa*m)*(mean-mu)*(mean-mu)/kappa_m;
	likelihood = -0.5*m*log(2*PI) + 0.5*log(kappa/kappa_m) + gsl_sf_lngamma(nu_m*0.5) - gsl_sf_lngamma(nu*0.5) + nu*0.5*log(lamda*0.5) - nu_m*0.5*log(lamda_m*0.5);
	
	return 0;
}
/*
int Block::metro_hast_run( long n_iter, double initial_temperature )
{
	grouping_sample.resize( num_pos );
	
	time_t t;
	t = time( NULL );
	gsl_rng_default_seed = t;
	gsl_rng * rng = gsl_rng_alloc( gsl_rng_taus );
	
	grouping_sample = initial_grouping;
	int ttt=grouping_sample.size();
	str_sample = initial_struct;
	
	vector< int > temp_grouping;
	int temp_str;
	
	
	current_post = 0;
	double propose_post = 0;
	double temperature = 0;
	
	long i=0;
	
	//debug
	int kkk5=calc_post( grouping_sample, str_sample, current_post );
	if (kkk5)
	//if( calc_post( grouping_sample, str_sample, current_post ) )
	{
		gsl_rng_free( rng );
		return 1;
	}
	
	for( i = 0; i < n_iter + 1; i ++ )
	{
		temp_grouping = grouping_sample;
		temp_str = str_sample;

		if( i < n_iter/2 )
			temperature = pow( initial_temperature, 1.0-2*i/n_iter );
		else
			temperature = 1;
		
		if( i %2 )
		{
			if( gsl_rng_uniform( rng ) < 0.5 )
			{
				if( propose_mutate( temp_grouping, temp_str, propose_post, rng ) )
				{
					gsl_rng_free( rng );
					return 1;
				}
			}else
			{
				if( propose_swap( temp_grouping, temp_str, propose_post, rng ) )
				{
					gsl_rng_free( rng );
					return 1;
				}
			}
		}else
		{
			if( propose_str( temp_grouping, temp_str, propose_post, rng ) )
			{
				gsl_rng_free( rng );
				return 1;
			}
		}
		
		if( gsl_rng_uniform( rng ) < exp((propose_post - current_post)/temperature) )
		{
			current_post = propose_post;
			grouping_sample = temp_grouping;
			str_sample = temp_str;
		}
	}
	
	gsl_rng_free( rng );
	
	return 0;
}
*/
long Block::get_num_sample()
{
	return num_sample;
}

long Block::get_start()
{
	return start;
}
/*
int Block::initialize()
{
	time_t t;
	t = time( NULL );
	gsl_rng_default_seed = t;
	gsl_rng * rng = gsl_rng_alloc( gsl_rng_taus );
	
	int i;
	
	if( grouping_sample.size() == num_pos )
	{
		initial_grouping = grouping_sample;
		initial_struct = str_sample;
		grouping_sample.clear();
		current_post = 0;
	}else
	{
		initial_struct = gsl_rng_uniform_int( rng, 2 );
		do
		{
			for( i = 0; i < num_pos; i ++ )
			{
				initial_grouping[i] = gsl_rng_uniform_int( rng, 3 );
			}
		}while( !is_valid_grouping( initial_grouping ) );
	}
	
	gsl_rng_free( rng );
	
	return 0;
}
*/
double Block::get_current_post()
{
	return current_post;
}
/*
int Block::get_last_grouping( vector< int > & grouping )
{
	grouping.clear();
	grouping = grouping_sample;
	return 0;
}

int Block::get_last_struct( int & str )
{
	str = str_sample;
	return 0;
}
*/