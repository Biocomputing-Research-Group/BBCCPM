/*
*  BlockRMS.cpp
*  
*  Created by Cong Li on 7/8/11.
*  Revised by Zhichao Lian on Apr/10/2013.
*  Copyright 2013 Yale University. All rights reserved.
*
*/

#include "BlockRMS.h"

RMS::RMS()
{
	num_sample = 0;
	num_pos = 0;

	num_iter_1 = 0;

	change_point_prior = -1;

	group_prior = -1;

	data = NULL;
}

RMS::~RMS()
{
	if( NULL != data )
	{
		gsl_matrix_free( data );
	}
}

int RMS::propose_divide_group( double & new_post, vector<vector< Block > > & vv_new_group_blocks, gsl_rng * rng, vector<int> & v_group_info){
	// choose the block to divide
	long block_num = gsl_rng_uniform_int( rng, vv_new_group_blocks.at(0).size() );

	long block_start = vv_new_group_blocks.at(0)[block_num].get_start();
	long block_num_sample = vv_new_group_blocks.at(0)[block_num].get_num_sample();

	if( block_num_sample <= 1 ){
		new_post = GSL_NEGINF;
		return 0;
	}

	long change_num = gsl_rng_uniform_int( rng, block_num_sample - 1 );
	change_num ++;

	int i_group_size = vv_new_group_blocks.size();

	new_post = 0;

	for(int i=0;i<i_group_size;i++){
		Block block_1( num_pos, hyper_par_group, hyper_par_single );
		Block block_2( num_pos, hyper_par_group, hyper_par_single );
		block_1.set_data( block_start, change_num, v_data.at(v_group_info.at(i)) );
		block_2.set_data( block_start + change_num, block_num_sample - change_num, v_data.at(v_group_info.at(i)) );
		long num_sample_1 =block_1.get_num_sample();
		long num_sample_2 =block_2.get_num_sample();
		if( num_sample_1 < 2 || num_sample_2 < 2 ){
			new_post = GSL_NEGINF;
			return 0;	
		}
		block_1.calc_post();
		block_2.calc_post();

		vector< Block > temp_blocks;
		temp_blocks.clear();

		for( long j = 0; j < vv_new_group_blocks.at(i).size(); j++ ){
			if( j == block_num ){
				temp_blocks.push_back( block_1 );
				temp_blocks.push_back( block_2 );
				new_post += block_1.get_current_post();
				new_post += block_2.get_current_post();
			}else{
				temp_blocks.push_back( vv_new_group_blocks.at(i).at(j));
				new_post += vv_new_group_blocks.at(i).at(j).get_current_post();
			}
		}
		vv_new_group_blocks.at(i) = temp_blocks;
	}

	for( long j = 0; j < vv_new_group_blocks.at(0).size(); j ++ ){
		new_post += change_point_prior;
	}

	return 0;
}

int RMS::propose_divide( double & new_post, vector< Block > & new_blocks, gsl_rng * rng )
{
	long block_num = gsl_rng_uniform_int( rng, new_blocks.size() );

	long block_start = new_blocks[block_num].get_start();
	long block_num_sample = new_blocks[block_num].get_num_sample();

	if( block_num_sample <= 1 )
	{
		new_post = GSL_NEGINF;
		return 0;
	}
	long change_num = gsl_rng_uniform_int( rng, block_num_sample - 1 );
	change_num ++;

	Block block_1( num_pos, hyper_par_group, hyper_par_single );
	Block block_2( num_pos, hyper_par_group, hyper_par_single );

	block_1.set_data( block_start, change_num, data );
	block_2.set_data( block_start + change_num, block_num_sample - change_num, data );
	long num_sample_1 =block_1.get_num_sample();
	long num_sample_2 =block_2.get_num_sample();

	if( num_sample_1 < 2 || num_sample_2 < 2 )
	{
		new_post = GSL_NEGINF;
		return 0;	
	}
	block_1.calc_post();
	block_2.calc_post();
	//	block_1.initialize();
	//	block_2.initialize();

	//	block_1.metro_hast_run( num_iter_2, 10 );
	//	block_2.metro_hast_run( num_iter_2, 10 );

	vector< Block > temp_blocks;
	temp_blocks.clear();
	new_post = 0;

	for( long i = 0; i < new_blocks.size(); i ++ )
	{
		if( i == block_num )
		{
			temp_blocks.push_back( block_1 );
			temp_blocks.push_back( block_2 );
			new_post += block_1.get_current_post() + change_point_prior;
			new_post += block_2.get_current_post() + change_point_prior;
		}else
		{
			temp_blocks.push_back( new_blocks[i] );
			new_post += new_blocks[i].get_current_post() + change_point_prior;
		}
	}

	new_blocks = temp_blocks;
	return 0;
}

int RMS::propose_combine_group( double & new_post, vector<vector< Block > > & vv_new_group_blocks, gsl_rng * rng, vector<int> & v_group_info){
	if( vv_new_group_blocks.at(0).size() <= 1 ){
		new_post = GSL_NEGINF;
		return 0;
	}

	long block_num = gsl_rng_uniform_int( rng, vv_new_group_blocks.at(0).size() - 1 );
	long block_start = vv_new_group_blocks.at(0)[block_num].get_start();
	long block_num_sample = vv_new_group_blocks.at(0)[block_num].get_num_sample() + vv_new_group_blocks.at(0)[block_num+1].get_num_sample();

	int i_group_size = vv_new_group_blocks.size();

	new_post = 0;

	for(int i=0;i<i_group_size;i++){
		Block combined_block( num_pos, hyper_par_group, hyper_par_single );
		combined_block.set_data( block_start, block_num_sample, v_data.at(v_group_info.at(i)) );
		combined_block.calc_post();

		vector< Block > temp_blocks;
		temp_blocks.clear();

		for( long j = 0; j < vv_new_group_blocks.at(i).size(); j ++ ){
			if( j == block_num ){
				temp_blocks.push_back( combined_block );
				new_post += combined_block.get_current_post();
				j++;
			}else{
				temp_blocks.push_back( vv_new_group_blocks.at(i)[j] );
				new_post += vv_new_group_blocks.at(i)[j].get_current_post();
			}
		}

		vv_new_group_blocks.at(i) = temp_blocks;
	}

	for( long j = 0; j < vv_new_group_blocks.at(0).size(); j ++ ){
		new_post += change_point_prior;
	}

	return 0;
}

int RMS::propose_combine( double & new_post, vector< Block > & new_blocks, gsl_rng * rng )
{
	if( new_blocks.size() <= 1 )
	{
		new_post = GSL_NEGINF;
		return 0;
	}

	long block_num = gsl_rng_uniform_int( rng, new_blocks.size() - 1 );
	long block_start = new_blocks[block_num].get_start();
	long block_num_sample = new_blocks[block_num].get_num_sample() + new_blocks[block_num+1].get_num_sample();
	Block combined_block( num_pos, hyper_par_group, hyper_par_single );
	combined_block.set_data( block_start, block_num_sample, data );
	combined_block.calc_post();

	//	combined_block.initialize();
	//	combined_block.metro_hast_run( num_iter_2, 10 );

	vector< Block > temp_blocks;
	temp_blocks.clear();
	new_post = 0;
	for( long i = 0; i < new_blocks.size(); i ++ )
	{
		if( i == block_num )
		{
			temp_blocks.push_back( combined_block );
			new_post += combined_block.get_current_post() + change_point_prior;
			i ++;
		}else
		{
			temp_blocks.push_back( new_blocks[i] );
			new_post += new_blocks[i].get_current_post() + change_point_prior;
		}
	}
	new_blocks = temp_blocks;
	return 0;
}

int RMS::propose_shift_group( double & new_post, vector<vector< Block > > & vv_new_group_blocks, gsl_rng * rng, vector<int> & v_group_info){
	if( vv_new_group_blocks.at(0).size() <= 1 ){
		new_post = GSL_NEGINF;
		return 0;
	}

	long block_num = gsl_rng_uniform_int( rng, vv_new_group_blocks.at(0).size() - 1 );

	long block_start_1 = vv_new_group_blocks.at(0)[block_num].get_start();
	long block_start_2 = vv_new_group_blocks.at(0)[block_num+1].get_start();
	long num_sample_1 = vv_new_group_blocks.at(0)[block_num].get_num_sample();
	long num_sample_2 = vv_new_group_blocks.at(0)[block_num+1].get_num_sample();

	int right = gsl_rng_uniform_int( rng, 2 );
	if( right ){
		block_start_2 += 1;
		num_sample_1 += 1;
		num_sample_2 -= 1;
	}else{
		block_start_2 -= 1;
		num_sample_1 -= 1;
		num_sample_2 += 1;
	}

	if( num_sample_1 < 2 || num_sample_2 < 2 ){
		new_post = GSL_NEGINF;
		return 0;	
	}

	new_post = 0;

	int i_group_size = vv_new_group_blocks.size();

	for(int i=0;i<i_group_size;i++){
		Block block_1( num_pos, hyper_par_group, hyper_par_single );
		Block block_2( num_pos, hyper_par_group, hyper_par_single );

		block_1.set_data( block_start_1, num_sample_1, v_data.at(v_group_info.at(i)) );
		block_2.set_data( block_start_2, num_sample_2, v_data.at(v_group_info.at(i)) );

		block_1.calc_post();
		block_2.calc_post();

		vector< Block > temp_blocks;
		temp_blocks.clear();

		for( long j = 0; j < vv_new_group_blocks.at(i).size(); j ++ ){
			if( j == block_num ){
				temp_blocks.push_back( block_1 );
				temp_blocks.push_back( block_2 );

				new_post += block_1.get_current_post();
				new_post += block_2.get_current_post();

				j++;
			}else{
				temp_blocks.push_back( vv_new_group_blocks.at(i)[j] );
				new_post += vv_new_group_blocks.at(i)[j].get_current_post();
			}
		}

		vv_new_group_blocks.at(i) = temp_blocks;
	}

	for( long j = 0; j < vv_new_group_blocks.at(0).size(); j ++ ){
		new_post += change_point_prior;
	}

	return 0;
}

int RMS::propose_shift( double & new_post, vector< Block > & new_blocks, gsl_rng * rng )
{
	if( new_blocks.size() <= 1 )
	{
		new_post = GSL_NEGINF;
		return 0;
	}

	long block_num = gsl_rng_uniform_int( rng, new_blocks.size() - 1 );

	Block block_1( num_pos, hyper_par_group, hyper_par_single );
	Block block_2( num_pos, hyper_par_group, hyper_par_single );

	long block_start_1 = new_blocks[block_num].get_start();
	long block_start_2 = new_blocks[block_num+1].get_start();
	long num_sample_1 = new_blocks[block_num].get_num_sample();
	long num_sample_2 = new_blocks[block_num+1].get_num_sample();

	int right = gsl_rng_uniform_int( rng, 2 );

	if( right )
	{
		block_start_2 += 1;
		num_sample_1 += 1;
		num_sample_2 -= 1;
	}else
	{
		block_start_2 -= 1;
		num_sample_1 -= 1;
		num_sample_2 += 1;
	}

	if( num_sample_1 < 2 || num_sample_2 < 2 )
	{
		new_post = GSL_NEGINF;
		return 0;	
	}

	block_1.set_data( block_start_1, num_sample_1, data );
	block_2.set_data( block_start_2, num_sample_2, data );

	block_1.calc_post();
	block_2.calc_post();

	//	block_1.initialize();
	//	block_2.initialize();

	//	block_1.metro_hast_run( num_iter_2, 10 );
	//	block_2.metro_hast_run( num_iter_2, 10 );

	vector< Block > temp_blocks;
	temp_blocks.clear();
	new_post = 0;
	for( long i = 0; i < new_blocks.size(); i ++ )
	{
		if( i == block_num )
		{
			temp_blocks.push_back( block_1 );
			temp_blocks.push_back( block_2 );

			new_post += block_1.get_current_post() + change_point_prior;
			new_post += block_2.get_current_post() + change_point_prior;

			i ++;
		}else
		{
			temp_blocks.push_back( new_blocks[i] );
			new_post += blocks[i].get_current_post() + change_point_prior;
		}

	}

	new_blocks = temp_blocks;
	return 0;
}

int RMS::group_init(gsl_rng * _rng){
	int i_num_indivduals = v_data.size();
	vi_group_assignment.resize(i_num_indivduals);
	vector<int> vi_group;
	int i_num_groups = gsl_rng_uniform_int(_rng, i_num_indivduals) + 1;
	for(int i=0;i<i_num_indivduals;i++){
		vi_group.push_back(gsl_rng_uniform_int(_rng, i_num_groups));
	}
	int i_min = i_num_indivduals;

	vector<int> vi_temp;
	for(int i=0;i<i_num_groups;i++){
		vi_temp.clear();
		i_min = i_num_indivduals;
		for(int j=0;j<i_num_indivduals;j++){
			if(vi_group.at(j) == i){
				vi_temp.push_back(j);
				if(j<i_min){
					i_min = j;
				}
			}
		}

		for(int j=0;j<vi_temp.size();j++){
			this->vi_group_assignment.at(vi_temp.at(j)) = i_min;
		}
	}
	vvi_group_info.clear();
	vvi_group_info.resize(i_num_indivduals, vector<int>());
	return 0;
}

int RMS::get_group_info(vector<vector<int> > & vvi_group_info, vector<int> & vi_group_assignment){
	int i_num_indivduals = v_data.size();
	for(int i=0;i<i_num_indivduals;i++){
		vvi_group_info.at(i).clear();
		vvi_group_info.at(vi_group_assignment.at(i)).push_back(i);
	}
	return 0;
}

int RMS::group_divide(vector<int> & _vi_new_hyper_group, vector<vector<int> > & vvi_hyper_group_info, gsl_rng * _rng){

	vector<int> vi_group_temp;
	for(int i=0;i<vvi_hyper_group_info.size();i++){
		if(vvi_hyper_group_info.at(i).size() != 0){
			vi_group_temp.push_back(i);
		}
	}

	int i_idx = gsl_rng_uniform_int( _rng, vi_group_temp.size());
	int i_group_1 = vi_group_temp.at(i_idx);

	if(vvi_hyper_group_info.at(i_group_1).size()==1){
		return 1;
	}

	vi_group_temp.clear();
	vector<int> vi_group_temp_2;
	vi_group_temp_2.clear();
	for(int i=0;i<vvi_hyper_group_info.at(i_group_1).size();i++){
		if(gsl_rng_uniform_int( _rng, 2)){
			vi_group_temp.push_back(vvi_hyper_group_info.at(i_group_1).at(i));
		}else{
			vi_group_temp_2.push_back(vvi_hyper_group_info.at(i_group_1).at(i));
		}
	}

	for(int i=0;i<vi_group_temp.size();i++){
		_vi_new_hyper_group.at(vi_group_temp.at(i)) = vi_group_temp.at(0);
	}

	for(int i=0;i<vi_group_temp_2.size();i++){
		_vi_new_hyper_group.at(vi_group_temp_2.at(i)) = vi_group_temp_2.at(0);
	}

	return 0;
}

int RMS::group_shift(vector<int> & _vi_new_hyper_group, vector<vector<int> > & vvi_hyper_group_info, gsl_rng * _rng){

	int i_idx = gsl_rng_uniform_int( _rng, _vi_new_hyper_group.size());

	if(gsl_rng_uniform_int(_rng, 2)){// become the new group
		if(_vi_new_hyper_group.at(i_idx) == i_idx && vvi_hyper_group_info.at(i_idx).size() == 1){
			return 1;
		}else{
			if(_vi_new_hyper_group.at(i_idx) != i_idx){
				_vi_new_hyper_group.at(i_idx) = i_idx;
			}else{
				for(int i=1;i<vvi_hyper_group_info.at(i_idx).size();i++){
					_vi_new_hyper_group.at(vvi_hyper_group_info.at(i_idx).at(i)) = vvi_hyper_group_info.at(i_idx).at(1);
				}
			}
		}
	}else{// add to exist group
		vector<int> vi_group_temp;
		for(int i=0;i<vvi_hyper_group_info.size();i++){
			if(i==_vi_new_hyper_group.at(i_idx)){
				continue;
			}
			if(vvi_hyper_group_info.at(i).size() != 0){
				vi_group_temp.push_back(i);
			}
		}

		if(vi_group_temp.size() == 0){
			return 1;
		}
		int i_group_idx = vi_group_temp.at(gsl_rng_uniform_int( _rng, vi_group_temp.size()));
		if(i_idx > i_group_idx){
			_vi_new_hyper_group.at(i_idx) = i_group_idx;
			for(int i=1;i<vvi_hyper_group_info.at(i_idx).size();i++){
				_vi_new_hyper_group.at(vvi_hyper_group_info.at(i_idx).at(i)) = vvi_hyper_group_info.at(i_idx).at(1);
			}
		}else{
			_vi_new_hyper_group.at(i_idx) = i_idx;
			for(int i=0;i<vvi_hyper_group_info.at(i_group_idx).size();i++){
				_vi_new_hyper_group.at(vvi_hyper_group_info.at(i_group_idx).at(i)) = i_idx;
			}
		}


	}

	return 0;
}

int RMS::group_swap(vector<int> & _vi_new_hyper_group, vector<vector<int> > & vvi_hyper_group_info, gsl_rng * _rng){
	vector<int> vi_group_temp;
	for(int i=0;i<vvi_hyper_group_info.size();i++){
		if(vvi_hyper_group_info.at(i).size() != 0){
			vi_group_temp.push_back(i);
		}
	}

	if(vi_group_temp.size() == 1){
		return 1;
	}

	int i_idx = gsl_rng_uniform_int( _rng, vi_group_temp.size());
	int i_group_1 = vi_group_temp.at(i_idx);
	vi_group_temp.erase(vi_group_temp.begin() + i_idx);
	i_idx = gsl_rng_uniform_int( _rng, vi_group_temp.size());
	int i_group_2 = vi_group_temp.at(i_idx);

	int i_ind_1 = 0;
	int i_ind_2 = 0;

	vector<int> vi_group_temp_2;
	vi_group_temp_2.clear();
	vi_group_temp.clear();
	for(int i=0;i<_vi_new_hyper_group.size();i++){
		if(_vi_new_hyper_group.at(i) == i_group_1){
			vi_group_temp.push_back(i);
		}else if(_vi_new_hyper_group.at(i) == i_group_2){
			vi_group_temp_2.push_back(i);
		}
	}

	i_ind_1 = vi_group_temp.at(gsl_rng_uniform_int( _rng, vi_group_temp.size()));
	i_ind_2 = vi_group_temp_2.at(gsl_rng_uniform_int( _rng, vi_group_temp_2.size()));

	_vi_new_hyper_group.at(i_ind_1) = i_group_2;
	_vi_new_hyper_group.at(i_ind_2) = i_group_1;

	if(i_ind_2 < i_group_1){
		for(int i=0;i<vi_group_temp.size();i++){
			if(vi_group_temp.at(i) == i_ind_1){
				continue;
			}
			_vi_new_hyper_group.at(vi_group_temp.at(i)) = i_ind_2;
		}
		_vi_new_hyper_group.at(i_ind_2) = i_ind_2;
	}

	if(i_ind_1 < i_group_2){
		for(int i=0;i<vi_group_temp_2.size();i++){
			if(vi_group_temp_2.at(i) == i_ind_2){
				continue;
			}
			_vi_new_hyper_group.at(vi_group_temp_2.at(i)) = i_ind_1;
		}
		_vi_new_hyper_group.at(i_ind_1) = i_ind_1;
	}

	if(vi_group_temp.size() == 1){
		_vi_new_hyper_group.at(i_ind_2) = i_ind_2;
	}

	if(vi_group_temp_2.size() == 1){
		_vi_new_hyper_group.at(i_ind_1) = i_ind_1;
	}

	return 0;
}

int RMS::group_combine(vector<int> & _vi_new_hyper_group, vector<vector<int> > & vvi_hyper_group_info, gsl_rng * _rng){
	vector<int> vi_group_temp;
	for(int i=0;i<vvi_hyper_group_info.size();i++){
		if(vvi_hyper_group_info.at(i).size() != 0){
			vi_group_temp.push_back(i);
		}
	}

	if(vi_group_temp.size() == 1){
		return 1;
	}

	int i_idx = gsl_rng_uniform_int( _rng, vi_group_temp.size());
	int i_group_1 = vi_group_temp.at(i_idx);
	vi_group_temp.erase(vi_group_temp.begin() + i_idx);
	i_idx = gsl_rng_uniform_int( _rng, vi_group_temp.size());
	int i_group_2 = vi_group_temp.at(i_idx);

	if(i_group_2 < i_group_1){
		int temp = i_group_2;
		i_group_2 = i_group_1;
		i_group_1 = temp;
	}

	for(int i=0;i<_vi_new_hyper_group.size();i++){
		if(_vi_new_hyper_group.at(i) == i_group_2){
			_vi_new_hyper_group.at(i) = i_group_1;
		}
	}

	return 0;
}

bool is_file_exist (const std::string& name) {
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	} else {
		return false;
	}   
}

int RMS::metro_hast_run_2_levels(){
	int keep_n = 100;
	vector<int> * v_group_info = new vector<int>();
	v_group_info->push_back(0);
	v_group_info->push_back(1);

	time_t t = time( NULL );
	gsl_rng_default_seed = t;
	gsl_rng * rng = gsl_rng_alloc( gsl_rng_taus );

	this->vi_group_assignment.clear();
	group_init(rng);
	get_group_info(this->vvi_group_info, this->vi_group_assignment);

	int i_operation = 0;
	vector<vector<int> > vvi_group_info_new;
	int i_num_individuals = v_data.size();
	vvi_group_info_new.clear();
	vvi_group_info_new.resize(v_data.size(), vector<int>());

	vector<double> old_post;
	old_post.clear();
	old_post.resize(i_num_individuals, 0);

	vector<double> new_post;
	new_post.clear();
	new_post.resize(i_num_individuals, 0);

	double group_post = 0;

	// get the first post
	for(long j = 0; j < i_num_individuals; j++){
		if(vvi_group_info.at(j).size() == 0){// empty group
			old_post.at(j) = 0;
		}else{
			metro_hast_inner_run(vvi_group_info.at(j));
			string trace_file = trace_file_prefix + "_ite_" + to_string(0) + "_group_" + to_string(j) + ".txt";
			string block_label_file = block_label_file_prefix + "_ite_" + to_string(0) + "_group_" + to_string(j) + ".txt";
			if(num_iter_outer <= keep_n){
				output_results(trace_file, block_label_file, thin );
			}
			old_post.at(j) = traces.back();
			group_post += group_prior;
		}
	}
	double sum_old_post = accumulate(old_post.begin(), old_post.end(), 0) + group_post;

	vd_group_traces.clear();
	vd_group_traces.resize(num_iter_outer);
	vd_group_traces.at(0) = sum_old_post;
	vvi_group_assignment.clear();
	vvi_group_assignment.resize(num_iter_outer);
	vvi_group_assignment.at(0) = vi_group_assignment;

	// MH run
	for(long i = 1; i < num_iter_outer; i++ ){// outer iteration
		cout << "Iteration " << i << endl;
		group_post = 0;
		i_operation = gsl_rng_uniform_int(rng, 3);
		vector<int> vi_group_assignment_new;
		vi_group_assignment_new = vi_group_assignment;

		switch (i_operation)
		{
		case 0:
			group_divide(vi_group_assignment_new, vvi_group_info, rng);
			break;
		case 1:
			group_combine(vi_group_assignment_new, vvi_group_info, rng);
			break;
		default:
			// group_swap(vi_group_assignment_new, vvi_group_info, rng);
			group_shift(vi_group_assignment_new, vvi_group_info, rng);
			break;
		}

		get_group_info(vvi_group_info_new, vi_group_assignment_new);

		long i_hyper_group_size = i_num_individuals;

		for(long j = 0; j < i_hyper_group_size; j++){// inner mh on each group of subjects
			if(vvi_group_info_new.at(j).size() == 0){// empty group
				new_post.at(j) = 0;
				continue;
			}
			// group prior
			group_post += group_prior;

			if(vvi_group_info_new.at(j) == vvi_group_info.at(j)){// the same, copy the old results
				new_post.at(j) = old_post.at(j);
			}else{
				// mh run
				metro_hast_inner_run(vvi_group_info_new.at(j));
				new_post.at(j) = traces.back();
			}
			// output the block info and trace
			if(	num_iter_outer - i <= keep_n){
				string trace_file = trace_file_prefix + "_ite_" + to_string(i) + "_group_" + to_string(j) + ".txt";
				string block_label_file = block_label_file_prefix + "_ite_" + to_string(i) + "_group_" + to_string(j) + ".txt";
				if(vvi_group_info_new.at(j) == vvi_group_info.at(j)){
					// check if the file exist of not
					string trace_file_old = trace_file_prefix + "_ite_" + to_string(i-1) + "_group_" + to_string(j) + ".txt";
					string block_label_file_old = block_label_file_prefix + "_ite_" + to_string(i-1) + "_group_" + to_string(j) + ".txt";
					if(is_file_exist(trace_file_old)&&is_file_exist(block_label_file_old)){
						ifstream src(trace_file_old, ios::binary);
						ofstream dst(trace_file, ios::binary);
						dst << src.rdbuf();
						dst.close();
						src.close();
						ifstream src_b(block_label_file_old, ios::binary);
						ofstream dst_b(block_label_file, ios::binary);
						dst_b << src_b.rdbuf();
						dst_b.close();
						src_b.close();
					}else{
						// mh run
						metro_hast_inner_run(vvi_group_info_new.at(j));
						new_post.at(j) = traces.back();
						output_results(trace_file, block_label_file, thin );
					}
				}else{
					output_results(trace_file, block_label_file, thin );
				}
			}
		}
		double sum_new_post = accumulate(new_post.begin(), new_post.end(), 0) + group_post;

		if( gsl_rng_uniform( rng ) < exp( sum_new_post - sum_old_post ) ){
			sum_old_post = sum_new_post;
			old_post = new_post;
			vvi_group_info = vvi_group_info_new;
			vi_group_assignment = vi_group_assignment_new;
		}else{
			if(	num_iter_outer - i <= keep_n){
				for(long j = 0; j < i_num_individuals; j++){
					if(vvi_group_info_new.at(j).size() != 0){// empty group
						string trace_file = trace_file_prefix + "_ite_" + to_string(i) + "_group_" + to_string(j) + ".txt";
						string block_label_file = block_label_file_prefix + "_ite_" + to_string(i) + "_group_" + to_string(j) + ".txt";
						remove(trace_file.c_str());
						remove(block_label_file.c_str());
					}

					if(vvi_group_info.at(j).size() == 0){// empty group
						continue;
					}else{
						string trace_file = trace_file_prefix + "_ite_" + to_string(i) + "_group_" + to_string(j) + ".txt";
						string block_label_file = block_label_file_prefix + "_ite_" + to_string(i) + "_group_" + to_string(j) + ".txt";
						string trace_file_old = trace_file_prefix + "_ite_" + to_string(i-1) + "_group_" + to_string(j) + ".txt";
						string block_label_file_old = block_label_file_prefix + "_ite_" + to_string(i-1) + "_group_" + to_string(j) + ".txt";
						if(is_file_exist(trace_file_old)&&is_file_exist(block_label_file_old)){
							ifstream src(trace_file_old, ios::binary);
							ofstream dst(trace_file, ios::binary);
							dst << src.rdbuf();
							dst.close();
							src.close();
							ifstream src_b(block_label_file_old, ios::binary);
							ofstream dst_b(block_label_file, ios::binary);
							dst_b << src_b.rdbuf();
							dst_b.close();
							src_b.close();
						}else{
							// mh run
							metro_hast_inner_run(vvi_group_info.at(j));
							new_post.at(j) = traces.back();
							output_results(trace_file, block_label_file, thin );
						}
					}
				}

			}
		}


		vd_group_traces.at(i) = sum_old_post;
		vvi_group_assignment.at(i) = vi_group_assignment;

	}

	string trace_file = trace_file_prefix + ".txt";
	string block_label_file = block_label_file_prefix + ".txt";
	output_results_group(trace_file, block_label_file);

	delete v_group_info;
	return 0;
}

int RMS::metro_hast_inner_run(vector<int> & v_group_info){
	long i, j, k;

	samples.clear();
	traces.clear();

	samples.resize( num_iter_1 );
	traces.resize( num_iter_1 );

	Block * block_ptr;
	long length = 0;

	long l_group_size = v_group_info.size();
	vector< vector<Block> > * p_vv_blocks = new vector< vector<Block> >();

	for(k=0;k<l_group_size;k++){
		length = 0;
		p_vv_blocks->push_back(vector<Block>());
		for( i = 0; i< init_block_lengths.size(); i ++ ){
			block_ptr = new Block(num_pos, hyper_par_group, hyper_par_single);
			block_ptr->set_data(length, init_block_lengths[i], v_data.at(v_group_info.at(k)));
			length += init_block_lengths[i];
			p_vv_blocks->at(k).push_back(*block_ptr);
			delete(block_ptr);
		}
	}

	time_t t;
	t = time( NULL );
	gsl_rng_default_seed = t;
	gsl_rng * rng = gsl_rng_alloc( gsl_rng_taus );

	int propose;
	vector< vector<Block> > * p_vv_new_blocks = new vector< vector<Block> >();
	double old_post = 0;
	double new_post = 0;

	for(k=0;k<l_group_size;k++){
		for( i = 0; i < p_vv_blocks->at(k).size(); i ++ ){
			p_vv_blocks->at(k)[i].calc_post();
			old_post += p_vv_blocks->at(k)[i].get_current_post();
		}
	}
	for( i = 0; i < p_vv_blocks->at(0).size(); i ++ ){
		old_post += change_point_prior;
	}

	samples[0] = p_vv_blocks->at(0);
	traces[0] = old_post;

	for( i = 1; i < num_iter_1; i ++ ){
		cout << i << "\r";

		(*p_vv_new_blocks) = (*p_vv_blocks);

		propose = gsl_rng_uniform_int( rng, 3 );

		switch( propose ){
		case 0:
			propose_divide_group(new_post, (*p_vv_new_blocks), rng, v_group_info);
			break;
		case 1:
			propose_combine_group(new_post, (*p_vv_new_blocks), rng, v_group_info);
			break;
		case 2:
			propose_shift_group(new_post, (*p_vv_new_blocks), rng, v_group_info);
			break;
		default:
			propose_shift_group(new_post, (*p_vv_new_blocks), rng, v_group_info);
			break;

		}

		if( new_post == GSL_NEGINF ){
			i--;
			continue;
		}

		if( gsl_rng_uniform( rng ) < exp( new_post - old_post ) ){
			(*p_vv_blocks) = (*p_vv_new_blocks);
			old_post= new_post;
		}

		old_post = 0;

		for(k=0;k<l_group_size;k++){
			for( j = 0; j < p_vv_blocks->at(k).size(); j ++ ){
				p_vv_blocks->at(k).at(j).calc_post();
				old_post += p_vv_blocks->at(k).at(j).get_current_post();
			}
		}
		for( j = 0; j < p_vv_blocks->at(0).size(); j ++ ){
			old_post += change_point_prior;
		}

		samples[i] = p_vv_blocks->at(0);
		traces[i] = old_post;
	}

	convert_to_block_labels();

	gsl_rng_free( rng );

	p_vv_blocks->clear();
	p_vv_new_blocks->clear();
	delete p_vv_blocks;
	delete p_vv_new_blocks;
	return 0;
}

int RMS::metro_hast_run()
{
	long i, j;

	samples.clear();
	traces.clear();

	samples.resize( num_iter_1 );
	traces.resize( num_iter_1 );

	Block * block_ptr;
	long length = 0;
	blocks.clear();
	for( i = 0; i< init_block_lengths.size(); i ++ )
	{
		block_ptr = new Block(num_pos, hyper_par_group, hyper_par_single);
		block_ptr->set_data(length, init_block_lengths[i], data);
		length += init_block_lengths[i];
		blocks.push_back(*block_ptr);
	}

	time_t t;
	t = time( NULL );
	gsl_rng_default_seed = t;
	gsl_rng * rng = gsl_rng_alloc( gsl_rng_taus );
	double  temp_prior=change_point_prior;

	int propose;
	vector< Block > new_blocks;
	double old_post = 0;
	double new_post = 0;

	for( i = 0; i < blocks.size(); i ++ )
	{
		blocks[i].calc_post();
		//		blocks[i].metro_hast_run( num_iter_2, 10 );
		old_post += blocks[i].get_current_post() + change_point_prior;
	}
	samples[0] = blocks;
	traces[0] = old_post;


	for( i = 1; i < num_iter_1; i ++ )
	{
		cout << i << endl;

		new_blocks = blocks;	
		propose = gsl_rng_uniform_int( rng, 3 );

		switch( propose )
		{
		case 0:
			propose_divide( new_post, new_blocks, rng );
			break;
		case 1:
			propose_combine( new_post, new_blocks, rng );
			break;
		case 2:
			propose_shift( new_post, new_blocks, rng );
			break;
		default:
			propose_shift( new_post, new_blocks, rng );
			break;
		}
		if( new_post == GSL_NEGINF )
		{
			i --;
			continue;
		}
		//cout << old_post << '\t' << new_post << '\t' << blocks.size() << '\t';

		if( gsl_rng_uniform( rng ) < exp( new_post - old_post ) )
		{
			blocks = new_blocks;
			old_post= new_post;
			//	cout << "Accepted\t";
		}

		old_post = 0;

		for( j = 0; j < blocks.size(); j ++ )
		{
			//blocks[j].initialize();

			//debug
			//blocks[j].metro_hast_run( num_iter_2, 10 );
			blocks[j].calc_post();
			//	blocks[j].metro_hast_run( 10, 1 );

			old_post += blocks[j].get_current_post() + change_point_prior;
		}
		//cout << old_post << '\t' << blocks.size() << endl;
		samples[i] = blocks;
		traces[i] = old_post;


	}
	convert_to_block_labels();

	gsl_rng_free( rng );

	return 0;
}

int RMS::load_multiple_data( long num_sample, long num_pos, string data_file, string par_file, string init_block_file)
{
	this -> num_sample = num_sample;
	this -> num_pos = num_pos;

	if( num_pos < 2 )
	{
		cout << "Data must contain at least 2 columns." << endl;
		return 1;
	}
	// get the number of data in data_file
	ifstream input_data_stream;
	input_data_stream.open(data_file.c_str());
	if( !input_data_stream.good() ){
		cout << "File " << data_file << " doesn't exist." << endl;
		return 1;
	}
	vector< string > v_filenames;
	string s_filename;
	while( input_data_stream >> s_filename ){
		v_filenames.push_back( s_filename );  //add the value of data_name to names.
	}
	input_data_stream.close();
	num_groups = v_filenames.size();

	// read the data one by one
	for(int i=0;i<num_groups;i++){
		gsl_matrix * p_one_data = gsl_matrix_calloc( num_sample, num_pos );
		FILE * data_stream;
		s_filename = v_filenames.at(i);
		data_stream = fopen(s_filename.c_str(), "r" );
		if( NULL == data_stream ){
			cout << "File " << s_filename << " doesn't exist." << endl;
			return 1;
		}

		if( GSL_EFAILED == gsl_matrix_fscanf( data_stream, p_one_data ) ){
			cout << "Wrong data file." << endl;
			return 1;
		}

		fclose( data_stream );

		this->v_data.push_back(p_one_data);
	}

	// bccpm original initial block
	ifstream init_block_stream;
	init_block_lengths.clear();
	init_block_stream.open(init_block_file.c_str());
	string length;
	if( init_block_stream.good())
	{
		while(init_block_stream >> length)
		{
			init_block_lengths.push_back(atoi(length.c_str()));
		}
	}
	init_block_stream.close();
	int sum=accumulate(init_block_lengths.begin(),init_block_lengths.end(),0);
	if (sum!=num_sample){
		init_block_lengths.clear();
		sum=0;
		time_t t;
		t = time( NULL );
		gsl_rng_default_seed = t;
		gsl_rng * rng = gsl_rng_alloc( gsl_rng_taus );
		for (int i=1;i<num_sample-1;i++){

			if (gsl_rng_uniform( rng)>0.9){
				if ((i-sum)>1){
					init_block_lengths.push_back(i-sum);
					sum=i;
				}
			}
		}
		if ((num_sample-sum)<2){
			int temp_back=init_block_lengths.back();
			init_block_lengths.pop_back();
			init_block_lengths.push_back(temp_back+num_sample-sum);
		}
		else
			init_block_lengths.push_back(num_sample-sum);		
	}

	// bccpm original hyper parameters
	ifstream par_stream;
	par_stream.open( par_file.c_str() );

	if( !par_stream.good() )
	{
		cout << "File doesn't exist." << endl;
		return 1;
	}
	vector< double > pars;
	double temp_par;
	string par_name;
	par_stream >> par_name;
	while( par_stream >> temp_par )
	{
		pars.push_back( temp_par );
		par_stream >> par_name;
	}
	if( pars.size() != 9 )
	{
		cout << "Wrong parameter file." << endl;
		return 1;
	}
	par_stream.close();
	hyper_par_single.mu = pars[0];
	hyper_par_single.nu = pars[1];
	hyper_par_single.kappa = pars[2];
	hyper_par_single.lamda = pars[3];

	hyper_par_group.mu = pars[4];
	hyper_par_group.nu = pars[5];
	hyper_par_group.kappa = pars[6];
	hyper_par_group.lamda_diag = pars[7];
	hyper_par_group.lamda_off_diag = pars[8];

	return 0;

}

int RMS::load_data( long num_sample, long num_pos, string data_file, string par_file, string init_block_file)
{
	this -> num_sample = num_sample;
	this -> num_pos = num_pos;

	if( num_pos < 2 )
	{
		cout << "Data must contain at least 2 columns." << endl;
		return 1;
	}


	data = gsl_matrix_calloc( num_sample, num_pos );
	FILE * data_stream;
	data_stream = fopen(data_file.c_str(), "r" );
	if( NULL == data_stream )
	{
		cout << "File doesn't exist." << endl;
		return 1;
	}

	if( GSL_EFAILED == gsl_matrix_fscanf( data_stream, data ) )
	{
		cout << "Wrong data file." << endl;
		return 1;
	}

	fclose( data_stream );

	ifstream init_block_stream;
	init_block_lengths.clear();
	init_block_stream.open(init_block_file.c_str());
	string length;
	if( init_block_stream.good())
	{
		while(init_block_stream >> length)
		{
			init_block_lengths.push_back(atoi(length.c_str()));
		}
	}
	init_block_stream.close();
	int sum=accumulate(init_block_lengths.begin(),init_block_lengths.end(),0);
	if (sum!=num_sample){
		init_block_lengths.clear();
		sum=0;
		time_t t;
		t = time( NULL );
		gsl_rng_default_seed = t;
		gsl_rng * rng = gsl_rng_alloc( gsl_rng_taus );
		for (int i=1;i<num_sample-1;i++){

			if (gsl_rng_uniform( rng)>0.9){
				if ((i-sum)>1){
					init_block_lengths.push_back(i-sum);
					sum=i;
				}
			}
		}
		if ((num_sample-sum)<2){
			int temp_back=init_block_lengths.back();
			init_block_lengths.pop_back();
			init_block_lengths.push_back(temp_back+num_sample-sum);
		}
		else
			init_block_lengths.push_back(num_sample-sum);		
	}

	ifstream par_stream;
	par_stream.open( par_file.c_str() );

	if( !par_stream.good() )
	{
		cout << "File doesn't exist." << endl;
		return 1;
	}
	vector< double > pars;
	double temp_par;
	string par_name;
	par_stream >> par_name;
	while( par_stream >> temp_par )
	{
		pars.push_back( temp_par );
		par_stream >> par_name;
	}
	if( pars.size() != 9 )
	{
		cout << "Wrong parameter file." << endl;
		return 1;
	}
	par_stream.close();
	hyper_par_single.mu = pars[0];
	hyper_par_single.nu = pars[1];
	hyper_par_single.kappa = pars[2];
	hyper_par_single.lamda = pars[3];

	hyper_par_group.mu = pars[4];
	hyper_par_group.nu = pars[5];
	hyper_par_group.kappa = pars[6];
	hyper_par_group.lamda_diag = pars[7];
	hyper_par_group.lamda_off_diag = pars[8];

	return 0;

}

int RMS::convert_to_block_labels()
{
	long num_iter = samples.size();
	block_labels.resize( num_iter );
	long i, j;
	for( i = 0; i < num_iter; i ++ )
	{
		block_labels[i].resize( num_sample );
		for( j = 0; j < num_sample; j ++ )
		{
			block_labels[i][j] = 0;
		}
		for( j = 0; j < samples[i].size(); j ++ )
		{
			block_labels[i][samples[i][j].get_start()] = 1;
		}
	}
	return 0;
}

int RMS::output_results_group( string trace_file, string block_label_file){
	int i, j;

	ofstream trace_stream;
	trace_stream.open( trace_file.c_str() );

	for( i = 0; i < vd_group_traces.size(); i ++ ){
		trace_stream << vd_group_traces[i] << endl;
	}
	trace_stream.close();

	ofstream block_label_stream;
	block_label_stream.open( block_label_file.c_str() );
	for( i = 0; i < vvi_group_assignment.size(); i ++ ){
		for(j=0;j<vvi_group_assignment.at(i).size();j++){
			block_label_stream << vvi_group_assignment.at(i).at(j) << "\t";
		}
		block_label_stream << endl;
	}
	block_label_stream.close();

	return 0;
}

int RMS::output_results( string trace_file, string block_label_file, int thin )
{
	set_thin( thin );

	int i, j;

	ofstream trace_stream;
	trace_stream.open( trace_file.c_str() );

	for( i = 0; i < thinned_size; i ++ )
	{
		trace_stream << traces[i*thin] << endl;
	}
	trace_stream.close();

	ofstream block_label_stream;
	block_label_stream.open( block_label_file.c_str() );
	for( i = 0; i < thinned_size; i ++ )
	{
		for( j = 0; j < block_labels[i*thin].size(); j ++ )
		{
			block_label_stream << block_labels[i*thin][j] << '\t';
		}
		block_label_stream << endl;
	}
	block_label_stream.close();

	return 0;
}




double RMS::validate_blocks( vector< int > & block_labels, long num_iter )
{
	Block temp_block( num_pos, hyper_par_group, hyper_par_single );
	long i;
	long start = 0;
	long num_sample = 0;
	for( i = 0; i < block_labels.size(); i ++ )
	{
		if( block_labels[i] == 1 )
		{
			if( num_sample )
			{
				temp_block.set_data( start, num_sample, data );
				blocks.push_back( temp_block );
			}
			start = i;
			num_sample = 0;
		}
		num_sample ++;
	}
	if( num_sample )
	{
		temp_block.set_data( start, num_sample, data );
		blocks.push_back( temp_block );
	}

	double post = 0;

	for( i = 0; i < blocks.size(); i ++ )
	{
		//		blocks[i].initialize();
		//		blocks[i].metro_hast_run( num_iter, 10 );
		blocks[i].calc_post();
		post += blocks[i].get_current_post() + change_point_prior;
	}
	return post;
}

int RMS::set_num_iter( long num_iter_1)
{
	this->num_iter_1 = num_iter_1;
	return 0;
}

int RMS::set_num_iter_outer(long num_iter_outer){
	this->num_iter_outer = num_iter_outer;
	return 0;
}

int RMS::set_output_info(string trace_file_prefix, string block_label_file_prefix, int thin){
	this->trace_file_prefix = trace_file_prefix;
	this->block_label_file_prefix = block_label_file_prefix;
	set_thin(thin);

	return 0;
}

int RMS::set_thin( int thin )
{
	this->thin = thin;
	this->thinned_size = num_iter_1/thin;

	return 0;
}


int RMS::set_change_point_prior( double change_point_prior )
{
	this->change_point_prior = change_point_prior;
	return 0;
}

int RMS::set_group_prior( double group_prior ){
	this->group_prior = group_prior;
	return 0;
}
