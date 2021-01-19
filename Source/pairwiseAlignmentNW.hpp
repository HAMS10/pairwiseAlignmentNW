/*
 * pairwiseAlignmentNW.hpp
 *
 *  Created on: 13 ene. 2021
 *      Author: medin
 */

#ifndef PAIRWISEALIGNMENTNW_HPP_
#define PAIRWISEALIGNMENTNW_HPP_

#include <vector>
#include <limits>

#include "../Include/scoreMatrix/scoreMatrix.hpp"
#include "../Include/sequence/sequence.hpp"

using namespace std;

class pairwiseAlignmentNW: public scoreMatrix, public sequence {

private:
	//horizontal city array E
	//vertical city array F
	//diagonal city array H
	vector<vector<float>> F, E, H;
	vector<vector<unsigned short int>> traceback_paths;
	float gap_opening, gap_extending;

	const float minus_infinity = -numeric_limits<float>::infinity();
	const string gap = "-";

	void fill_matrices(unsigned int sequence0_index, unsigned int sequence1_index);
	void fill_when_index_equal_zero(unsigned int &row_index, unsigned int &column_index);
	void fill_when_index_non_zero(unsigned int &sequence0_index, unsigned int &sequence1_index, unsigned int &row_index, unsigned int &column_index);
	float find_max_value_of(vector<float> compare_values);
	unsigned short int find_traceback_paths(vector<float> &compare_values, float &max);
	unsigned short int find_traceback_direction(unsigned short int &traceback_value);
	void traceback(unsigned int &sequence0_index, unsigned int &sequence1_index);
	void insert_gap_in_sequence_when_any_index_is_zero(unsigned int &sequence0_index, unsigned int &sequence1_index, unsigned int &row_index, unsigned int &column_index);
	void insert_gap_in_sequence_index(unsigned int &sequence_index, unsigned int &row_or_column_index, unsigned int &row_or_column_size);
	void print_vector_float(vector<vector<float>> &data_matrix);
	void print_vector_uint(vector<vector<unsigned short int>> &data_matrix);

public:
	void gaps_load(float gop, float gep);
	void needleman_wunsch(unsigned int sequence0_index, unsigned int sequence1_index);
	void print_matrices();

};

#endif /* PAIRWISEALIGNMENTNW_HPP_ */
