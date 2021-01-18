/*
 * pairwiseAlignmentNW.hpp
 *
 *  Created on: 13 ene. 2021
 *      Author: medin
 */

#ifndef PAIRWISEALIGNMENTNW_HPP_
#define PAIRWISEALIGNMENTNW_HPP_

#include "scoreMatrix.hpp"
#include "sequence.hpp"

#include <vector>

using namespace std;

class pairwiseAlignmentNW : public scoreMatrix, public sequence {

private:
  //horizontal city array E
  //vertical city array F
  //diagonal city array H
  vector<vector<float>> F, E, H;
  vector<vector<unsigned short int>> traceback_paths;
  float gap_opening, gap_extending;

  void fill_matrices( unsigned int sequence0_index, unsigned int sequence1_index );
  void fill_when_index_equal_zero( unsigned int& row_index, unsigned int& column_index);
  float find_max_value_of( vector<float> compare_values );
  unsigned short int find_traceback_paths( vector<float> compare_values, float& max );
  unsigned short int find_traceback_direction( unsigned short int& traceback_value);
  void traceback(unsigned int& sequence0_index, unsigned int& sequence1_index);
  void print_vector_float(vector<vector<float>>& data_matrix);
  void print_vector_uint(vector<vector<unsigned short int>>& data_matrix);

public:
  void gaps_load( float gop, float gep );
  void needleman_wunsch( unsigned int sequence0_index, unsigned int sequence1_index );
  void print_matrices();

};

#endif /* PAIRWISEALIGNMENTNW_HPP_ */
