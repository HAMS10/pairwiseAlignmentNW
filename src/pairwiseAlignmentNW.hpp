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
  vector<unsigned short int> traceback_paths;
  float gap_opening, gap_extending;

  void fill_when_index_equal_zero( unsigned int& row_index, unsigned int& column_index);

public:
  void gaps_load( float gop, float gep );
  void needleman_wunsch( unsigned int sequence0_index, unsigned int sequence1_index );
  
};

#endif /* PAIRWISEALIGNMENTNW_HPP_ */
