//============================================================================
// Name        : pairwiseAlignmentNW.cpp
// Author      : Harold Alexander Medina Santacruz
// Version     :
// Copyright   : Your copyright notice
// Description : Pairwise alignment Needleman and Wunch in C++, Ansi-style
//============================================================================

#include "pairwiseAlignmentNW.hpp"

#include <limits>

const float minus_infinity = - numeric_limits<float>::infinity();

void pairwiseAlignmentNW::gaps_load( float gop, float gep )
{
	gap_opening = gop;
	gap_extending = gep;
}

void pairwiseAlignmentNW::needleman_wunsch( unsigned int sequence0_index, unsigned int sequence1_index )
{
  float diagnal_value;
  float extending_value;
  bool index_equal_zero;

  //initialize the matrix with size of sequences and filling with zeros.
  E = vector<vector<float>> (sequences[sequence0_index].size(), vector<float>(sequences[sequence1_index].size(), 0.0));
  F = E;
  H = E;

  for (size_t row_index = 0; row_index < sequences[sequence0_index].size(); row_index++)
  {

    for (size_t column_index = 0; column_index < sequences[sequence1_index].size(); column_index++)
    {
      index_equal_zero = row_index == 0 || column_index == 0;

      if (index_equal_zero)
        fill_when_index_equal_zero( row_index, column_index);
      else
      {
        diagnal_value = H[row_index - 1][column_index] - gap_opening - gap_extending;
        extending_value = E[row_index - 1][column_index] - gap_extending;
        E[row_index][column_index] = max( diagnal_value, extending_value );

        diagnal_value = H[row_index][column_index - 1] - gap_opening - gap_extending;
        extending_value = F[row_index][column_index - 1] - gap_extending;
        F[row_index][column_index] = max( diagnal_value, extending_value );

				diagonal = comparison_score[row_index - 1][column_index - 1]  +  comparison_score( sequences[sequence0_index][row_index], sequences[sequence1_index][column_index] );

				H[row_index][column_index] = max( diagonal , E[row_index][column_index], F[row_index][column_index] );
								
      }
    }
  }
}

void pairwiseAlignmentNW::fill_when_index_equal_zero( unsigned int& row_index, unsigned int& column_index){
  if (row_index == 0 && column_index == 0) {
    E[row_index][column_index] = minus_infinity;
    F[row_index][column_index] = minus_infinity;
    H[row_index][column_index] = 0;
  }
  else if(!(row_index == 0) && column_index == 0)
  {
    F[row_index][column_index] = minus_infinity;
    E[row_index][column_index] = - gap_opening - row_index * gap_extending;
    H[row_index][column_index] = E[row_index][column_index];
  }
  else
  {
    E[row_index][column_index] = minus_infinity;
    F[row_index][column_index] = - gap_opening - column_index * gap_extending;
    H[row_index][column_index] = F[row_index][column_index];
  }
}
