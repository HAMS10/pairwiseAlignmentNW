/*
 * method_fill_matrices.cpp
 *
 *  Created on: 19 ene. 2021
 *      Author: medin
 */


//-------------------------------------------------Private || #Potected  methods

#include <algorithm>

#include "pairwiseAlignmentNW.hpp"
//#include <iomanip>

void pairwiseAlignmentNW::fill_matrices()
{
  bool index_equal_zero;

  //initialize the matrices with size of sequences and filling with zeros.
  unsigned int sequence_0_size = pairwise_output_sequences[0].size() + 1;
  unsigned int sequence_1_size = pairwise_output_sequences[1].size() + 1;

  E = vector<vector<float>> (sequence_0_size, vector<float>(sequence_1_size, 0.0));
  F = E;
  H = E;
  traceback_paths = vector<vector<unsigned short int>> (sequence_0_size, vector<unsigned short int>(sequence_1_size, 0.0));

  for (size_t row_index = 0; row_index < sequence_0_size; row_index++)
  {

    for (size_t column_index = 0; column_index < sequence_1_size; column_index++)
    {
      index_equal_zero = row_index == 0 || column_index == 0;

      if (index_equal_zero)
        fill_when_index_equal_zero( row_index, column_index);
      else
      {
        fill_when_index_non_zero( row_index, column_index );
      }
    }
  }
}

void pairwiseAlignmentNW::fill_when_index_equal_zero( unsigned int& row_index, unsigned int& column_index){
  if (row_index == 0 && column_index == 0)
  {
    E[row_index][column_index] = minus_infinity;
    F[row_index][column_index] = minus_infinity;
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
    traceback_paths[row_index][column_index] = 1;
  }
}

void pairwiseAlignmentNW::fill_when_index_non_zero(unsigned int& row_index, unsigned int& column_index)
{
  vector<float> matrix_values_2_compare;

  matrix_values_2_compare =
  {
    H[row_index - 1][column_index] - gap_opening - gap_extending,
    E[row_index - 1][column_index] - gap_extending
  };

  E[row_index][column_index] = find_max_value_of( matrix_values_2_compare );

  matrix_values_2_compare =
  {
    H[row_index][column_index - 1] - gap_opening - gap_extending,
    F[row_index][column_index - 1] - gap_extending
  };

  F[row_index][column_index] = find_max_value_of(matrix_values_2_compare);

  matrix_values_2_compare =
  {
    F[row_index][column_index],
    E[row_index][column_index],
    H[row_index - 1][column_index - 1]  +  comparison_score( pairwise_output_sequences[0][row_index - 1], pairwise_output_sequences[1][column_index - 1])
  };

  H[row_index][column_index] = find_max_value_of( matrix_values_2_compare );

  traceback_paths[row_index][column_index] = find_traceback_paths( matrix_values_2_compare, H[row_index][column_index] );
}

float pairwiseAlignmentNW::find_max_value_of( vector<float> compare_values )
{
	sort(compare_values.begin(), compare_values.end(), greater<float>() );
	return compare_values[ 0 ];
}

unsigned short int pairwiseAlignmentNW::find_traceback_paths( vector<float>& compare_values, float& max )
{
	unsigned short int traceback_values = 0;

  for (size_t value_index = 0; value_index < compare_values.size(); value_index++) {

    if (compare_values[value_index] == max) {

      if (value_index == 0 )
        traceback_values = 1;
      else
        traceback_values = traceback_values | (value_index << 1);
    }
  }
  traceback_values = find_traceback_direction( traceback_values );

	return traceback_values;
}

unsigned short int pairwiseAlignmentNW::find_traceback_direction( unsigned short int& traceback_value)
{
  // value 0 left | value 1 up | value 2 diagonal
  switch (traceback_value) {
    case 1:
    case 3:
    case 5:
    case 7:
      return 0;
    case 2:
    case 6:
      return 1;
    case 0:
    case 4:
      return 2;
    default:
      return 2;
  }
}
