/*
 * traceback.cpp
 *
 *  Created on: 19 ene. 2021
 *      Author: medin
 */
#include "pairwiseAlignmentNW.hpp"

void pairwiseAlignmentNW::traceback(unsigned int& sequence0_index, unsigned int& sequence1_index)
{
  unsigned int row_size = traceback_paths.size() - 1;
  unsigned int column_size = traceback_paths[0].size() - 1;

  unsigned int path;

  unsigned int row_index = row_size;
  unsigned int column_index = column_size;

  bool is_row_or_column_index_greater_zero = row_index > 0 || column_index > 0;
  bool is_row_or_column_index_equal_zero = row_index == 0 || column_index == 0;

  while ( is_row_or_column_index_greater_zero )
  {
    is_row_or_column_index_equal_zero = row_index == 0||column_index == 0;

    path = traceback_paths[row_index][column_index];

    if (is_row_or_column_index_equal_zero)
    	insert_gap_in_sequence_when_any_index_is_zero(sequence0_index, sequence1_index, row_index, column_index);
    else
    {
      if (path == 0)
      {
        insert_gap_in_sequence_index( sequence1_index, column_index, column_size);
        row_index--;
      }
      else if (path == 1) {
        insert_gap_in_sequence_index( sequence0_index, row_index, row_size);
        column_index--;
      }
      else
      {
        row_index--;
        column_index--;
      }
    }
    is_row_or_column_index_greater_zero = row_index > 0 || column_index > 0;
  }
}

void pairwiseAlignmentNW::insert_gap_in_sequence_when_any_index_is_zero(unsigned int& sequence0_index, unsigned int& sequence1_index, unsigned int& row_index, unsigned int& column_index)
{
  if (column_index == 0)
    {
      sequences[sequence1_index].insert(column_index, gap);
      row_index--;
    }
    else if (row_index == 0)
    {
      sequences[sequence0_index].insert(row_index, gap);
      column_index--;
    }
}

void pairwiseAlignmentNW::insert_gap_in_sequence_index(unsigned int& sequence_index, unsigned int& row_or_column_index, unsigned int& row_or_column_size)
{
    if (row_or_column_index != row_or_column_size)
      sequences[sequence_index].insert(row_or_column_index, gap);
    else
      sequences[sequence_index].append(gap);
}

void pairwiseAlignmentNW::print_vector_float(vector<vector<float>>& data_matrix)
{
  for (auto row_index : data_matrix)
  {
    for(auto element : row_index){
      cout << element << "  ";
    }
    cout << '\n';
  }
}

void pairwiseAlignmentNW::print_vector_uint(vector<vector<unsigned short int>>& data_matrix)
{
  for (auto row_index : data_matrix)
  {
    for(auto element : row_index){
      cout << element << "  ";
    }
    cout << '\n';
  }
}


