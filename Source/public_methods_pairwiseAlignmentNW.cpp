//============================================================================
// Name        : pairwiseAlignmentNW.cpp
// Author      : Harold Alexander Medina Santacruz
// Version     :
// Copyright   : Your copyright notice
// Description : Pairwise alignment Needleman and Wunch in C++, Ansi-style
//============================================================================

#include <algorithm>

#include <string>

#include "pairwiseAlignmentNW.hpp"

//---------------------------------------------------------------Public methods

void pairwiseAlignmentNW::gaps_load( float gop, float gep )
{
	gap_opening = gop;
	gap_extending = gep;
}

void pairwiseAlignmentNW::needleman_wunsch( unsigned int sequence0_index, unsigned int sequence1_index )
{
  fill_matrices(sequence0_index, sequence1_index);
  traceback(sequence0_index, sequence1_index);
}

void pairwiseAlignmentNW::print_matrices()
{
  cout << "Matrix E" << '\n';
  print_vector_float(E);
  cout << "Matrix F" << '\n';
  print_vector_float(F);
  cout << "Matrix H" << '\n';
  print_vector_float(H);
  cout << "Matrix Traceback" << '\n';
  print_vector_uint(traceback_paths);
}
