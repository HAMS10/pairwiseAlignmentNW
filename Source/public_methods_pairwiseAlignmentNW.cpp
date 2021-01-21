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

void pairwiseAlignmentNW::needleman_wunsch(unsigned int sequence_0_index, unsigned int sequence_1_index)
{
	pairwise_output_sequences = {sequences[sequence_0_index], sequences[sequence_1_index]};
	process_sequences_index ={sequence_0_index, sequence_1_index};
	fill_matrices();
	traceback();
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

void pairwiseAlignmentNW::print_output_pairwise()
{
	for( unsigned int sequence_index = 0 ; sequence_index < pairwise_output_sequences.size(); sequence_index++ ){
		cout << sequences_name[process_sequences_index[sequence_index]] << '\n';
		cout << pairwise_output_sequences[sequence_index] << '\n';
	}
}
