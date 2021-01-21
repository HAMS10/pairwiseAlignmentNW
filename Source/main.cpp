/*
 * main.cpp
 *
 *  Created on: 13 ene. 2021
 *      Author: medin
 */

#include "pairwiseAlignmentNW.hpp"

int main()
{
  pairwiseAlignmentNW align;
  align.import_score_matrix("Data/substitution_matrices/NUC44");
  align.import_sequences("Data/query/msa.fasta");
  align.print_scoring_matrix();
  align.gaps_load(5, 2);
  unsigned int a = 3, b = 4;
  align.needleman_wunsch(a, b);
  align.print_matrices();
  align.print_input_sequences();
  align.print_output_pairwise();


  return 0;
}
