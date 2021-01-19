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
  align.import_score_matrix("data\\substitution_matrices\\NUC44");
  align.import_sequences("data\\query\\msa.fasta");
  align.print_scoring_matrix();
  align.print_sequences();
  align.gaps_load(10, 0.1);
  align.needleman_wunsch(0, 4);
  align.print_matrices();
  align.print_sequences();

  return 0;
}
