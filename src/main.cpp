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
  align.print_scoring_matrix();
  cout << align.comparison_score('A', 'A') ;
  return 0;
}
