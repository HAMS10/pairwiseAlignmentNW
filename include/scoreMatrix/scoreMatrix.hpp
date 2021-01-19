/*
 * scoreMatrix.hpp
 *
 *  Created on: 10 ene. 2021
 *      Author: medin
 */

#ifndef SCOREMATRIX_HPP_
#define SCOREMATRIX_HPP_

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <string>

#include <map>

using namespace std;

class scoreMatrix {
protected:
  float values[15][15];
  string bases;
  map<char, int> score_matrix_dictionary;

  void skip_comments( ifstream& file, string& actual_line );
  void generate_score_matrix_index();
  unsigned int score_matrix_index( char base );

public:
  void import_score_matrix ( string file_path );
  int comparison_score( char first_base, char second_base);
  void print_scoring_matrix();


};

#endif /* SCOREMATRIX_HPP_ */
