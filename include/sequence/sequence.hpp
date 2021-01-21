/*
 * import_sequences.hpp
 *
 *  Created on: 12 ene. 2021
 *      Author: medin
 */

#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <string>

#include <map>

using namespace std;

class sequence {
protected:
  vector<string> sequences;
  vector<string> sequences_name;
public:
  void import_sequences( string file_path );
  vector<string>& share_sequences();
  void print_input_sequences();
};

#endif /* SEQUENCE_HPP_ */
