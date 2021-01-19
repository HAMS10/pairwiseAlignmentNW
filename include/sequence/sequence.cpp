//============================================================================
// Name        : importSequences.cpp
// Author      : Harold Alexander Medina Santacruz
// Version     :
// Copyright   : Your copyright notice
// Description : import dna/aminoacid sequences in C++, Ansi-style
//============================================================================
#include "sequence.hpp"

void sequence::import_sequences( string file_path )
{
  string actual_line, word;
  int sequence_index = 0;
  bool is_sequence_name, is_first_line = true, is_first_sequence = true;
  ifstream my_file( file_path );

  if (my_file.is_open()) {

    while (getline(my_file,actual_line)) {
      stringstream separate_words(actual_line);
      is_first_sequence = true;

      while(separate_words >> word){
        is_sequence_name = word[0] == '>';

        if (is_sequence_name)
        {
          if (is_first_sequence)
            is_first_sequence = false;
          else
          {
            is_first_line = true;
            sequence_index++;
          }
          sequence_name.push_back(word);
          break;
        }
        else
        {
          if (is_first_line) {
            sequences.push_back(word);
          }
          else
          {
            sequences[sequence_index].append(word);
          }
        }

      }

    }

  }
  else
    cout << "Unable to open file" << '\n';
}

vector<string>& sequence::share_sequences()
{
  return sequences;
}

void sequence::print_sequences() {
  for( unsigned int sequence_index = 0 ; sequence_index < sequence_name.size(); sequence_index++ ){
    cout << sequence_name[sequence_index] << '\n';
    cout << sequences[sequence_index] << '\n';
  }
}
