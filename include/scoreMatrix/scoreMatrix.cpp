
#include "scoreMatrix.hpp"

//private || Protected methods
void scoreMatrix::skip_comments( ifstream& file, string& actual_line ){
	bool is_line_commented;

	while (getline( file, actual_line )) {
		is_line_commented = (actual_line[0] == '#');

		if (!is_line_commented)
			break;
	}
}

void scoreMatrix::generate_score_matrix_index()
{
  for (size_t base_index = 0; base_index < bases.length(); base_index++) {
	  score_matrix_dictionary[bases[base_index]] = base_index;
  }
}

unsigned int scoreMatrix::score_matrix_index(char base)
{
  return (base == '-') ? score_matrix_dictionary[ 'N' ] : score_matrix_dictionary[ base ];
}

//Public methods

void scoreMatrix::import_score_matrix( string file_path )
{
	string actual_line, word;
	ifstream my_file(file_path);
	bool is_a_character;
	int row_index, column_index;

	if (my_file.is_open()) {
		skip_comments( my_file, actual_line );
		row_index = 0;

		while (getline(my_file, actual_line )) {
			stringstream separate_words(actual_line);
			is_a_character = true;
			column_index = 0;

			while (separate_words >> word) {

				if (!is_a_character) {
					values[row_index][column_index] = stof(word);
					column_index++;
				}
				else
				{
				  bases.append(word);
				  is_a_character = false;
				}
			}
			row_index++;
		}
		my_file.close();
		generate_score_matrix_index();
	}
	else
		cout << "Unable to open file" << '\n';
}

int scoreMatrix::comparison_score( char first_base, char second_base)
{
	return values[ score_matrix_index(first_base) ][ score_matrix_index(second_base) ];
}

void scoreMatrix::print_scoring_matrix(){
  std::cout << ' ' << ' ';

  for (size_t row_index = 0; row_index < bases.length(); row_index++) {
    cout << bases.at(row_index) << "  ";
  }
  cout << '\n';

  for (size_t row_index = 0; row_index < bases.length(); row_index++)
  {
    cout << bases.at(row_index) << ' ';

    for (size_t column_index = 0; column_index < bases.length(); column_index++)
    {
      cout << values[row_index][column_index] << ' ';
    }
    std::cout << '\n';
  }
}
