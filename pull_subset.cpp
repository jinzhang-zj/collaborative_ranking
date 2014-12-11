// partition the given comparison dataset into training and testing data based on given parameters
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

void partition(char* input_filename, double per, char* train_filename, char* test_filename) {
	ifstream input_file(input_filename);
	ofstream train_file(train_filename);
	ofstream test_file(test_filename);

	string line;
	if (input_file) {
		while (getline(input_file, line) ) {
			if ((double) rand() / RAND_MAX < per) {
				train_file << line << endl;
			} else {
				test_file << line << endl;
			}
		}
	} else {
		cout << "input file is not open" << endl;
	}

	input_file.close();
	train_file.close();
	test_file.close();
}

int main(int argc, char** argv) {
	if (argc < 5) {
		cout << "Usage: " << argv[0] << " [data_file] [percentage of training] [training_filename] [testing_filename]" << endl;
		cout << "if you specified percentage x of the dataset as training, the rest 1-x will be testing data" << endl;
	}

	char* input_filename = argv[1];
	double per = atof(argv[2]);
	char* train_filename = argv[3];
	char* test_filename = argv[4];

	partition(input_filename, per, train_filename, test_filename);

	return 0;
}
