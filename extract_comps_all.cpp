#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "collrank.h"

using namespace std;

bool makeComparisons(vector<rating> &ratings, ofstream &output_file, int new_user_id) {

	int n_ratings_current_user = ratings.size();

	if (!n_ratings_current_user) {
		return false;
	}

	vector<comparison> comp_list(0);

	// Construct the whole comparison list for user i
	for(int j1=0; j1<n_ratings_current_user; j1++) {
		for(int j2=j1+1; j2<n_ratings_current_user; j2++) {
					
			if (ratings[j1].score > ratings[j2].score)
				comp_list.push_back(comparison(0, ratings[j1].item_id, ratings[j2].item_id));
			else if (ratings[j2].score < ratings[j1].score)
				comp_list.push_back(comparison(0, ratings[j2].item_id, ratings[j1].item_id));

		}
	}

	vector<int> order;
	for (int j=0; j<comp_list.size(); j++) {
		order.push_back(j);
	}
	random_shuffle(order.begin(), order.end() );

	// Sample comparisons for training 
	for(int j=0; j<comp_list.size() ; j++) {
		int idx = order[j];
		output_file << new_user_id << ' ' << comp_list[idx].item1_id << ' ' << comp_list[idx].item2_id << endl;
	}

	return true;
}

int main(int argc, char **argv) {

	char *input_filename = nullptr, *output_filename = nullptr;

	if (argc < 2) {
		cout << "Extracting all comparions dataset from a rating dataset" << endl;
		cout << "Usage   : ./extract_comps_all [options]" << endl;
		cout << "           -i (input file name for rating dataset) " << endl; 
		cout << "           -o (output file name for traning comparisons) " << endl; 

		return 0;
	}

	for(int i=1; i<argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				case 'i':	// filename
					input_filename = argv[++i];
					cout << input_filename << endl;
					break;

				case 'o':
					output_filename = argv[++i];
					cout << output_filename << endl;
					break;
			}
		}
	}

	if (input_filename == nullptr) { cerr << "Input file required!" << endl; exit(11); }
	
	// Read the dataset
	ifstream input_file;
	input_file.open(input_filename);
	if (!input_file.is_open()) { cerr << "File not opened!" << endl; exit(11); }

	int n_users, n_items, n_ratings;
	input_file >> n_users >> n_items >> n_ratings;
	cout << n_users << " users, " << n_items << " items, " << n_ratings << " ratings" << endl;

	ofstream output_file;
	if (output_filename == nullptr) output_file.open("comp_data.txt"); else output_file.open(output_filename);

	string line;
	rating r;
	vector<rating> ratings(0);
	srand(time(NULL));

	int current_user_id = 1, new_user_id = 1;

	for(int i=0; i<n_ratings; i++) {
		input_file >> r.user_id >> r.item_id >> r.score;
		getline(input_file, line);

		if (current_user_id < r.user_id) {
			if (makeComparisons(ratings, output_file, new_user_id) ) new_user_id++;
			ratings.clear();
			current_user_id = r.user_id;
		}
		
		ratings.push_back(r);
	}
	if (makeComparisons(ratings, output_file, new_user_id) ) new_user_id++;

	new_user_id--;

	input_file.close();
	output_file.close();

	cout << "Comparisons for " << new_user_id << " users extracted" << endl;
}
