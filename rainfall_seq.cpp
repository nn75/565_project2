#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;


int main(int argc, char* argv[]) {
  if(argc != 6) {
    cout << "usage:./rainfall <P> <M> <A> <N> <elevation_file>" << endl;
    exit(EXIT_FAILURE);
  }
  int P = atoi(argv[1]);   // number of parallel threads
  int M = atoi(argv[2]);   // time stamp
  float A = atof(argv[3]); // absorption rate
  int N = atoi(argv[4]);   // dimension of landscape
  string elevation_file_name = argv[5];
  vector<vector<int> >  elevation_array;

  
  // read elevation file
  elevation_array.assign(N, vector<int>(N, 0));
  ifstream infile;
  infile.open(elevation_file_name); 
  for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			infile >> elevation_array[i][j];
		}
  }
  infile.close();

  // print 
  for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
      cout << setw(6) <<elevation_array[i][j];
		}
    cout << endl;
  }
  

  exit(EXIT_SUCCESS);
}
