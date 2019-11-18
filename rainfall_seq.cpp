#include <fstream>
#include <iomanip>
#include "RainfallSimulator.h"

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
  vector<vector<int> >  elevation_mat(N, vector<int>(N));

  // read elevation file
  ifstream infile;
  infile.open(elevation_file_name); 
  for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			infile >> elevation_mat[i][j];
		}
  }
  infile.close();

  // print 
  for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
      cout << setw(6) <<elevation_mat[i][j];
		}
    cout << endl;
  }

  struct timespec start_time, end_time;
  clock_gettime(CLOCK_MONOTONIC, &start_time);
  RainfallSimulator simulator;
  simulator.simulate_seq(elevation_mat, M, A, N);
  clock_gettime(CLOCK_MONOTONIC, &end_time);
  double simulation_time = simulator.calculateTimeInSecond(start_time, end_time);
  cout << simulation_time << endl;
  
  exit(EXIT_SUCCESS);
}
