#ifndef __RAINFALLSIMULATOR_H__
#define __RAINFALLSIMULATOR_H__
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <time.h>

using namespace std;

class RainfallSimulator {
private:
  bool has_lower_neigh(int curRow, int curCol, vector<vector<int>> &lower_neigh_dir);
  int simulate_seq_helper(vector<vector<double>> &res);
  int endAbsorption(vector<vector<double>> &curWater, vector<vector<double>> &res);
  bool flowToLowerNeighbour(vector<vector<double>> &curWater);
public:
   vector<vector<int>> elevation_mat;
  int M;
  double A;
  int N;
  vector<vector<double>> simulate_seq();
  double calculateTimeInSecond(struct timespec start_time, struct timespec end_time);
  void setSimulationParameters(const vector<vector<int>> &elevation_mat_in, int M_in, double A_in, int N_in); 
};

#endif
