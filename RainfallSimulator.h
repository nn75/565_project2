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
  bool has_lower_neigh(int curRow, int curCol, const vector<vector<int>> &elevation_mat, vector<vector<int>> &lower_neigh_dir);
  int simulate_seq_helper(const vector<vector<int>> &elevation_mat, vector<vector<double>> &res, int M, double A);
  int endAbsorption(vector<vector<double>> &curWater, vector<vector<double>> &res, double A);
  bool flowToLowerNeighbour(const vector<vector<int>> &elevation_mat, vector<vector<double>> &curWater);
public:
  vector<vector<double>> simulate_seq(const vector<vector<int>> &elevation_mat, int M, double A, int N);
  double calculateTimeInSecond(struct timespec start_time, struct timespec end_time);
};

#endif
