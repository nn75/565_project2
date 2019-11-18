#include <math.h>
#include "RainfallSimulator.h"
#define SECUNIT 1000000000.0

void printMat(const vector<vector<double>> &matrix) {
  int M = matrix.size();
  if (M == 0) {
    cout << "matrix is empty" << endl;
    return;
  }
  int N = matrix[0].size();
  cout << "matrix: " << endl;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      cout << matrix[i][j] << " ";
    }
    cout << endl;
  }
}
void printMatInt(const vector<vector<int>> &matrix) {
  int M = matrix.size();
  if (M == 0) {
    cout << "matrixInt is empty" << endl;
    return;
  }
  int N = matrix[0].size();
  cout << "matrixInt: " << endl;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      cout << matrix[i][j] << " ";
    }
    cout << endl;
  }
}

bool RainfallSimulator::has_lower_neigh(int curRow, int curCol, const vector<vector<int>> &elevation_mat, vector<vector<int>> &lower_neigh_dir) {
  int N = elevation_mat.size();
  int min_elevation = elevation_mat[curRow][curCol];
  vector<vector<int>> directions{{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
  bool has_lower_neigh = false;
  // find the min elevation among its neighbours
  for (int i = 0; i < directions.size(); i++) {
    int neighRow = curRow + directions[i][0];
    int neighCol = curCol + directions[i][1];
    if (neighRow < 0 || neighRow >= N
        || neighCol < 0 || neighCol >= N) {
      continue;
    }
    if (elevation_mat[neighRow][neighCol] < min_elevation) {
      has_lower_neigh = true;
      min_elevation = elevation_mat[neighRow][neighCol];
    }
  }
  if (!has_lower_neigh) {
    return false;
  }
  for (int i = 0; i < directions.size(); i++) {
    int neighRow = curRow + directions[i][0];
    int neighCol = curCol + directions[i][1];
    if (neighRow < 0 || neighRow >= N
        || neighCol < 0 || neighCol >= N) {
      continue;
    }
    if (elevation_mat[neighRow][neighCol] == min_elevation) {
      vector<int> lower_neigh{neighRow, neighCol};
      lower_neigh_dir.push_back(lower_neigh);
    }
  }
  return true;
}

int RainfallSimulator::endAbsorption(vector<vector<double>> &curWater, vector<vector<double>> &res, double A) {
  int max_time_step = 0;
  int N = curWater.size();
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int cur_time_step = (int) ceil(curWater[i][j] / A);
      max_time_step = max(cur_time_step, max_time_step);
      res[i][j] += curWater[i][j];
      curWater[i][j] = 0.0;
    }
  }
  return max_time_step;
}

bool RainfallSimulator::flowToLowerNeighbour(const vector<vector<int>> &elevation_mat, vector<vector<double>> &curWater) {
  int N = elevation_mat.size();
  vector<vector<double>> flowWaterMat(N, vector<double>(N, 0.0));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      vector<vector<int>> lower_neigh_dir;
      if (curWater[i][j] > 0 && RainfallSimulator::has_lower_neigh(i, j, elevation_mat, lower_neigh_dir)) {
	int neigh_num = lower_neigh_dir.size();
	double flow_volumn = min(curWater[i][j], 1.0);
	double water_portion = flow_volumn / (double) neigh_num;
	curWater[i][j] -= flow_volumn;
	for (int k = 0; k < lower_neigh_dir.size(); k++) {
	  int neighRow = lower_neigh_dir[k][0];
	  int neighCol = lower_neigh_dir[k][1];
	  flowWaterMat[neighRow][neighCol] += water_portion;
	}
      }
    }
  }
  bool no_flow = true;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (flowWaterMat[i][j] != 0) {
	no_flow = false;
	curWater[i][j] += flowWaterMat[i][j];
      }
    }
  }
  return no_flow;
}

int RainfallSimulator::simulate_seq_helper(const vector<vector<int>> &elevation_mat, vector<vector<double>> &res, int M, double A) {
  int time_step = 0;
  int N = elevation_mat.size();
  vector<vector<double>> curWater(N, vector<double>(N, 0.0));
  bool no_flow = false;
  // stop when there is no water flow for each spot
  while (!no_flow) {
    // initialize no_flow at the beginning of each iteration
    no_flow = true;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
	// let water drop on the spot
	if (time_step < M) {
	  curWater[i][j] += 1.0;
	}
	// skip dry spots
	if (curWater[i][j] == 0) {
	  continue;
	}
	// absorb water
	double cur_absorption = min(curWater[i][j], A);
        res[i][j] += cur_absorption;
        curWater[i][j] -= cur_absorption;
      }
    }

    // flow the remaining water to lower neighbour
    no_flow = RainfallSimulator::flowToLowerNeighbour(elevation_mat, curWater);
    
    time_step++;
  }
  int remaining_time_step = RainfallSimulator::endAbsorption(curWater, res, A);
  time_step += remaining_time_step;
  return time_step;
}

vector<vector<double>> RainfallSimulator::simulate_seq(const vector<vector<int>> &elevation_mat, int M, double A, int N) {
  vector<vector<double>> res(N, vector<double>(N, 0.0));
  int absorption_cycle = RainfallSimulator::simulate_seq_helper(elevation_mat, res, M, A);
  cout << absorption_cycle << endl;
  printMat(res);
  return res;
}

double RainfallSimulator::calculateTimeInSecond(struct timespec start_time, struct timespec end_time) {
  double start_sec = (double)start_time.tv_sec * SECUNIT + (double)start_time.tv_nsec;
  double end_sec = (double)end_time.tv_sec * SECUNIT + (double)end_time.tv_nsec;
  return start_sec < end_sec ? (end_sec - start_sec) / SECUNIT : 0;
}
