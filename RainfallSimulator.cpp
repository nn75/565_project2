#include <math.h>
#include "RainfallSimulator.h"

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

int RainfallSimulator::simulate_one_drop(const vector<vector<int>> &elevation_mat, vector<vector<double>> &res, double A) {
  int time_step = 0;
  int N = elevation_mat.size();
  vector<vector<double>> curWater(N, vector<double>(N, 1.0));
  bool no_flow = true;
  // stop when there is no water flow for each spot
  while (no_flow) {
    // initialize no_flow at the beginning of each iteration
    no_flow = true;
    for (int i = 0; i < elevation_mat.size(); i++) {
      for (int j = 0; j < elevation_mat[0].size(); j++) {
	// skip dry spots
	if (curWater[i][j] == 0) {
	  continue;
	}
	// absorb water
	double cur_absorption = min(curWater[i][j], (double) A);
        res[i][j] += cur_absorption;
        curWater[i][j] -= cur_absorption;
        // flow to lower neighbour
        vector<vector<int>> lower_neigh_dir;
        if (curWater[i][j] > 0 && RainfallSimulator::has_lower_neigh(i, j, elevation_mat, lower_neigh_dir)) {
	  no_flow = false;
	  int neigh_num = lower_neigh_dir.size();
	  //cout << "neigh_num: " << neigh_num << endl;
	  //printMatInt(lower_neigh_dir);
	  double water_portion = curWater[i][j] / neigh_num;
	  curWater[i][j] = 0.0;
	  for (int k = 0; k < lower_neigh_dir.size(); k++) {
	    int neighRow = lower_neigh_dir[k][0];
            int neighCol = lower_neigh_dir[k][1];
            curWater[neighRow][neighCol] += water_portion;
          }
        }
      }
    }
    time_step++;
  }
  int remaining_time_step = RainfallSimulator::endAbsorption(curWater, res, A);
  time_step += remaining_time_step;
  return time_step;
}

vector<vector<double>> RainfallSimulator::simulate_seq(const vector<vector<int>> &elevation_mat, int M, double A, int N) {
  vector<vector<double>> res(N, vector<double>(N, 0.0));
  int absorption_cycle = RainfallSimulator::simulate_one_drop(elevation_mat, res, A);
  if (absorption_cycle > M) {
    cout << "Invalid M because water cannot be absorbed in a cycle." << endl;
  }
  cout << absorption_cycle << endl;
  printMat(res);
  return res;
}
