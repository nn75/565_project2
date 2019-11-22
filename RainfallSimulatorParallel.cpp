#include <math.h>
#include "RainfallSimulatorParallel.h"
#include <pthread.h>
#define SECUNIT 1000000000.0

vector<vector<int>> elevation_mat;
int M;
double A;
int N;
int P;
  
int thread_num;
pthread_mutex_t** mutexes;
pthread_barrier_t barrier;
vector<vector<double>> absorptionMat;
vector<vector<double>> inflowMat;
vector<vector<double>> curWaterMat;
int time_step;
bool no_flow;

void printMat(vector<vector<double>> &matrix) {
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
void printMatInt(vector<vector<int>> &matrix) {
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

bool has_lower_neigh(int curRow, int curCol, vector<vector<int>> &lower_neigh_dir) {
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

int RainfallSimulatorParallel::endAbsorption(vector<vector<double>> &curWater, vector<vector<double>> &res) {
  int max_time_step = 0;
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

bool RainfallSimulatorParallel::flowToLowerNeighbour(vector<vector<double>> &curWater) {
  vector<vector<double>> flowWaterMat(N, vector<double>(N, 0.0));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      vector<vector<int>> lower_neigh_dir;
      if (curWater[i][j] > 0 && has_lower_neigh(i, j, lower_neigh_dir)) {
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

int RainfallSimulatorParallel::simulate_seq_helper(vector<vector<double>> &res) {
  int time_step = 0;
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
    no_flow = RainfallSimulatorParallel::flowToLowerNeighbour(curWater);
    
    time_step++;
  }
  int remaining_time_step = RainfallSimulatorParallel::endAbsorption(curWater, res);
  time_step += remaining_time_step;
  return time_step;
}

void RainfallSimulatorParallel::setSimulationParameters(const vector<vector<int>> &elevation_mat_in,
						int P_in, int M_in, double A_in, int N_in) {
  for (int i = 0; i < N_in; i++) {
    vector<int> temp;
    for (int j = 0; j < N_in; j++) {
      temp.push_back(elevation_mat_in[i][j]);
    }
    elevation_mat.push_back(temp);
  }
  P = P_in;
  M = M_in;
  A = A_in;
  N = N_in;
}

vector<vector<double>> RainfallSimulatorParallel::simulate_seq() {
  vector<vector<double>> res(N, vector<double>(N, 0.0));
  int absorption_cycle = RainfallSimulatorParallel::simulate_seq_helper(res);
  cout << absorption_cycle << endl;
  printMat(res);
  return res;
}

vector<vector<double>> RainfallSimulatorParallel::simulate_pt() {
  vector<vector<double>> res(N, vector<double>(N, 0.0));
  //int absorption_cycle = RainfallSimulatorParallel::simulate_seq_helper(res);
  //cout << absorption_cycle << endl;
  //printMat(res);
  RainfallSimulatorParallel::initializePthreads();
  int absorption_cycle = RainfallSimulatorParallel::simulate_pt_helper(res);
  cout << absorption_cycle << endl;
  printMat(absorptionMat);
  RainfallSimulatorParallel::freeMutexes();
  return res;
}


void RainfallSimulatorParallel::initializePthreads() {
  // to make sure each row will be only assigned a thread
  thread_num = min(N, P);
  // we will have thread_num threads. The matrix will be divided into
  // thread_num parts. We only need to care about the elements on the boundaries.
  // Then we only need (thread_num - 1) * 2 * N mutexes
  mutexes = (pthread_mutex_t**) malloc(((thread_num - 1) * 2) * sizeof(pthread_mutex_t*));
  for (int i = 0; i < (thread_num - 1) * 2; i++) {
    mutexes[i] = (pthread_mutex_t*) malloc(N * sizeof(pthread_mutex_t));
  }
  for (int i = 0; i < (thread_num - 1) * 2; i++) {
    for (int j = 0; j < N; j++) {
      mutexes[i][j] = PTHREAD_MUTEX_INITIALIZER;
    }
  }
  pthread_barrier_init(&barrier, NULL, thread_num);
  
}

void RainfallSimulatorParallel::freeMutexes() {
  for (int i = 0; i < (thread_num - 1) * 2; i++) {
    free(mutexes[i]);
  }
  free(mutexes);
}

void runOneCyclePT(int id) {
  int startRow = id * N / thread_num;
  int endRow = (id + 1) * (N / thread_num) - 1;
  int i, j, k;
  for (i = startRow; i <= endRow; i++) {
    for (j = 0; j < N; j++) {
      if (time_step < M) {
	curWaterMat[i][j] += 1.0;
      }
      if (curWaterMat[i][j] == 0) {
	continue;
      }
      double cur_absorption = min(curWaterMat[i][j], A);
      absorptionMat[i][j] += cur_absorption;
      curWaterMat[i][j] -= cur_absorption;

      // reset the corresponding spot at inflow matrix to 0
      inflowMat[i][j] = 0.0;
    }
  }
  // wait all threads for updating current water matrix and initialize
  // the inflow matrix  
  pthread_barrier_wait(&barrier);
  
  for (i = startRow; i <= endRow; i++) {
    for (j = 0; j < N; j++) {
      vector<vector<int>> lower_neigh_dir;
      if (curWaterMat[i][j] > 0 && has_lower_neigh(i, j, lower_neigh_dir)) {
	int neigh_num = lower_neigh_dir.size();
	double flow_volumn = min(curWaterMat[i][j], 1.0);
	double water_portion = flow_volumn / (double) neigh_num;
	curWaterMat[i][j] -= flow_volumn;
	for (k = 0; k < lower_neigh_dir.size(); k++) {
	  int neighRow = lower_neigh_dir[k][0];
	  int neighCol = lower_neigh_dir[k][1];
	  if (neighRow != 0 && neighRow != N - 1 &&
	     (neighRow == startRow || neighRow == startRow - 1
	      || neighRow == endRow || neighRow == endRow + 1)) {
	    // need to lock to make sure there is only one thread working in the critical section
	    int require_mutex_index;
	    if (neighRow == startRow - 1) {
	      require_mutex_index = id == thread_num - 1 ? 2 * (thread_num - 1) - 2 : (id - 1) * 2;
	    } else if (neighRow == startRow) {
	      require_mutex_index = id == thread_num - 1 ? 2 * (thread_num - 1) - 1 : (id - 1) * 2 + 1;
	    } else if (neighRow == endRow) {
	      require_mutex_index = id == 0 ? 0 : (id - 1) * 2 + 2;
	    } else { // neighRow == endRow + 1
	      require_mutex_index = id == 0 ? 1 : (id - 1) * 2 + 3;
	    }
	    pthread_mutex_lock(&mutexes[require_mutex_index][neighCol]);
	    inflowMat[neighRow][neighCol] += water_portion;
	    pthread_mutex_unlock(&mutexes[require_mutex_index][neighCol]);
	  } else {
	    inflowMat[neighRow][neighCol] += water_portion;
	  }
	}
      }
    }
  }
  pthread_barrier_wait(&barrier);
  for (i = startRow; i <= endRow; i++) {
    for (j = 0; j < N; j++) {
      if (inflowMat[i][j] != 0) {
	no_flow = false;
	curWaterMat[i][j] += inflowMat[i][j];
      }
    }
  }
}

void* worker(void *arg) {
  int id = *((int *) arg);
  runOneCyclePT(id);
  return NULL;
}

int RainfallSimulatorParallel::simulate_pt_helper(vector<vector<double>> &res) {
  int i, j;
  int *p;
  time_step = 0;
  for (i = 0; i < N; i++) {
    vector<double> tempAbsorptionVec;
    vector<double> tempInflowVec;
    vector<double> tempCurWaterVec;
    for (j = 0; j < N; j++) {
      tempAbsorptionVec.push_back(0.0);
      tempInflowVec.push_back(0.0);
      tempCurWaterVec.push_back(0.0);
    }
    absorptionMat.push_back(tempAbsorptionVec);
    inflowMat.push_back(tempInflowVec);
    curWaterMat.push_back(tempCurWaterVec);
  }
  no_flow = false;
  pthread_t *threads = (pthread_t *) malloc(thread_num * sizeof(pthread_t));
  // stop when there is no water flow for each spot
  while (!no_flow) {
    // initialize no_flow at the beginning of each iteration
    no_flow = true;
    for (i = 0; i < thread_num; i++) {
      p = (int *) malloc(sizeof(int));
      *p = i;
      pthread_create(&threads[i], NULL, worker, (void *)(p)); 
    }

    for (i = 0; i < thread_num; i++) {
      pthread_join(threads[i], NULL);
    }
    time_step++;
  }
  int remaining_time_step = RainfallSimulatorParallel::endAbsorption(curWaterMat, absorptionMat);
  time_step += remaining_time_step;
  free(p);
  free(threads);
  return time_step;
}

double RainfallSimulatorParallel::calculateTimeInSecond(struct timespec start_time, struct timespec end_time) {
  double start_sec = (double)start_time.tv_sec * SECUNIT + (double)start_time.tv_nsec;
  double end_sec = (double)end_time.tv_sec * SECUNIT + (double)end_time.tv_nsec;
  return start_sec < end_sec ? (end_sec - start_sec) / SECUNIT : 0;
}
