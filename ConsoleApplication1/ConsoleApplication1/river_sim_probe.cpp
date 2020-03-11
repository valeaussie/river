
#include <random>
#include <numeric>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_river.h"

using namespace std;

//DEFINITIONS

const double theta = 0.3;
const double phi = 0.3;
const int N = 10;

vector < vector < int > > X;
vector < vector < int > > Obs;
size_t trials{ 0 };
int first_observed{ 0 };



//void print_matrix(vector < vector < int > > m);
//void print_vector(vector < int > v);


/* this is the code to simulate the river model and the matix of the observations.
Simulate the invasion: each time, at each side of the invasion, simulate with proability theta
that the next cell has been invaded (1 for invaded 0 for not invaded).
Simulate the observation: The observations are done by the probe. Each time, at each side
of the invasion, we simulate with probability phi that the next cell has been observed as
invaded, if the invasion has been observed on that cell simulate the next cell. Continue until
the cell is simulated as not observed (1 for invaded 0 for not invaded). */



int simprobe() {

	random_device rd;
	mt19937 generator(rd());

	//sampling from a bernoulli distribution simulate the invasion model. Put the values in a matrix "X"

	//simulate the first invaded cell 
	vector < int > x(N, 0);
	int end_of_range = N - 1;
	int first_invaded = rd() % end_of_range;
	x[first_invaded] = 1;
	X.push_back(x);

	//simulate invasion for all other cells
	vector < int > tempx(N, 0);
	for (size_t i = 0; i < X.size(); i++) {
		for (size_t j = 0; j < N - 1; j++) {
			if (X[i][j] == 1 && X[i][j + 1] == 0) {
				bernoulli_distribution berd(theta);
				tempx[j] = x[j];
				tempx[j + 1] = berd(generator);
			}
		}
		for (size_t j = 1; j < N; j++) {
			if (X[i][j] == 1 && X[i][j - 1] == 0) {
				bernoulli_distribution berd(theta);
				tempx[j] = x[j];
				tempx[j - 1] = berd(generator);
			}
		}
		for (size_t j = 0; j < N - 1; j++) {
			x[j] = tempx[j];
		}
		X.push_back(tempx);
		if (tempx[N - 1] == 1 && tempx[0] == 1) { break; }
	}

	//Creating a csv file with the values of X
	//calling it "river_invasion.csv""
	ofstream outFile001("./river_invasion.csv");
	outFile001 << endl;
	for (size_t lin = 0; lin < X.size(); lin++) {
		for (size_t col = 0; col < N - 1; col++) {
			outFile001 << X[lin][col] << ",";
		}
		outFile001 << X[lin][N - 1];
		outFile001 << endl;
	}
	outFile001.close();

	//Sampling from a bernoulli distribution with probabiltity phi simulate the observations model. 
	//Put the values in a matrix "Obs".
	vector < int > obs(N, 0);
	//calcualte when the first observation will happen with a geometric distribution
	//note that we can observe only in the range where there invasion is
	vector < int > range;
	geometric_distribution<int> geo(phi);
	do {
		trials = geo(generator) + 1;
	} while (trials > X.size());
	for (int i = 0; i < N; i++) {
		if (X[trials - 1][i] == 1) {
			range.push_back(i);
		}
	}
	for (size_t i = 0; i < trials - 1; i++) {
		Obs.push_back(obs);
	}
	if (range[range.size() - 1] == range[0]) {
		first_observed = range[0];
	}
	else {
		first_observed = rd() % (range[range.size() - 1] - range[0]) + range[0];
	}
	obs[first_observed] = 1;
	//Obs.push_back(obs);
	//simulate now the rest of the observations
	
	//calculate left and right which are the first and last invaded cells
	int left{ 0 };
	int right{ 0 };
	vector < int > tempobs(N, 0);
	for (size_t i = trials - 1; i < X.size(); i++) {
		for (size_t j = 0; j < N - 1; j++) {
			if (X[i][j] == 1 && X[i][j + 1] == 0) {
				right = j;
			}
			else if (X[i][N - 1] == 1) { right = N - 1; }
		}
		for (size_t j = 1; j < N; j++) {
			if (X[i][j] == 1 && X[i][j - 1] == 0) {
				left = j;
			}
			else if (X[i][0] == 1) { left = 0; }
		}
		//sample between left and right
		for (size_t j = 0; j < N; j++) {
			tempobs[j] = obs[j];
		}
		if (left == right) {
			Obs.push_back(obs);
		}
		else {
			for (int j = left; j < right; j++) {
				int last_inv_left{ 0 };
				if (obs[j] == 0 && obs[j + 1] == 1) {
					last_inv_left = j + 1;
					for (int j = last_inv_left; j >= left; j--) {
						bernoulli_distribution berd(phi);
						int l = berd(generator);
						if (l == 1) {
							tempobs[j + 1] = obs[j + 1];
							tempobs[j] = l;
							obs[j] = l;
						}
						else { break; }
					}
				}
				int last_inv_right{ 0 };
				if (obs[j] == 1 && obs[j + 1] == 0) {
					last_inv_right = j;
					for (int j = last_inv_right; j < right; j++) {
						bernoulli_distribution berd(phi);
						int l = berd(generator);
						if (l == 1) {
							tempobs[j] = obs[j];
							tempobs[j + 1] = l;
							obs[j + 1] = l;
						}
						else { break; }
					}
				}
			}
			for (size_t j = 0; j < N; j++) {
				obs[j] = tempobs[j];
			}
			Obs.push_back(tempobs);
		}
	}

	//Creating a csv file with the values of Obs
	//calling it "river_observations.csv""
	ofstream outFile002("./river_observations.csv");
	outFile002 << endl;
	for (size_t lin = 0; lin < Obs.size(); lin++) {
		for (size_t col = 0; col < N - 1; col++) {
			outFile002 << Obs[lin][col] << ",";
		}
		outFile002 << Obs[lin][N - 1];
		outFile002 << endl;
	}
	outFile002.close();

	return 0;
}
