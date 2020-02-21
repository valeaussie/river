#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_river.h"

using namespace std;

void print_matrix(vector < vector < int > > m);
void print_vector(vector < int > v);


//This is the code for the method:
//Firstly I calculate the lower triangular 3-dimentional matrices called for the sampled events
//and for the samples with substitutions for every particle.
//I then calculate the lower triangular 3-dimensional matrices for the sampled observations
//and for the sampled observations with substitutions
//Finally I calculate the weights.

int main() {

	random_device rd;
	mt19937 generator(rd());

	simprobe();

	print_matrix(X);
	print_matrix(Obs);


	


	/*
	//DEFINITIONS

	//define the container for the sampled events and the sampled observations (0s and 1s)
	vector < vector < vector < double > > > sample;
	vector < vector < vector < size_t > > > sam_obs;
	//define the container for the new sampled events and the new sampled observations (0s and 1s)
	vector < vector < vector < size_t > > > new_sam_obs;
	vector < vector < vector < double > > > new_sample;
	vector < vector < vector < double > > > resampled;
	//define the containter for the unnormalised weights
	vector < vector < vector < double > > > un_weights;
	//define the container for the normalised weights
	vector < vector < vector < double > > > weights;
	//define the number of particles
	double n = 1000;

	//Sampling from a normal distribution with mean 0
	//and variance sigma^2/(1-phi^2) for every particle
	//and store this in a vector clalled "vector_y0"
	//This vector will be used as the starting point to sample all other vectors of events
	//that will populate the matrix "sample"
	vector < double > vector_y0;
	for (unsigned j = 0; j < n; j++) {
		normal_distribution < double > normalDist(0, sigmasq / (1 - phi * phi));
		vector_y0.push_back(normalDist(generator));
	}

	//Sampling for every particle from a normal distribution centred in the previous event times phi
	//and with variance sigma^2, filling the container "sample".
	//Making the substitiution every time I have an observation in real life,
	//filling the container for the new updated events "new_sample".
	for (size_t j = 0; j < n; j++) {
		vector < vector < double > > matrix_sample;
		vector < vector < double > > matrix_new_sample;
		vector < double > row_matrix_sample;
		vector < double > row_matrix_new_sample;
		double y;
		y = vector_y0[j];
		row_matrix_sample.push_back(y);
		if (obs[0][0] == 1) {
			row_matrix_new_sample.push_back(X[0]);
		}
		else { row_matrix_new_sample.push_back(y); }
		matrix_sample.push_back(row_matrix_sample);
		matrix_new_sample.push_back(row_matrix_new_sample);
		for (size_t i = 1; i < N; i++) {
			vector < size_t > row_obs;
			row_obs = obs[i];
			normal_distribution < double > normalDist(phi * row_matrix_new_sample[i - 1], sigmasq);
			double gen = normalDist(generator);
			row_matrix_new_sample.push_back(gen);
			row_matrix_sample.push_back(gen);
			for (size_t k = 0; k < i; k++) {
				if (row_obs[k] == 1) {
					row_matrix_sample[k] = row_matrix_new_sample[k];
				}
			}
			for (size_t k = 0; k < i + 1; k++) {
				if (row_obs[k] == 1) {
					row_matrix_new_sample[k] = X[k];
				}
			}
			matrix_sample.push_back(row_matrix_sample);
			matrix_new_sample.push_back(row_matrix_new_sample);
		}
		row_matrix_sample.clear();
		row_matrix_new_sample.clear();
		sample.push_back(matrix_sample);
		new_sample.push_back(matrix_new_sample);
		matrix_sample.clear();
		matrix_new_sample.clear();
	}

	//Sampling for every particle from a bernoulli distribution with probability p
	//filling the container of the sampled observations "sam_obs".
	//Substituting a0 with a 1 every time I have an observation in real life,
	//filling the matrix of updated sampled observations "new_sam_obs"
	for (size_t j = 0; j < n; j++) {
		vector < vector < size_t > > matrix_obs;
		vector < vector < size_t > > matrix_new_obs;
		vector < size_t > row_matrix_obs;
		vector < size_t > row_matrix_new_obs;
		for (size_t i = 0; i < N; i++) {
			vector < size_t > row_obs;
			row_obs = obs[i];
			bernoulli_distribution BerDist(p);
			double gen = BerDist(generator);
			row_matrix_obs.push_back(0);
			row_matrix_new_obs.push_back(0);
			for (size_t k = 0; k < i + 1; k++) {
				if (row_obs[k] == 1) {
					row_matrix_new_obs[k] = 1;
				}
				else if (row_obs[k] == 0) {
					row_matrix_obs[k] = gen;
				}
			}
			for (size_t k = 0; k < i + 1; k++) {
				if (row_obs[k] == 0 && row_matrix_obs[k] == 1) {
					row_matrix_new_obs[k] = 0;
				}
				else if (row_obs[k] == 0 && row_matrix_obs[k] == 0) {
					row_matrix_new_obs[k] = 0;
				}
			}
			matrix_obs.push_back(row_matrix_obs);
			matrix_new_obs.push_back(row_matrix_new_obs);
		}
		for (size_t i = 0; i < N - 1; i++) {
			for (size_t k = 0; k < i + 1; k++) {
				if (matrix_new_obs[i][k] == 1) {
					matrix_obs[i + 1][k] = matrix_new_obs[i][k];
				}
			}
		}
		row_matrix_obs.clear();
		row_matrix_new_obs.clear();
		sam_obs.push_back(matrix_obs);
		new_sam_obs.push_back(matrix_new_obs);
		matrix_obs.clear();
		matrix_new_obs.clear();
	}

	//Finding the unnormalised weights (using log then exponentiating)
	//filling the container "un_weights".
	//This is an important part of the code, should be always sure it is correct.

	vector < vector < vector < double > > > tresampled;
	for (size_t j = 0; j < N; j++) {
		vector < vector < double > > matrix_un_weights;
		vector < double > vector_w;
		for (size_t i = 0; i < n; i++) {
			vector < double > vector_log_num;
			vector < double > vector_log_den;
			vector < double > row_sample;
			row_sample = sample[i][j];
			vector < size_t > row_obs;
			row_obs = sam_obs[i][j];
			vector < double > row_new_sample;
			row_new_sample = new_sample[i][j];
			vector < size_t > row_new_obs;
			row_new_obs = new_sam_obs[i][j];
			for (size_t k = 0; k < j + 1; k++) {
				double w{ 1 };
				double num{ 0 };
				double den{ 0 };
				double exp_sum_num{ 1 };
				double exp_sum_den{ 1 };
				if (k == 0) { w = 1; }
				else if (((row_new_sample[k] == row_sample[k]) && (row_new_sample[k - 1] == row_sample[k - 1]))) {}
				else {
					num = ((row_new_sample[k] - phi * row_new_sample[k - 1]) * (row_new_sample[k] - phi * row_new_sample[k - 1]));
					den = ((row_sample[k] - phi * row_sample[k - 1]) * (row_sample[k] - phi * row_sample[k - 1]));
				}
				if (k == 0) {}
				else if (row_new_obs[k] == row_obs[k]) {}
				else if (row_new_obs[k] == 1 && row_obs[k] == 0) { num = num + log(p); den = den + log(1 - p); }
				else { num = num + log(1 - p); den = den + log((p)); }
				vector_log_num.push_back(num);
				vector_log_den.push_back(den);
				double sum_num = accumulate(vector_log_num.begin(), vector_log_num.end(), 0.0);
				double sum_den = accumulate(vector_log_den.begin(), vector_log_den.end(), 0.0);
				exp_sum_num = exp(sum_num);
				exp_sum_den = exp(sum_den);
				if (exp_sum_den != 0) { w = exp_sum_num / exp_sum_den; }
				vector_w.push_back(w);
			}
			matrix_un_weights.push_back(vector_w);
			vector_w.clear();
		}

		// begin resampling step (resampling every time)
		vector < vector < double > > temp_matrix_x;
		vector < double > temp_x;
		for (size_t k = 0; k < j + 1; k++) {
			if (obs[j][k] == 1) {
				for (size_t l = 0; l < n; l++) {
					temp_x.push_back(new_sample[l][j][k]);
				}
				temp_matrix_x.push_back(temp_x);
				temp_x.clear();
			}
			else {
				vector < double > temp_w;
				for (size_t l = 0; l < n; l++) {
					temp_w.push_back(matrix_un_weights[l][k]);
				}
				for (size_t l = 0; l < n; l++) {
					discrete_distribution < int > discrete(temp_w.begin(), temp_w.end());
					temp_x.push_back(new_sample[discrete(generator)][j][k]);
				}
				temp_matrix_x.push_back(temp_x);
				temp_x.clear();
			}
		}
		vector < vector < double > > ttmatrix;
		vector < double > ttvector;
		for (size_t k = 0; k < n; k++) {
			for (size_t l = 0; l < j + 1; l++) {
				ttvector.push_back(temp_matrix_x[l][k]);
			}
			ttmatrix.push_back(ttvector);
			ttvector.clear();
		}
		resampled.push_back(ttmatrix);
		ttmatrix.clear();
	}

	//Create .csv files for the plots

	//Create a .csv file with the resampled particles
	ofstream outFile3("./resampled.csv");
	outFile3 << endl;
	for (size_t lin = 0; lin < n; lin++) {
		for (size_t col = 0; col < N; col++) {
			outFile3 << resampled[N - 1][lin][col] << ",";
		}
		outFile3 << endl;
	}
	outFile3.close();

	cout << "another test" << endl;
	*/
	return 0;


}



//functions definitions

//this function prints a vector of integers
void print_vector(vector < int > v) {
	for (const int x : v) cout << x << ' ';
	cout << endl;
}

//this function prints a matrix of integers
void print_matrix(vector < vector < int > > m) {
	for (const vector < int > v : m) {
		for (int x : v) cout << x << ' ';
		cout << endl;
	}
}
