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
//Firstly I calculate the 3-dimentional matrix called for the sampled events
//making the substitutions for every particle.
//I then calculate the 3-dimensional matrix for the sampled observations
//and for the sampled observations with substitutions
//Finally I calculate the weights and resample.

int main() {

	random_device rd;
	mt19937 generator(rd());

	simprobe();

	//DEFINITIONS

	//number of particles
	int n = 1;
	//define the container for the sampled events and the sampled observations
	vector < vector < vector < int > > > sample;
	int m = X.size();
	vector < vector < vector < int > > > sam_obs(n, vector < vector < int > >(m, vector < int >(N)));
	//define the container for the new sampled events and the new sampled observations
	vector < vector < vector < int > > > new_sample;
	vector < vector < vector < int > > > resampled;
	//define the containter for the unnormalised weights
	vector < vector < vector < int > > > un_weights;
	//define the container for the normalised weights
	vector < vector < vector < int > > > weights;
	//define the number of particles

	print_matrix(X);
	print_matrix(Obs);

	//simulate the first invaded cell for all particles
	//notice that there is a limitation on which cell can be infected
	//which will depend on when the first observation happened and
	//on which cell was first invaded
	vector < vector < int > > samplem;
	vector < int > first_sim_vect;
	int start_of_range = first_observed - (trials - 1);
	int end_of_range = first_observed + (trials + 1);
	if (start_of_range < 0) { start_of_range = 0; }
	if (end_of_range > N - 1) { end_of_range = N - 1; }
	for (int i = 0; i < n; i++) {
		vector < int > samplev(N, 0);
		int first_simulated = start_of_range + rd() % (end_of_range - start_of_range-1);
		samplev[first_simulated] = 1;
		samplem.push_back(samplev);
		first_sim_vect.push_back(first_simulated);
	}
	//Sampling the invasion for every particle from a bernoulli distribution
	//with probability theta, filling the container "sample".
	//Making a substitiution every time I have an observation in real life,
	//filling the container for the new updated events "new_sample".
	for (int j = 0; j < n; j++) {
		vector < vector < int > > matrix_sample;
		matrix_sample.push_back(samplem[j]);
		vector < vector < int > > matrix_new_sample;
		matrix_new_sample.push_back(samplem[j]);
		vector < int > row_matrix_sample(N, 0);
		row_matrix_sample[first_sim_vect[j]] = 1;
		vector < int > row_matrix_new_sample(N, 0);
		row_matrix_new_sample[first_sim_vect[j]] = 1;
		vector < int > temp(N, 0);
		//make all the simulations for all the times before we have an observations
		for (size_t i = 0; i < trials-1; i++) {
			for (int k = 0; k < N - 1; k++) {
				if (row_matrix_sample[k] == 1 && row_matrix_sample[k + 1] == 0) {
					bernoulli_distribution berd(theta);
					temp[k] = row_matrix_sample[k];
					if (temp[N - 1] != 1) {
						temp[k + 1] = berd(generator);
					}
					else { temp[k + 1] = 1; }
				}
			}
			for (int k = 1; k < N; k++) {
				if (row_matrix_sample[k] == 1 && row_matrix_sample[k - 1] == 0) {
					bernoulli_distribution berd(theta);
					temp[k] = row_matrix_sample[k];
					temp[k - 1] = berd(generator);
					if (temp[0] == 1) { break; }
				}
			}
			for (int k = 0; k < N - 1; k++) {
				row_matrix_sample[k] = temp[k];
				row_matrix_new_sample[k] = temp[k];
			}
			matrix_sample.push_back(temp);
			matrix_new_sample.push_back(temp);
		}
		//at the time of the first observation make the substitutions for the first
		//observation and for the previous times if needed
		for (int k = 0; k < N; k++) {
			if (Obs[trials - 1][k] == 1)
			{
				matrix_new_sample[trials - 1][k] = 1;
			}
		}
		if (first_observed != 0) {
			if (first_observed <= (trials - 1)){
				for (size_t i = 0; i < first_observed; i++) {
					for (size_t k = 0; k < first_observed - i +1; k++) {
						if (matrix_new_sample[trials - 1 - i][k] == 1) {
							if (matrix_new_sample[trials - 1 - i][k + 1] == 0) {
								matrix_new_sample[trials - 1 - i][k + 1] = 1;
							}
						}
					}
				}
			}
			else {
				for (size_t i = 0; i < trials - 1; i++) {
					for (size_t k = 0; k < trials - i; k++) {
						if (matrix_new_sample[trials - 1 - i][k] == 1) {
							if (matrix_new_sample[trials - 1 - i][k + 1] == 0) {
								matrix_new_sample[trials - 1 - i][k + 1] = 1;
							}
						}
					}
				}
			}
		}
		if (first_observed != N) {
			for (size_t i = 0; i < trials - 1; i++) {
				for (size_t k = N - 1; k >= first_observed + 1 + i; k--) {
					if (matrix_new_sample[trials - 1 - i][k] == 1) {
						if (matrix_new_sample[trials - 1 - i][k - 1] == 0) {
							matrix_new_sample[trials - 1 - i][k - 1] = 1;
						}
					}
				}
			}
		}
		//for all the subsequent times sample as usual
		//then make a correction
		vector < int > temp2(N, 0);
		vector < int > temp3(N, 0);
		for (int k = 0; k < N; k++){
			temp2[k] = matrix_new_sample[trials - 1][k];
			temp3[k] = matrix_new_sample[trials - 1][k];
		}
		for (size_t i = trials; i < X.size(); i++) {
			if (first_observed != N) {
				for (int k = 0; k < N - 1; k++) {
					if (matrix_sample[i - 1][k] == 1 && matrix_sample[i - 1][k + 1] == 0) {
						bernoulli_distribution berd(theta);
						temp2[k] = matrix_sample[i - 1][k];
						temp2[k + 1] = berd(generator);
					}
					if (Obs[i][k + 1] == 1) {
						temp3[k + 1] = 1;
					}
					else {
						temp3[k + 1] = temp2[k + 1];
					}
				}
			}
			if (first_observed != 0) {
				for (int k = 1; k < N; k++) {
					if (matrix_sample[i - 1][k] == 1 && matrix_sample[i - 1][k - 1] == 0) {
						bernoulli_distribution berd(theta);
						temp2[k] = matrix_sample[i - 1][k];
						temp2[k - 1] = berd(generator);
					}
					if (Obs[i][k - 1] == 1) {
						temp3[k - 1] = 1;
					}
					else {
						temp3[k - 1] = temp2[k - 1];
					}
				}
			}
			matrix_sample.push_back(temp2);
			matrix_new_sample.push_back(temp3);
		}
		sample.push_back(matrix_sample);
		new_sample.push_back(matrix_new_sample);
	}
	//simulate the first observation for each particle with a geometric distribution
    //note that we can observe only in the range where the invasion is
	//if the first observation in real life happens after or at the same time
	//as the simulated first observation, do nothing
	//otherwise simulate observations up to the first observation in real life

	//simulate first observation
	vector < int > vec_sim_trials;
	for (int j = 0; j < n; j++) {
		vector < int > obs(N, 0);
		vector < int > range;
		int sim_first_obs{ 0 };
		geometric_distribution < int > geo(phi);
		int sim_trials{ 0 };
		do {
		sim_trials = geo(generator) + 1;
		} while (sim_trials > X.size());
		vec_sim_trials.push_back(sim_trials);
		if (sim_trials < trials) {
			for (int i = 0; i < N; i++) {
				if (sample[j][sim_trials - 1][i] == 1) {
					range.push_back(i);
				}
			}
			if (range[range.size() - 1] == range[0]) {
				sim_first_obs = range[0];
			}
			else {
				sim_first_obs = rd() % (range[range.size() - 1] - range[0]) + range[0];
			}
			sam_obs[j][sim_trials - 1][sim_first_obs] = 1;
		}
	}

	//simulate the rest up to the first observation in real life
	for (int k = 0; k < n; k++) {
		//calculate left and right which are the first and last invaded cells
		if (vec_sim_trials[k] < trials) {
			for (size_t i = vec_sim_trials[k] - 1; i < trials - 1; i++) {
				int left{ 0 };
				int right{ 0 };
				for (size_t j = 0; j < N - 1; j++) {
					if (new_sample[k][i + 1][j] == 1 && new_sample[k][i + 1][j + 1] == 0) {
						right = j;
					}
					else if (new_sample[k][i + 1][N - 1] == 1) { right = N - 1; }
				}
				for (size_t j = 1; j < N; j++) {
					if (new_sample[k][i + 1][j] == 1 && new_sample[k][i + 1][j - 1] == 0) {
						left = j;
					}
					else if (new_sample[k][i + 1][0] == 1) { left = 0; }
				}
				//sample between left and right
				for (int j = 0; j < N; j++) {
					sam_obs[k][i + 1][j] = sam_obs[k][i][j];
				}
				if (left == right) {}
				else {
					for (int j = left; j < right; j++) {
						int last_inv_left{ 0 };
						if (sam_obs[k][i][j] == 0 && sam_obs[k][i][j + 1] == 1) {
							last_inv_left = j + 1;
							for (int j = last_inv_left; j >= left; j--) {
								bernoulli_distribution berd(phi);
								int l = berd(generator);
								if (l == 1) {
									sam_obs[k][i + 1][j] = l;
								}
								else { break; }
							}
						}
						int last_inv_right{ 0 };
						if (sam_obs[k][i][j] == 1 && sam_obs[k][i][j + 1] == 0) {
							last_inv_right = j;
							for (int j = last_inv_right; j < right; j++) {
								bernoulli_distribution berd(phi);
								int l = berd(generator);
								if (l == 1) {
									sam_obs[k][i + 1][j + 1] = l;
								}
								else { break; }
							}
						}
					}
				}
			}
		}
		else {
			for (int j = 0; j < vec_sim_trials[k]; j++) {
				for (int i = 0; i < N; i++) {
					sam_obs[k][j][i] = Obs[j][i];
				}
			}
		}
	}


	//Sampling the observations for every particle from a bernoulli distribution
	//with probability phi, filling the container "sam_obs".
	//Making a substitiution every time I have an observation in real life,
	//the matrix for the corrected observations is equivalent to the real life observation matrix "Obs".
	if (trials == 1) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < N; i++) {
				sam_obs[j][trials - 1][i] = Obs[trials - 1][i];
			}
		}
	}
	vector < int > temp(N, 0);
	for (int j = 0; j < N; j++) {
		temp[j] = Obs[trials - 1][j];
	}
	for (int k = 0; k < n; k++) {
		//from the time of the first observation sample next step
		//then make the substitutions
		if (vec_sim_trials[k] < trials) {
			for (int j = 0; j < N; j++) {
				sam_obs[k][trials - 1][j] = Obs[trials - 1][j];
			}
			for (int i = trials - 1; i < X.size(); i++) {
				//calculate left and right which are the first and last invaded cells
				int left{ 0 };
				int right{ 0 };
				for (size_t j = 0; j < N - 1; j++) {
					if (new_sample[k][i][j] == 1 && new_sample[k][i][j + 1] == 0) {
						right = j;
					}
					else if (new_sample[k][i][N - 1] == 1) { right = N - 1; }
				}
				for (size_t j = 1; j < N; j++) {
					if (new_sample[k][i][j] == 1 && new_sample[k][i][j - 1] == 0) {
						left = j;
					}
					else if (new_sample[k][i][0] == 1) { left = 0; }
				}
				//sample between left and right
				if (left == right) {
					for (int j = 0; j < N; j++) {
						sam_obs[k][i][j] = temp[j];
					}
				}
				else {
					for (int j = 0; j < N; j++) {
						temp[j] = Obs[i][j];
					}
					for (int j = left; j < right; j++) {
						int last_inv_left{ 0 };
						if (Obs[i][j] == 0 && Obs[i][j + 1] == 1) {
							last_inv_left = j + 1;
							for (int j = last_inv_left; j >= left; j--) {
								bernoulli_distribution berd(phi);
								int l = berd(generator);
								if (l == 1) {
									temp[j] = l;
								}
								else { break; }
							}
						}
						int last_inv_right{ 0 };
						if (Obs[i][j] == 1 && Obs[i][j + 1] == 0) {
							last_inv_right = j;
							for (int j = last_inv_right; j < right; j++) {
								bernoulli_distribution berd(phi);
								int l = berd(generator);
								if (l == 1) {
									temp[j + 1] = l;
								}
								else { break; }
							}
						}
					}
					for (int j = 0; j < N; j++) {
						sam_obs[k][i][j] = temp[j];
					}
				}
			}
		}
		else {
			for (int j = 0; j < N; j++) {
				temp[j] = Obs[vec_sim_trials[k] - 1][j];
			}
			for (int i = vec_sim_trials[k] - 1; i < X.size(); i++) {
				//calculate left and right which are the first and last invaded cells
				int left{ 0 };
				int right{ 0 };
				for (size_t j = 0; j < N - 1; j++) {
					if (new_sample[k][i][j] == 1 && new_sample[k][i][j + 1] == 0) {
						right = j;
					}
					else if (new_sample[k][i][N - 1] == 1) { right = N - 1; }
				}
				for (size_t j = 1; j < N; j++) {
					if (new_sample[k][i][j] == 1 && new_sample[k][i][j - 1] == 0) {
						left = j;
					}
					else if (new_sample[k][i][0] == 1) { left = 0; }
				}
				//sample between left and right
				for (int j = 0; j < N; j++) {
					temp[j] = Obs[i][j];
				}
				for (int j = left; j < right; j++) {
					sam_obs[k][i][j] = temp[j];
					temp[j] = Obs[i][j];
					int last_inv_left{ 0 };
					if (Obs[i][j] == 0 && Obs[i][j + 1] == 1) {
						last_inv_left = j + 1;
						for (int j = last_inv_left; j >= left; j--) {
							bernoulli_distribution berd(phi);
							int l = berd(generator);
							if (l == 1) {
								temp[j] = l;
							}
							else { break; }
						}
					}
					int last_inv_right{ 0 };
					if (Obs[i][j] == 1 && Obs[i][j + 1] == 0) {
						last_inv_right = j;
						for (int j = last_inv_right; j < right; j++) {
							bernoulli_distribution berd(phi);
							int l = berd(generator);
							if (l == 1) {
								temp[j + 1] = l;
							}
							else { break; }
						}
					}
				}
				for (int j = 0; j < N; j++) {
					sam_obs[k][i][j] = temp[j];
				}
			}
		}
		cout << "new sample[j]" << endl;
		print_matrix(new_sample[k]);
		cout << "sam_obs[j]" << endl;
		print_matrix(sam_obs[k]);
	}

	/*
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
