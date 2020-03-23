#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_river.h"


using namespace std;

//This is the code for the method:
//Firstly I calculate the 3-dimentional matrix called for the sampled events
//making the substitutions for every particle.
//I then calculate the 3-dimensional matrix for the sampled observations
//and for the sampled observations with substitutions
//Finally I calculate the weights and resample.

//FUNCTIONS

void print_matrix(vector < vector < int > > m);
void print_vector(vector < int > v);
int left(vector < int > v);
int right(vector < int > v);
int function_h(int k, int l, int L);
//bool is_zero(vector < int > v);



int main() {

	random_device rd;
	mt19937 generator(rd());

	simprobe();

	//DEFINITIONS

	//number of particles
	int n = 4;
	int m = X.size();
	//define the container for the samples 
	vector < vector < vector < int > > > sample;
	//define the container for the corrected samples
	vector < vector < vector < int > > > new_sample;
	//define the empty matrix for the sampled observations
	vector < vector < vector < int > > > sam_obs(n, vector < vector < int > >(m, vector < int >(N)));
	//define the container for the resampling
	vector < vector < vector < int > > > resampled;
	//define the containter for the unnormalised weights
	vector < vector < vector < int > > > un_weights;
	//define the container for the normalised weights
	vector < vector < vector < int > > > weights;
	//define the number of particles

	print_matrix(X);
	print_matrix(Obs);
	
	//simulate the first invaded cell for all particles
	vector < vector < int > > samplem;
	vector < int > first_sim_vect;
	//for each particle, define the cells where there is an invasion 
	//and randomly select the first observation in that range
	//int end_of_range = right_first_observed + (trials - 1);
	//int start_of_range = left_first_observed - (trials - 1);
	//if (start_of_range < 0) { start_of_range = 0; }
	//else if (end_of_range > N - 1) { end_of_range = N - 1; }
	//else {}
	for (int j = 0; j < n; j++) {
		vector < int > samplev(N, 0);
		uniform_int_distribution<> unif(0, N - 1);
		int first_simulated = unif(generator);
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
		//make all the simulations for all the times before we have observations
		for (size_t i = 0; i < trials - 1; i++) {
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
		int left_first_observed = left(Obs[trials - 1]);
		int right_first_observed = Obs[trials - 1].size() - right(Obs[trials - 1]) - 1;
		if (first_sim_vect[j] < right_first_observed) {
			if (right_first_observed != 0) {
				for (size_t i = 0; i < trials; i++) {
					for (size_t k = first_sim_vect[j] + 1; k < right_first_observed - i + 1; k++) {
						matrix_new_sample[trials - 1 - i][k] = 1;
					}
				}
			}
		}
		if (first_sim_vect[j] > left_first_observed) {
			if (left_first_observed != N) {
				for (size_t i = 0; i < trials; i++) {
					for (size_t k = left_first_observed + i; k < first_sim_vect[j] + 1; k++) {
						matrix_new_sample[trials - 1 - i][k] = 1;
					}
				}
			}
		}
		//for all the subsequent times sample as usual
		//then make a correction and put the corrected terms in "new_sample"
		vector < int > temp2(N, 0);
		vector < int > temp3(N, 0);
		for (int k = 0; k < N; k++){
			temp2[k] = matrix_new_sample[trials - 1][k];
			temp3[k] = matrix_new_sample[trials - 1][k];
		}
		for (size_t i = trials; i < X.size(); i++) {
			if (first_sim_vect[j] != N) {
				for (int k = 0; k < N - 1; k++) {
					if (matrix_new_sample[i - 1][k] == 1 && matrix_new_sample[i - 1][k + 1] == 0) {
						bernoulli_distribution berd(theta);
						temp2[k] = 1;
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
			if (first_sim_vect[j] != 0) {
				for (int k = 1; k < N; k++) {
					if (matrix_new_sample[i - 1][k] == 1 && matrix_new_sample[i - 1][k - 1] == 0) {
						bernoulli_distribution berd(theta);
						temp2[k] = 1;
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
			//this ensures that there are no 0s between 1s
			for (vector < int >::iterator i = temp3.begin(); i != temp3.end() - 1; i++) {
				if (*i == 1 && *(i + 1) == 0){
					for (vector < int >::iterator j = temp3.end() - 1; j != i; j--) {
						if (*j == 1 && *(j - 1) == 0) {
							*(j - 1) = 1;
						}
					}
				}
			}
			for (int k = 0; k < N; k++) {
				temp2[k] = temp3[k];
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
		int lf = left(sample[j][sim_trials - 1]);
		int rg = N - right(sample[j][sim_trials - 1]) - 1;
		if (lf == rg) {
			sim_first_obs = lf;
		}
		else {
			uniform_int_distribution<> unif(lf, rg);
			sim_first_obs = unif(generator);
		}
		sam_obs[j][sim_trials - 1][sim_first_obs] = 1;
		//simulate all other observations making a correction each time
		for (int k = sim_trials - 1; k < X.size(); k++) {
			vector < int > corrected(N, 0);
			for (int i = 0; i < N; i++) {
				corrected[i] = Obs[k][i];
			}
			int l = left(sample[j][k]);
			int r = N - right(sample[j][k]) - 1;
			bool all_0 = none_of(corrected.begin(), corrected.end(), [](int x) {return x == 1; });
			if (all_0 == 1) {
				int sim_another_1{ 0 };
				if (l == r) {
					sim_another_1 = l;
				}
				else {
					uniform_int_distribution<> unif(l, r);
					sim_another_1 = unif(generator);
				}
				sam_obs[j][k][sim_another_1] = 1;
			}
			else {
				int last_observed_left = left(corrected);
				for (int i = last_observed_left; i >= l; i--) {
					bernoulli_distribution berd(phi);
					int next_left = berd(generator);
					if (next_left == 1) {
						corrected[i] = 1;
					}
					else { break; }
				}
				int last_observed_right = N - right(corrected) - 1;
				for (int i = last_observed_right; i < r + 1; i++) {
					bernoulli_distribution berd(phi);
					int next_right = berd(generator);
					if (next_right == 1) {
						corrected[i] = 1;
					}
					else { break; }
				}
				for (int i = 0; i < N; i++) {
					sam_obs[j][k][i] = corrected[i];
				}
			}
		}
	} 

	//Finding the unnormalised weights (using log then exponentiating)
	//filling the container "unnorm_weights".
	/*
	vector < int > unnorm_weights{1/n};
	for (int k = 0; k < n; k++) {
		for (int i = 1; i < X.size(); i++) {
			//these are all the variable used to calculate h(k_t) primed
			int num_prev_invaded = count(sample[0][i-1].begin(), sample[0][i-1].end(), 1);
			int num_now_invaded = count(sample[0][i].begin(), sample[0][i].end(), 1);
			int num_prev_observed = count(sam_obs[k][i-1].begin(), sam_obs[k][i-1].end(), 1);
			int num_now_observed = count(sam_obs[k][i].begin(), sam_obs[k][i].end(), 1) - num_prev_observed;
			int alpha = left_unobserved(sam_obs[k][i-1]);
			int beta = right_unobserved(sam_obs[k][i-1]);
			int minimum = min(alpha, beta);
			int maximum = max(alpha, beta);
			//these are all the variable used to calculate h(k_t) primed
			int num_prev_invaded_primed = count(new_sample[k][i-1].begin(), new_sample[0][i-1].end(), 1);
			int num_now_invaded_primed = count(new_sample[k][i].begin(), new_sample[0][i].end(), 1);
			int num_prev_observed_primed = count(Obs[i-1].begin(), Obs[i-1].end(), 1);
			int num_now_observed_primed = count(Obs[i].begin(), Obs[i].end(), 1) - num_prev_observed_primed;
			int alpha_primed = left_unobserved(Obs[i-1]);
			int beta_primed = right_unobserved(Obs[i-1]);
			int minimum_primed = min(alpha_primed, beta_primed);
			int maximum_primed = max(alpha_primed, beta_primed);
			//calculate h(k_t)
			int hk{ 0 };
			if (i <= vec_sim_trials[k] - 2) {}
			else {
				hk = function_h(num_now_observed, minimum, maximum);
			}
			//cout << "alpha " << alpha << endl;
			//cout << "hk " << hk << endl;
		}
	}
	*/
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

//this function calculates the number of 0 on the left before the first 1
int left(vector < int > v) {
	vector < int >::iterator it = find(v.begin(), v.end(), 1);
	int result = distance(v.begin(), it);
	return result;
}

//this function calculates the number of 0 on the right before the first 1
int right(vector < int > v) {
	vector < int >::reverse_iterator it = find(v.rbegin(), v.rend(), 1);
	int result = distance(v.rbegin(), it);
	return result;
}

//this function calculates the discrete function h(k_t)
int function_h(int k, int l, int L) {
	int h{ 0 };
	if (k <= l - 1) {
		h = k + 1;
	}
	else if (k >= l || k <= L) {
		h = l + 1;
	}
	else {
		h = l - k + L + 1;
	}
	return h;
}


