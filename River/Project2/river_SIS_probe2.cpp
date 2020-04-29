#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_river2.h"


using namespace std;

//This is the code for the method:
//Firstly I calculate the 3-dimentional matrix called for the sampled events
//making the substitutions for every particle.
//I then calculate the 3-dimensional matrix for the sampled observations
//and for the sampled observations with substitutions
//Finally I calculate the weights and resample.

//FUNCTIONS

void print_matrix(vector < vector < int > > m);
void print_matrix(vector < vector < double > > m);
void print_vector(vector < int > v);
void print_vector(vector < double > v);
int left(vector < int > v);
int right(vector < int > v);


int main() {

	random_device rd;
	mt19937 generator(rd());

	simprobe();

	//DEFINITIONS

	//number of particles
	int n = 1000;
	//define the container for the samples 
	vector < vector < vector < int > > > sample;
	//define the container for the corrected samples
	vector < vector < vector < int > > > new_sample;
	//define the empty matrix for the sampled observations
	vector < vector < vector < int > > > sam_obs(n, vector < vector < int > >((X.size()), vector < int >(N)));
	//define the container for the weights
	vector < vector < double > > weights((X.size()), vector < double >(n, 1));
	//define the container for the resampling
	vector < vector < vector < int > > > resampled(n, vector < vector < int > >((X.size()), vector < int >(N)));

	//first observation
	vector < vector < int > > samplem;
	samplem.push_back(X[0]);

	cout << "invasion" << endl;
	print_matrix(X);
	cout << "observations" << endl;
	print_matrix(Obs);

	//Sampling the invasion for every particle from a Bernoulli distribution
	//with probability theta, filling the container "sample".
	//Making a substitiution every time I have an observation in real life,
	//filling the container for the new updated events "new_sample".
	for (int j = 0; j < n; j++) {
		vector < vector < int > > matrix_sample;
		matrix_sample.push_back(samplem[0]);
		vector < vector < int > > matrix_new_sample;
		matrix_new_sample.push_back(samplem[0]);
		vector < int > temp2(N, 0);
		vector < int > temp3(N, 0);
		for (int k = 0; k < N; k++) {
			temp2[k] = matrix_new_sample[0][k];
			temp3[k] = matrix_new_sample[0][k];
		}
		for (size_t i = 1; i < X.size(); i++) {
			if (first_invaded != N-1) {
				for (int k = 0; k < N - 1; k++) {
					if (matrix_new_sample[i-1][k] == 1 && matrix_new_sample[i-1][k + 1] == 0) {
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
			if (first_invaded != 0) {
				for (int k = 1; k < N; k++) {
					if (matrix_new_sample[i-1][k] == 1 && matrix_new_sample[i-1][k - 1] == 0) {
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
				if (*i == 1 && *(i + 1) == 0) {
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

	//simulate the observations
	vector < int > vec_sim_trials;
	for (int j = 0; j < n; j++) {
		sam_obs[j][0][first_invaded] = 1;
		//simulate all other observations making a correction each time
		for (size_t k = 0; k < X.size(); k++) {
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
				for (int i = last_observed_left - 1; i >= l; i--) {
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

	for (int k = 0; k < n; k++) {
		resampled[k][0] = new_sample[k][0];
	}
	//Finding the unnormalised weights and resampling
	vector < double > vector_weights;
	for (int k = 0; k < n; k++) {
		vector < double > sum_logs((X.size()), 0);
		vector < double > log_obs_left((X.size()), 0);
		vector < double > log_obs_right((X.size()), 0);
		vector < double > log_obs_left_corr((X.size()), 0);
		vector < double > log_obs_right_corr((X.size()), 0);
		vector < double > log_invasion_left((X.size()), 0);
		vector < double > log_invasion_right((X.size()), 0);
		vector < double > log_invasion_left_corr((X.size()), 0);
		vector < double > log_invasion_right_corr((X.size()), 0);
		vector < double > log_weights((X.size()), 0);
		vector_weights.push_back(1 / N);
		for (int i = 1; i < X.size(); i++) {;
			//log of detections on the left
			int left_exp = left(Obs[i - 1]) - left(sam_obs[k][i]);
			int left_exp_corr = left(Obs[i - 1]) - left(Obs[i]);
			if (left_exp == left_exp_corr) {
				log_obs_left[i] = 1;	
			}
			else {
				double phiexpleft{ 1 };
				double phiexpleft_corr{ 1 };
				for (int j = 0; j < left_exp; j++) {
					phiexpleft = phiexpleft * phi;
				}
				for (int j = 0; j < left_exp_corr; j++) {
					phiexpleft_corr = phiexpleft_corr * phi;
				}
				double log_phiexpleft = log(phiexpleft);
				double log_phiexpleft_corr = log(phiexpleft_corr);
				if (left(sam_obs[k][i]) != 0) {
					log_obs_left[i] = log_phiexpleft + log(1 - phi);
				}
				else {
					log_obs_left[i] = log_phiexpleft;
				}
				if (left(Obs[i]) != 0) {
					log_obs_left_corr[i] = log_phiexpleft_corr + log(1 - phi);
				}
				else {
					log_obs_left_corr[i] = log_phiexpleft_corr;
				}
			}
			//log of detections on the right
			int right_exp = right(Obs[i - 1]) - right(sam_obs[k][i]);
			int right_exp_corr = right(Obs[i - 1]) - right(Obs[i]);
			if (left_exp == right_exp_corr) {
				log_obs_right[i] = 1;
			}
			else {
				double phiexpright{ 1 };
				double phiexpright_corr{ 1 };
				for (int j = 0; j < right_exp; j++) {
					phiexpright = phiexpright * phi;
				}
				for (int j = 0; j < right_exp_corr; j++) {
					phiexpright_corr = phiexpright_corr * phi;
				}
				double log_phiexpright = log(phiexpright);
				double log_phiexpright_corr = log(phiexpright_corr);
				if (right(sam_obs[k][i]) != N - 1) {
					log_obs_right[i] = log_phiexpright + log(1 - phi);
				}
				else {
					log_obs_right[i] = log_phiexpright;
				}
				if (right(Obs[i]) != N - 1) {
					log_obs_right_corr[i] = log_phiexpright_corr + log(1 - phi);
				}
				else {
					log_obs_right_corr[i] = log_phiexpright_corr;
				}
			}
			double log_obs = log_obs_left[i] - log_obs_left_corr[i] + log_obs_right[i] - log_obs_right_corr[i];
			//log of expansion on the left
			int expk = left(sample[k][i - 1]) - left(sample[k][i]);
			int expk_corr = left(new_sample[k][i - 1]) - left(new_sample[k][i]);
			if (expk == expk_corr) {
				log_invasion_left[i] = 1;
			}
			else {
				if (left(sample[k][i - 1]) != 0 && left(sample[k][i]) != 0) {
					if (expk == 0) {
						log_invasion_left[i] = log(1 - theta);
					}
					else {
						log_invasion_left[i] = log(theta);
					}
				}
				else {}
				if (left(new_sample[k][i - 1]) != 0 && left(new_sample[k][i]) != 0) {
					if (expk == 0) {
						log_invasion_left_corr[i] = log(1 - theta);
					}
					else {
						log_invasion_left_corr[i] = log(theta);
					}
				}
				else {}
			}
			//log of expansion on the right
			int exph = right(sample[k][i]) - right(sample[k][i - 1]);
			int exph_corr = right(new_sample[k][i - 1]) - right(new_sample[k][i]);
			if (exph == exph_corr) {
				log_invasion_left[i] = 1;
			}
			else {
				if (right(sample[k][i - 1]) != N - 1 && left(sample[k][i]) != N - 1) {
					if (exph == 0) {
						log_invasion_right[i] = log(1 - theta);
					}
					else {
						log_invasion_right[i] = log(theta);
					}
				}
				else {}
				if (right(new_sample[k][i - 1]) != N - 1 && left(new_sample[k][i]) != N - 1) {
					if (exph == 0) {
						log_invasion_right_corr[i] = log(1 - theta);
					}
					else {
						log_invasion_right_corr[i] = log(theta);
					}
				}
				else {}
			}
			double log_inv = log_invasion_left[i] - log_invasion_left_corr[i] + log_invasion_right[i] - log_invasion_right_corr[i];
			sum_logs[i] = log_obs + log_inv;
			double w = exp(sum_logs[i]);
			weights[i][k] = w;
		}
	}
	for (int i = 0; i < X.size(); i++) {
		for (int k = 0; k < n; k++) {
			discrete_distribution < int > discrete(weights[i].begin(), weights[i].end());
			resampled[k][i] = new_sample[discrete(generator)][i];
		}
	}

	//calculate the expectations for each cell at each time
	vector < vector < double > > expectations((X.size()), vector < double >(N, 0));
	for (int k = 0; k < X.size(); k++) {
		for (int i = 0; i < N; i++) {
			if (Obs[k][i] == 0) {
				for (int j = 0; j < n; j++) {
					if (resampled[j][k][i] == 1) {
						expectations[k][i] = expectations[k][i] + 1;
					}
				}
				expectations[k][i] = expectations[k][i] / n;
			}
			else { expectations[k][i] = 1; }
		}
	}
	print_matrix(expectations);

	//Create a .csv file with the resampled particles
	ofstream outFile3("./resampled.csv");
	outFile3 << endl;
	int half = round(X.size() / 2);
	for (size_t k = 0; k < n; k++) {
		for (size_t col = 0; col < N - 1; col++) {
			outFile3 << resampled[k][half][col] << ",";
		}
		outFile3 << resampled[k][half][N-1];
		outFile3 << endl;
	}
	outFile3.close();

	//Create a .csv file with the resampled particles
	ofstream outFile4("./expectations.csv");
	outFile4 << endl;
	for (size_t k = 0; k < X.size(); k++) {
		for (size_t col = 0; col < N - 1; col++) {
			outFile4 << expectations[k][col] << ",";
		}
		outFile4 << expectations[k][N - 1];
		outFile4 << endl;
	}
	outFile4.close();
	
	return 0;


}

//functions definitions

//this function prints a vector of integers
void print_vector(vector < int > v) {
	for (const int x : v) cout << x << ' ';
	cout << endl;
}

//this function prints a vector of double
void print_vector(vector < double > v) {
	for (const double x : v) cout << x << ' ';
	cout << endl;
}

//this function prints a matrix of integers
void print_matrix(vector < vector < int > > m) {
	for (const vector < int > v : m) {
		for (int x : v) cout << x << ' ';
		cout << endl;
	}
}

//this function prints a matrix of double
void print_matrix(vector < vector < double > > m) {
	for (const vector < double > v : m) {
		for (double x : v) cout << x << ' ';
		cout << endl;
	}
}

//this function calculates the position of the first 1 on the left
int left(vector < int > v) {
	vector < int >::iterator it = find(v.begin(), v.end(), 1);
	int result = distance(v.begin(), it);
	return result;
}

//this function calculates the position of first 1 on the right
int right(vector < int > v) {
	vector < int >::reverse_iterator it = find(v.rbegin(), v.rend(), 1);
	int result = distance(v.rbegin(), it);
	return result;
}



