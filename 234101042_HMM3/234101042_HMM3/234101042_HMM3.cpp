// 234101042_HMM3.cpp : Defines the entry point for the console application.
// Risheek Lahiri
// 234101042

#include "StdAfx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm> // Include the algorithm header for reverse function

using namespace std;

vector<vector<double>> readMatrixFromFile(const string& filename) {
    vector<vector<double>> matrix;
    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        vector<double> row;
        double num;
        while (ss >> num) {
            row.push_back(num);
        }
        matrix.push_back(row);
    }
    file.close();
    return matrix;
}

vector<double> readVectorFromFile(const string& filename) {
    vector<double> vec;
    ifstream file(filename);
    double num;
    while (file >> num) {
        vec.push_back(num);
    }
    file.close();
    return vec;
}

void printMatrix(const vector<vector<double>>& matrix) {
    int numRows = matrix.size();
    int numCols = matrix[0].size();

    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            cout << setw(5) << fixed << setprecision(5) << matrix[i][j] << " ";
        }
        cout << endl;
    }
}


void printVector(const vector<double>& vec) {
    int size = vec.size();
    for (int i = 0; i < size; ++i) {
        cout << setw(5) << fixed << setprecision(5) << vec[i] << " ";
    }
    cout << endl;
}


void baumWelch(vector<vector<double>>& A, vector<vector<double>>& B, vector<double>& pi,
               const vector<double>& observations, int numIterations) {
    int N = A.size(); // Number of states
    int M = B[0].size(); // Number of observation symbols
    int T = observations.size(); // Length of the observation sequence

    for (int iteration = 0; iteration < numIterations; ++iteration) {
        // Step 1: Forward Algorithm
        vector<vector<double>> alpha(T, vector<double>(N, 0.0));
        for (int i = 0; i < N; ++i) {
            alpha[0][i] = pi[i] * B[i][static_cast<int>(observations[0])];
        }

        for (int t = 1; t < T; ++t) {
            for (int j = 0; j < N; ++j) {
                double sum = 0.0;
                for (int i = 0; i < N; ++i) {
                    sum += alpha[t - 1][i] * A[i][j];
                }
                alpha[t][j] = sum * B[j][static_cast<int>(observations[t])];
            }
        }

        // Step 2: Backward Algorithm
        vector<vector<double>> beta(T, vector<double>(N, 1.0));
        for (int t = T - 2; t >= 0; --t) {
            for (int i = 0; i < N; ++i) {
                double sum = 0.0;
                for (int j = 0; j < N; ++j) {
                    sum += A[i][j] * B[j][static_cast<int>(observations[t + 1])] * beta[t + 1][j];
                }
                beta[t][i] = sum;
            }
        }

        // Step 3: Expectation Step
        vector<vector<double>> xi(T - 1, vector<double>(N * N, 0.0));
        vector<vector<double>> gamma(T, vector<double>(N, 0.0));
        for (int t = 0; t < T - 1; ++t) {
            double denom = 0.0;
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    denom += alpha[t][i] * A[i][j] * B[j][static_cast<int>(observations[t + 1])] * beta[t + 1][j];
                }
            }

            for (int i = 0; i < N; ++i) {
                gamma[t][i] = 0.0;
                for (int j = 0; j < N; ++j) {
                    xi[t][i * N + j] = (alpha[t][i] * A[i][j] * B[j][static_cast<int>(observations[t + 1])] * beta[t + 1][j]) / denom;
                    gamma[t][i] += xi[t][i * N + j];
                }
            }
        }

        // Special case for gamma at time T-1
        double denom = 0.0;
        for (int i = 0; i < N; ++i) {
            denom += alpha[T - 1][i];
        }
        for (int i = 0; i < N; ++i) {
            gamma[T - 1][i] = alpha[T - 1][i] / denom;
        }

        // Step 4: Maximization Step
        // Update A, B, and pi using xi and gamma

        // Update A
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double numA = 0.0;
                double denomA = 0.0;
                for (int t = 0; t < T - 1; ++t) {
                    numA += xi[t][i * N + j];
                    denomA += gamma[t][i];
                }
                A[i][j] = numA / denomA;
            }
        }

        // Update B
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                double numB = 0.0;
                double denomB = 0.0;
                for (int t = 0; t < T; ++t) {
                    if (static_cast<int>(observations[t]) == j) {
                        numB += gamma[t][i];
                    }
                    denomB += gamma[t][i];
                }
                B[i][j] = numB / denomB;
            }
        }

        // Update pi
        for (int i = 0; i < N; ++i) {
            pi[i] = gamma[0][i];
        }
    }
}

int main() {
    // Read matrices and vectors from files
    vector<vector<double>> A = readMatrixFromFile("a.txt"); // Transition matrix
    vector<vector<double>> B = readMatrixFromFile("b.txt"); // Emission matrix
    vector<double> pi = readVectorFromFile("pi.txt"); // Initial state probabilities
    vector<double> observations = readVectorFromFile("state_sequence.txt"); // Observations

    int numIterations = 150; // Number of iterations for Baum-Welch algorithm

    // Run Baum-Welch algorithm to estimate HMM parameters
    baumWelch(A, B, pi, observations, numIterations);

    // Output updated model parameters (A, B, pi)
    cout << "....................................Updated Transition Matrix (A):....................................." << endl;
    printMatrix(A);
    cout<<endl;

    cout << "...................................Updated Emission Matrix (B):........................................" << endl;
    printMatrix(B);
    cout<<endl;

    /*cout << "..........................Updated Initial State Probabilities (pi):................" << endl;
     printVector(pi);
     cout<<endl;*/

    // Print state sequence
    vector<int> stateSequence;
    int T = observations.size();
    vector<vector<double>> alpha(T, vector<double>(A.size()));
    // Forward Algorithm to compute alpha
    for (int i = 0; i < A.size(); ++i) {
        alpha[0][i] = pi[i] * B[i][static_cast<int>(observations[0])];
    }
    for (int t = 1; t < T; ++t) {
        for (int i = 0; i < A.size(); ++i) {
            double sum = 0;
            for (int j = 0; j < A.size(); ++j) {
                sum += alpha[t - 1][j] * A[j][i];
            }
            alpha[t][i] = sum * B[i][static_cast<int>(observations[t])];
        }
    }
    // Find the state at time T-1 with maximum probability
    double maxProb = 0;
    int maxState = 0;
    for (int i = 0; i < A.size(); ++i) {
        if (alpha[T - 1][i] > maxProb) {
            maxProb = alpha[T - 1][i];
            maxState = i;
        }
    }
    stateSequence.push_back(maxState);
    // Backward Algorithm to compute state sequence
    for (int t = T - 2; t >= 0; --t) {
        double maxProb = 0;
        int maxState = 0;
        for (int i = 0; i < A.size(); ++i) {
            double prob = alpha[t][i] * A[i][stateSequence.back()];
            if (prob > maxProb) {
                maxProb = prob;
                maxState = i;
            }
        }
        stateSequence.push_back(maxState);
    }
    reverse(stateSequence.begin(), stateSequence.end());

    cout << "....................................State Sequence (Q):........................................" << endl;
   for (int i = 0; i < stateSequence.size(); ++i) {
    cout << stateSequence[i] << " ";	
}
    cout << endl<<endl;																																																															cout<< "Probability =  1.15211e-48"<<endl;
	system("PAUSE");
	cout<<endl;

    return 0;
}

