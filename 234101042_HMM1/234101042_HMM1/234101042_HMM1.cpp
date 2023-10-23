// 234101042_HMM1.cpp : Defines the entry point for the console application.
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
using namespace std;

// Function to read matrix from a file
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

// Function to read vector from a file
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


// Output state sequence logic using forward and backward variables
vector<int> getStateSequence(const vector<vector<double>>& forward, const vector<vector<double>>& backward) {
    int T = forward.size();
    int N = forward[0].size();

    // Compute state probabilities at each time step
    vector<vector<double>> stateProbabilities(T, vector<double>(N));
    for (int t = 0; t < T; ++t) {
        for (int i = 0; i < N; ++i) {
            stateProbabilities[t][i] = forward[t][i] * backward[t][i];
        }
    }

    // Determine the most probable state at each time step
    vector<int> stateSequence(T);
    for (int t = 0; t < T; ++t) {
        double maxProbability = stateProbabilities[t][0];
        int maxState = 0;
        for (int i = 1; i < N; ++i) {
            if (stateProbabilities[t][i] > maxProbability) {
                maxProbability = stateProbabilities[t][i];
                maxState = i;
            }
        }
        stateSequence[t] = maxState;
    }

    return stateSequence;
}


int main() {
    // Read matrices and vectors from files
    vector<vector<double>> A = readMatrixFromFile("a.txt"); // Transition matrix
    vector<vector<double>> B = readMatrixFromFile("b.txt"); // Emission matrix
    vector<double> pi = readVectorFromFile("pi.txt"); // Initial state probabilities
    vector<double> observations = readVectorFromFile("state_sequence.txt"); // Observations

    // Number of states and observation symbols
    int N = A.size();
    int T = observations.size();

    // Forward algorithm
    vector<vector<double>> forward(T, vector<double>(N));
    for (int i = 0; i < N; ++i) {
        forward[0][i] = pi[i] * B[i][observations[0]];
    }

    for (int t = 1; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            double sum = 0.0;
            for (int i = 0; i < N; ++i) {
                sum += forward[t - 1][i] * A[i][j];
            }
            forward[t][j] = sum * B[j][observations[t]];
        }
    }

    // Backward algorithm
    vector<vector<double>> backward(T, vector<double>(N, 1.0));
    for (int t = T - 2; t >= 0; --t) {
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j) {
                sum += A[i][j] * B[j][observations[t + 1]] * backward[t + 1][j];
            }
            backward[t][i] = sum;
        }
    }

    // Calculate the probability of the observation sequence using backward variables
    double probability = 0.0;
    for (int i = 0; i < N; ++i) {
        probability += pi[i] * B[i][observations[0]] * backward[0][i];
    }

    // Output results
    cout << "...............................State Sequence:........................................... " << endl;
    // Output state sequence logic can be added here using forward and backward variables
     vector<int> stateSequence = getStateSequence(forward, backward);
    // cout << "State Sequence: ";
    for (int t = 0; t < T; ++t) {
        cout << stateSequence[t] << " ";
    }
    cout << endl<<endl;

    cout << "...................................Final Aij Matrix:......................................" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout<<endl;

    cout << "...................................Final Bij Matrix:....................................\n" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < B[i].size(); j++) {
            cout << B[i][j] << " ";
        }
        
    }
    cout<<endl<<endl;

    cout << "Probability =  " << probability << endl<<endl;
	system("PAUSE");
    return 0;
}


