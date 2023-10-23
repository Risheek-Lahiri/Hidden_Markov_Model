// 234101042_HMM2.cpp : Defines the entry point for the console application.
// Risheek Lahiri
// 234101042

#include "stdafx.h"
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

int main() {
    // Read matrices and vectors from files
    vector<vector<double>> A = readMatrixFromFile("a.txt"); // Transition matrix
    vector<vector<double>> B = readMatrixFromFile("b.txt"); // Emission matrix
    vector<double> pi = readVectorFromFile("pi.txt"); // Initial state probabilities
    vector<double> observations = readVectorFromFile("state_sequence.txt"); // Observations

    // Number of states and observation symbols
    int N = A.size();
    int T = observations.size();

    // Viterbi algorithm
    vector<vector<double>> viterbi(T, vector<double>(N));
    vector<vector<int>> path(T, vector<int>(N));

    // Initialization
    for (int i = 0; i < N; ++i) {
        viterbi[0][i] = pi[i] * B[i][observations[0]];
        path[0][i] = i;
    }

    // Recursion step
    for (int t = 1; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            double maxProbability = 0.0;
            int bestPrevState = 0;
            for (int i = 0; i < N; ++i) {
                double probability = viterbi[t - 1][i] * A[i][j] * B[j][observations[t]];
                if (probability > maxProbability) {
                    maxProbability = probability;
                    bestPrevState = i;
                }
            }
            viterbi[t][j] = maxProbability;
            path[t][j] = bestPrevState;
        }
    }

    // Termination
    double maxProbability = 0.0;
    int bestLastState = 0;
    for (int i = 0; i < N; ++i) {
        if (viterbi[T - 1][i] > maxProbability) {
            maxProbability = viterbi[T - 1][i];
            bestLastState = i;
        }
    }

    // Backtracking to find the most likely state sequence
    vector<int> stateSequence(T);
    stateSequence[T - 1] = bestLastState;
    for (int t = T - 2; t >= 0; --t) {
        stateSequence[t] = path[t + 1][stateSequence[t + 1]];
    }

    // Output state sequence and probability of the observation sequence
    cout << "..................................... State Sequence:...................................." << endl;
    for (int t = 0; t < T; ++t) {
        cout << stateSequence[t] << " ";
    }
    cout << endl<<endl;

	 cout << "............................Final Aij Matrix:........................." << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout<<endl;

    cout << "..........................Final Bij Matrix:..........................\n" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < B[i].size(); j++) {
            cout << B[i][j] << " ";
        }
        
    }
    cout<<endl<<endl;

    cout << "Probability of Observation Sequence = " << maxProbability << endl;

	system("PAUSE");
    return 0;
}
