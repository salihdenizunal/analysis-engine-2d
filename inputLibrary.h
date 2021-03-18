/**
    Input Library: Reads and changes the input file to be used in the program, creates input matrices.
    @author Salih Deniz Unal
    @version 2.0 24/08/2020
*/

#ifndef INPUTLIBRARY_H_INCLUDED
#define INPUTLIBRARY_H_INCLUDED

#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
using namespace Eigen;

/**
    Reads the input file "input.txt", creates a numerical input file "numericInput.txt". Constructs a input vector.
    @param inputVector input vector of the data provided in "input.txt".
*/
int readFile(vector <double>& inputVector) {

    double input = 0.0;
    string line;

    ifstream inputFile("input.txt", ios::in);

    if (inputFile.is_open()) {

        ofstream numericInputFile("numericInput.txt");


        while (!inputFile.eof()) {

            getline(inputFile, line);

            if (line[0] != '#') {
                numericInputFile << line << "\n";
            }
        }

        inputFile.close();
        numericInputFile.close();

        ifstream data("numericInput.txt", ios::in);


        while (data >> input) {

            inputVector.push_back(input);
        }

        data.close();
    }
    else {

        cout << "Unable to open input file!" << endl;
    }

    return 0;
}

/**
    Create input matrices.
    @param inputIndex start index of the matrix in the input vector.
    @param numRow number of rows of the matrix.
    @param numCol number of columns of the matrix.
    @param m created input matrix.
    @param inputVector input vector of the data provided in "input.txt".
*/
int createInputMatrix(int& inputIndex, int numRow, int numCol, MatrixXd& m, vector <double> inputVector) {

    int i, j;

    m.resize(numRow, numCol);


    for (i = 0; i < numRow; i++) {

        for (j = 0; j < numCol; j++) {

            m(i, j) = inputVector[inputIndex++];
        }
    }

    return 0;
}


#endif // INPUTLIBRARY_H_INCLUDED
