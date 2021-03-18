/**
    Frame analysis program: Writes the global stiffness matrix, global force vector, structural displacements
    and member end forces into "output.txt" using the data provided in "input.txt". Can be used for 2D or 3D,
    frame or trusses.
    @author Salih Deniz Unal
    @version 2.0 24/08/2020
*/

#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <string>

#include "inputLibrary.h"
#include "nodeClass.h"
#include "elementClass.h"

using namespace std;
using namespace Eigen;

namespace globalMatrices {

    /**
        Constructs the global stiffness matrix from element stiffness matrices.
        @param numDof number of degrees of freedom of each node of the structure.
        @param elementGlobalDof id of global degree of freedom of the element's node.
        @param k element stiffness matrix.
        @param K global striffnes matrix.
    */
    int constructGlobalStiff(int numDof, MatrixXd elementGlobalDof, MatrixXd k, MatrixXd& K) {

        int p, q;

        for (p = 0; p < numDof * 2; p++) {

            for (q = 0; q < numDof * 2; q++) {

                if (elementGlobalDof(p) != 0 && elementGlobalDof(q) != 0) {

                    K(elementGlobalDof(p) - 1, elementGlobalDof(q) - 1) = K(elementGlobalDof(p) - 1, elementGlobalDof(q) - 1) + k(p, q);
                }
            }
        }

        return 0;
    }

    /**
        Constructs the global load vector.
        @param i row number of loaded joint in the input load matrix L.
        @param numDof number of degrees of freedom of each node of the structure.
        @param numLoadJoint number of loaded joints in the structure.
        @param F global load vector.
        @param L input load matrix.
        @param nodeGlobalDof id of global degree of freedom of the loaded node.
    */
    int constructGlobalLoad(int i, int numDof, int numLoadJoint, VectorXd& F, MatrixXd L, MatrixXd nodeGlobalDof) {


        for (int j = 0; j < numDof; j++) {

            if (nodeGlobalDof(j) != 0) {

                F(nodeGlobalDof(j) - 1) = F(nodeGlobalDof(j) - 1) + L(i, j + 1);
            }
        }

        return 0;
    }
}

// Controls operation of the program.
int main() {

    int i, j;

    vector <double> inputVector;
    readFile(inputVector);


    // get the number of dimensions, nodes, degree of freedoms, elements, supports, loaded joints, material properties
    int inputIndex = 0;

    int numDimension = inputVector[inputIndex++];
    int numNode = inputVector[inputIndex++];
    int numDof = inputVector[inputIndex++];
    int numElem = inputVector[inputIndex++];
    int numSupport = inputVector[inputIndex++];
    int numLoadJoint = inputVector[inputIndex++];
    int numDistLoad = inputVector[inputIndex++];
    int numPointLoad = inputVector[inputIndex++];
    int numConcMoment = inputVector[inputIndex++];
    int numMaterialProp = inputVector[inputIndex++];


    // create the matrices from input file, send current input index, number of rows, number of columns, matrix and input vector using "inputLibrary.h"

    // coordinates
    MatrixXd XY;
    createInputMatrix(inputIndex, numNode, numDimension, XY, inputVector);

    // material properties
    MatrixXd M;
    createInputMatrix(inputIndex, numMaterialProp, 3, M, inputVector);

    // connectivity
    MatrixXd C;
    createInputMatrix(inputIndex, numElem, 3, C, inputVector);

    // supports
    MatrixXd S;
    createInputMatrix(inputIndex, numSupport, numDof + 1, S, inputVector);

    // loads
    MatrixXd L;
    createInputMatrix(inputIndex, numLoadJoint, numDof + 1, L, inputVector);

    // distributed load information
    MatrixXd distLoad;
    createInputMatrix(inputIndex, numDistLoad, 2, distLoad, inputVector);

    // point load information
    MatrixXd pointLoad;
    createInputMatrix(inputIndex, numPointLoad, 2, pointLoad, inputVector);

    // point load location
    MatrixXd locationPointLoad;
    createInputMatrix(inputIndex, numPointLoad, 2, locationPointLoad, inputVector);

    // concentrated moment information
    MatrixXd concMoment;
    createInputMatrix(inputIndex, numConcMoment, 2, concMoment, inputVector);

    // concentrated moment location
    MatrixXd locationConcMoment;
    createInputMatrix(inputIndex, numConcMoment, 2, locationConcMoment, inputVector);


    // equation numbers

    MatrixXd E;
    E = MatrixXd::Constant(numNode, numDof, 0);

    //Node support[numSupport];
    vector <Node> supportVector;

    // insert the restrained degrees of freedoms of supports to E

    for (i = 0; i < numSupport; i++) {
        
        Node support;

        support.getID(i, S);
        support.getIndex();


        for (j = 0; j < numDof; j++) {

            E(support.index, j) = S(i, j + 1);
        }

        supportVector.push_back(support);
    }


    // assign consecutive numbers to active degrees of freedoms to find number of equations

    int counter = 1;
    int row, col;

    for (i = 0; i < E.size(); i++) {

        row = floor(i / numDof);
        col = i % numDof;

        if (E(row, col) == 0) {

            E(row, col) = counter++;
        }
        else {

            E(row, col) = 0;
        }

    }

    int numEq = counter - 1;  // number of equations


    // elements

    //Element element[numElem];
    vector <Element> elementVector;

    for (i = 0; i < numElem; i++) {

        Element element;

        element.getMaterialType(i, C, M);
        element.getCoordinates(i, XY, C);
        element.getkPrime(numDof);
        element.getRotation(numDof);
        element.getStiffnes();
        element.getGlobalDof(numDof, E);

        element.getFEofConcMoment(numDof, i, numConcMoment, concMoment, locationConcMoment);
        element.getFEofPointLoad(numDof, i, numPointLoad, pointLoad, locationPointLoad);
        element.getFEofDistLoad(numDof, i, numDistLoad, distLoad);

        element.getGlobalFE(numEq, numDof);

        elementVector.push_back(element);
    }




    // global stiffness matrix

    MatrixXd K;
    K = MatrixXd::Constant(numEq, numEq, 0);

    for (i = 0; i < numElem; i++) {

        globalMatrices::constructGlobalStiff(numDof, elementVector[i].globalDof, elementVector[i].k, K);
    }


    // global load vector

    VectorXd F;
    F = MatrixXd::Constant(numEq, 1, 0);

    //Node loadJoint[numLoadJoint];
    vector <Node> loadJointVector;

    for (i = 0; i < numLoadJoint; i++) {
        
        Node loadJoint;

        loadJoint.getID(i, L);
        loadJoint.getIndex();
        loadJoint.getGlobalDof(numDof, E);

        loadJointVector.push_back(loadJoint);

        globalMatrices::constructGlobalLoad(i, numDof, numLoadJoint, F, L, loadJointVector[i].globalDof);
    }


    for (i = 0; i < numElem; i++) {

        F -= elementVector[i].globalFixedEndForces;

    }



    // structural displacements

    VectorXd D;
    D = K.colPivHouseholderQr().solve(F);


    // getting the output file

    ofstream outputFile;
    outputFile.open("output.txt");
    outputFile << "Global Stiffness Matrix:\n" << K << endl << "Global Force Vector:\n" << F << endl << "Structural Displacements:\n" << D << endl;


    // member end forces

    for (i = 0; i < numElem; i++) {

        elementVector[i].getEndForces(numDof, D);

        outputFile << i + 1 << ". member end forces in local coordinates:\n" << elementVector[i].f << endl;
    }


    outputFile.close();

    return 0;
}
