/**
    Node Class: Constructs a class of nodes to be used in a structure.
    @author Salih Deniz Unal
    @version 2.0 24/08/2020
*/

#ifndef NODECLASS_H_INCLUDED
#define NODECLASS_H_INCLUDED

#include <iostream>
#include <Eigen/Dense>
#include <vector>

class Node {

public:

    int id, index;

    MatrixXd globalDof;


    /**
        Gets the id.
        @param i row number of node in the input matrix.
        @param nodeInfo the input matrix for node information to be found.
    */
    void getID(int i, MatrixXd nodeInfo) {

        id = nodeInfo(i, 0);
    }

    /**
        Calculates index(row number) of node in the matrices from id.
    */
    void getIndex() {

        index = id - 1;
    }


    /**
        Constructs id of global degree of freedom of the node for truss and frame.
        @param numDof number of degrees of freedom of each node of the structure.
        @param E equation number matrix.
    */
    void getGlobalDof(int numDof, MatrixXd E) {

        globalDof = MatrixXd::Constant(numDof, 1, 0);

        for (int i = 0; i < numDof; i++) {

            globalDof(i) = E(index, i);
        }
    }
};

#endif // NODECLASS_H_INCLUDED
