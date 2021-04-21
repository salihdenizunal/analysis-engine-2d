#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <vector>

class MembraneElement {

public:

	int materialType;
	int node1, node2, node3, node4;

	double x1, y1, x2, y2, x3, y3, x4, y4;
	double len, depth, t;

	double E,mu;

	MatrixXd kPrime;
	MatrixXd R;
	MatrixXd k;
	MatrixXd globalDof;

	VectorXd d, dPrime, f;

	MatrixXd fixedEndForces;
	MatrixXd rotatedFixedEndForces;
	MatrixXd globalFixedEndForces;


	/**
		Gets the material type, calculates EA and EI.
		@param i row number of element in the input connectivity and material type matrix C.
		@param C the input connectivity and material type matrix.
		@param M the input material parameters matrix.
	*/
	void getMaterialType(int i, MatrixXd C, MatrixXd M) {

		materialType = C(i, 4);

		t = M(materialType - 1, 0);
		E = M(materialType - 1, 1);
		mu = M(materialType - 1, 2);
	}


	/**
		Gets the coordinates, calculates length.
		@param i row number of element in the input connectivity and material type matrix C.
		@param XY the input coordinate matrix of nodes.
		@param C the input connectivity and material type matrix.
	*/
	void getCoordinates(int i, MatrixXd XY, MatrixXd C) {

		node1 = C(i, 0);
		node2 = C(i, 1);
		node3 = C(i, 2);
		node4 = C(i, 3);

		x1 = XY(node1 - 1, 0);
		y1 = XY(node1- 1, 1);

		x2 = XY(node2 - 1, 0);
		y2 = XY(node2 - 1, 1);

		x3 = XY(node3 - 1, 0);
		y3 = XY(node3 - 1, 1);

		x4 = XY(node4 - 1, 0);
		y4 = XY(node4 - 1, 1);

		len = x2 - x1;
		depth = y3 - y2;
	}


	/**
		Constructs k' for truss and frame.
		@param numDof number of degrees of freedom of each node of the structure.
	*/
	void getkPrime(int numDof) {

		kPrime = MatrixXd::Constant(4 * numDof, 4 * numDof, 0);

		switch (numDof) {

		case 2:

			kPrime << -1.0*t*(0.25*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), t*(0.5*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), 1.0*depth*len*t*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))), t*(0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.5*depth*len*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))),
				-1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.5*depth*len*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))), t*(0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), 1.0*depth*len*t*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))), -1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), t*(0.5*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))),
				t*(0.5*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.5*depth*len*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), 1.0*depth*len*t*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))), -1.0*t*(0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))),
				-1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.5*depth*len*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), t*(0.5*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))), -1.0*t*(0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), 1.0*depth*len*t*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))),
				1.0*depth*len*t*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))), t*(0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.5*depth*len*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), t*(0.5*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))),
				t*(0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), 1.0*depth*len*t*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))), -1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), t*(0.5*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.5*depth*len*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))),
				-1.0*t*(0.5*depth*len*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), 1.0*depth*len*t*((0.16666666667*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))), -1.0*t*(0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), t*(0.5*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(len,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(depth,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))),
				t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), t*(0.5*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))), -1.0*t*(0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), 1.0*depth*len*t*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))), -1.0*t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) + (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.5*depth*len*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) + (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0)))), t*(0.25*depth*len*((0.62200846793*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.5*depth*len*((0.16666666667*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.16666666667*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E*mu) / (depth*len*(pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (depth*len*(pow(mu,2) - 1.0)))), -1.0*t*(0.25*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.62200846793*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.62200846793*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))) + 0.25*depth*len*((0.044658198739*E) / (pow(depth,2) * (pow(mu,2) - 1.0)) - (0.044658198739*E*(0.5*mu - 0.5)) / (pow(len,2) * (pow(mu,2) - 1.0))));


			break;

		case 3:

			break;
		}
	}


	/**
		Constructs rotation matrix R for truss and frame.
		@param numDof number of degrees of freedom of each node of the structure.
	*/
	void getRotation(int numDof) {

		R = MatrixXd::Constant(2 * numDof, 2 * numDof, 0);

	}


	/**
		Constructs element stiffness matrix for truss and frame.
	*/
	void getStiffnes() {

		//k = R.transpose() * kPrime * R;
		k = kPrime;
	}


	/**
		Constructs id of global degree of freedom of the element's node for truss and frame.
		@param numDof number of degrees of freedom of each node of the structure.
		@param E equation number matrix.
	*/
	void getGlobalDof(int numDof, MatrixXd E) {

		globalDof = MatrixXd::Constant(2 * 2 * numDof, 1, 0);
			
		globalDof(0) = E(node1 - 1, 0);
		globalDof(1) = E(node1 - 1, 1);

		globalDof(2) = E(node2 - 1, 0);
		globalDof(3) = E(node2 - 1, 1);

		globalDof(4) = E(node3 - 1, 0);
		globalDof(5) = E(node3 - 1, 1);

		globalDof(6) = E(node4 - 1, 0);
		globalDof(7) = E(node4 - 1, 1);

	}


	/**
		Calculates end forces of the element for truss and frame.
		@param numDof number of degrees of freedom of each node of the structure.
		@param D global displacement vector.
	*/
	void getEndForces(int numDof, VectorXd D) {


	}


	/**
		Calculates fixed end forces due to the concentrated moment on the element of a frame.
		@param numDof number of degrees of freedom of each node of the structure.
		@param elementIndex index of the element (element id - 1)
		@param numConcMoment number of concentrated moments on the structure
		@param concMoment concentrated moments' information (element id and magnitude)
		@param locationConcMoment concentrated moments' location
	*/
	void getFEofConcMoment(int numDof, int elementIndex, int numConcMoment, MatrixXd concMoment, MatrixXd locationConcMoment) {


	}


	/**
		Calculates fixed end forces due to the point load on the element of a frame.
		@param numDof number of degrees of freedom of each node of the structure.
		@param elementIndex index of the element (element id - 1)
		@param numPointLoad number of point loads on the structure
		@param pointLoad point loads' information (element id and magnitude)
		@param locationPointLoad point loads' location
	*/
	void getFEofPointLoad(int numDof, int elementIndex, int numPointLoad, MatrixXd pointLoad, MatrixXd locationPointLoad) {

		
	}


	/**
		Calculates fixed end forces due to the point load on the element of a frame.
		@param numDof number of degrees of freedom of each node of the structure.
		@param elementIndex index of the element (element id - 1)
		@param numDistLoad
		@param pointLoad point loads' information (element id and magnitude)
	*/
	void getFEofDistLoad(int numDof, int elementIndex, int numDistLoad, MatrixXd distLoad) {

	}


	/**
		Changes fixed end forces of the frame element to the global degrees of freedom locations.
		@param numDof number of degrees of freedom of each node of the structure.
		@param numEq number of active degree of freedoms in the structure (number of equations)
	*/
	void getGlobalFE(int numEq, int numDof) {

	}
};


