/**
    Element Class: Constructs a class of elements to be used in a structure.
    @author Salih Deniz Unal
    @version 2.2 06/09/2020
*/

#ifndef ELEMENTCLASS_H_INCLUDED
#define ELEMENTCLASS_H_INCLUDED

#include <iostream>
#include <Eigen/Dense>
#include <vector>

class Element {

public:

    int materialType;
    int startNode, endNode;

    double xStart, yStart, xEnd, yEnd;
    double length;
    double cosTeta, sinTeta;

    double EA, EI;

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

        materialType = C(i, 2);

        EA = M(materialType - 1, 2) * M(materialType - 1, 0);
        EI = M(materialType - 1, 2) * M(materialType - 1, 1);
    }


    /**
        Gets the coordinates, calculates length.
        @param i row number of element in the input connectivity and material type matrix C.
        @param XY the input coordinate matrix of nodes.
        @param C the input connectivity and material type matrix.
    */
    void getCoordinates(int i, MatrixXd XY, MatrixXd C) {

        startNode = C(i, 0);
        endNode = C(i, 1);

        xStart = XY(startNode - 1, 0);
        xEnd = XY(endNode - 1, 0);
        yStart = XY(startNode - 1, 1);
        yEnd = XY(endNode - 1, 1);

        length = sqrt(pow(xEnd - xStart, 2) + pow(yEnd - yStart, 2));
    }


    /**
        Constructs k' for truss and frame.
        @param numDof number of degrees of freedom of each node of the structure.
    */
    void getkPrime(int numDof) {

        kPrime = MatrixXd::Constant(2 * numDof, 2 * numDof, 0);

        switch (numDof) {

        case 2:

            kPrime << EA / length, 0, -EA / length, 0,
                0, 0, 0, 0,
                -EA / length, 0, EA / length, 0,
                0, 0, 0, 0;

            break;

        case 3:

            kPrime << EA / length, 0, 0, -EA / length, 0, 0,
                0, (12 * EI) / pow(length, 3), (6 * EI) / pow(length, 2), 0, -(12 * EI) / pow(length, 3), (6 * EI) / pow(length, 2),
                0, (6 * EI) / pow(length, 2), (4 * EI) / length, 0, -(6 * EI) / pow(length, 2), (2 * EI) / length,
                -EA / length, 0, 0, EA / length, 0, 0,
                0, -(12 * EI) / pow(length, 3), -(6 * EI) / pow(length, 2), 0, (12 * EI) / pow(length, 3), -(6 * EI) / pow(length, 2),
                0, (6 * EI) / pow(length, 2), (2 * EI) / length, 0, -(6 * EI) / pow(length, 2), (4 * EI) / length;

            break;
        }
    }


    /**
        Constructs rotation matrix R for truss and frame.
        @param numDof number of degrees of freedom of each node of the structure.
    */
    void getRotation(int numDof) {

        R = MatrixXd::Constant(2 * numDof, 2 * numDof, 0);

        cosTeta = (xEnd - xStart) / length;
        sinTeta = (yEnd - yStart) / length;

        switch (numDof) {

        case 2:

            R << cosTeta, sinTeta, 0, 0,
                -sinTeta, cosTeta, 0, 0,
                0, 0, cosTeta, sinTeta,
                0, 0, -sinTeta, cosTeta;

            break;

        case 3:

            R << cosTeta, sinTeta, 0, 0, 0, 0,
                -sinTeta, cosTeta, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0,
                0, 0, 0, cosTeta, sinTeta, 0,
                0, 0, 0, -sinTeta, cosTeta, 0,
                0, 0, 0, 0, 0, 1;

            break;
        }
    }


    /**
        Constructs element stiffness matrix for truss and frame.
    */
    void getStiffnes() {

        k = R.transpose() * kPrime * R;
    }


    /**
        Constructs id of global degree of freedom of the element's node for truss and frame.
        @param numDof number of degrees of freedom of each node of the structure.
        @param E equation number matrix.
    */
    void getGlobalDof(int numDof, MatrixXd E) {

        globalDof = MatrixXd::Constant(2 * numDof, 1, 0);

        for (int i = 0; i < numDof * 2; i++) {

            if (i < numDof) {

                globalDof(i) = E(startNode - 1, i);
            }
            else {

                int j = i % numDof;
                globalDof(i) = E(endNode - 1, j);
            }
        }
    }


    /**
        Calculates end forces of the element for truss and frame.
        @param numDof number of degrees of freedom of each node of the structure.
        @param D global displacement vector.
    */
    void getEndForces(int numDof, VectorXd D) {

        d = MatrixXd::Constant(numDof * 2, 1, 0);
        dPrime = MatrixXd::Constant(numDof * 2, 1, 0);
        f = MatrixXd::Constant(numDof * 2, 1, 0);

        for (int i = 0; i < numDof * 2; i++) {

            if (globalDof(i) != 0) {

                d(i) = D(globalDof(i) - 1);
            }
        }

        dPrime = R * d;
        f = kPrime * dPrime + fixedEndForces;

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

        fixedEndForces = MatrixXd::Constant(2 * numDof, 1, 0);

        if (numDof == 3) {

            for (int i = 0; i < numConcMoment; i++) {

                if (elementIndex + 1 == concMoment(i, 0)) {

                    double M = concMoment(i, 1);
                    double a = locationConcMoment(i, 1);
                    double b = length - a;

                    /*
                                            if ((globalDof(0) == 0 || globalDof(1) == 0) && globalDof(2)){

                                                fixedEndForces(5) = (M * (pow(length,2) - 3 * pow(a,2))) / (2 * pow(length,2));
                                                fixedEndForces(1) = (M + fixedEndForces(5)) / length;
                                                fixedEndForces(4) = -fixedEndForces(1);
                                            }
                                            else if ((globalDof(3) == 0 || globalDof(4) == 0) && globalDof(5)){

                                                fixedEndForces(2) = -(M * (pow(length,2) - 3 * pow(a,2))) / (2 * pow(length,2));
                                                fixedEndForces(1) = (M + fixedEndForces(2)) / length;
                                                fixedEndForces(4) = -fixedEndForces(1);
                                            }
                                            else {
                    */
                    fixedEndForces(1) = (6 * M * a * b) / (pow(length, 3));
                    fixedEndForces(2) = (M * b * (2 * a - b)) / pow(length, 2);

                    fixedEndForces(4) = -fixedEndForces(1);
                    fixedEndForces(5) = (M * a * (2 * b - a)) / pow(length, 2);
                }
            }
        }
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

        //fixedEndForces = MatrixXd::Constant(2*numDof,1,0);

        if (numDof == 3) {

            for (int i = 0; i < numPointLoad; i++) {

                if (elementIndex + 1 == pointLoad(i, 0)) {

                    double P = pointLoad(i, 1);
                    double a = locationPointLoad(i, 1);
                    double b = length - a;

                    /*
                                            if (0){
                                            //if ((globalDof(0) == 0 || globalDof(1) == 0) && globalDof(2)){

                                                fixedEndForces(5) += (P * a * b * (length + a)) / (2 * pow(length,2));
                                                fixedEndForces(4) += -(P * a + fixedEndForces(5)) / length;

                                                fixedEndForces(1) += -(P + fixedEndForces(4));
                                            }
                                            //else if ((globalDof(3) == 0 || globalDof(4) == 0) && globalDof(5)){
                                            else if (0){

                                                fixedEndForces(2) += -(P * a * b * (length + b)) / (2 * pow(length,2));
                                                fixedEndForces(4) += -(P * a + fixedEndForces(2)) / length;

                                                fixedEndForces(1) += -(P + fixedEndForces(4));
                                            }
                                            else {
                    */
                    fixedEndForces(1) += -(P * pow(b, 2) * (3 * a + b)) / (pow(length, 3));
                    fixedEndForces(2) += -(P * a * pow(b, 2)) / pow(length, 2);

                    fixedEndForces(4) += -(P * pow(a, 2) * (3 * b + a)) / (pow(length, 3));
                    fixedEndForces(5) += (P * b * pow(a, 2)) / pow(length, 2);
                }
            }
        }
    }


    /**
        Calculates fixed end forces due to the point load on the element of a frame.
        @param numDof number of degrees of freedom of each node of the structure.
        @param elementIndex index of the element (element id - 1)
        @param numDistLoad
        @param pointLoad point loads' information (element id and magnitude)
    */
    void getFEofDistLoad(int numDof, int elementIndex, int numDistLoad, MatrixXd distLoad) {

        //fixedEndForces = MatrixXd::Constant(2*numDof,1,0);

        if (numDof == 3) {

            for (int i = 0; i < numDistLoad; i++) {

                if (elementIndex + 1 == distLoad(i, 0)) {

                    double q = distLoad(i, 1);

                    fixedEndForces(1) += -(q * length) / 2;
                    fixedEndForces(2) += -(q * pow(length, 2)) / 12;

                    fixedEndForces(4) += fixedEndForces(1);
                    fixedEndForces(5) += -fixedEndForces(2);
                }
            }
        }
    }


    /**
        Changes fixed end forces of the frame element to the global degrees of freedom locations.
        @param numDof number of degrees of freedom of each node of the structure.
        @param numEq number of active degree of freedoms in the structure (number of equations)
    */
    void getGlobalFE(int numEq, int numDof) {

        globalFixedEndForces = MatrixXd::Constant(numEq, 1, 0);
        rotatedFixedEndForces = MatrixXd::Constant(2 * numDof, 1, 0);

        rotatedFixedEndForces = R.transpose() * fixedEndForces;


        for (int i = 0; i < numDof * 2; i++) {

            if (globalDof(i) != 0) {

                globalFixedEndForces(globalDof(i) - 1) = rotatedFixedEndForces(i);
            }
        }
    }
};


#endif // ELEMENTCLASS_H_INCLUDED
