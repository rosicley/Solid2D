#pragma once
#include <iostream>
#include <vector>
#include "Node.h"
#include "Material.h"
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

class Element
{
public:
    Element();

    Element(const int &index,
            const std::vector<Node *> &connection,
            Material *material,
            const double &thickness);

    ~Element();

    int getIndex();

    std::vector<Node *> getConnection();

    Material *getMaterial();

    vector<double> domainShapeFunction(const double &xsi1, const double &xsi2);

    matrix<double> domainDerivativeShapeFunction(const double &xsi1, const double &xsi2);

    matrix<double> hammerQuadrature();

    bounded_matrix<double, 2, 2> referenceJacobianMatrix(const double &xsi1, const double &xsi2);

    bounded_matrix<double, 2, 2> currentJacobianMatrix(const double &xsi1, const double &xsi2);

    double jacobianDeterminant(const bounded_matrix<double, 2, 2> &jacobianMatrix);

    std::pair<vector<double>, matrix<double>> elementContributions();

    // std::vector<double> InternalForce();

    // bounded_matrix<double, 6, 6> localHessian();

    // bounded_matrix<double, 6, 6> localMassMatrix();

    // void setArea(const double &area);

private:
    int index_;

    std::vector<Node *> connection_;

    Material *material_;

    double thickness_;

    double gravity_;

    double beta_, deltat_;  
};