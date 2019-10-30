#pragma once
#include <iostream>
#include <vector>
#include <string>
#include "Node.h"
#include "Material.h"
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>


using namespace boost::numeric::ublas;

class Element
{
public:
    Element();

    Element(const int &index,
            const std::vector<Node *> &connection,
            Material *material,
            const double &thickness,
            const std::string &elementType);

    ~Element();

    int getIndex();

    std::vector<Node *> getConnection();

    Material *getMaterial();

    vector<double> domainShapeFunction(const double &xsi1, const double &xsi2);

    matrix<double> domainDerivativeShapeFunction(const double &xsi1, const double &xsi2);

    matrix<double> hammerQuadrature(const int& nh);

    bounded_matrix<double, 2, 2> referenceJacobianMatrix(const double &xsi1, const double &xsi2);

    bounded_matrix<double, 2, 2> currentJacobianMatrix(const double &xsi1, const double &xsi2);

    double jacobianDeterminant(const bounded_matrix<double, 2, 2> &jacobianMatrix);

    std::pair<vector<double>, matrix<double>> elementContributions(const std::string &ep, const std::string &typeAnalyze, const int &step, const int &numberOfStep);

    void setShapeForce(const std::vector<double> &shapeForce);

    void setAnalysisParameters(const double &numberOfDomainIntegrationPoints, const double &deltat, const double &beta, const double &gamma);

    void StressCalculate(const std::string &ep);

    // std::vector<double> InternalForce();

    // bounded_matrix<double, 6, 6> localHessian();

    // bounded_matrix<double, 6, 6> localMassMatrix();

    // void setArea(const double &area);

private:
    int index_;

    int order_;

    int numberOfDomainIntegrationPoints_;

    double thickness_;

    double beta_;

    double gamma_;

    double deltat_;

    std::string elementType_;

    bounded_vector<double, 2> shapeForce_;

    std::vector<Node *> connection_;

    Material *material_;

    //double gravity_;
};