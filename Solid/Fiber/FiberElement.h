#pragma once
#include <iostream>
#include <vector>
#include "FiberNode.h"
#include "FiberMaterial.h"
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

class FiberElement
{
public:
    FiberElement();

    FiberElement(const int &index,
                 const std::vector<FiberNode *> &connection,
                 FiberMaterial *material,
                 const double &area);

    ~FiberElement();

    int getIndex();

    std::vector<FiberNode *> getConnection();

    std::pair<vector<double>, matrix<double>> fiberLocalContributions(const std::string &typeAnalyze, const double &deltat, const double &beta);

    bounded_matrix<double, 4, 4> localMassMatrix();

    void updateNormalForce();

private:
    int index_;

    std::vector<FiberNode *> connection_;

    FiberMaterial *material_;

    double area_;

    double initialLength_;
};