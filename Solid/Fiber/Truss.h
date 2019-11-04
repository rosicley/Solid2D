#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "Node.h"
#include "Element.h"
#include "Material.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <boost/numeric/bindings/lapack/driver/ggev.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getri.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>

using namespace boost::numeric::ublas;

class Truss
{
public:
    std::vector<Node *> getNodes();

    std::vector<Element *> getElements();

    std::vector<Material *> getMaterial();

    vector<double> InternalForces();

    vector<double> TemperatureForces(const int &numberOfSteps, const int &currentStep);

    vector<double> InertialForces(const double &beta, const double &gamma, const double &deltat);

    vector<double> ExternalForces();

    vector<int> BoundaryCondition();

    matrix<double> Hessian();

    matrix<double> TemperatureHessian(const int &numberOfSteps, const int &currentStep);

    matrix<double> MassMatrix();

    int solveStaticProblem(const int &numberOfSteps, const double &tolerance);

    int solveDynamicProblem(const int &numberOfTimes, const double &tolerance);

    void exportToParaview(const int &loadstep);

    void DisplacementVesursLoad(const double &force);

    void addNode(const int &index,
                 const std::vector<double> &initialCoordinate);

    void addElement(const int &index,
                    const std::vector<int> &connection,
                    const int &material,
                    const double &area);

    void addMaterial(const int &index,
                     const double &young,
                     const double &plastStrain,
                     const double &hardeningModulus,
                     const double &density,
                     const double &expansionCoef);

    void readInput(const std::string &read,
                   const std::string &typeAnalyze);

private:
    std::vector<Node *> nodes_;

    std::vector<Element *> elements_;

    std::vector<Material *> materials_;

    std::vector<int> boundaryConditions_;

    std::vector<double> externalForces_;

    // int numberOfSteps_;

    // double tolerance_;

    std::string name_;
};