#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "Node.h"
#include "Element.h"
#include "Material.h"
#include "DirichletCondition.h"
#include "NeumannCondition.h"
#include "Fiber/FiberElement.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/bindings/lapack/driver/posv.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getri.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>

using namespace boost::numeric::ublas;

class Solid
{
public:
    Solid();

    ~Solid();

    void addMaterial(const int &index,
                     const double &young,
                     const double &poisson,
                     const double &density);

    void addFiberMaterial(const int &index,
                          const double &young,
                          const double &plastStrain,
                          const double &hardeningModulus,
                          const double &density);

    void addNode(const int &index,
                 const bounded_vector<double, 2> &initialCoordinate);

    void addFiberNode(const int &index,
                      const bounded_vector<double, 2> &initialCoordinate);

    void addElement(const int &index,
                    const std::vector<int> &nodesIndex,
                    const int &materialIndex,
                    const double &thickness,
                    const std::string &elementType);

    void addFiberElement(const int &index,
                         const std::vector<int> &nodesIndex,
                         const int &materialIndex,
                         const double &area);

    void addDirichletCondition(const int &index,
                               const int &direction,
                               const double &value);

    void addNeumannCondition(const int &index,
                             const int &direction,
                             const double &value);

    void setAnalysisParameters(const std::string &planeState, const std::string &elementType, const int &numberOfHammer);

    void setDynamicAnalysisParameters(const double &deltat, const double &beta, const double &gama);

    vector<double> domainShapeFunction(const double &xsi1, const double &xsi2);

    std::pair<vector<double>, matrix<double, column_major>> globalSolid(const std::string &typeAnalyze, const int &step, const int &numberOfStep);

    // std::vector<Node *> getNodes();

    // std::vector<Element *> getElements();

    // std::vector<Material *> getMaterial();

    vector<double> ExternalForces();

    int solveStaticProblem(const int &numberOfSteps, const int &maximumOfInteration, const double &tolerance);

    bounded_matrix<double, 2 ,2> inverseMatrix(const bounded_matrix<double, 2, 2> &matrix);

    //int solveDynamicProblem(const int &numberOfTimes, const double &tolerance);

    void exportToParaview(const int &loadstep);

    void readAnsysInput(const std::string &read);

    void readFibersInput(const std::string &read);

    void incidenceOfFibers();

private:
    std::vector<Node *> nodes_;

    std::vector<Element *> elements_;

    std::vector<Material *> materials_;

    std::vector<DirichletCondition *> dirichletConditions_;

    std::vector<NeumannCondition *> neumannConditions_;

    std::vector<FiberNode *> fiberNodes_;

    std::vector<FiberElement *> fiberElements_;

    std::vector<FiberElement *> fiberInsideSolid_;

    std::vector<FiberMaterial *> fiberMaterials_;

    int numberOfSteps_;

    double tolerance_;

    //std::string name_;

    std::string elementType_;

    std::string planeState_;

    double deltat_;

    double gamma_;

    double beta_;

    bounded_vector<double, 2> shapeForces_;

    int numberOfHammer_;

    int order_;
};