#include "Element.h"

Element::Element() {}

Element::Element(const int &index,
                 const std::vector<Node *> &connection,
                 Material *material,
                 const double &thickness,
                 const std::string &elementType)
{
    index_ = index;
    connection_ = connection;
    material_ = material;
    thickness_ = thickness;
    elementType_ = elementType;
    //shapeForce_ = 0.0;

    if (elementType == "T3")
    {
        order_ = 1;
    }
    else if (elementType == "T6")
    {
        order_ = 2;
    }
    else if (elementType == "T10")
    {
        order_ = 3;
    }

    int n = (order_ + 1) * (order_ + 2) / 2.0;

    connection_.reserve(n);
}

Element::~Element() {}

int Element::getIndex()
{
    return index_;
}

std::vector<Node *> Element::getConnection()
{
    return connection_;
}

Material *Element::getMaterial()
{
    return material_;
}

void Element::setShapeForce(const std::vector<double> &shapeForce)
{
    bounded_vector<double, 2> shapeee;
    shapeee(0) = shapeForce[0];
    shapeee(1) = shapeForce[1];

    shapeForce_ = shapeee;
}

void Element::setAnalysisParameters(const double &numberOfDomainIntegrationPoints, const double &deltat, const double &beta, const double &gamma)
{
    numberOfDomainIntegrationPoints_ = numberOfDomainIntegrationPoints;
    deltat_ = deltat;
    beta_ = beta;
    gamma_ = gamma;
}

vector<double> Element::domainShapeFunction(const double &xsi1, const double &xsi2)
{
    int n = (order_ + 1) * (order_ + 2) / 2.0; //number of nodes per element
    vector<double> phi(n, 0.0);

    if (order_ == 1)
    {
        phi(0) = xsi1;
        phi(1) = xsi2;
        phi(2) = 1.0 - xsi1 - xsi2;
    }
    else if (order_ == 2)
    {
        phi(0) = xsi1 * (2.0 * xsi1 - 1.0);
        phi(1) = xsi2 * (2.0 * xsi2 - 1.0);
        phi(2) = (xsi2 + xsi1 - 1.0) * (2.0 * xsi2 + 2.0 * xsi1 - 1.0);
        phi(3) = 4.0 * xsi1 * xsi2;
        phi(4) = -4.0 * xsi2 * (xsi2 + xsi1 - 1.0);
        phi(5) = -4.0 * xsi1 * (xsi2 + xsi1 - 1.0);
    }
    else if (order_ == 3)
    {
        phi(0) = (xsi1 * (3.0 * xsi1 - 2.0) * (3.0 * xsi1 - 1.0)) / 2.0;
        phi(1) = (xsi2 * (3.0 * xsi2 - 2.0) * (3.0 * xsi2 - 1.0)) / 2.0;
        phi(2) = -((xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0) * (3.0 * xsi2 + 3.0 * xsi1 - 1.0)) / 2.0;
        phi(3) = (9.0 * xsi1 * xsi2 * (3.0 * xsi1 - 1.0)) / 2.0;
        phi(4) = (9.0 * xsi1 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        phi(5) = -(9.0 * xsi2 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 - 1.0)) / 2.0;
        phi(6) = (9.0 * xsi2 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0)) / 2.0;
        phi(7) = (9.0 * xsi1 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0)) / 2.0;
        phi(8) = -(9.0 * xsi1 * (3.0 * xsi1 - 1.0) * (xsi2 + xsi1 - 1.0)) / 2.0;
        phi(9) = -27.0 * xsi1 * xsi2 * (xsi2 + xsi1 - 1.0);
    }

    return phi;
}

matrix<double> Element::domainDerivativeShapeFunction(const double &xsi1, const double &xsi2)
{
    int n = (order_ + 1) * (order_ + 2) / 2.0; //number of nodes per element
    matrix<double> dphi_dxsi(2, n, 0.0);

    if (order_ == 1)
    {
        dphi_dxsi(0, 0) = 1.0;
        dphi_dxsi(0, 1) = 0.0;
        dphi_dxsi(0, 2) = -1.0;

        dphi_dxsi(1, 0) = 0.0;
        dphi_dxsi(1, 1) = 1.0;
        dphi_dxsi(1, 2) = -1.0;
    }
    else if (order_ == 2)
    {
        dphi_dxsi(0, 0) = 4.0 * xsi1 - 1.0;
        dphi_dxsi(0, 1) = 0.0;
        dphi_dxsi(0, 2) = 4.0 * xsi2 + 4.0 * xsi1 - 3.0;
        dphi_dxsi(0, 3) = 4.0 * xsi2;
        dphi_dxsi(0, 4) = -4.0 * xsi2;
        dphi_dxsi(0, 5) = -4.0 * (xsi2 + 2.0 * xsi1 - 1.0);

        dphi_dxsi(1, 0) = 0.0;
        dphi_dxsi(1, 1) = 4.0 * xsi2 - 1.0;
        dphi_dxsi(1, 2) = 4.0 * xsi2 + 4.0 * xsi1 - 3.0;
        dphi_dxsi(1, 3) = 4.0 * xsi1;
        dphi_dxsi(1, 4) = -4.0 * (2.0 * xsi2 + xsi1 - 1.0);
        dphi_dxsi(1, 5) = -4.0 * xsi1;
    }
    else if (order_ == 3)
    {
        dphi_dxsi(0, 0) = (27.0 * xsi1 * xsi1 - 18.0 * xsi1 + 2.0) / 2.0;
        dphi_dxsi(0, 1) = 0.0;
        dphi_dxsi(0, 2) = -(27.0 * xsi2 * xsi2 + 54.0 * xsi1 * xsi2 - 36.0 * xsi2 + 27.0 * xsi1 * xsi1 - 36.0 * xsi1 + 11.0) / 2.0;
        dphi_dxsi(0, 3) = (9.0 * xsi2 * (6.0 * xsi1 - 1.0)) / 2.0;
        dphi_dxsi(0, 4) = (9.0 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        dphi_dxsi(0, 5) = -(9.0 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        dphi_dxsi(0, 6) = (9.0 * xsi2 * (6.0 * xsi2 + 6.0 * xsi1 - 5.0)) / 2.0;
        dphi_dxsi(0, 7) = (9.0 * (3.0 * xsi2 * xsi2 + 12.0 * xsi1 * xsi2 - 5.0 * xsi2 + 9.0 * xsi1 * xsi1 - 10.0 * xsi1 + 2.0)) / 2.0;
        dphi_dxsi(0, 8) = -(9.0 * (6.0 * xsi1 * xsi2 - xsi2 + 9.0 * xsi1 * xsi1 - 8.0 * xsi1 + 1.0)) / 2.0;
        dphi_dxsi(0, 9) = -27.0 * xsi2 * (xsi2 + 2.0 * xsi1 - 1.0);

        dphi_dxsi(1, 0) = 0.0;
        dphi_dxsi(1, 1) = (27.0 * xsi2 * xsi2 - 18.0 * xsi2 + 2) / 2.0;
        dphi_dxsi(1, 2) = -(27.0 * xsi2 * xsi2 + 54.0 * xsi1 * xsi2 - 36.0 * xsi2 + 27.0 * xsi1 * xsi1 - 36.0 * xsi1 + 11.0) / 2.0;
        dphi_dxsi(1, 3) = (9.0 * xsi1 * (3.0 * xsi1 - 1.0)) / 2.0;
        dphi_dxsi(1, 4) = (9.0 * xsi1 * (6.0 * xsi2 - 1.0)) / 2.0;
        dphi_dxsi(1, 5) = -(9.0 * (9.0 * xsi2 * xsi2 + 6.0 * xsi1 * xsi2 - 8.0 * xsi2 - xsi1 + 1.0)) / 2.0;
        dphi_dxsi(1, 6) = (9.0 * (9.0 * xsi2 * xsi2 + 12.0 * xsi1 * xsi2 - 10.0 * xsi2 + 3.0 * xsi1 * xsi1 - 5.0 * xsi1 + 2.0)) / 2.0;
        dphi_dxsi(1, 7) = (9.0 * xsi1 * (6.0 * xsi2 + 6.0 * xsi1 - 5.0)) / 2.0;
        dphi_dxsi(1, 8) = -(9.0 * xsi1 * (3.0 * xsi1 - 1.0)) / 2.0;
        dphi_dxsi(1, 9) = -27.0 * xsi1 * (2.0 * xsi2 + xsi1 - 1.0);
    }

    return dphi_dxsi;
}

matrix<double> Element::hammerQuadrature(const int &nh)
{
    matrix<double> hammer(nh, 3, 0.0);

    if (nh == 7)
    {
        hammer(0, 0) = 1.0 / 3.0;
        hammer(0, 1) = 1.0 / 3.0;
        hammer(0, 2) = 9.0 / 80.0;

        hammer(1, 0) = (9.0 + 2.0 * sqrt(15.0)) / 21.0;
        hammer(1, 1) = (6.0 - sqrt(15.0)) / 21.0;
        hammer(1, 2) = (155.0 - sqrt(15.0)) / 2400.0;

        hammer(2, 0) = (6.0 - sqrt(15.0)) / 21.0;
        hammer(2, 1) = (9.0 + 2.0 * sqrt(15.0)) / 21.0;
        hammer(2, 2) = (155.0 - sqrt(15.0)) / 2400.0;

        hammer(3, 0) = (6.0 - sqrt(15.0)) / 21.0;
        hammer(3, 1) = (6.0 - sqrt(15.0)) / 21.0;
        hammer(3, 2) = (155.0 - sqrt(15.0)) / 2400.0;

        hammer(4, 0) = (6.0 + sqrt(15.0)) / 21.0;
        hammer(4, 1) = (6.0 + sqrt(15.0)) / 21.0;
        hammer(4, 2) = (155.0 + sqrt(15.0)) / 2400.0;

        hammer(5, 0) = (9.0 - 2.0 * sqrt(15.0)) / 21.0;
        hammer(5, 1) = (6.0 + sqrt(15.0)) / 21.0;
        hammer(5, 2) = (155.0 + sqrt(15.0)) / 2400.0;

        hammer(6, 0) = (6.0 + sqrt(15.0)) / 21.0;
        hammer(6, 1) = (9.0 - 2.0 * sqrt(15.0)) / 21.0;
        hammer(6, 2) = (155.0 + sqrt(15.0)) / 2400.0;
    }

    else if (nh == 12)
    {
        hammer(0, 0) = 0.501426509658179;
        hammer(0, 1) = 0.249286745170910;
        hammer(0, 2) = 0.05839313786319;

        hammer(1, 0) = 0.249286745170910;
        hammer(1, 1) = 0.249286745170910;
        hammer(1, 2) = 0.05839313786319;

        hammer(2, 0) = 0.249286745170910;
        hammer(2, 1) = 0.501426509658179;
        hammer(2, 2) = 0.05839313786319;

        hammer(3, 0) = 0.873821971016996;
        hammer(3, 1) = 0.063089014491502;
        hammer(3, 2) = 0.025422453185104;

        hammer(4, 0) = 0.063089014491502;
        hammer(4, 1) = 0.063089014491502;
        hammer(4, 2) = 0.025422453185104;

        hammer(5, 0) = 0.063089014491502;
        hammer(5, 1) = 0.873821971016996;
        hammer(5, 2) = 0.025422453185104;

        hammer(6, 0) = 0.053145049844816;
        hammer(6, 1) = 0.310352451033785;
        hammer(6, 2) = 0.041425537809187;

        hammer(7, 0) = 0.310352451033785;
        hammer(7, 1) = 0.636502499121399;
        hammer(7, 2) = 0.041425537809187;

        hammer(8, 0) = 0.636502499121399;
        hammer(8, 1) = 0.053145049844816;
        hammer(8, 2) = 0.041425537809187;

        hammer(9, 0) = 0.310352451033785;
        hammer(9, 1) = 0.053145049844816;
        hammer(9, 2) = 0.041425537809187;

        hammer(10, 0) = 0.636502499121399;
        hammer(10, 1) = 0.310352451033785;
        hammer(10, 2) = 0.041425537809187;

        hammer(11, 0) = 0.053145049844816;
        hammer(11, 1) = 0.636502499121399;
        hammer(11, 2) = 0.041425537809187;
    }
    return hammer;
}

bounded_matrix<double, 2, 2> Element::referenceJacobianMatrix(const double &xsi1, const double &xsi2)
{
    matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);
    double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

    for (int i = 0; i < connection_.size(); i++)
    {
        bounded_vector<double, 2> initialCoord = connection_[i]->getInitialCoordinate();
        dx1_dxsi1 += initialCoord(0) * dphi_dxsi(0, i);
        dx1_dxsi2 += initialCoord(0) * dphi_dxsi(1, i);
        dx2_dxsi1 += initialCoord(1) * dphi_dxsi(0, i);
        dx2_dxsi2 += initialCoord(1) * dphi_dxsi(1, i);
    }

    bounded_matrix<double, 2, 2> referenceJacobianMatrix;
    referenceJacobianMatrix(0, 0) = dx1_dxsi1;
    referenceJacobianMatrix(1, 0) = dx2_dxsi1;
    referenceJacobianMatrix(0, 1) = dx1_dxsi2;
    referenceJacobianMatrix(1, 1) = dx2_dxsi2;

    return referenceJacobianMatrix;
}

bounded_matrix<double, 2, 2> Element::currentJacobianMatrix(const double &xsi1, const double &xsi2)
{
    matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);
    double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

    for (int i = 0; i < connection_.size(); i++)
    {
        bounded_vector<double, 2> currentCoord = connection_[i]->getCurrentCoordinate();
        dx1_dxsi1 += currentCoord(0) * dphi_dxsi(0, i);
        dx1_dxsi2 += currentCoord(0) * dphi_dxsi(1, i);
        dx2_dxsi1 += currentCoord(1) * dphi_dxsi(0, i);
        dx2_dxsi2 += currentCoord(1) * dphi_dxsi(1, i);
    }

    bounded_matrix<double, 2, 2> currentJacobianMatrix;
    currentJacobianMatrix(0, 0) = dx1_dxsi1;
    currentJacobianMatrix(1, 0) = dx2_dxsi1;
    currentJacobianMatrix(0, 1) = dx1_dxsi2;
    currentJacobianMatrix(1, 1) = dx2_dxsi2;

    return currentJacobianMatrix;
}

double Element::jacobianDeterminant(const bounded_matrix<double, 2, 2> &jacobianMatrix)
{
    return (jacobianMatrix(0, 0) * jacobianMatrix(1, 1) - jacobianMatrix(0, 1) * jacobianMatrix(1, 0));
}

std::pair<vector<double>, matrix<double>> Element::elementContributions(const std::string &ep, const std::string &typeAnalyze, const int &step, const int &numberOfStep)
{
    vector<double> rhs(2 * connection_.size(), 0.0);
    matrix<double> tangent(2 * connection_.size(), 2 * connection_.size(), 0.0);

    matrix<double> domainIntegrationPoints_ = hammerQuadrature(numberOfDomainIntegrationPoints_);

    for (int ih = 0; ih < numberOfDomainIntegrationPoints_; ih++)
    {
        double xsi1 = domainIntegrationPoints_(ih, 0);
        double xsi2 = domainIntegrationPoints_(ih, 1);
        double weight = domainIntegrationPoints_(ih, 2);

        vector<double> phi = domainShapeFunction(xsi1, xsi2);
        matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);  //row = direction, column = node
        bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
        double j0 = jacobianDeterminant(A0);
        bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
        A0I(0, 0) = A0(1, 1) / j0;
        A0I(1, 1) = A0(0, 0) / j0;
        A0I(0, 1) = -A0(0, 1) / j0;
        A0I(1, 0) = -A0(1, 0) / j0;
        bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2); //current configuration map
        double j1 = jacobianDeterminant(A1);
        bounded_matrix<double, 2, 2> A1I; //inverse current configuration map
        A1I(0, 0) = A1(1, 1) / j1;
        A1I(1, 1) = A1(0, 0) / j1;
        A1I(1, 0) = -A1(1, 0) / j1;
        A1I(0, 1) = -A1(0, 1) / j1;
        bounded_matrix<double, 2, 2> Ac = prod(A1, A0I); //current deformation gradient
        double jac = jacobianDeterminant(Ac);
        bounded_matrix<double, 2, 2> AcI; //inverse current deformation gradient
        AcI(0, 0) = Ac(1, 1) / jac;
        AcI(1, 1) = Ac(0, 0) / jac;
        AcI(1, 0) = -Ac(1, 0) / jac;
        AcI(0, 1) = -Ac(0, 1) / jac;
        identity_matrix<double> I(2);                                      //identity matrix
        bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor

        bounded_matrix<double, 2, 2> S; //second piola kirchhoff stress tensor
        double young = material_->getYoung();
        double poisson = material_->getPoisson();
        double density;

        if (typeAnalyze == "STATIC")
        {
            density = 0.0;
        }
        else
        {
            density = material_->getDensity();
        }

        if (ep == "EPD")
        {
            S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }
        else
        {
            S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }

        bounded_matrix<double, 2, 2> dA_dy;
        bounded_matrix<double, 2, 2> dA_dy2;

        for (int i = 0; i < connection_.size(); i++)
        {
            for (int j = 0; j < 2; j++)
            {
                if (j == 0)
                {
                    dA_dy(0, 0) = dphi_dxsi(0, i);
                    dA_dy(0, 1) = dphi_dxsi(1, i);
                    dA_dy(1, 0) = 0.0;
                    dA_dy(1, 1) = 0.0;
                }
                else
                {
                    dA_dy(1, 0) = dphi_dxsi(0, i);
                    dA_dy(1, 1) = dphi_dxsi(1, i);
                    dA_dy(0, 0) = 0.0;
                    dA_dy(0, 1) = 0.0;
                }

                bounded_matrix<double, 2, 2> mat1 = prod(trans(A0I), trans(dA_dy));
                bounded_matrix<double, 2, 2> mat2 = prod(dA_dy, A0I);

                bounded_matrix<double, 2, 2> dE_dy = 0.5 * (prod(mat1, Ac) + prod(trans(Ac), mat2)); //first derivative of E regarding i,j

                double r = dE_dy(0, 0) * S(0, 0) + dE_dy(1, 1) * S(1, 1) + dE_dy(0, 1) * S(0, 1) + dE_dy(1, 0) * S(1, 0); //internal force

                double accel = 0.0;
                double shape = 0.0;

                shape = phi(i) * shape;

                if (typeAnalyze == "STATIC")
                {
                    for (int m = 0; m < connection_.size(); m++)
                    {
                        shape += phi(m) * (1.0 * shapeForce_(j) * step / (1.0 * numberOfStep));
                    }
                }
                else
                {
                    for (int m = 0; m < connection_.size(); m++)
                    {
                        accel += phi(m) * connection_[m]->getCurrentAcceleration()(j);
                        shape += phi(m) * shapeForce_(j);
                    }
                }

                double m = density * phi(i) * accel; //inertial force

                shape = phi(i) * shape;

                // // bounded_vector<double, 2> b;
                // // b(0) = 0.0;
                // // b(1) = phi(i) * density * gravity_; //domain force

                rhs(2 * i + j) -= (r + m - shape) * weight * j0 * thickness_;

                // rhs(2 * i + j) -= (r)*weight * j0 * thickness_;

                for (int k = 0; k < connection_.size(); k++)
                {
                    for (int l = 0; l < 2; l++)
                    {
                        if (l == 0)
                        {
                            dA_dy2(0, 0) = dphi_dxsi(0, k);
                            dA_dy2(0, 1) = dphi_dxsi(1, k);
                            dA_dy2(1, 0) = 0.0;
                            dA_dy2(1, 1) = 0.0;
                        }
                        else
                        {
                            dA_dy2(1, 0) = dphi_dxsi(0, k);
                            dA_dy2(1, 1) = dphi_dxsi(1, k);
                            dA_dy2(0, 0) = 0.0;
                            dA_dy2(0, 1) = 0.0;
                        }

                        bounded_matrix<double, 2, 2> mat3 = prod(trans(A0I), trans(dA_dy2));
                        bounded_matrix<double, 2, 2> mat4 = prod(dA_dy2, A0I);

                        bounded_matrix<double, 2, 2> dE_dy2 = 0.5 * (prod(mat3, Ac) + prod(trans(Ac), mat4)); //first derivative of E regarding k,l
                        bounded_matrix<double, 2, 2> dE_dy3 = 0.5 * (prod(mat1, mat4) + prod(mat3, mat2));    //second derivative of E regarding i,j,k,l

                        bounded_matrix<double, 2, 2> dS_dy;
                        dS_dy(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(0, 0) + poisson * dE_dy(1, 1)));
                        dS_dy(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(1, 1) + poisson * dE_dy(0, 0)));
                        dS_dy(1, 0) = (young / (1.0 + poisson)) * dE_dy(1, 0);
                        dS_dy(0, 1) = (young / (1.0 + poisson)) * dE_dy(0, 1);

                        double v = dE_dy3(0, 0) * S(0, 0) + dE_dy3(1, 1) * S(1, 1) + dE_dy3(0, 1) * S(0, 1) + dE_dy3(1, 0) * S(1, 0) +                //second part of equation 5.88
                                   dE_dy2(0, 0) * dS_dy(0, 0) + dE_dy2(1, 1) * dS_dy(1, 1) + dE_dy2(0, 1) * dS_dy(0, 1) + dE_dy2(1, 0) * dS_dy(1, 0); //viscosity and pressure contribution

                        // dS_dy(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(0, 0) + poisson * dE_dy(1, 1)));
                        // dS_dy(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(1, 1) + poisson * dE_dy(0, 0)));
                        // dS_dy(1, 0) = 2.0 * (young / (2.0 * (1.0 + poisson))) * dE_dy(1, 0);
                        // dS_dy(0, 1) = 2.0 * (young / (2.0 * (1.0 + poisson))) * dE_dy(0, 1);

                        // double v = dE_dy3(0, 0) * S(0, 0) + dE_dy3(1, 1) * S(1, 1) + dE_dy3(0, 1) * S(0, 1) + dE_dy3(1, 0) * S(1, 0) +                //second part of equation 5.88
                        //            dE_dy2(0, 0) * dS_dy(0, 0) + dE_dy2(1, 1) * dS_dy(1, 1) + dE_dy2(0, 1) * dS_dy(0, 1) + dE_dy2(1, 0) * dS_dy(1, 0); //viscosity and pressure contribution

                        double m = density * phi(i) * phi(k); //mass contribution

                        tangent(2 * i + j, 2 * k + l) += v * weight * j0 * thickness_;
                        if (j == l)
                            tangent(2 * i + j, 2 * k + l) += (m / (beta_ * deltat_ * deltat_)) * weight * j0 * thickness_;
                    }
                }
            }
        }
    }
    return std::make_pair(rhs, tangent);
}

void Element::StressCalculate(const std::string &ep)
{
    if (elementType_ == "T6")
    {
        for (int i = 0; i < 6; i++)
        {
            bounded_matrix<double, 6, 2> xsi;
            xsi(0, 0) = 1.0;
            xsi(0, 1) = 0.0;

            xsi(1, 0) = 0.0;
            xsi(1, 1) = 1.0;

            xsi(2, 0) = 0.0;
            xsi(2, 1) = 0.0;

            xsi(3, 0) = 0.5;
            xsi(3, 1) = 0.5;

            xsi(4, 0) = 0.0;
            xsi(4, 1) = 0.5;

            xsi(5, 0) = 0.5;
            xsi(5, 1) = 0.0;

            bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi(i, 0), xsi(i, 1)); //initial configuration map
            double j0 = jacobianDeterminant(A0);
            bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
            A0I(0, 0) = A0(1, 1) / j0;
            A0I(1, 1) = A0(0, 0) / j0;
            A0I(0, 1) = -A0(0, 1) / j0;
            A0I(1, 0) = -A0(1, 0) / j0;
            bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi(i, 0), xsi(i, 1));
            bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                   //current deformation gradient
            identity_matrix<double> I(2);                                      //identity matrix
            bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor
            bounded_matrix<double, 2, 2> S;                                    //second piola kirchhoff stress tensor
            double young = material_->getYoung();
            double poisson = material_->getPoisson();

            if (ep == "EPD")
            {
                S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
                S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
                S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
                S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
            }
            else
            {
                S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
                S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
                S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
                S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
            }

            bounded_matrix<double, 2, 2> sigma;
            double jac = jacobianDeterminant(Ac);
            bounded_matrix<double, 2, 2> mat1;
            mat1 = prod(Ac, S);
            sigma = (1.0 / jac) * prod(mat1, trans(Ac));

            bounded_vector<double, 3> cauchy;

            cauchy(0) = sigma(0, 0);
            cauchy(1) = sigma(1, 1);
            cauchy(2) = sigma(0, 1);

            getConnection()[i]->setStressState(cauchy);
        }
    }

    else if (elementType_ == "T10")
    {
        for (int i = 0; i < 10; i++)
        {
            bounded_matrix<double, 10, 2> xsi;
            xsi(0, 0) = 1.0;
            xsi(0, 1) = 0.0;

            xsi(1, 0) = 0.0;
            xsi(1, 1) = 1.0;

            xsi(2, 0) = 0.0;
            xsi(2, 1) = 0.0;

            xsi(3, 0) = 2.0 / 3.0;
            xsi(3, 1) = 1.0 / 3.0;

            xsi(4, 0) = 1.0 / 3.0;
            xsi(4, 1) = 2.0 / 3.0;

            xsi(5, 0) = 0.0;
            xsi(5, 1) = 2.0 / 3.0;

            xsi(6, 0) = 0.0;
            xsi(6, 1) = 1.0 / 3.0;

            xsi(7, 0) = 1.0 / 3.0;
            xsi(7, 1) = 0.0;

            xsi(8, 0) = 2.0 / 3.0;
            xsi(8, 1) = 0.0;

            xsi(9, 0) = 1.0 / 3.0;
            xsi(9, 1) = 1.0 / 3.0;

            bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi(i, 0), xsi(i, 1)); //initial configuration map
            double j0 = jacobianDeterminant(A0);
            bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
            A0I(0, 0) = A0(1, 1) / j0;
            A0I(1, 1) = A0(0, 0) / j0;
            A0I(0, 1) = -A0(0, 1) / j0;
            A0I(1, 0) = -A0(1, 0) / j0;
            bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi(i, 0), xsi(i, 1));
            bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                   //current deformation gradient
            identity_matrix<double> I(2);                                      //identity matrix
            bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor
            bounded_matrix<double, 2, 2> S;                                    //second piola kirchhoff stress tensor
            double young = material_->getYoung();
            double poisson = material_->getPoisson();

            if (ep == "EPD")
            {
                S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
                S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
                S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
                S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
            }
            else
            {
                S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
                S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
                S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
                S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
            }

            bounded_matrix<double, 2, 2> sigma;
            double jac = jacobianDeterminant(Ac);
            bounded_matrix<double, 2, 2> mat1;
            mat1 = prod(Ac, S);
            sigma = (1.0 / jac) * prod(mat1, trans(Ac));

            bounded_vector<double, 3> cauchy;

            cauchy(0) = sigma(0, 0);
            cauchy(1) = sigma(1, 1);
            cauchy(2) = sigma(0, 1);

            getConnection()[i]->setStressState(cauchy);
        }
    }
}