#include "FiberElement.h"

FiberElement::FiberElement() {}

FiberElement::FiberElement(const int &index,
                           const std::vector<FiberNode *> &connection,
                           FiberMaterial *material,
                           const double &area)
{
    index_ = index;
    connection_ = connection;
    material_ = material;
    area_ = area;
    initialLength_ = norm_2(connection_[0]->getInitialCoordinate() - connection_[1]->getInitialCoordinate());
}

FiberElement::~FiberElement() {}

int FiberElement::getIndex()
{
    return index_;
}

std::vector<FiberNode *> FiberElement::getConnection()
{
    return connection_;
}

std::pair<vector<double>, matrix<double>> FiberElement::fiberContributions(const std::string &typeAnalyze, const double &deltat, const double &beta)
{
    vector<double> first(4, 0.0);
    matrix<double> second(4, 4, 0.0);

    bounded_vector<double, 2> initialNode = connection_[0]->getCurrentCoordinate();
    bounded_vector<double, 2> endNode = connection_[1]->getCurrentCoordinate();

    double currentLength = norm_2(endNode - initialNode);
    double green = 0.5 * ((currentLength * currentLength) / (initialLength_ * initialLength_) - 1.0);
    double s, young;

    if (green <= material_->getPlastStrain() && green >= material_->getPlastStrain() * -1.0)
    {
        s = green * (material_->getYoung());
        young = material_->getYoung();
    }
    else if (green > material_->getPlastStrain())
    {
        s = material_->getPlastStrain() * (material_->getYoung()) + (green - material_->getPlastStrain()) * material_->getHardeningModulus();
        young = material_->getHardeningModulus();
    }
    else if (green < material_->getPlastStrain() * -1.0)
    {
        s = -(material_->getPlastStrain()) * (material_->getYoung()) + (green + material_->getPlastStrain()) * material_->getHardeningModulus();
        young = material_->getHardeningModulus();
    }

    first(0) = -1.0 * area_ * s * (endNode(0) - initialNode(0)) / initialLength_;
    first(1) = -1.0 * area_ * s * (endNode(1) - initialNode(1)) / initialLength_;
    first(2) = area_ * s * (endNode(0) - initialNode(0)) / initialLength_;
    first(3) = area_ * s * (endNode(1) - initialNode(1)) / initialLength_;

    for (int beta = 0; beta < 2; beta++)
    {
        for (int alfa = 0; alfa < 2; alfa++)
        {
            for (int i = 0; i < 2; i++)
            {
                for (int k = 0; k < 2; k++)
                {
                    int l = 2 * beta + i;
                    int m = 2 * alfa + k;
                    double de_dy = pow(-1, beta + 1) * pow(-1, alfa + 1) * area_ * young * (endNode(i) - initialNode(i)) * (endNode(k) - initialNode(k)) / (initialLength_ * initialLength_ * initialLength_);
                    double d2edy2 = 0.0;

                    if (i == k)
                    {
                        d2edy2 = pow(-1, beta + 1) * pow(-1, alfa + 1) * area_ / initialLength_ * s;
                    }
                    second(l, m) = de_dy + d2edy2;
                }
            }
        }
    }

    if (typeAnalyze == "DYNAMIC") //NÃO ESTÁ CONSIDERANDO O AMORTECIMENTO!
    {
        identity_matrix<double> identity(4);
        bounded_matrix<double, 4, 4> mass;
        mass = 0.5 * area_ * (material_->getDensity()) * initialLength_ * identity;

        vector<double> acceleration(4, 0.0);
        acceleration(0) = connection_[0]->getCurrentAcceleration()(0);
        acceleration(1) = connection_[0]->getCurrentAcceleration()(1);
        acceleration(2) = connection_[1]->getCurrentAcceleration()(0);
        acceleration(3) = connection_[1]->getCurrentAcceleration()(1);
        first = first + prod(mass, acceleration);

        second = second + mass / (beta * deltat * deltat);
    }

    return std::make_pair(first, second);
}

void FiberElement::updateNormalForce()
{
    bounded_vector<double, 2> initialNode = connection_[0]->getCurrentCoordinate();
    bounded_vector<double, 2> endNode = connection_[1]->getCurrentCoordinate();

    double currentLength = norm_2(endNode - initialNode);
    double green = 0.5 * ((currentLength * currentLength) / (initialLength_ * initialLength_) - 1.0);
    double s;

    if (green <= material_->getPlastStrain() && green >= material_->getPlastStrain() * -1.0)
    {
        s = green * (material_->getYoung());
    }
    else if (green > material_->getPlastStrain())
    {
        s = material_->getPlastStrain() * (material_->getYoung()) + (green - material_->getPlastStrain()) * material_->getHardeningModulus();
    }
    else if (green < material_->getPlastStrain() * -1.0)
    {
        s = -(material_->getPlastStrain()) * (material_->getYoung()) + (green + material_->getPlastStrain()) * material_->getHardeningModulus();
    }

    connection_[0]->setNormalForce(s*area_);
    connection_[1]->setNormalForce(s*area_);
}