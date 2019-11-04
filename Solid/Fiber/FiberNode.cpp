#include "FiberNode.h"

FiberNode::FiberNode() {}

FiberNode::FiberNode(const int &index,
                     const bounded_vector<double, 2> &initialCoordinate)
{
    index_ = index;
    initialCoordinate_ = initialCoordinate;
    pastCoordinate_ = initialCoordinate;
    currentCoordinate_ = initialCoordinate;
    for (size_t i = 0; i < 2; i++)
    {
        currentVelocity_(i) = 0.0;
        currentAcceleration_(i) = 0.0;
        pastVelocity_(i) = 0.0;
        pastAcceleration_(i) = 0.0;
        incidence_(i) = 0.0;
    }
    indexElement_ = -1;
    normalForce_=0.0;
}

FiberNode::~FiberNode() {}

int FiberNode::getIndex()
{
    return index_;
}

bounded_vector<double, 2> FiberNode::getInitialCoordinate()
{
    return initialCoordinate_;
}

bounded_vector<double, 2> FiberNode::getPastCoordinate()
{
    return pastCoordinate_;
}

bounded_vector<double, 2> FiberNode::getPastVelocity()
{
    return pastVelocity_;
}

bounded_vector<double, 2> FiberNode::getPastAcceleration()
{
    return pastAcceleration_;
}

bounded_vector<double, 2> FiberNode::getCurrentCoordinate()
{
    return currentCoordinate_;
}

bounded_vector<double, 2> FiberNode::getCurrentVelocity()
{
    return currentVelocity_;
}

bounded_vector<double, 2> FiberNode::getCurrentAcceleration()
{
    return currentAcceleration_;
}

double FiberNode::getNormalForce()
{
    return normalForce_;
}

bounded_vector<double, 2> FiberNode::getDimensionlessCoordinates()
{
    return incidence_;
}

int FiberNode::getIndexSolidElement()
{
    return indexElement_;
}

void FiberNode::setPastCoordinate(const bounded_vector<double, 2> &pastCoordinate)
{
    pastCoordinate_ = pastCoordinate;
}

void FiberNode::setPastVelocity(const bounded_vector<double, 2> &pastVelocity)
{
    pastVelocity_ = pastVelocity;
}

void FiberNode::setPastAcceleration(const bounded_vector<double, 2> &pastAcceleration)
{
    pastAcceleration_ = pastAcceleration;
}

void FiberNode::setCurrentCoordinate(const bounded_vector<double, 2> &currentCoordinate)
{
    currentCoordinate_ = currentCoordinate;
}

void FiberNode::setCurrentVelocity(const bounded_vector<double, 2> &currentVelocity)
{
    currentVelocity_ = currentVelocity;
}

void FiberNode::setCurrentAcceleration(const bounded_vector<double, 2> &currentAcceleration)
{
    currentAcceleration_ = currentAcceleration;
}

void FiberNode::setIncidenceInSolid(const int &indexElement, const bounded_vector<double, 2> &incidence)
{
    indexElement_=indexElement;
    incidence_ = incidence;
}

void FiberNode::setNormalForce(const double &force)
{
    normalForce_=force;
}