#include "Node.h"

Node::Node() {}

Node::Node(const int &index,
           const bounded_vector<double, 2> &initialCoordinate)
{
    index_ = index;
    initialCoordinate_ = initialCoordinate;
    pastCoordinate_ = initialCoordinate;
    currentCoordinate_ = initialCoordinate;
    pastVelocity_ = 0.0;
    pastAcceleration_ = 0.0;
    currentVelocity_ = 0.0;
    currentAcceleration_ = 0.0;
}

Node::~Node() {}

int Node::getIndex()
{
    return index_;
}

bounded_vector<double, 2> Node::getInitialCoordinate()
{
    return initialCoordinate_;
}

bounded_vector<double, 2> Node::getPastCoordinate()
{
    return pastCoordinate_;
}

bounded_vector<double, 2> Node::getPastVelocity()
{
    return pastVelocity_;
}

bounded_vector<double, 2> Node::getPastAcceleration()
{
    return pastAcceleration_;
}

bounded_vector<double, 2> Node::getCurrentCoordinate()
{
    return currentCoordinate_;
}

bounded_vector<double, 2> Node::getCurrentVelocity()
{
    return currentVelocity_;
}

bounded_vector<double, 2> Node::getCurrentAcceleration()
{
    return currentAcceleration_;
}

void Node::setPastCoordinate(const bounded_vector<double, 2> &pastCoordinate)
{
    pastCoordinate_ = pastCoordinate;
}

void Node::setPastVelocity(const bounded_vector<double, 2> &pastVelocity)
{
    pastVelocity_ = pastVelocity;
}

void Node::setPastAcceleration(const bounded_vector<double, 2> &pastAcceleration)
{
    pastAcceleration_ = pastAcceleration;
}

void Node::setCurrentCoordinate(const bounded_vector<double, 2> &currentCoordinate)
{
    currentCoordinate_ = currentCoordinate;
}

void Node::setCurrentVelocity(const bounded_vector<double, 2> &currentVelocity)
{
    currentVelocity_ = currentVelocity;
}

void Node::setCurrentAcceleration(const bounded_vector<double, 2> &currentAcceleration)
{
    currentAcceleration_ = currentAcceleration;
}