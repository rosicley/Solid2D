#include "DirichletCondition.h"

DirichletCondition::DirichletCondition() {}

DirichletCondition::DirichletCondition(Node *node,
                                       const int &direction,
                                       const double &value)
{
    node_ = node;
    direction_ = direction;
    value_ = value;
}

DirichletCondition::~DirichletCondition() {}

Node *DirichletCondition::getNode()
{
    return node_;
}

int DirichletCondition::getDirection()
{
    return direction_;
}

double DirichletCondition::getValue()
{
    return value_;
}