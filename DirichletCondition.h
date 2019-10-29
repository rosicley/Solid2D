#pragma once
#include "Node.h"

class DirichletCondition
{
public:
    DirichletCondition();

    DirichletCondition(Node *node,
                       const int &direction,
                       const double &value);

    ~DirichletCondition();

    Node *getNode();

    int getDirection();

    double getValue();

private:
    Node *node_;

    int direction_;

    double value_;
};