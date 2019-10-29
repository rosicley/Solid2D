#pragma once
#include "Node.h"

class NeumannCondition
{
public:
    NeumannCondition();

    NeumannCondition(Node *node,
                       const int &direction,
                       const double &value);

    ~NeumannCondition();

    Node *getNode();

    int getDirection();

    double getValue();

private:
    Node *node_;

    int direction_;

    double value_;
};