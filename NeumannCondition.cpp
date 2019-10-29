#include "NeumannCondition.h"

NeumannCondition::NeumannCondition() {}

NeumannCondition::NeumannCondition(Node *node,
                                   const int &direction,
                                   const double &value)
{
    node_ = node;
    direction_ = direction;
    value_ = value;
}

NeumannCondition::~NeumannCondition() {}

Node *NeumannCondition::getNode()
{
    return node_;
}

int NeumannCondition::getDirection()
{
    return direction_;
}

double NeumannCondition::getValue()
{
    return value_;
}