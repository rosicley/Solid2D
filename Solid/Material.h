#pragma once
#include <iostream>

class Material
{
public:
    Material();

    Material(const int &index,
             const double &young,
             const double &poisson,
             const double &density);

    ~Material();

    int getIndex();

    double getYoung();

    double getPoisson();

    double getDensity();

private:
    int index_;

    double young_;

    double poisson_;

    double density_;
};
