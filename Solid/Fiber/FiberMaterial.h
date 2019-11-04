#pragma once
#include <iostream>

class FiberMaterial
{
public:
    FiberMaterial();

    FiberMaterial(const int &index,
             const double &young,
             const double &plastStrain,
             const double &hardeningModulus,
             const double &density);

    ~FiberMaterial();

    int getIndex();

    double getYoung();

    double getDensity();

    double getPlastStrain();

    double getHardeningModulus();


private:
    int index_;

    double young_;

    double density_;

    double plastStrain_;

    double hardeningModulus_;
};
