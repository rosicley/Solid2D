#include "FiberMaterial.h"

FiberMaterial::FiberMaterial() {}

FiberMaterial::FiberMaterial(const int &index,
                   const double &young,
                   const double &plastStrain,
                   const double &hardeningModulus,
                   const double &density)
{
    index_ = index;
    young_ = young;
    plastStrain_=plastStrain;
    hardeningModulus_=hardeningModulus;
    density_ = density;
}

FiberMaterial::~FiberMaterial() {}

int FiberMaterial::getIndex()
{
    return index_;
}

double FiberMaterial::getYoung()
{
    return young_;
}

double FiberMaterial::getDensity()
{
    return density_;
}

double FiberMaterial::getPlastStrain()
{
    return plastStrain_;
}

double FiberMaterial::getHardeningModulus()
{
    return hardeningModulus_;
}
