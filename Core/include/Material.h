#pragma once
#include "XYData.h"
#include <memory>

class Material
{
public:
    double E;
    double G;
    double PoissonsRatio;
    double UnitWeight;
    double InitialModulus;
    std::shared_ptr<XYData> StrainStressCurve;

    Material(double e, double v, double rho);
    Material(double e, double v, double rho, std::shared_ptr<XYData> sigmaEps);
    Material();
    ~Material();

    double getStressAt(double eps);
    double getTangentModulusAt(double eps);
    double getSecantModulusAt(double eps);
};