#pragma once
#include "XYData.h"

class Material
{
public:
    double E;
    double G;
    double PoissonsRatio;
    double UnitWeight;
    XYData StrainStressCurve;

    Material(double e, double v, double rho);
    Material(double e, double v, double rho, XYData sigmaEps);
    Material();
    ~Material();

    double getStressAt(double eps);
    double getTangentStiffnessAt(double eps);
    double getSecantStiffnessAt(double eps);
};