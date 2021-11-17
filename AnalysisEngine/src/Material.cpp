#include "Material.h"

Material::Material()
    : E(0), G(0), PoissonsRatio(0), UnitWeight(0) { }

Material::Material(double e, double v, double rho)
    : E(e), PoissonsRatio(v), UnitWeight(rho) {
    G = E / (2 * (1 + PoissonsRatio));
}

Material::Material(double e, double v, double rho, XYData sigmaEps) 
    : E(e), PoissonsRatio(v), UnitWeight(rho), StrainStressCurve(sigmaEps) {
    G = E / (2 * (1 + PoissonsRatio));
}

Material::~Material() {}

double Material::getStressAt(double eps)
{
    return this->StrainStressCurve.getY(eps);
}

double Material::getTangentStiffnessAt(double eps)
{
    auto epsBefore = 0.99 * eps;
    auto epsAfter = 1.01 * eps;
    auto stressBefore = StrainStressCurve.getY(epsBefore);
    auto stressAfter = StrainStressCurve.getY(epsAfter);
    return (stressAfter - stressBefore) / (epsAfter - epsBefore);
}

double Material::getSecantStiffnessAt(double eps)
{
    auto initialStrain = StrainStressCurve.getXData()[0];
    auto initialStress = StrainStressCurve.getYData()[0];
    auto stressAtEps = StrainStressCurve.getY(eps);
    return (stressAtEps - initialStress) / (eps - initialStrain);
}

