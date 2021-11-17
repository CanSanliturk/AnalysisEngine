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

double Material::getTangentModulusAt(double eps)
{
    if (abs(eps) < 1e-8)
        return this->E;

    auto epsBefore = 0.99 * eps;
    auto epsAfter = 1.01 * eps;
    auto stressBefore = StrainStressCurve.getY(epsBefore);
    auto stressAfter = StrainStressCurve.getY(epsAfter);
    return (stressAfter - stressBefore) / (epsAfter - epsBefore);
}

double Material::getSecantModulusAt(double eps)
{
    if (abs(eps) < 1.0e-8)
        return this->E;
    auto stressAtEps = StrainStressCurve.getY(eps);
    return stressAtEps / eps;
}

