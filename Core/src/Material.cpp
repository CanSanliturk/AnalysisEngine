#include "Material.h"
#include <iostream>

Material::Material() {
        E = 0;
        G = 0;
        PoissonsRatio = 0;
        UnitWeight = 0;
        InitialModulus = 0;
        auto ex = ExtrapolationMethod::Zero;
        StrainStressCurve = std::make_shared<XYData>(ex);
}

Material::Material(double e, double v, double rho)
    : E(e), PoissonsRatio(v), UnitWeight(rho), InitialModulus(e) {
    G = E / (2 * (1 + PoissonsRatio));
    auto ex = ExtrapolationMethod::Zero;
    StrainStressCurve = std::make_shared<XYData>(ex);
    
}

Material::Material(double e, double v, double rho, std::shared_ptr<XYData> sigmaEps) 
    : E(e), PoissonsRatio(v), UnitWeight(rho), InitialModulus(e), StrainStressCurve(sigmaEps) {
    G = E / (2 * (1 + PoissonsRatio));
}

Material::~Material() {}

double Material::getStressAt(double eps)
{
    return this->StrainStressCurve->getY(eps);
}

double Material::getTangentModulusAt(double eps)
{
    auto maxNeg = StrainStressCurve->getMaxNegative();
    auto minPos = StrainStressCurve->getMinPositive();

    if ((maxNeg < eps) &&  (eps < minPos))
        return this->E;

    auto epsBefore = 0.99 * eps;
    auto epsAfter = 1.01 * eps;
    auto stressBefore = StrainStressCurve->getY(epsBefore);
    auto stressAfter = StrainStressCurve->getY(epsAfter);
    return (stressAfter - stressBefore) / (epsAfter - epsBefore);
}

double Material::getSecantModulusAt(double eps)
{
    auto maxNeg = StrainStressCurve->getMaxNegative();
    auto minPos = StrainStressCurve->getMinPositive();

    if ((maxNeg < eps) &&  (eps < minPos))
        return this->E;

    auto stressAtEps = StrainStressCurve->getY(eps);
    return stressAtEps / eps;
}
