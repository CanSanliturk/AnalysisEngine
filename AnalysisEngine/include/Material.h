#pragma once
class Material
{
public:
    double E;
    double G;
    double PoissonsRatio;
    double UnitWeight;

    Material(double e, double v, double rho)
        : E(e), PoissonsRatio(v), UnitWeight(rho)
    { 
        G = E / (2 * (1 + PoissonsRatio));
    };

    Material()
        : E(0.0), G(0.0),
        PoissonsRatio(0.0), UnitWeight(0.0)
    { };

    ~Material()
    { };
};