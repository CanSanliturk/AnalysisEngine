#pragma once

#include <vector>
#include "Node.h"

/// <summary>
/// Boundary condition for fixed forces and moments.
/// </summary>
class Hinge
{
public:
    double ForceX;
    double ForceY;
    double ForceZ;
    double MomentX;
    double MomentY;
    double MomentZ;

    bool IsReleasedForceX = false;
    bool IsReleasedForceY = false;
    bool IsReleasedForceZ = false;
    bool IsReleasedMomentX = false;
    bool IsReleasedMomentY = false;
    bool IsReleasedMomentZ = false;

    Hinge(std::vector<bool> isReleased, std::vector<double> release)
        : IsReleasedForceX(isReleased[0]), IsReleasedForceY(isReleased[1]), IsReleasedForceZ(isReleased[2]),
        IsReleasedMomentX(isReleased[3]), IsReleasedMomentY(isReleased[4]), IsReleasedMomentZ(isReleased[5]),
        ForceX(release[0]), ForceY(release[1]), ForceZ(release[2]),
        MomentX(release[3]), MomentY(release[4]), MomentZ(release[5])
    { };

    Hinge()
        : IsReleasedForceX(false), IsReleasedForceY(false), IsReleasedForceZ(false),
        IsReleasedMomentX(false), IsReleasedMomentY(false), IsReleasedMomentZ(false),
        ForceX(0.0), ForceY(0.0), ForceZ(0.0),
        MomentX(0.0), MomentY(0.0), MomentZ(0.0)
    { };

    ~Hinge()
    { };

};
