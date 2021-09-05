#pragma once

#include <math.h>
#include "UtilMethods.h"

struct XYZPoint
{
    double X;
    double Y;
    double Z;

    XYZPoint(double x, double y, double z)
        : X(x), Y(y), Z(z)
    { };

    XYZPoint()
        : X(0.0), Y(0.0), Z(0.0)
    { };

    ~XYZPoint()
    { };

    double DistanceTo(XYZPoint pt)
    {
        return sqrt(((X - pt.X) * (X - pt.X)) + ((Y - pt.Y) * (Y - pt.Y)) + ((Z - pt.Z) * (Z - pt.Z)));
    };

    bool operator==(const XYZPoint& that)
    {
        return Utils::AreEqual(this->X, that.X) && 
            Utils::AreEqual(this->Y, that.Y) && 
            Utils::AreEqual(this->Z, that.Z);
    }
};
