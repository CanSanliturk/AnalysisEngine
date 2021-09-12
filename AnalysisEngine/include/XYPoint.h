#pragma once
#include "UtilMethods.h"

struct XYPoint
{
    double X;
    double Y;

    XYPoint(double x, double y)
        : X(x), Y(y)
    { };

    XYPoint()
        : X(0.0), Y(0.0)
    { };

    ~XYPoint()
    { };

    double DistanceTo(XYPoint pt)
    {
        return sqrt(((X - pt.X) * (X - pt.X)) + ((Y - pt.Y) * (Y - pt.Y)));
    };

    bool operator==(const XYPoint& that)
    {
        return Utils::AreEqual(this->X, that.X) &&
            Utils::AreEqual(this->Y, that.Y);
    }
};