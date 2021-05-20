#pragma once

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
};