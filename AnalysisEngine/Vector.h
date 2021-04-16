#pragma once

#include "XYZPoint.h"

class Vector
{
public:
    double X;
    double Y;
    double Z;
    double Length;

    double AngleTo(Vector v);

    /// <summary>
    /// Creates a vector starting from origin
    /// </summary>
    Vector(XYZPoint Point);

    Vector(double x, double y, double z);

    /// <summary>
    /// Creates a vector starting from I-End to J-End
    /// </summary>
    Vector(XYZPoint IEnd, XYZPoint JEnd);

    Vector();

    ~Vector();

    Vector getUnitVector();

    Vector operator+(const Vector& v);
    Vector operator-(const Vector& v);
    Vector operator*(const Vector& v);
    Vector operator*(const double mult);
    bool operator==(const Vector& v);

};

