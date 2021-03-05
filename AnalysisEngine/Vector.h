#pragma once

#include "XYZPoint.h"

struct Vector
{
    double X;
    double Y;
    double Z;
    double Length;

    double AngleTo(Vector v);

    /// <summary>
    /// Creates a vector starting from origin
    /// </summary>
    Vector(XYZPoint Point);

    template <class T>
    Vector(T x, T y, T z);

    /// <summary>
    /// Creates a vector starting from I-End to J-End
    /// </summary>
    Vector(XYZPoint IEnd, XYZPoint JEnd);

    Vector();

    ~Vector();

    Vector operator+(const Vector& v);
    Vector operator-(const Vector& v);
};

