#pragma once

#include <vector>
#include "XYPoint.h"
#include "Piece.h"

/// <summary>
/// Axes are 1 and 2
/// </summary>
///				2 (Local Y-Axis)
///   __________|___________
///  |			|			|
///	 |			|			|
///  |			|			|
///  |			|			|
///	 |			|			|
///  |----------|-----------|----1 (Local Z-Axis)
///  |			|			|
///	 |			|			|
///  |			|			|
///  |			|			|
///	 |			|			|
///  |__________|___________|

class Section
{
public:
    std::vector<Piece> Pieces;
    std::vector<XYPoint> Vertices;
    double Area;
    double Inertia11;
    double Inertia22;
    double Inertia12;

    Section(double area, double inertia11, double inertia22, double inertia12)
        : Area(area), Inertia11(inertia11), Inertia22(inertia22), Inertia12(inertia12)
    { };

    Section()
        : Area(0.0), Inertia11(0.0), Inertia22(0.0), Inertia12(0.0)
    { };

    ~Section()
    { };

private:

};