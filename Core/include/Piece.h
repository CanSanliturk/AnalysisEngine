#pragma once

#include <vector>
#include "XYPoint.h"

/// <summary>
/// Axes are 1 and 2
/// </summary>
///				2
///   __________|___________
///  |			|			|
///	 |			|			|
///  |			|			|
///  |			|			|
///	 |			|			|
///  |----------|-----------|----1
///  |			|			|
///	 |			|			|
///  |			|			|
///  |			|			|
///	 |			|			|
///  |__________|___________|

class Piece
{
public:
	std::vector<XYPoint> Vertices;
	std::vector<XYPoint> Holes;
	double Area;
	double Inertia11;
	double Inertia22;
	double Inertia12;

private:

};
