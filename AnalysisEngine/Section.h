#pragma once

#include <vector>
#include "XYPoint.h"
#include "Piece.h"

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
	{
		this->Area = area;
		this->Inertia11 = inertia11;
		this->Inertia22 = inertia22;
		this->Inertia12 = inertia12;
	};

	Section()
	{
		this->Area = -123.456678;
		this->Inertia11 = -123.456678;
		this->Inertia22 = -123.456678;
		this->Inertia12 = -123.456678;
	};

	~Section()
	{

	};

private:

};