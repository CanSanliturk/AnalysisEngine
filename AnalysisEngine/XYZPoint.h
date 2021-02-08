#pragma once

#include <math.h>

struct XYZPoint
{
	double X;
	double Y;
	double Z;

	XYZPoint(double x, double y, double z)
	{
		this->X = x;
		this->Y = y;
		this->Z = z;
	};

	XYZPoint()
	{
		this->X = 0;
		this->Y = 0;
		this->Z = 0;
	};

	~XYZPoint()
	{

	};

	double DistanceTo(XYZPoint pt)
	{
		return sqrt(((this->X - pt.X) * (this->X - pt.X)) + ((this->Y - pt.Y) * (this->Y - pt.Y)) + ((this->Z - pt.Z) * (this->Z - pt.Z)));
	};
};
