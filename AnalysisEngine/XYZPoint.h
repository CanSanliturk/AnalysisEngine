#pragma once

#include <math.h>

class XYZPoint
{
public:
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
		this->X = -123.456789;
		this->Y = -123.456789;
		this->Z = -123.456789;
	};

	~XYZPoint()
	{

	};

	double DistanceTo(XYZPoint pt)
	{		
		return sqrt(((this->X - pt.X) * (this->X - pt.X)) + ((this->Y - pt.Y) * (this->Y - pt.Y)) + ((this->Z - pt.Z) * (this->Z - pt.Z)));
	};

private:
};
