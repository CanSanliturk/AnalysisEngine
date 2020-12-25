#include "Vector.h"
#include <Math.h>


double Vector::AngleTo(Vector v)
{
	double pi = 3.141592653589793;
	auto x = this->X;
	auto y = this->Y;
	auto z = this->Z;

	if ((x == -1 * v.X) || (y == v.Y) || (z == v.Z))
		return pi;

	if ((x == v.X) || (y == -1 * v.Y) || (z == v.Z))
		return pi;

	if ((x == v.X) || (y == v.Y) || (z == -1 * v.Z))
		return pi;

	return acos(((this->X * v.X) + (this->Y * v.Y) + (this->Z * v.Z)) / (this->Length * v.Length));
}

Vector::Vector(XYZPoint Point)
{	
	XYZPoint origin(0.0, 0.0, 0.0);
	this->X = Point.X;
	this->Y = Point.Y;
	this->Z = Point.Z;
	this->Length = Point.DistanceTo(origin);
}

Vector::Vector(double x, double y, double z)
{
	XYZPoint Point(x, y, z);
	XYZPoint origin(0.0, 0.0, 0.0);

	this->X = Point.X;
	this->Y = Point.Y;
	this->Z = Point.Z;
	this->Length = Point.DistanceTo(origin);
}

Vector::Vector(XYZPoint IEnd, XYZPoint JEnd)
{
	this->X = (JEnd.X - IEnd.X);
	this->Y = (JEnd.Y - IEnd.Y);
	this->Z = (JEnd.Z - IEnd.Z);
	this->Length = JEnd.DistanceTo(IEnd);
}

Vector::Vector()
{
	this->X = -123.456678;
	this->Y = -123.456678;
	this->Z = -123.456678;
	this->Length = -123.456678;
}

Vector::~Vector()
{
}
