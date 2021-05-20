#include "Vector.h"
#include "UtilMethods.h"
#include <Math.h>

Vector::Vector(XYZPoint Point)
    : X(Point.X), Y(Point.Y), Z(Point.Z),
    Length(sqrt((X * X) + (Y * Y) + (Z * Z)))
{ }

Vector::Vector(double x, double y, double z)
    : X(x), Y(y), Z(z),
    Length(sqrt((X* X) + (Y * Y) + (Z * Z)))
{ }

Vector::Vector(XYZPoint IEnd, XYZPoint JEnd)
    : X((JEnd.X - IEnd.X)), Y((JEnd.Y - IEnd.Y)),
    Z((JEnd.Z - IEnd.Z)), Length(JEnd.DistanceTo(IEnd))
{ }

Vector::Vector()
    : X(0.0), Y(0.0), Z(0.0), Length(0.0)
{ }

Vector::~Vector()
{ }

Vector Vector::getUnitVector()
{
    Vector v;
    v.X = this->X / this->Length;
    v.Y = this->Y / this->Length;
    v.Z = this->Z / this->Length;
    v.Length = 1;
    return v;
}

double Vector::AngleTo(Vector v)
{
    double pi = 3.141592653589793;
    auto x = this->X;
    auto y = this->Y;
    auto z = this->Z;

    return acos(((this->X * v.X) + (this->Y * v.Y) + (this->Z * v.Z)) / (this->Length * v.Length));
}

Vector Vector::operator+(const Vector& v)
{
    Vector vec;
    XYZPoint Point(this->X + v.X, this->Y + v.Y, this->Z + v.Z);
    XYZPoint origin(0.0, 0.0, 0.0);

    vec.X = Point.X;
    vec.Y = Point.Y;
    vec.Z = Point.Z;
    vec.Length = Point.DistanceTo(origin);

    return vec;
}

Vector Vector::operator-(const Vector& v)
{
    Vector vec;
    XYZPoint Point(this->X - v.X, this->Y - v.Y, this->Z - v.Z);
    XYZPoint origin(0.0, 0.0, 0.0);

    vec.X = Point.X;
    vec.Y = Point.Y;
    vec.Z = Point.Z;
    vec.Length = Point.DistanceTo(origin);

    return vec;
}

Vector Vector::operator*(const Vector& v)
{
    Vector vec;

    // Multiply
    vec.X = (this->Y * v.Z) - (this->Z * v.Y);
    vec.Y = (this->Z * v.X) - (this->X * v.Z);
    vec.Z = (this->X * v.Y) - (this->Y * v.X);
    XYZPoint p(vec.X, vec.Y, vec.Z);
    XYZPoint origin(0.0, 0.0, 0.0);
    auto length = p.DistanceTo(origin);

    // Make it unit vector
    vec.X /= length;
    vec.Y /= length;
    vec.Z /= length;
    XYZPoint p2(vec.X, vec.Y, vec.Z);
    vec.Length = p2.DistanceTo(origin);

    return vec;
}

Vector Vector::operator*(const double mult)
{
    // Multiply
    Vector vec;
    vec.X = mult * this->X;
    vec.Y = mult * this->Y;
    vec.Z = mult * this->Z;
    XYZPoint p(vec.X, vec.Y, vec.Z);
    XYZPoint origin(0.0, 0.0, 0.0);
    vec.Length = p.DistanceTo(origin);
    return vec;
}

bool Vector::operator==(const Vector& v)
{
    return (Utils::AreEqual(this->X, v.X) && Utils::AreEqual(this->Y, v.Y) && Utils::AreEqual(this->Z, v.Z));
}
