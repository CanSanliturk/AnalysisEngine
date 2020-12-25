#pragma once

class XYPoint
{
public:
	double X;
	double Y;

	XYPoint(double x, double y)
	{
		this->X = x;
		this->Y = y;
	};

	XYPoint()
	{
		this->X = -123.456789;
		this->Y = -123.456789;
	};

	~XYPoint()
	{

	};

private:
};