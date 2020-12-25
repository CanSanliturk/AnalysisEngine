#pragma once
class Material 
{
public:
	double E;
	double G;
	double PoissonsRatio;
	double UnitWeight;

	Material(double e, double v, double rho)
	{
		this->E = e;
		this->PoissonsRatio = v;
		this->G = (this->E / (2 * (1 + this->PoissonsRatio)));
		this->UnitWeight = rho;
	};

	Material()
	{
		this->E = -123.456789;
		this->G = -123.456789;
		this->PoissonsRatio = -123.456789;
		this->UnitWeight = -123.456789;
	};

	~Material()
	{

	};

};