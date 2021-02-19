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
		this->E = 0;
		this->G = 0;
		this->PoissonsRatio = 0;
		this->UnitWeight = 0;
	};

	~Material()
	{

	};
};