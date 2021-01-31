#pragma once

#include <vector>
#include "Node.h"

/// <summary>
/// Boundary condition for fixed forces and moments.
/// </summary>
class Hinge
{
public:
	double ForceX;
	double ForceY;
	double ForceZ;
	double MomentX;
	double MomentY;
	double MomentZ;

	bool IsReleasedForceX = false;
	bool IsReleasedForceY = false;
	bool IsReleasedForceZ = false;
	bool IsReleasedMomentX = false;
	bool IsReleasedMomentY = false;
	bool IsReleasedMomentZ = false;

	Hinge(std::vector<bool> isReleased, std::vector<double> release)
	{
		this->IsReleasedForceX = isReleased.at(0);
		this->IsReleasedForceY = isReleased.at(1);
		this->IsReleasedForceZ = isReleased.at(2);
		this->IsReleasedMomentX = isReleased.at(3);
		this->IsReleasedMomentY = isReleased.at(4);
		this->IsReleasedMomentZ = isReleased.at(5);

		this->ForceX = release.at(0);
		this->ForceY = release.at(1);
		this->ForceZ = release.at(2);
		this->MomentX = release.at(3);
		this->MomentY = release.at(4);
		this->MomentZ = release.at(5);
	};

	Hinge()
	{
		this->IsReleasedForceX = false;
		this->IsReleasedForceY = false;
		this->IsReleasedForceZ = false;
		this->IsReleasedMomentX = false;
		this->IsReleasedMomentY = false;
		this->IsReleasedMomentZ = false;

		this->ForceX = 0.0;
		this->ForceY = 0.0;
		this->ForceZ = 0.0;
		this->MomentX = 0.0;
		this->MomentY = 0.0;
		this->MomentZ = 0.0;
	};

	~Hinge()
	{

	};

};