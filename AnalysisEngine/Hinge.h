#pragma once

#include <vector>
#include "Node.h"

/// <summary>
/// Boundary condition for fixed forces and moments.
/// </summary>
class Hinge
{
public:
	Node* RestrainedNode;

	double ForceX;
	double ForceY;
	double ForceZ;
	double MomentX;
	double MomentY;
	double MomentZ;

	bool IsRestraintForceX = false;
	bool IsRestraintForceY = false;
	bool IsRestraintForceZ = false;
	bool IsRestraintMomentX = false;
	bool IsRestraintMomentY = false;
	bool IsRestraintMomentZ = false;

	Hinge(Node* releasedNode, std::vector<bool> isReleased, std::vector<double> release)
	{
		this->RestrainedNode = releasedNode;

		this->IsRestraintForceX = isReleased.at(0);
		this->IsRestraintForceY = isReleased.at(1);
		this->IsRestraintForceZ = isReleased.at(2);
		this->IsRestraintMomentX = isReleased.at(3);
		this->IsRestraintMomentY = isReleased.at(4);
		this->IsRestraintMomentZ = isReleased.at(5);

		this->ForceX = release.at(0);
		this->ForceY = release.at(1);
		this->ForceZ = release.at(2);
		this->MomentX = release.at(3);
		this->MomentY = release.at(4);
		this->MomentZ = release.at(5);
	};

	Hinge()
	{
		Node node;
		this->RestrainedNode = &node;

		this->IsRestraintForceX = false;
		this->IsRestraintForceY = false;
		this->IsRestraintForceZ = false;
		this->IsRestraintMomentX = false;
		this->IsRestraintMomentY = false;
		this->IsRestraintMomentZ = false;

		this->ForceX = -123.456789;
		this->ForceY = -123.456789;
		this->ForceZ = -123.456789;
		this->MomentX = -123.456789;
		this->MomentY = -123.456789;
		this->MomentZ = -123.456789;
	};

	~Hinge()
	{

	};

};