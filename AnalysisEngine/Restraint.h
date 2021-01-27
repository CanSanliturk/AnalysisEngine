#pragma once

#include <vector>
#include "Node.h"

class Restraint
{
public:
	Node* RestrainedNode;

	double TranslationX;
	double TranslationY;
	double TranslationZ;
	double RotationX;
	double RotationY;
	double RotationZ;

	bool IsRestraintTranslationX = false;
	bool IsRestraintTranslationY = false;
	bool IsRestraintTranslationZ = false;
	bool IsRestraintRotationX = false;
	bool IsRestraintRotationY = false;
	bool IsRestraintRotationZ = false;

	Restraint(Node* resNode, std::vector<bool> isRest, std::vector<double> rest)
	{
		this->RestrainedNode = resNode;

		this->IsRestraintTranslationX = isRest.at(0);
		this->IsRestraintTranslationY = isRest.at(1);
		this->IsRestraintTranslationZ = isRest.at(2);
		this->IsRestraintRotationX = isRest.at(3);
		this->IsRestraintRotationY = isRest.at(4);
		this->IsRestraintRotationZ = isRest.at(5);

		this->TranslationX = rest.at(0);
		this->TranslationY = rest.at(1);
		this->TranslationZ = rest.at(2);
		this->RotationX = rest.at(3);
		this->RotationY = rest.at(4);
		this->RotationZ = rest.at(5);
	};

	Restraint()
	{
		Node node;
		this->RestrainedNode = &node;

		this->IsRestraintTranslationX = false;
		this->IsRestraintTranslationY = false;
		this->IsRestraintTranslationZ = false;
		this->IsRestraintRotationX = false;
		this->IsRestraintRotationY = false;
		this->IsRestraintRotationZ = false;

		this->TranslationX = -123.456789;
		this->TranslationY = -123.456789;
		this->TranslationZ = -123.456789;
		this->RotationX = -123.456789;
		this->RotationY = -123.456789;
		this->RotationZ = -123.456789;
	};

	~Restraint()
	{

	};

};