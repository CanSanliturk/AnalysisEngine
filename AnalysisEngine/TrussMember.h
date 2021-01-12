#pragma once
#include "Node.h"
#include "Section.h"
#include "Material.h"
#include "Element.h"
#include <vector>

class TrussMember : public Element
{
public:
	unsigned int ElementIndex;
	Node* Nodes[2];
	Section* TrussSection;
	Material* TrussMaterial;
	ElmType::ElementType Type;
	double Length;

	double LocalCoordinateStiffnessMatrix[6][6];
	double LocalCoordinateMassMatrix[6][6];
	double RotationMatrix[6][6];
	double GlobalCoordinateStiffnessMatrix[6][6];
	double GlobalCoordinateMassMatrix[6][6];

	TrussMember(unsigned int elmIndex, Node* iNode, Node* jNode, Section* section, Material* material);
	TrussMember();
	~TrussMember();

	unsigned int GetElementIndex() override { return this->ElementIndex; };
	unsigned int GetNumberOfDoF() override { return 12; };
	ElmType::ElementType GetElementType() override { return this->Type; };
	void* GetLocalCoordinateStiffnessMatrix() override { return &this->LocalCoordinateStiffnessMatrix; };
	void* GetLocalCoordinateMassMatrix() override { return &this->LocalCoordinateMassMatrix; };
	void* GetGlobalCoordinateStiffnessMatrix() override { return &this->GlobalCoordinateStiffnessMatrix; };
	void* GetGlobalCoordinateMassMatrix() override { return &this->GlobalCoordinateMassMatrix; };
	void* GetRotationMatrix() override { return &this->RotationMatrix; };
	std::vector<Node*> GelElementNodes() override
	{
		std::vector<Node*> retVal;
		retVal.push_back(this->Nodes[0]);
		retVal.push_back(this->Nodes[1]);
		return retVal;
	};

private:
	void AssembleElementLocalStiffnessMatrix();
	void AssembleElementLocalMassMatrix();
	void AssembleElementRotationMatrix();
	void AssembleElementGlobalStiffnessMatrix();
	void AssembleElementGlobalMassMatrix();

};