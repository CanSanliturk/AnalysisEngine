#pragma once
#include <vector>
#include "Hinge.h"
#include "Node.h"
#include "Section.h"
#include "Material.h"
#include "Element.h"

class FrameMember : public Element
{
public:
    unsigned int ElementIndex;
    Node* Nodes[2];
    Hinge* Hinges[2];
    Section* FrameSection;
    Material* FrameMaterial;
    ElmType::ElementType Type;
    double Length;

    double LocalCoordinateStiffnessMatrix[12][12];
    double LocalCoordinateMassMatrix[12][12];
    double RotationMatrix[12][12];
    double GlobalCoordinateStiffnessMatrix[12][12];
    double GlobalCoordinateMassMatrix[12][12];

    FrameMember(unsigned int elmIndex, Node* iNode, Node* jNode, Section* section, Material* material, Hinge* iEndHinge, Hinge* jEndHinge);
    FrameMember();
    ~FrameMember();

    unsigned int GetElementIndex() { return this->ElementIndex; };
    unsigned int GetNumberOfDoF() { return 12; };
    ElmType::ElementType GetElementType() { return this->Type; };
    void* GetLocalCoordinateStiffnessMatrix() { return &this->LocalCoordinateStiffnessMatrix; };
    void* GetLocalCoordinateMassMatrix() { return &this->LocalCoordinateMassMatrix; };
    void* GetGlobalCoordinateStiffnessMatrix() { return &this->GlobalCoordinateStiffnessMatrix; };
    void* GetGlobalCoordinateMassMatrix() { return &this->GlobalCoordinateMassMatrix; };
    void* GetRotationMatrix() { return &this->RotationMatrix; };
    std::vector<Node*> GelElementNodes()
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