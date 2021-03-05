#pragma once
#include <vector>
#include <memory>
#include "Hinge.h"
#include "Node.h"
#include "Section.h"
#include "Material.h"
#include "Element.h"

class FrameMember : public Element
{
public:
    unsigned int ElementIndex;
    std::vector<std::shared_ptr<Node>> Nodes;
    std::vector<std::shared_ptr<Hinge>> Hinges;
    std::shared_ptr<Section> FrameSection;
    std::shared_ptr<Material> FrameMaterial;
    ElmType::ElementType Type;
    double Length;

    double LocalCoordinateStiffnessMatrix[12][12];
    double LocalCoordinateMassMatrix[12][12];
    double RotationMatrix[12][12];
    double GlobalCoordinateStiffnessMatrix[12][12];
    double GlobalCoordinateMassMatrix[12][12];

    FrameMember(unsigned int elmIndex, std::shared_ptr<Node> iNode, std::shared_ptr<Node> jNode, 
        std::shared_ptr<Section> section, std::shared_ptr<Material> material, std::shared_ptr<Hinge> iEndHinge, std::shared_ptr<Hinge> jEndHinge);
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
    std::vector<std::shared_ptr<Node>> GelElementNodes()
    {
        std::vector<std::shared_ptr<Node>> retVal;
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