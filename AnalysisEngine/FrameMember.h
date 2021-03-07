#pragma once
#include <vector>
#include <memory>
#include "Material.h"
#include "Section.h"
#include "Element.h"
#include "Matrix.h"
#include "Hinge.h"
#include "Node.h"

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
    bool isLumpedMassMatrix = true;

    std::shared_ptr<Matrix<double>> LocalCoordinateStiffnessMatrix;
    std::shared_ptr<Matrix<double>> LocalCoordinateMassMatrix;
    std::shared_ptr<Matrix<double>> RotationMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateStiffnessMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateMassMatrix;

    FrameMember(unsigned int elmIndex, std::shared_ptr<Node> iNode, std::shared_ptr<Node> jNode, 
        std::shared_ptr<Section> section, std::shared_ptr<Material> material, std::shared_ptr<Hinge> iEndHinge, std::shared_ptr<Hinge> jEndHinge, bool isLumpedMassMatrix);
    FrameMember();
    ~FrameMember();

    unsigned int GetElementIndex() { return this->ElementIndex; };
    unsigned int GetNumberOfDoF() { return 12; };
    ElmType::ElementType GetElementType() { return this->Type; };
    std::shared_ptr<Matrix<double>> GetLocalCoordinateStiffnessMatrix() { return this->LocalCoordinateStiffnessMatrix; };
    std::shared_ptr<Matrix<double>> GetLocalCoordinateMassMatrix() { return this->LocalCoordinateMassMatrix; };
    std::shared_ptr<Matrix<double>> GetGlobalCoordinateStiffnessMatrix() { return this->GlobalCoordinateStiffnessMatrix; };
    std::shared_ptr<Matrix<double>> GetGlobalCoordinateMassMatrix() { return this->GlobalCoordinateMassMatrix; };
    std::shared_ptr<Matrix<double>> GetRotationMatrix() { return this->RotationMatrix; };
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