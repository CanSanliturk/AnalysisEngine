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

    std::shared_ptr<Matrix<double>> LocalCoordinateMassMatrix;
    std::shared_ptr<Matrix<double>> LocalCoordinateStiffnessMatrix;
    std::shared_ptr<Matrix<double>> LocalCoordinateDampingMatrix;
    std::shared_ptr<Matrix<double>> RotationMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateMassMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateStiffnessMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateDampingMatrix;
    std::shared_ptr<Matrix<double>> ElementLoads;

    FrameMember(unsigned int elmIndex, std::shared_ptr<Node> iNode, std::shared_ptr<Node> jNode, 
        std::shared_ptr<Section> section, std::shared_ptr<Material> material, bool isLumpedMassMatrix,
        std::shared_ptr<Matrix<double>> elementLoads = nullptr,
        std::shared_ptr<Hinge> iEndHinge = nullptr, std::shared_ptr<Hinge> jEndHinge = nullptr,
        double rayleighDampingMassMultiplier = 0, double rayleighDampingStiffnessMultiplier = 0);
    FrameMember();
    ~FrameMember();

    unsigned int GetElementIndex();
    unsigned int GetNumberOfDoF();
    ElmType::ElementType GetElementType();
    std::shared_ptr<Matrix<double>> GetLocalCoordinateMassMatrix();
    std::shared_ptr<Matrix<double>> GetLocalCoordinateStiffnessMatrix();
    std::shared_ptr<Matrix<double>> GetLocalCoordinateDampingMatrix();
    std::shared_ptr<Matrix<double>> GetGlobalCoordinateMassMatrix();
    std::shared_ptr<Matrix<double>> GetGlobalCoordinateStiffnessMatrix();
    std::shared_ptr<Matrix<double>> GetGlobalCoordinateDampingMatrix();
    std::shared_ptr<Matrix<double>> GetRotationMatrix();
    std::shared_ptr<Matrix<double>> GetElementLoads();
    std::vector<std::shared_ptr<Node>> GelElementNodes();

private:
    void AssembleElementLocalMassMatrix();
    void AssembleElementLocalStiffnessMatrix();
    void AssembleElementLocalDampingMatrix(double mult1, double mult2);
    void AssembleElementRotationMatrix();
    void AssembleElementGlobalMassMatrix();
    void AssembleElementGlobalStiffnessMatrix();
    void AssembleElementGlobalDampingMatrix(double mult1, double mult2);
};
