#pragma once
#include "Material.h"
#include "Element.h"

enum class MembraneType
{
    NONE = 0,
    Bilinear = 1,
    Incompatible = 2,
    Drilling = 3
};

enum class PlateType
{
    NONE = 0,
    Thin = 1,
    Thick = 2
};

class ShellMember : public Element
{
public:
    unsigned int ElementIndex;
    std::vector<std::shared_ptr<Node>> Nodes;
    std::shared_ptr<Material> ShellMaterial;
    double Thickness = 0;
    ElmType::ElementType Type;
    bool isLumpedMassMatrix = true;
    MembraneType membraneType = MembraneType::NONE;
    PlateType plateType = PlateType::NONE;

    std::shared_ptr<Matrix<double>> LocalCoordinateMassMatrix;
    std::shared_ptr<Matrix<double>> LocalCoordinateStiffnessMatrix;
    std::shared_ptr<Matrix<double>> LocalCoordinateDampingMatrix;
    std::shared_ptr<Matrix<double>> RotationMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateMassMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateStiffnessMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateDampingMatrix;

    ShellMember(unsigned int elmIndex, 
        std::shared_ptr<Node> iNode, std::shared_ptr<Node> jNode, std::shared_ptr<Node> kNode, std::shared_ptr<Node> lNode,
        std::shared_ptr<Material> material, double thickness, MembraneType memType, PlateType pltType, bool isLumpedMassMatrix = true, 
        double rayleighDampingMassMultiplier = 0, double rayleighDampingStiffnessMultiplier = 0);
    ShellMember();
    ~ShellMember();

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

