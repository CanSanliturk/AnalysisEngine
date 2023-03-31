#pragma once
#pragma once
#include "Material.h"
#include "Element.h"

class SerendipityShell : public Element
{
public:
    unsigned int ElementIndex;
    std::vector<std::shared_ptr<Node>> Nodes;
    std::shared_ptr<Material> SerendipityMaterial;
    double Thickness = 0;
    ElmType::ElementType Type;
    bool isLumpedMassMatrix = true;

    std::shared_ptr<Matrix<double>> LocalCoordinateMassMatrix;
    std::shared_ptr<Matrix<double>> LocalCoordinateStiffnessMatrix;
    std::shared_ptr<Matrix<double>> LocalCoordinateDampingMatrix;
    std::shared_ptr<Matrix<double>> RotationMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateMassMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateStiffnessMatrix;
    std::shared_ptr<Matrix<double>> GlobalCoordinateDampingMatrix;
    std::shared_ptr<Matrix<double>> K1212;

    SerendipityShell(unsigned int elmIndex,
        std::shared_ptr<Node> iNode, std::shared_ptr<Node> jNode, std::shared_ptr<Node> kNode, std::shared_ptr<Node> lNode,
        std::shared_ptr<Node> ijNode, std::shared_ptr<Node> jkNode, std::shared_ptr<Node> klNode, std::shared_ptr<Node> liNode,
        std::shared_ptr<Material> material, double thickness, bool isLumpedMassMatrix = true,
        double rayleighDampingMassMultiplier = 0, double rayleighDampingStiffnessMultiplier = 0);
    SerendipityShell();
    ~SerendipityShell();

    unsigned int GetElementIndex();
    unsigned int GetNumberOfDoF();
    ElmType::ElementType GetElementTdsaype();
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


