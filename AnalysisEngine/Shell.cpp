#include "Shell.h"

ShellMember::ShellMember(unsigned int elmIndex,
    std::shared_ptr<Node> iNode, std::shared_ptr<Node> jNode, std::shared_ptr<Node> kNode, std::shared_ptr<Node> lNode,
    std::shared_ptr<Material> material, double thickness, bool isMembrane, bool isPlate,
    bool isLumpedMassMatrix, double rayleighDampingMassMultiplier, double rayleighDampingStiffnessMultiplier)
{
    ElementIndex = elmIndex;
    
    Nodes.resize(4);
    Nodes[0] = iNode;
    Nodes[1] = jNode;
    Nodes[2] = kNode;
    Nodes[3] = lNode;
    
    ShellMaterial = material;

    Thickness = thickness;

    Type = (ElmType::ElementType::Shell);

    isLumpedMassMatrix = isLumpedMassMatrix;
    isMembraneAction = isMembrane;
    isPlateAction = isPlate;
    
    LocalCoordinateMassMatrix = std::make_shared<Matrix<double>>(24);
    LocalCoordinateStiffnessMatrix = std::make_shared<Matrix<double>>(24);
    LocalCoordinateDampingMatrix = std::make_shared<Matrix<double>>(24);
    RotationMatrix = std::make_shared<Matrix<double>>(24);
    GlobalCoordinateStiffnessMatrix = std::make_shared<Matrix<double>>(24);
    GlobalCoordinateMassMatrix = std::make_shared<Matrix<double>>(24);
    GlobalCoordinateDampingMatrix = std::make_shared<Matrix<double>>(24);

    AssembleElementLocalMassMatrix();
    AssembleElementLocalStiffnessMatrix();
    AssembleElementLocalDampingMatrix(rayleighDampingMassMultiplier, rayleighDampingStiffnessMultiplier);
    AssembleElementRotationMatrix();
    AssembleElementGlobalMassMatrix();
    AssembleElementGlobalStiffnessMatrix();
    AssembleElementGlobalDampingMatrix(rayleighDampingMassMultiplier, rayleighDampingStiffnessMultiplier);

    iNode->ConnectedElements.push_back(elmIndex);
    jNode->ConnectedElements.push_back(elmIndex);
    kNode->ConnectedElements.push_back(elmIndex);
    lNode->ConnectedElements.push_back(elmIndex);
}

ShellMember::ShellMember()
{
}

ShellMember::~ShellMember()
{
}

unsigned int ShellMember::GetElementIndex()
{
    return 0;
}

unsigned int ShellMember::GetNumberOfDoF()
{
    return 0;
}

ElmType::ElementType ShellMember::GetElementType()
{
    return ElmType::ElementType();
}

std::shared_ptr<Matrix<double>> ShellMember::GetLocalCoordinateMassMatrix()
{
    return std::shared_ptr<Matrix<double>>();
}

std::shared_ptr<Matrix<double>> ShellMember::GetLocalCoordinateStiffnessMatrix()
{
    return std::shared_ptr<Matrix<double>>();
}

std::shared_ptr<Matrix<double>> ShellMember::GetLocalCoordinateDampingMatrix()
{
    return std::shared_ptr<Matrix<double>>();
}

std::shared_ptr<Matrix<double>> ShellMember::GetGlobalCoordinateMassMatrix()
{
    return std::shared_ptr<Matrix<double>>();
}

std::shared_ptr<Matrix<double>> ShellMember::GetGlobalCoordinateStiffnessMatrix()
{
    return std::shared_ptr<Matrix<double>>();
}

std::shared_ptr<Matrix<double>> ShellMember::GetGlobalCoordinateDampingMatrix()
{
    return std::shared_ptr<Matrix<double>>();
}

std::shared_ptr<Matrix<double>> ShellMember::GetRotationMatrix()
{
    return std::shared_ptr<Matrix<double>>();
}


std::shared_ptr<Matrix<double>> ShellMember::GetElementLoads()
{
    return std::shared_ptr<Matrix<double>>();
}

std::vector<std::shared_ptr<Node>> ShellMember::GelElementNodes()
{
    return std::vector<std::shared_ptr<Node>>();
}

void ShellMember::AssembleElementLocalMassMatrix()
{
}

void ShellMember::AssembleElementLocalStiffnessMatrix()
{

}

void ShellMember::AssembleElementLocalDampingMatrix(double mult1, double mult2)
{
}

void ShellMember::AssembleElementRotationMatrix()
{
}

void ShellMember::AssembleElementGlobalMassMatrix()
{
}

void ShellMember::AssembleElementGlobalStiffnessMatrix()
{
}

void ShellMember::AssembleElementGlobalDampingMatrix(double mult1, double mult2)
{
}
