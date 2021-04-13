#include "TrussMember.h"
#include "XYZPoint.h"
#include "Vector.h"
#include "GeometryHelper.h"
#include <math.h>

constexpr double G = 9.807;

TrussMember::TrussMember(unsigned int elmIndex, std::shared_ptr<Node> iNode, std::shared_ptr<Node> jNode,
    std::shared_ptr<Section> section, std::shared_ptr<Material> material,
    double rayleighDampingMassMultiplier, double rayleighDampingStiffnessMultiplier)
{
    Nodes.resize(2);
    ElementIndex = elmIndex;
    Nodes[0] = iNode;
    Nodes[1] = jNode;
    TrussSection = section;
    TrussMaterial = material;
    Type = (ElmType::ElementType::Truss);
    Length = jNode->Coordinate.DistanceTo(iNode->Coordinate);
    LocalCoordinateStiffnessMatrix = std::make_shared<Matrix<double>>(12);
    LocalCoordinateMassMatrix = std::make_shared<Matrix<double>>(12);
    RotationMatrix = std::make_shared<Matrix<double>>(12);
    GlobalCoordinateStiffnessMatrix = std::make_shared<Matrix<double>>(12);
    GlobalCoordinateMassMatrix = std::make_shared<Matrix<double>>(12);
    AssembleElementLocalMassMatrix();
    AssembleElementLocalStiffnessMatrix();
    AssembleElementLocalDampingMatrix(rayleighDampingMassMultiplier, rayleighDampingStiffnessMultiplier);
    AssembleElementRotationMatrix();
    AssembleElementGlobalMassMatrix();
    AssembleElementGlobalStiffnessMatrix();
    AssembleElementGlobalDampingMatrix(rayleighDampingMassMultiplier, rayleighDampingStiffnessMultiplier);
    iNode->ConnectedElements.push_back(elmIndex);
    jNode->ConnectedElements.push_back(elmIndex);
}

TrussMember::TrussMember()
{
    std::shared_ptr<Node> iN;
    std::shared_ptr<Node> jN;
    std::shared_ptr<Section> s;
    std::shared_ptr<Material> m;
    this->Nodes[0] = iN;
    this->Nodes[1] = jN;
    this->TrussSection = s;
    this->TrussMaterial = m;
    this->Length = 0;
}

TrussMember::~TrussMember()
{
}

void TrussMember::AssembleElementLocalMassMatrix()
{
    auto A = this->TrussSection->Area;
    auto rho = this->TrussMaterial->UnitWeight * A / G;
    auto rX2 = this->TrussSection->Inertia11 / A;
    double L = this->Length;
    double m = rho * L / 420.0;
    double mElm[12][12];

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            mElm[i][j] = 0;

    mElm[0][0] = m;
    mElm[6][6] = m;

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            (*this->LocalCoordinateMassMatrix)(i, j) = mElm[i][j];
}

void TrussMember::AssembleElementLocalStiffnessMatrix()
{
    auto E = this->TrussMaterial->E;
    auto A = this->TrussSection->Area;
    auto L = this->Length;
    double kElm[12][12];

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            kElm[i][j] = 0;

    kElm[0][0] = E * A / L;
    kElm[0][6] = -1 * E * A / L;
    kElm[6][0] = -1 * E * A / L;
    kElm[6][6] = E * A / L;

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            (*this->LocalCoordinateStiffnessMatrix)(i, j) = kElm[i][j];
}

void TrussMember::AssembleElementLocalDampingMatrix(double mult1, double mult2)
{
    this->LocalCoordinateDampingMatrix =
        std::make_shared<Matrix<double>>((*LocalCoordinateMassMatrix * mult1) + (*LocalCoordinateStiffnessMatrix * mult2));
}

void TrussMember::AssembleElementRotationMatrix()
{
    double pi = 3.141592653589793;

    Vector elmVector(this->Nodes[0]->Coordinate, this->Nodes[1]->Coordinate);
    auto minorRotMat = GeometryHelper::GetTranslationalRotationMatrix(elmVector, 0);

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            (*this->RotationMatrix)(i, j) = 0.0;

    for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i, j);

    for (unsigned int i = 3; i < 6; i++)
        for (unsigned int j = 3; j < 6; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 3, j - 3);

    for (unsigned int i = 6; i < 9; i++)
        for (unsigned int j = 6; j < 9; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 6, j - 6);

    for (unsigned int i = 9; i < 12; i++)
        for (unsigned int j = 9; j < 12; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 9, j - 9);
}

void TrussMember::AssembleElementGlobalMassMatrix()
{
    auto rotTrans = (*this->RotationMatrix).transpose();
    *this->GlobalCoordinateMassMatrix =
        (rotTrans * (*this->LocalCoordinateMassMatrix)) * (*this->RotationMatrix);
}

void TrussMember::AssembleElementGlobalStiffnessMatrix()
{
    auto rotTrans = (*this->RotationMatrix).transpose();
    *this->GlobalCoordinateStiffnessMatrix =
        (rotTrans * (*this->LocalCoordinateStiffnessMatrix)) * (*this->RotationMatrix);
}

void TrussMember::AssembleElementGlobalDampingMatrix(double mult1, double mult2)
{
    this->GlobalCoordinateDampingMatrix =
        std::make_shared<Matrix<double>>((*GlobalCoordinateMassMatrix * mult1) + (*GlobalCoordinateStiffnessMatrix * mult2));
}

unsigned int TrussMember::GetElementIndex() { return ElementIndex; };
unsigned int TrussMember::GetNumberOfDoF() { return 12; };
ElmType::ElementType TrussMember::GetElementType() { return Type; };
std::shared_ptr<Matrix<double>> TrussMember::GetLocalCoordinateMassMatrix() { return LocalCoordinateMassMatrix; };
std::shared_ptr<Matrix<double>> TrussMember::GetLocalCoordinateStiffnessMatrix() { return LocalCoordinateStiffnessMatrix; };
std::shared_ptr<Matrix<double>> TrussMember::GetLocalCoordinateDampingMatrix() { return LocalCoordinateDampingMatrix; };
std::shared_ptr<Matrix<double>> TrussMember::GetGlobalCoordinateMassMatrix() { return GlobalCoordinateMassMatrix; };
std::shared_ptr<Matrix<double>> TrussMember::GetGlobalCoordinateStiffnessMatrix() { return GlobalCoordinateStiffnessMatrix; };
std::shared_ptr<Matrix<double>> TrussMember::GetGlobalCoordinateDampingMatrix() { return GlobalCoordinateDampingMatrix; };
std::shared_ptr<Matrix<double>> TrussMember::GetRotationMatrix() { return RotationMatrix; };
std::shared_ptr<Matrix<double>> TrussMember::GetElementLoads() { return ElementLoads; };
std::vector<std::shared_ptr<Node>> TrussMember::GelElementNodes()
{
    std::vector<std::shared_ptr<Node>> retVal = { Nodes[0], Nodes[1] };
    return retVal;
};