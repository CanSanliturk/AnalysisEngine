#include "TrussMember.h"
#include "XYZPoint.h"
#include "Vector.h"
#include "MatrixHelper.h"
#include <math.h>

TrussMember::TrussMember(unsigned int ElmIndex, Node* iNode, Node* jNode, Section* section, Material* material)
    : ElementIndex(ElmIndex), TrussSection(section), TrussMaterial(material), Type(ElmType::ElementType::Truss),
    Length(jNode->Coordinate.DistanceTo(iNode->Coordinate))
{
    this->Nodes[0] = iNode;
    this->Nodes[1] = jNode;
    AssembleElementLocalStiffnessMatrix();
    AssembleElementLocalMassMatrix();
    AssembleElementRotationMatrix();
    AssembleElementGlobalStiffnessMatrix();
    AssembleElementGlobalMassMatrix();
    iNode->ConnectedElements.push_back(ElmIndex);
}

TrussMember::TrussMember()
{
    Node iN;
    Node jN;
    Section s;
    Material m;
    this->ElementIndex = 0;
    this->Nodes[0] = &iN;
    this->Nodes[1] = &jN;
    this->TrussSection = &s;
    this->TrussMaterial = &m;
    this->Length = -123.456789;
}

TrussMember::~TrussMember()
{
}

void TrussMember::AssembleElementLocalStiffnessMatrix()
{
    auto A = this->TrussSection->Area;
    auto E = this->TrussMaterial->E;
    auto L = this->Length;
    double kElm[12][12];

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            kElm[i][j] = 0;

    kElm[0][0] = E * A / L;
    kElm[0][6] = -1 * E * A / L;
    kElm[6][0] = -1 * E * A / L;
    kElm[6][6] = E * A / L;

    for (unsigned int i = 0; i < 6; i++)
        for (unsigned int j = 0; j < 6; j++)
            this->LocalCoordinateStiffnessMatrix[i][j] = kElm[i][j];
}

void TrussMember::AssembleElementLocalMassMatrix()
{
    auto A = this->TrussSection->Area;
    auto rho = this->TrussMaterial->UnitWeight;
    auto rX2 = this->TrussSection->Inertia11 / A;
    auto a = this->Length / 2.0;
    auto a2 = a * a;
    auto mult = rho * A * a / 210.0;
    double mElm[12][12];

    for (unsigned int i = 0; i < 6; i++)
        for (unsigned int j = 0; j < 6; j++)
            mElm[i][j] = 0;

    for (unsigned int i = 0; i < 6; i++)
        for (unsigned int j = 0; j < 6; j++)
            this->LocalCoordinateMassMatrix[i][j] = mElm[i][j];
}

void TrussMember::AssembleElementRotationMatrix()
{
    double pi = 3.141592653589793;

    Vector elmVector(this->Nodes[0]->Coordinate, this->Nodes[1]->Coordinate);
    auto minorRotMat = MatrixHelper::GetTranslationalRotationMatrix(elmVector, 0);

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            this->RotationMatrix[i][j] = 0.0;

    for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
            this->RotationMatrix[i][j] = minorRotMat.at(i).at(j);

    for (unsigned int i = 3; i < 6; i++)
        for (unsigned int j = 3; j < 6; j++)
            this->RotationMatrix[i][j] = minorRotMat.at(i - 3).at(j - 3);

    for (unsigned int i = 6; i < 9; i++)
        for (unsigned int j = 6; j < 9; j++)
            this->RotationMatrix[i][j] = minorRotMat.at(i - 6).at(j - 6);

    for (unsigned int i = 9; i < 12; i++)
        for (unsigned int j = 9; j < 12; j++)
            this->RotationMatrix[i][j] = minorRotMat.at(i - 9).at(j - 9);
}

void TrussMember::AssembleElementGlobalStiffnessMatrix()
{
    std::vector<std::vector<double>> kElm;
    std::vector<std::vector<double>> rot;

    for (unsigned int i = 0; i < 12; i++)
    {
        std::vector<double> kElmRow;
        std::vector<double> rotRow;

        for (unsigned int j = 0; j < 12; j++)
        {
            kElmRow.push_back(this->LocalCoordinateStiffnessMatrix[i][j]);
            rotRow.push_back(this->RotationMatrix[i][j]);
        }

        kElm.push_back(kElmRow);
        rot.push_back(rotRow);
    }

    auto rotTrans = MatrixHelper::GetTranspose(rot);

    auto firstStep = MatrixHelper::MultiplyMatrices(rotTrans, kElm);
    auto kElmGlob = MatrixHelper::MultiplyMatrices(firstStep, rot);

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            this->GlobalCoordinateStiffnessMatrix[i][j] = kElmGlob.at(i).at(j);
}

void TrussMember::AssembleElementGlobalMassMatrix()
{
    std::vector<std::vector<double>> mElm;
    std::vector<std::vector<double>> rot;

    for (unsigned int i = 0; i < 6; i++)
    {
        std::vector<double> mElmRow;
        std::vector<double> rotRow;

        for (unsigned int j = 0; j < 6; j++)
        {
            mElmRow.push_back(this->LocalCoordinateMassMatrix[i][j]);
            rotRow.push_back(this->RotationMatrix[i][j]);
        }

        mElm.push_back(mElmRow);
        rot.push_back(rotRow);
    }

    auto rotTrans = MatrixHelper::GetTranspose(rot);

    auto firstStep = MatrixHelper::MultiplyMatrices(rotTrans, mElm);
    auto mElmGlob = MatrixHelper::MultiplyMatrices(firstStep, rot);

    for (unsigned int i = 0; i < 6; i++)
        for (unsigned int j = 0; j < 6; j++)
            this->GlobalCoordinateMassMatrix[i][j] = mElmGlob.at(i).at(j);
}
