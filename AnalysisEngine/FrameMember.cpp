#include "FrameMember.h"
#include "XYZPoint.h"
#include "Vector.h"
#include "MatrixHelper.h"
#include <math.h>

constexpr double G = 9.807;

FrameMember::FrameMember(unsigned int elmIndex, std::shared_ptr<Node> iNode, std::shared_ptr<Node> jNode,
    std::shared_ptr<Section> section, std::shared_ptr<Material> material, std::shared_ptr<Hinge> iEndHinge, std::shared_ptr<Hinge> jEndHinge)
{
    this->Nodes.resize(2);
    this->Hinges.resize(2);

    this->ElementIndex = elmIndex;
    this->Nodes[0] = iNode;
    this->Nodes[1] = jNode;
    this->Hinges[0] = iEndHinge;
    this->Hinges[1] = jEndHinge;
    this->FrameSection = section;
    this->FrameMaterial = material;
    this->Type = (ElmType::ElementType::Frame);
    this->Length = jNode->Coordinate.DistanceTo(iNode->Coordinate);
    AssembleElementLocalStiffnessMatrix();
    AssembleElementLocalMassMatrix();
    AssembleElementRotationMatrix();
    AssembleElementGlobalStiffnessMatrix();
    AssembleElementGlobalMassMatrix();
    iNode->ConnectedElements.push_back(elmIndex);
    jNode->ConnectedElements.push_back(elmIndex);
}

FrameMember::FrameMember()
{
    std::shared_ptr<Node> iN;
    std::shared_ptr<Node> jN;
    std::shared_ptr<Section> s;
    std::shared_ptr<Material> m;
    this->Nodes[0] = iN;
    this->Nodes[1] = jN;
    this->FrameSection = s;
    this->FrameMaterial = m;
    this->Length = 0;
}

FrameMember::~FrameMember()
{
}

void FrameMember::AssembleElementLocalStiffnessMatrix()
{
    auto A = this->FrameSection->Area;
    auto E = this->FrameMaterial->E;
    auto G = this->FrameMaterial->G;
    auto I11 = this->FrameSection->Inertia11;
    auto I22 = this->FrameSection->Inertia22;
    auto J = this->FrameSection->Inertia12;
    auto L = this->Length;
    auto L2 = L * L;
    auto L3 = L * L * L;
    double kElm[12][12];

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            kElm[i][j] = 0;

    kElm[0][0] = E * A / L;
    kElm[0][6] = -1 * E * A / L;

    kElm[1][1] = 12 * E * I11 / L3;
    kElm[1][5] = 6 * E * I11 / L2;
    kElm[1][7] = -12 * E * I11 / L3;
    kElm[1][11] = 6 * E * I11 / L2;

    kElm[2][2] = 12 * E * I22 / L3;
    kElm[2][4] = -6 * E * I22 / L2;
    kElm[2][8] = -12 * E * I22 / L3;
    kElm[2][10] = -6 * E * I22 / L2;

    kElm[3][3] = G * J / L;
    kElm[3][9] = -1 * G * J / L;

    kElm[4][2] = kElm[2][4];
    kElm[4][4] = 4 * E * I11 / L;
    kElm[4][8] = 6 * E * I11 / L2;
    kElm[4][10] = 2 * E * I11 / L;

    kElm[5][1] = kElm[1][5];
    kElm[5][5] = 4 * E * I22 / L;
    kElm[5][7] = -6 * E * I22 / L2;
    kElm[5][11] = 2 * E * I22 / L;

    kElm[6][0] = kElm[0][6];
    kElm[6][6] = E * A / L;

    kElm[7][1] = kElm[1][7];
    kElm[7][5] = kElm[5][7];
    kElm[7][7] = 12 * E * I11 / L3;
    kElm[7][11] = -6 * E * I11 / L2;

    kElm[8][2] = kElm[2][8];
    kElm[8][4] = kElm[4][8];
    kElm[8][8] = 12 * E * I22 / L3;
    kElm[8][10] = 6 * E * I22 / L2;

    kElm[9][3] = kElm[3][9];
    kElm[9][9] = G * J / L;

    kElm[10][2] = kElm[2][10];
    kElm[10][4] = kElm[4][10];
    kElm[10][8] = kElm[8][10];
    kElm[10][10] = 4 * E * I11 / L;

    kElm[11][1] = kElm[1][11];
    kElm[11][5] = kElm[5][11];
    kElm[11][7] = kElm[7][11];
    kElm[11][11] = 4 * E * I22 / L;

    // Check end releases
    for (int i = 0; i < 2; i++)
    {
        auto end = this->Hinges[i];

        if (!end)
            continue;

        if (end->IsReleasedForceX)
        {

            for (int j = 0; j < 12; j++) kElm[(6 * i) + 0][j] = 0.0;
            for (int j = 0; j < 12; j++) kElm[j][(6 * i) + 0] = 0.0;
        }

        if (end->IsReleasedForceY)
        {
            for (int j = 0; j < 12; j++) kElm[(6 * i) + 1][j] = 0.0;
            for (int j = 0; j < 12; j++) kElm[j][(6 * i) + 1] = 0.0;
        }

        if (end->IsReleasedForceZ)
        {
            for (int j = 0; j < 12; j++) kElm[(6 * i) + 2][j] = 0.0;
            for (int j = 0; j < 12; j++) kElm[j][(6 * i) + 2] = 0.0;
        }

        if (end->IsReleasedMomentX)
        {
            for (int j = 0; j < 12; j++) kElm[(6 * i) + 3][j] = 0.0;
            for (int j = 0; j < 12; j++) kElm[j][(6 * i) + 3] = 0.0;
        }

        if (end->IsReleasedMomentY)
        {
            for (int j = 0; j < 12; j++) kElm[(6 * i) + 4][j] = 0.0;
            for (int j = 0; j < 12; j++) kElm[j][(6 * i) + 4] = 0.0;
        }

        if (end->IsReleasedMomentY)
        {
            for (int j = 0; j < 12; j++) kElm[(6 * i) + 5][j] = 0.0;
            for (int j = 0; j < 12; j++) kElm[j][(6 * i) + 5] = 0.0;
        }
    }

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            this->LocalCoordinateStiffnessMatrix[i][j] = kElm[i][j];
}

/// <summary>
/// Taken from
/// http://what-when-how.com/the-finite-element-method/fem-for-frames-finite-element-method-part-1/
/// </summary>
void FrameMember::AssembleElementLocalMassMatrix()
{
    auto A = this->FrameSection->Area;
    auto rho = this->FrameMaterial->UnitWeight * A / G;
    auto rX2 = this->FrameSection->Inertia11 / A;
    double L = this->Length;
    double m = rho * L / 420.0;
    double mElm[12][12];

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            mElm[i][j] = 0;

    m = 0.5 * rho * L;
    mElm[0][0] = m;
    mElm[1][1] = m;
    mElm[2][2] = m;
    mElm[3][3] = 0 * m / 2;
    mElm[4][4] = 0 * m / 2;
    mElm[5][5] = 0 * m / 2;
    mElm[6][6] = m;
    mElm[7][7] = m;
    mElm[8][8] = m;
    mElm[9][9] = 0 * m / 2;
    mElm[10][10] = 0 * m / 2;
    mElm[11][11] = 0 * m / 2;

    /*mElm[0][0] = mElm[6][6] = m * 140.0;
    mElm[0][6] = mElm[6][0] = m * 70.0;
    mElm[3][3] = mElm[9][9] = m * rX2 * 140.0;
    mElm[3][9] = mElm[9][3] = m * rX2 * 70.0;

    mElm[2][2] = mElm[8][8] = m * 156.0;
    mElm[2][8] = mElm[8][2] = m * 54.0;
    mElm[4][4] = mElm[10][10] = m * 4.0 * L * L;
    mElm[4][10] = mElm[10][4] = -m * 3.0 * L * L;
    mElm[2][4] = mElm[4][2] = -m * 22.0 * L;
    mElm[8][10] = mElm[10][8] = -mElm[2][4];
    mElm[2][10] = mElm[10][2] = m * 13.0 * L;
    mElm[4][8] = mElm[8][4] = -mElm[2][10];

    mElm[1][1] = mElm[7][7] = m * 156.0;
    mElm[1][7] = mElm[7][1] = m * 54.0;
    mElm[5][5] = mElm[11][11] = m * 4.0 * L * L;
    mElm[5][11] = mElm[11][5] = -m * 3.0 * L * L;
    mElm[1][5] = mElm[5][1] = m * 22.0 * L;
    mElm[7][11] = mElm[11][7] = -mElm[1][5];
    mElm[1][11] = mElm[11][1] = -m * 13.0 * L;
    mElm[5][7] = mElm[7][5] = -mElm[1][11];*/

    /*auto a = this->Length / 2;
    auto a2 = a * a;
    auto mult = rho * a / 105;

    mElm[0][0] = mult * 70;
    mElm[0][6] = mult * 25;

    mElm[1][1] = mult * 78;
    mElm[1][5] = mult * 22 * a;
    mElm[1][7] = mult * 27;
    mElm[1][11] = mult * -13 * a;

    mElm[2][2] = mult * 78;
    mElm[2][4] = mult * -22 * a;
    mElm[2][8] = mult * 27;
    mElm[2][10] = mult * 13 * a;

    mElm[3][3] = mult * 70 * rX2;
    mElm[3][9] = mult * -35 * rX2;

    mElm[4][2] = mElm[2][4];
    mElm[4][4] = mult * 8 * a2;
    mElm[4][8] = mult * -13 * a;
    mElm[4][10] = mult * -6 * a2;

    mElm[5][1] = mElm[1][5];
    mElm[5][5] = mult * 8 * a2;
    mElm[5][7] = mult * 12 * a;
    mElm[5][11] = mult * -6 * a2;

    mElm[6][0] = mElm[0][6];
    mElm[6][6] = mult * 70;

    mElm[7][1] = mElm[1][7];
    mElm[7][5] = mElm[5][7];
    mElm[7][7] = mult * 78;
    mElm[7][11] = mult * -22 * a;

    mElm[8][2] = mElm[2][8];
    mElm[8][4] = mElm[4][8];
    mElm[8][8] = mult * 78;
    mElm[8][10] = mult * 22 * a;

    mElm[9][3] = mElm[3][9];
    mElm[9][9] = mult * 70 * rX2;

    mElm[10][2] = mElm[2][10];
    mElm[10][4] = mElm[4][10];
    mElm[10][8] = mElm[8][10];
    mElm[10][10] = mult * 8 * a2;

    mElm[11][1] = mElm[1][11];
    mElm[11][5] = mElm[5][11];
    mElm[11][7] = mElm[7][11];
    mElm[11][11] = mult * 8 * a2;*/

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            this->LocalCoordinateMassMatrix[i][j] = mElm[i][j];
}

void FrameMember::AssembleElementRotationMatrix()
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

void FrameMember::AssembleElementGlobalStiffnessMatrix()
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

void FrameMember::AssembleElementGlobalMassMatrix()
{
    std::vector<std::vector<double>> mElm;
    std::vector<std::vector<double>> rot;

    for (unsigned int i = 0; i < 12; i++)
    {
        std::vector<double> mElmRow;
        std::vector<double> rotRow;

        for (unsigned int j = 0; j < 12; j++)
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

    for (unsigned int i = 0; i < 12; i++)
        for (unsigned int j = 0; j < 12; j++)
            this->GlobalCoordinateMassMatrix[i][j] = mElmGlob.at(i).at(j);
}
