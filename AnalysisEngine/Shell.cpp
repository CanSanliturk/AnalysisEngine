#include "Shell.h"
#include "Vector.h"
#include "GeometryHelper.h"

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

    this->isLumpedMassMatrix = isLumpedMassMatrix;
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
    return ElementIndex;
}

unsigned int ShellMember::GetNumberOfDoF()
{
    return 24;
}

ElmType::ElementType ShellMember::GetElementType()
{
    return ElmType::ElementType::Shell;
}

std::shared_ptr<Matrix<double>> ShellMember::GetLocalCoordinateMassMatrix()
{
    return this->LocalCoordinateMassMatrix;
}

std::shared_ptr<Matrix<double>> ShellMember::GetLocalCoordinateStiffnessMatrix()
{
    return this->LocalCoordinateStiffnessMatrix;
}

std::shared_ptr<Matrix<double>> ShellMember::GetLocalCoordinateDampingMatrix()
{
    return this->LocalCoordinateDampingMatrix;
}

std::shared_ptr<Matrix<double>> ShellMember::GetGlobalCoordinateMassMatrix()
{
    return this->GlobalCoordinateMassMatrix;
}

std::shared_ptr<Matrix<double>> ShellMember::GetGlobalCoordinateStiffnessMatrix()
{
    return this->GlobalCoordinateStiffnessMatrix;
}

std::shared_ptr<Matrix<double>> ShellMember::GetGlobalCoordinateDampingMatrix()
{
    return this->GlobalCoordinateDampingMatrix;
}

std::shared_ptr<Matrix<double>> ShellMember::GetRotationMatrix()
{
    return this->RotationMatrix;
}


std::shared_ptr<Matrix<double>> ShellMember::GetElementLoads()
{
    return std::shared_ptr<Matrix<double>>();
}

std::vector<std::shared_ptr<Node>> ShellMember::GelElementNodes()
{
    return this->Nodes;
}

void ShellMember::AssembleElementLocalMassMatrix()
{
    auto asd = std::make_shared<Matrix<double>>(24, 24);
    this->LocalCoordinateMassMatrix = asd;
}

void ShellMember::AssembleElementLocalStiffnessMatrix()
{
    // Stiffness matrix is calculated at four Gauss points using Gauss Quadrature
    // It is assumed that nodes are oriented in counter-clock wise direction

    // Material parameters
    auto e = this->ShellMaterial->E;
    auto v = this->ShellMaterial->PoissonsRatio;
    auto eMult = e / (1 - (v * v));
    Matrix<double> eMat(3, 3);
    eMat(0, 0) = eMult * 1; eMat(0, 1) = eMult * v;
    eMat(1, 0) = eMult * v; eMat(1, 1) = eMult * 2;
    eMat(2, 2) = eMult * (1 - v) / 2;

    // Thickness
    auto thickness = this->Thickness;

    // Map coordinates of flat plane to 2-D surface
    auto d1 = Nodes[0]->Coordinate.DistanceTo(Nodes[1]->Coordinate);
    auto d2 = Nodes[1]->Coordinate.DistanceTo(Nodes[2]->Coordinate);
    auto d3 = Nodes[2]->Coordinate.DistanceTo(Nodes[3]->Coordinate);

    Vector p1V(Nodes[0]->Coordinate);
    Vector p2V(Nodes[1]->Coordinate);
    Vector p3V(Nodes[2]->Coordinate);
    Vector p4V(Nodes[3]->Coordinate);

    // Angle between first line and second line
    auto firstVector1 = p1V - p2V;
    auto secondVector1 = p3V - p2V;
    auto alpha1 = firstVector1.AngleTo(secondVector1);

    // Angle between second line and third line
    auto firstVector2 = p2V - p3V;
    auto secondVector2 = p4V - p3V;
    auto alpha2 = firstVector2.AngleTo(secondVector2);

    // Map 3D coordinates to 2D plane using angles and length found above to be able to
    // use natural coordinates
    auto x1 = 0.0; auto y1 = 0.0;
    auto x2 = d1; auto y2 = 0.0;
    auto x3 = x2 - (d2 * cos(alpha1)); auto y3 = d2 * sin(alpha1);
    auto x4 = x3 - (d3 * sin(alpha2)); auto y4 = y3 + (d3 * cos(alpha2));

    Matrix<double> mappedCoords(4, 2);
    mappedCoords(0, 0) = x1; mappedCoords(0, 1) = y1;
    mappedCoords(1, 0) = x2; mappedCoords(1, 1) = y2;
    mappedCoords(2, 0) = x3; mappedCoords(2, 1) = y3;
    mappedCoords(3, 0) = x4; mappedCoords(3, 1) = y4;

    // Membrane action resists in-plane translational degrees of freedom and the plate 
    // action resists bending effects. 
    Matrix<double> elmK(24, 24);

    if (isMembraneAction)
    {
        Matrix<double> kMembrane(8, 8);

        auto gp = 1 / sqrt(3);
        double gaussPoints[4][2] = { {-gp, -gp}, {gp, -gp}, {gp, gp}, {-gp, gp} };

        auto rowCounter = 0;
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                // Get Gauss point
                auto ksi = gaussPoints[rowCounter][0]; auto eta = gaussPoints[rowCounter][1];

                // Calculate jacobi
                Matrix<double> j1(2, 4);
                j1(0, 0) = eta - 1; j1(0, 1) = 1 - eta; j1(0, 2) = eta + 1; j1(0, 3) = -eta - 1;
                j1(1, 0) = ksi - 1; j1(1, 1) = -ksi - 1; j1(1, 2) = ksi + 1; j1(1, 3) = 1 - ksi;
                auto j2 = j1 * mappedCoords;
                auto jacobi = j2 * 0.25;

                Matrix<double> inversejacobi(2, 2);
                auto detjacobi = (jacobi(0, 0) * jacobi(1, 1)) - (jacobi(0, 1) * jacobi(1, 0));
                inversejacobi(0, 0) = jacobi(1, 1) / detjacobi; inversejacobi(0, 1) = -1 * jacobi(0, 1) / detjacobi;
                inversejacobi(1, 0) = -1 * jacobi(1, 0) / detjacobi; inversejacobi(1, 1) = jacobi(0, 0) / detjacobi;

                // Calculate strain-displacement matrix (B)
                Matrix<double> mat1(3, 4);
                mat1(0, 0) = 1; mat1(1, 3) = 1; mat1(2, 1) = 1; mat1(2, 2) = 1;

                Matrix<double> mat2(4, 4);
                mat2(0, 0) = inversejacobi(0, 0); mat2(0, 1) = inversejacobi(0, 1); mat2(1, 0) = inversejacobi(1, 0); mat2(1, 1) = inversejacobi(1, 1);
                mat2(2, 2) = inversejacobi(0, 0); mat2(2, 3) = inversejacobi(0, 1); mat2(3, 2) = inversejacobi(1, 0); mat2(3, 3) = inversejacobi(1, 1);

                Matrix<double> mat3(4, 8);
                mat3(0, 0) = eta - 1; mat3(0, 2) = 1 - eta; mat3(0, 4) = eta + 1; mat3(0, 6) = -eta - 1;
                mat3(1, 0) = ksi - 1; mat3(1, 2) = -ksi - 1; mat3(1, 4) = ksi + 1; mat3(1, 6) = 1 - ksi;
                mat3(2, 1) = eta - 1; mat3(2, 3) = 1 - eta; mat3(2, 5) = eta + 1; mat3(2, 7) = -eta - 1;
                mat3(3, 1) = ksi - 1; mat3(3, 3) = -ksi - 1; mat3(3, 5) = ksi + 1; mat3(3, 7) = 1 - ksi;
                mat3 *= 0.25;

                auto b = mat1 * mat2 * mat3;
                auto bT = b.transpose();

                auto kPt = bT * eMat * b * thickness * detjacobi;

                for (size_t row = 0; row < 8; row++)
                    for (size_t col = 0; col < 8; col++)
                        kMembrane(row, col) += kPt(row, col);

                rowCounter++;
            }
        }

        elmK(0, 0) = kMembrane(0, 0); elmK(0, 1) = kMembrane(0, 1); elmK(0, 6) = kMembrane(0, 2); elmK(0, 7) = kMembrane(0, 3); elmK(0, 11) = kMembrane(0, 4); elmK(0, 12) = kMembrane(0, 5); elmK(0, 17) = kMembrane(0, 6); elmK(0, 18) = kMembrane(0, 7);
        elmK(1, 0) = kMembrane(1, 0); elmK(1, 1) = kMembrane(1, 1); elmK(1, 6) = kMembrane(1, 2); elmK(1, 7) = kMembrane(1, 3); elmK(1, 11) = kMembrane(1, 4); elmK(1, 12) = kMembrane(1, 5); elmK(1, 17) = kMembrane(1, 6); elmK(1, 18) = kMembrane(1, 7);
        elmK(6, 0) = kMembrane(2, 0); elmK(6, 1) = kMembrane(2, 1); elmK(6, 6) = kMembrane(2, 2); elmK(6, 7) = kMembrane(2, 3); elmK(6, 11) = kMembrane(2, 4); elmK(6, 12) = kMembrane(2, 5); elmK(6, 17) = kMembrane(2, 6); elmK(6, 18) = kMembrane(2, 7);
        elmK(7, 0) = kMembrane(3, 0); elmK(7, 1) = kMembrane(3, 1); elmK(7, 6) = kMembrane(3, 2); elmK(7, 7) = kMembrane(3, 3); elmK(7, 11) = kMembrane(3, 4); elmK(7, 12) = kMembrane(3, 5); elmK(7, 17) = kMembrane(3, 6); elmK(7, 18) = kMembrane(3, 7);

        elmK(12, 0) = kMembrane(4, 0); elmK(12, 1) = kMembrane(4, 1); elmK(12, 6) = kMembrane(4, 2); elmK(12, 7) = kMembrane(4, 3); elmK(12, 11) = kMembrane(4, 4); elmK(12, 12) = kMembrane(4, 5); elmK(12, 17) = kMembrane(4, 6); elmK(12, 18) = kMembrane(4, 7);
        elmK(13, 0) = kMembrane(5, 0); elmK(13, 1) = kMembrane(5, 1); elmK(13, 6) = kMembrane(5, 2); elmK(13, 7) = kMembrane(5, 3); elmK(13, 11) = kMembrane(5, 4); elmK(13, 12) = kMembrane(5, 5); elmK(13, 17) = kMembrane(5, 6); elmK(13, 18) = kMembrane(5, 7);
        elmK(18, 0) = kMembrane(6, 0); elmK(18, 1) = kMembrane(6, 1); elmK(18, 6) = kMembrane(6, 2); elmK(18, 7) = kMembrane(6, 3); elmK(18, 11) = kMembrane(6, 4); elmK(18, 12) = kMembrane(6, 5); elmK(18, 17) = kMembrane(6, 6); elmK(18, 18) = kMembrane(6, 7);
        elmK(19, 0) = kMembrane(7, 0); elmK(19, 1) = kMembrane(7, 1); elmK(19, 6) = kMembrane(7, 2); elmK(19, 7) = kMembrane(7, 3); elmK(19, 11) = kMembrane(7, 4); elmK(19, 12) = kMembrane(7, 5); elmK(19, 17) = kMembrane(7, 6); elmK(19, 18) = kMembrane(7, 7);
    }

    // Plate action is not implemented yet
    if (isPlateAction)
    {
        // Place holder
    }

    auto asd = std::make_shared<Matrix<double>>(elmK);
    this->LocalCoordinateStiffnessMatrix = asd;
}

void ShellMember::AssembleElementLocalDampingMatrix(double mult1, double mult2)
{
    auto mat1 = *LocalCoordinateMassMatrix * mult1;
    auto mat2 = *LocalCoordinateStiffnessMatrix * mult2;
    this->LocalCoordinateDampingMatrix =
        std::make_shared<Matrix<double>>(mat1 + mat2);
}

void ShellMember::AssembleElementRotationMatrix()
{
    // For 1-D elements, local and global matrices are same if their element vectors are in global positive X-Direction. For 2-D members,
    // local and global matrices are same if their normal vectors are in global positive X-Direction. Right now, element matrices are set
    // like it lies in XY-Plane. So, two rotations are required. First, rotation from XY-plane to YZ-plane to make like its stiffness matrix
    // is set for YZ-plane. Then, normal rotation.

    // Calculate element normal vector
    Vector p1V(Nodes[0]->Coordinate);
    Vector p2V(Nodes[1]->Coordinate);
    Vector p3V(Nodes[2]->Coordinate);
    auto planeVector1 = p2V - p1V;
    auto planeVector2 = p3V - p2V;
    auto planeNormal = planeVector1 * planeVector2;
    Vector corrector(0, 0, 1);
    
    Vector opt1(1, 0, 0);

    auto rotVec = (planeNormal == corrector) ? opt1 : corrector * planeNormal;

    auto minorRotMat = GeometryHelper::GetTranslationalRotationMatrix(rotVec, 0);

    /*for (unsigned int i = 0; i < 3; i++)
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

    for (unsigned int i = 12; i < 15; i++)
        for (unsigned int j = 12; j < 15; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 12, j - 12);

    for (unsigned int i = 15; i < 18; i++)
        for (unsigned int j = 15; j < 18; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 15, j - 15);

    for (unsigned int i = 18; i < 21; i++)
        for (unsigned int j = 18; j < 21; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 18, j - 18);

    for (unsigned int i = 21; i < 24; i++)
        for (unsigned int j = 21; j < 24; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 21, j - 21);*/
    for (size_t i = 0; i < 24; i++)
        for (size_t j = 0; j < 24; j++)
            if (i == j)
                (*this->RotationMatrix)(i, j) = 1;
}

void ShellMember::AssembleElementGlobalMassMatrix()
{
    auto rotTrans = (*this->RotationMatrix).transpose();
    *this->GlobalCoordinateMassMatrix =
        (rotTrans * (*this->LocalCoordinateMassMatrix)) * (*this->RotationMatrix);
}

void ShellMember::AssembleElementGlobalStiffnessMatrix()
{
    this->GlobalCoordinateStiffnessMatrix = this->LocalCoordinateStiffnessMatrix;
    //auto rotTrans = (*this->RotationMatrix).transpose();
    //*this->GlobalCoordinateStiffnessMatrix =
    //    (rotTrans * (*this->LocalCoordinateStiffnessMatrix)) * (*this->RotationMatrix);
}

void ShellMember::AssembleElementGlobalDampingMatrix(double mult1, double mult2)
{
    auto mat1 = *GlobalCoordinateMassMatrix * mult1;
    auto mat2 = *GlobalCoordinateStiffnessMatrix * mult2;
    this->GlobalCoordinateDampingMatrix =
        std::make_shared<Matrix<double>>(mat1 + mat2);
}
