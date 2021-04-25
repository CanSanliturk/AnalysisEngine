#include "Shell.h"
#include "Vector.h"
#include "GeometryHelper.h"


ShellMember::ShellMember(unsigned int elmIndex,
    std::shared_ptr<Node> iNode, std::shared_ptr<Node> jNode, std::shared_ptr<Node> kNode, std::shared_ptr<Node> lNode,
    std::shared_ptr<Material> material, double thickness, MembraneType memType, PlateType pltType,
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

    Type = ElmType::ElementType::Shell;

    this->isLumpedMassMatrix = isLumpedMassMatrix;
    membraneType = memType;
    plateType = pltType;

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
    eMat(1, 0) = eMult * v; eMat(1, 1) = eMult * 1;
    eMat(2, 2) = eMult * (1 - v) / 2;

    // Thickness
    auto thickness = this->Thickness;

    // Map coordinates of flat plane to 2-D surface
    auto d1 = Nodes[0]->Coordinate.DistanceTo(Nodes[1]->Coordinate);
    auto d2 = Nodes[1]->Coordinate.DistanceTo(Nodes[2]->Coordinate);
    auto d3 = Nodes[2]->Coordinate.DistanceTo(Nodes[3]->Coordinate);
    auto d4 = Nodes[3]->Coordinate.DistanceTo(Nodes[0]->Coordinate);

    Vector p1V(Nodes[0]->Coordinate);
    Vector p2V(Nodes[1]->Coordinate);
    Vector p3V(Nodes[2]->Coordinate);
    Vector p4V(Nodes[3]->Coordinate);

    // Angle between first line and fourth line
    auto firstVector0 = p2V - p1V;
    auto secondVector0 = p4V - p1V;
    auto alpha0 = firstVector0.AngleTo(secondVector0);

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
    auto x4 = d4 * cos(alpha0); auto y4 = d4 * sin(alpha0);

    Matrix<double> mappedCoords(4, 2);
    mappedCoords(0, 0) = x1; mappedCoords(0, 1) = y1;
    mappedCoords(1, 0) = x2; mappedCoords(1, 1) = y2;
    mappedCoords(2, 0) = x3; mappedCoords(2, 1) = y3;
    mappedCoords(3, 0) = x4; mappedCoords(3, 1) = y4;

    // Membrane action resists in-plane translational degrees of freedom and the plate 
    // action resists bending effects. 
    Matrix<double> elmK(24, 24);

    Matrix<double> kMembrane(8, 8);
    if (membraneType == MembraneType::Bilinear)
    {

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

                // Stiffness calculated at given Gauss point
                auto kPt = bT * eMat * b * thickness * detjacobi;

                // Update stifness
                for (size_t row = 0; row < 8; row++)
                    for (size_t col = 0; col < 8; col++)
                        kMembrane(row, col) += kPt(row, col);

                rowCounter++;
            }
        }
    }
    else if (membraneType == MembraneType::Incompatible)
    {
        // Populate 12 x 12 stiffness matrix including effects of generalized displacements. Then, condense out 
        // last 4 x 4 part of the stiffness matrix and equate it to kMembrane.

        auto gp = 1 / sqrt(3);
        double gaussPoints[4][2] = { {-gp, -gp}, {gp, -gp}, {gp, gp}, {-gp, gp} };
        auto rowCounter = 0;
        Matrix<double> k1212(12, 12);

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

                // Calculate Bu
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

                auto bu = mat1 * mat2 * mat3;

                // Insert effect of generalized displacement to strain-displacement matrix along with strain-displacement matrix of nodal displacements
                Matrix<double> b(3, 12);

                for (size_t bRowCounter = 0; bRowCounter < 3; bRowCounter++)
                    for (size_t bColCounter = 0; bColCounter < 8; bColCounter++)
                        b(bRowCounter, bColCounter) = bu(bRowCounter, bColCounter);

                auto invJ11 = inversejacobi(0, 0); auto invJ12 = inversejacobi(0, 1);
                auto invJ21 = inversejacobi(1, 0); auto invJ22 = inversejacobi(1, 1);

                b(0, 8) = invJ11 * -2 * ksi; b(0, 9) = invJ12 * -2 * eta;
                b(1, 10) = invJ21 * -2 * ksi; b(1, 11) = invJ22 * -2 * eta;
                b(2, 8) = invJ21 * -2 * ksi; b(2, 9) = invJ22 * -2 * eta; b(2, 10) = invJ11 * -2 * ksi; b(2, 11) = invJ12 * -2 * eta;

                auto bT = b.transpose();

                // Stiffness calculated at given Gauss point
                auto kPt = bT * eMat * b * thickness * detjacobi;

                // Update stifness
                for (size_t row = 0; row < 12; row++)
                    for (size_t col = 0; col < 12; col++)
                        k1212(row, col) += kPt(row, col);

                rowCounter++;
            }
        }

        this->K1212 = std::make_shared<Matrix<double>>(k1212);

        // After calculation 12 x 12 stiffness matrix, condense out last 4 x 4 part
        auto kRR = k1212.getSubmatrix(0, 7, 0, 7);
        auto kRC = k1212.getSubmatrix(0, 7, 8, 11);
        auto kCR = kRC.transpose();
        auto kCC = k1212.getSubmatrix(8, 11, 8, 11);
        auto invKCC = this->InvertMatrix4(kCC);

        auto midStep = (kRC * invKCC * kCR);
        kMembrane = kRR - midStep;
    }
    else if (membraneType == MembraneType::Drilling)
    {
        auto gp = 1 / sqrt(3);
        double gaussPoints[4][2] = { {-gp, -gp}, {gp, -gp}, {gp, gp}, {-gp, gp} };
        Matrix<double> k1212(12, 12);
        int MS[4][2] = { {7, 4}, {4, 5}, {5, 6}, {6, 7} };
        int ML[4] = { 3, 0, 1, 2 };
        int JK[4][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
        double vValue = 0.0;
        Matrix<double> initialKStorage(12, 12);
        // Find delta-Y's and minus delta-X's for edges
        Matrix<double> nL(4, 2);
        for (size_t edgeCounter = 0; edgeCounter < 4; edgeCounter++)
        {
            auto edgeStartNodeCoordX = this->Nodes[JK[edgeCounter][0]]->Coordinate.X;
            auto edgeStartNodeCoordY = this->Nodes[JK[edgeCounter][0]]->Coordinate.Y;
            auto edgeEndNodeCoordX = this->Nodes[JK[edgeCounter][1]]->Coordinate.X;
            auto edgeEndNodeCoordY = this->Nodes[JK[edgeCounter][1]]->Coordinate.Y;
            nL(edgeCounter, 0) = edgeEndNodeCoordY - edgeStartNodeCoordY;
            nL(edgeCounter, 1) = edgeStartNodeCoordX - edgeEndNodeCoordX;
        }

        for (int rowCounter = 0; rowCounter < 4; rowCounter++)
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

            // Calculate strain-displacement matrix
            // Initially, calculate it like bilinear strain-displacement matrix. Then, insert columns for rotational dofs.
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

            auto littleBMatrix = mat1 * mat2 * mat3;

            // Create a new B-Matrix and insert elements of little B matrix
            Matrix<double> bMatrix(3, 12);
            int mapper[12] = { 0, 1, -1, 2, 3, -1, 4, 5, -1, 6, 7, -1 };
            for (size_t colIdx = 0; colIdx < 12; colIdx++)
                if (mapper[colIdx] != -1)
                    for (size_t rowIdx = 0; rowIdx < 3; rowIdx++)
                        bMatrix(rowIdx, colIdx) = littleBMatrix(rowIdx, mapper[colIdx]);

            // Strain-displacement relation is established for translational dofs. Now, elements about strain-displacement relation for rotational dofs
            // should be inserted

            // Derivatives of quadratic shape functions in natural coordinate system
            Matrix<double> dN8Ksi(8, 1);
            dN8Ksi(4, 0) = 0.25 * 4 * ksi * (eta - 1);
            dN8Ksi(5, 0) = 0.25 * (2 - (2 * eta * eta));
            dN8Ksi(6, 0) = 0.25 * -4 * ksi * (eta + 1);
            dN8Ksi(7, 0) = 0.25 * (-2 + (2 * eta * eta));

            Matrix<double> dN8Eta(8, 1);
            dN8Eta(4, 0) = 0.25 * (-2 + (2 * ksi * ksi));
            dN8Eta(5, 0) = 0.25 * (-4 * (1 + ksi) * eta);
            dN8Eta(6, 0) = 0.25 * (2 - (2 * ksi * ksi));
            dN8Eta(7, 0) = 0.25 * 4 * (ksi - 1) * eta;

            // Derivatives of quadratic shape functions in cartesian coordinate system
            auto JI11 = inversejacobi(0, 0); auto JI12 = inversejacobi(0, 1);
            auto JI21 = inversejacobi(1, 0); auto JI22 = inversejacobi(1, 1);

            Matrix<double> dNS_dx(8, 1);
            dNS_dx(4, 0) = (JI11 * dN8Ksi(4, 0)) + (JI12 * dN8Eta(4, 0));
            dNS_dx(5, 0) = (JI11 * dN8Ksi(5, 0)) + (JI12 * dN8Eta(5, 0));
            dNS_dx(6, 0) = (JI11 * dN8Ksi(6, 0)) + (JI12 * dN8Eta(6, 0));
            dNS_dx(7, 0) = (JI11 * dN8Ksi(7, 0)) + (JI12 * dN8Eta(7, 0));

            Matrix<double> dNS_dy(8, 1);
            dNS_dx(4, 0) = (JI21 * dN8Ksi(4, 0)) + (JI22 * dN8Eta(4, 0));
            dNS_dx(5, 0) = (JI21 * dN8Ksi(5, 0)) + (JI22 * dN8Eta(5, 0));
            dNS_dx(6, 0) = (JI21 * dN8Ksi(6, 0)) + (JI22 * dN8Eta(6, 0));
            dNS_dx(7, 0) = (JI21 * dN8Ksi(7, 0)) + (JI22 * dN8Eta(7, 0));


            for (size_t i = 0; i < 4; i++)
            {
                auto m = (3 * i) + 2;

                auto MSi1 = MS[i][0];
                auto MSi2 = MS[i][1];
                auto MLi = ML[i];

                bMatrix(0, m) = ((dNS_dx(MSi1, 0) * nL(MLi, 0)) - (dNS_dx(MSi2, 0) * nL(i, 0))) / 8;
                bMatrix(1, m) = ((dNS_dy(MSi1, 0) * nL(MLi, 1)) - (dNS_dy(MSi2, 0) * nL(i, 1))) / 8;
                bMatrix(2, m) = (((dNS_dy(MSi1, 0) * nL(MLi, 0)) - (dNS_dy(MSi2, 0) * nL(i, 0))) +
                    ((dNS_dx(MSi1, 0) * nL(MLi, 1)) - (dNS_dx(MSi2, 0) * nL(i, 1)))) / 8;
            }

            double viValue = thickness * detjacobi;
            vValue += viValue;

            auto littleK = bMatrix.transpose() * eMat * bMatrix * viValue;
            initialKStorage += littleK;
        }

        // Stabilization matrix for Drilling DOF
        Matrix<double> kD(12, 12);
        auto kDia = (1 / 150000) * this->ShellMaterial->E * vValue;
        
        for (size_t i = 0; i < 4; i++)
            kD(3 * i, 3 * i) = 1.75;
        
        kD(2, 5) = -0.75; kD(2, 8) = -0.25; kD(2, 11) = -0.75;
        kD(5, 8) = -0.75; kD(5, 11) = -0.25;
        kD(8, 11) = -0.75;
        kD(5, 2) = kD(2, 5);
        kD(8, 2) = kD(2, 8);
        kD(11, 2) = kD(2, 11);
        kD(8, 5) = kD(5, 8);
        kD(11, 5) = kD(5, 11);
        kD(11, 8) = kD(8, 11);
        
        kD *= kDia;
        
        auto membK = initialKStorage + kD;

        elmK(0, 0) = membK(0, 0); elmK(0, 1) = membK(0, 1); elmK(0, 3) = membK(0, 2); elmK(0, 6) = membK(0, 3); elmK(0, 7) = membK(0, 4); elmK(0, 9) = membK(0, 5); elmK(0, 12) = membK(0, 6); elmK(0, 13) = membK(0, 7); elmK(0, 15) = membK(0, 8); elmK(0, 18) = membK(0, 9); elmK(0, 19) = membK(0, 10); elmK(0, 21) = membK(0, 11);
        elmK(1, 0) = membK(1, 0); elmK(1, 1) = membK(1, 1); elmK(1, 3) = membK(1, 2); elmK(1, 6) = membK(1, 3); elmK(1, 7) = membK(1, 4); elmK(1, 9) = membK(1, 5); elmK(1, 12) = membK(1, 6); elmK(1, 13) = membK(1, 7); elmK(1, 15) = membK(1, 8); elmK(1, 18) = membK(1, 9); elmK(1, 19) = membK(1, 10); elmK(1, 21) = membK(1, 11);
        elmK(3, 0) = membK(2, 0); elmK(3, 1) = membK(2, 1); elmK(3, 3) = membK(2, 2); elmK(3, 6) = membK(2, 3); elmK(3, 7) = membK(2, 4); elmK(3, 9) = membK(2, 5); elmK(3, 12) = membK(2, 6); elmK(3, 13) = membK(2, 7); elmK(3, 15) = membK(2, 8); elmK(3, 18) = membK(2, 9); elmK(3, 19) = membK(2, 10); elmK(3, 21) = membK(2, 11);
        elmK(6, 0) = membK(3, 0); elmK(6, 1) = membK(3, 1); elmK(6, 3) = membK(3, 2); elmK(6, 6) = membK(3, 3); elmK(6, 7) = membK(3, 4); elmK(6, 9) = membK(3, 5); elmK(6, 12) = membK(3, 6); elmK(6, 13) = membK(3, 7); elmK(6, 15) = membK(3, 8); elmK(6, 18) = membK(3, 9); elmK(6, 19) = membK(3, 10); elmK(6, 21) = membK(3, 11);
        elmK(7, 0) = membK(4, 0); elmK(7, 1) = membK(4, 1); elmK(7, 3) = membK(4, 2); elmK(7, 6) = membK(4, 3); elmK(7, 7) = membK(4, 4); elmK(7, 9) = membK(4, 5); elmK(7, 12) = membK(4, 6); elmK(7, 13) = membK(4, 7); elmK(7, 15) = membK(4, 8); elmK(7, 18) = membK(4, 9); elmK(7, 19) = membK(4, 10); elmK(7, 21) = membK(4, 11);
        elmK(9, 0) = membK(5, 0); elmK(9, 1) = membK(5, 1); elmK(9, 3) = membK(5, 2); elmK(9, 6) = membK(5, 3); elmK(9, 7) = membK(5, 4); elmK(9, 9) = membK(5, 5); elmK(9, 12) = membK(5, 6); elmK(9, 13) = membK(5, 7); elmK(9, 15) = membK(5, 8); elmK(9, 18) = membK(5, 9); elmK(9, 19) = membK(5, 10); elmK(9, 21) = membK(5, 11);
        elmK(12, 0) = membK(6, 0); elmK(12, 1) = membK(6, 1); elmK(12, 3) = membK(6, 2); elmK(12, 6) = membK(6, 3); elmK(12, 7) = membK(6, 4); elmK(12, 9) = membK(6, 5); elmK(12, 12) = membK(6, 6); elmK(12, 13) = membK(6, 7); elmK(12, 15) = membK(6, 8); elmK(12, 18) = membK(6, 9); elmK(12, 19) = membK(6, 10); elmK(12, 21) = membK(6, 11);
        elmK(13, 0) = membK(7, 0); elmK(13, 1) = membK(7, 1); elmK(13, 3) = membK(7, 2); elmK(13, 6) = membK(7, 3); elmK(13, 7) = membK(7, 4); elmK(13, 9) = membK(7, 5); elmK(13, 12) = membK(7, 6); elmK(13, 13) = membK(7, 7); elmK(13, 15) = membK(7, 8); elmK(13, 18) = membK(7, 9); elmK(13, 19) = membK(7, 10); elmK(13, 21) = membK(7, 11);
        elmK(15, 0) = membK(8, 0); elmK(15, 1) = membK(8, 1); elmK(15, 3) = membK(8, 2); elmK(15, 6) = membK(8, 3); elmK(15, 7) = membK(8, 4); elmK(15, 9) = membK(8, 5); elmK(15, 12) = membK(8, 6); elmK(15, 13) = membK(8, 7); elmK(15, 15) = membK(8, 8); elmK(15, 18) = membK(8, 9); elmK(15, 19) = membK(8, 10); elmK(15, 21) = membK(8, 11);
        elmK(18, 0) = membK(9, 0); elmK(18, 1) = membK(9, 1); elmK(18, 3) = membK(9, 2); elmK(18, 6) = membK(9, 3); elmK(18, 7) = membK(9, 4); elmK(18, 9) = membK(9, 5); elmK(18, 12) = membK(9, 6); elmK(18, 13) = membK(9, 7); elmK(18, 15) = membK(9, 8); elmK(18, 18) = membK(9, 9); elmK(18, 19) = membK(9, 10); elmK(18, 21) = membK(9, 11);
        elmK(19, 0) = membK(10, 0); elmK(19, 1) = membK(10, 1); elmK(19, 3) = membK(10, 2); elmK(19, 6) = membK(10, 3); elmK(19, 7) = membK(10, 4); elmK(19, 9) = membK(10, 5); elmK(19, 12) = membK(10, 6); elmK(19, 13) = membK(10, 7); elmK(19, 15) = membK(10, 8); elmK(19, 18) = membK(10, 9); elmK(19, 19) = membK(10, 10); elmK(19, 21) = membK(10, 11);
        elmK(21, 0) = membK(11, 0); elmK(21, 1) = membK(11, 1); elmK(21, 3) = membK(11, 2); elmK(21, 6) = membK(11, 3); elmK(21, 7) = membK(11, 4); elmK(21, 9) = membK(11, 5); elmK(21, 12) = membK(11, 6); elmK(21, 13) = membK(11, 7); elmK(21, 15) = membK(11, 8); elmK(21, 18) = membK(11, 9); elmK(21, 19) = membK(11, 10); elmK(21, 21) = membK(11, 11);
    }

    if ((this->membraneType != MembraneType::Drilling) && (this->membraneType != MembraneType::NONE))
    {
        elmK(0, 0) = kMembrane(0, 0); elmK(0, 1) = kMembrane(0, 1); elmK(0, 6) = kMembrane(0, 2); elmK(0, 7) = kMembrane(0, 3); elmK(0, 12) = kMembrane(0, 4); elmK(0, 13) = kMembrane(0, 5); elmK(0, 18) = kMembrane(0, 6); elmK(0, 19) = kMembrane(0, 7);
        elmK(1, 0) = kMembrane(1, 0); elmK(1, 1) = kMembrane(1, 1); elmK(1, 6) = kMembrane(1, 2); elmK(1, 7) = kMembrane(1, 3); elmK(1, 12) = kMembrane(1, 4); elmK(1, 13) = kMembrane(1, 5); elmK(1, 18) = kMembrane(1, 6); elmK(1, 19) = kMembrane(1, 7);
        elmK(6, 0) = kMembrane(2, 0); elmK(6, 1) = kMembrane(2, 1); elmK(6, 6) = kMembrane(2, 2); elmK(6, 7) = kMembrane(2, 3); elmK(6, 12) = kMembrane(2, 4); elmK(6, 13) = kMembrane(2, 5); elmK(6, 18) = kMembrane(2, 6); elmK(6, 19) = kMembrane(2, 7);
        elmK(7, 0) = kMembrane(3, 0); elmK(7, 1) = kMembrane(3, 1); elmK(7, 6) = kMembrane(3, 2); elmK(7, 7) = kMembrane(3, 3); elmK(7, 12) = kMembrane(3, 4); elmK(7, 13) = kMembrane(3, 5); elmK(7, 18) = kMembrane(3, 6); elmK(7, 19) = kMembrane(3, 7);

        elmK(12, 0) = kMembrane(4, 0); elmK(12, 1) = kMembrane(4, 1); elmK(12, 6) = kMembrane(4, 2); elmK(12, 7) = kMembrane(4, 3); elmK(12, 12) = kMembrane(4, 4); elmK(12, 13) = kMembrane(4, 5); elmK(12, 18) = kMembrane(4, 6); elmK(12, 19) = kMembrane(4, 7);
        elmK(13, 0) = kMembrane(5, 0); elmK(13, 1) = kMembrane(5, 1); elmK(13, 6) = kMembrane(5, 2); elmK(13, 7) = kMembrane(5, 3); elmK(13, 12) = kMembrane(5, 4); elmK(13, 13) = kMembrane(5, 5); elmK(13, 18) = kMembrane(5, 6); elmK(13, 19) = kMembrane(5, 7);
        elmK(18, 0) = kMembrane(6, 0); elmK(18, 1) = kMembrane(6, 1); elmK(18, 6) = kMembrane(6, 2); elmK(18, 7) = kMembrane(6, 3); elmK(18, 12) = kMembrane(6, 4); elmK(18, 13) = kMembrane(6, 5); elmK(18, 18) = kMembrane(6, 6); elmK(18, 19) = kMembrane(6, 7);
        elmK(19, 0) = kMembrane(7, 0); elmK(19, 1) = kMembrane(7, 1); elmK(19, 6) = kMembrane(7, 2); elmK(19, 7) = kMembrane(7, 3); elmK(19, 12) = kMembrane(7, 4); elmK(19, 13) = kMembrane(7, 5); elmK(19, 18) = kMembrane(7, 6); elmK(19, 19) = kMembrane(7, 7);
    }

    // Plate action is not implemented yet

    this->LocalCoordinateStiffnessMatrix = std::make_shared<Matrix<double>>(elmK);
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
    // local and global matrices are same if their normal vectors are in global positive Z-Direction and their local X-Directions coincide with 
    // global X-Direction. So, use elements local x-axis and skew angle to find rotation matrix. Skew angle is the angle between global-Z direction
    // and element normak vector, which is local-z vector. Note that nodes are assumed to be in counter-clockwise orientation.

    // Find element local-x vector
    Vector xVector(Nodes[0]->Coordinate, Nodes[1]->Coordinate);

    // Find skew angle
    //  Multiplication of two vector in a plane gives normal vector
    Vector secondPlaneVector(Nodes[1]->Coordinate, Nodes[2]->Coordinate);
    auto normalVector = xVector * secondPlaneVector;
    Vector globalZ(0.0, 0.0, 1.0);
    double skewAngle = globalZ.AngleTo(normalVector);

    auto minorRotMat = GeometryHelper::GetTranslationalRotationMatrix(xVector, skewAngle);

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
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 21, j - 21);

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
    auto rotTrans = (*this->RotationMatrix).transpose();
    *this->GlobalCoordinateStiffnessMatrix =
        (rotTrans * (*this->LocalCoordinateStiffnessMatrix)) * (*this->RotationMatrix);
}

void ShellMember::AssembleElementGlobalDampingMatrix(double mult1, double mult2)
{
    auto mat1 = *GlobalCoordinateMassMatrix * mult1;
    auto mat2 = *GlobalCoordinateStiffnessMatrix * mult2;
    this->GlobalCoordinateDampingMatrix =
        std::make_shared<Matrix<double>>(mat1 + mat2);
}

Matrix<double> ShellMember::InvertMatrix4(Matrix<double> m)
{
    double A2323 = m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2);
    double A1323 = m(2, 1) * m(3, 3) - m(2, 3) * m(3, 1);
    double A1223 = m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1);
    double A0323 = m(2, 0) * m(3, 3) - m(2, 3) * m(3, 0);
    double A0223 = m(2, 0) * m(3, 2) - m(2, 2) * m(3, 0);
    double A0123 = m(2, 0) * m(3, 1) - m(2, 1) * m(3, 0);
    double A2313 = m(1, 2) * m(3, 3) - m(1, 3) * m(3, 2);
    double A1313 = m(1, 1) * m(3, 3) - m(1, 3) * m(3, 1);
    double A1213 = m(1, 1) * m(3, 2) - m(1, 2) * m(3, 1);
    double A2312 = m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2);
    double A1312 = m(1, 1) * m(2, 3) - m(1, 3) * m(2, 1);
    double A1212 = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1);
    double A0313 = m(1, 0) * m(3, 3) - m(1, 3) * m(3, 0);
    double A0213 = m(1, 0) * m(3, 2) - m(1, 2) * m(3, 0);
    double A0312 = m(1, 0) * m(2, 3) - m(1, 3) * m(2, 0);
    double A0212 = m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0);
    double A0113 = m(1, 0) * m(3, 1) - m(1, 1) * m(3, 0);
    double A0112 = m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0);

    auto det = m(0, 0) * (m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223)
        - m(0, 1) * (m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223)
        + m(0, 2) * (m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123)
        - m(0, 3) * (m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123);
    det = 1 / det;

    Matrix<double> im(4, 4);

    im(0, 0) = det * (m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223);
    im(0, 1) = det * -(m(0, 1) * A2323 - m(0, 2) * A1323 + m(0, 3) * A1223);
    im(0, 2) = det * (m(0, 1) * A2313 - m(0, 2) * A1313 + m(0, 3) * A1213);
    im(0, 3) = det * -(m(0, 1) * A2312 - m(0, 2) * A1312 + m(0, 3) * A1212);
    im(1, 0) = det * -(m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223);
    im(1, 1) = det * (m(0, 0) * A2323 - m(0, 2) * A0323 + m(0, 3) * A0223);
    im(1, 2) = det * -(m(0, 0) * A2313 - m(0, 2) * A0313 + m(0, 3) * A0213);
    im(1, 3) = det * (m(0, 0) * A2312 - m(0, 2) * A0312 + m(0, 3) * A0212);
    im(2, 0) = det * (m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123);
    im(2, 1) = det * -(m(0, 0) * A1323 - m(0, 1) * A0323 + m(0, 3) * A0123);
    im(2, 2) = det * (m(0, 0) * A1313 - m(0, 1) * A0313 + m(0, 3) * A0113);
    im(2, 3) = det * -(m(0, 0) * A1312 - m(0, 1) * A0312 + m(0, 3) * A0112);
    im(3, 0) = det * -(m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123);
    im(3, 1) = det * (m(0, 0) * A1223 - m(0, 1) * A0223 + m(0, 2) * A0123);
    im(3, 2) = det * -(m(0, 0) * A1213 - m(0, 1) * A0213 + m(0, 2) * A0113);
    im(3, 3) = det * (m(0, 0) * A1212 - m(0, 1) * A0212 + m(0, 2) * A0112);

    return im;
}
