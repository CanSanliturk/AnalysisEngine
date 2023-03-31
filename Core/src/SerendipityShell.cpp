#include "SerendipityShell.h"
// #include <armadillo>
#include "Vector.h"
#include "GeometryHelper.h"
#include "UtilMethods.h"

SerendipityShell::SerendipityShell(unsigned int elmIndex, std::shared_ptr<Node> iNode, std::shared_ptr<Node> jNode, std::shared_ptr<Node> kNode, std::shared_ptr<Node> lNode, std::shared_ptr<Node> ijNode, std::shared_ptr<Node> jkNode, std::shared_ptr<Node> klNode, std::shared_ptr<Node> liNode, std::shared_ptr<Material> material, double thickness, bool isLumpedMassMatrix, double rayleighDampingMassMultiplier, double rayleighDampingStiffnessMultiplier)
{
    ElementIndex = elmIndex;

    Nodes.resize(8);
    Nodes[0] = iNode;
    Nodes[1] = jNode;
    Nodes[2] = kNode;
    Nodes[3] = lNode;
    Nodes[4] = ijNode;
    Nodes[5] = jkNode;
    Nodes[6] = klNode;
    Nodes[7] = liNode;

    SerendipityMaterial = material;

    Thickness = thickness;

    Type = ElmType::ElementType::SerendipityShell;

    this->isLumpedMassMatrix = isLumpedMassMatrix;

    LocalCoordinateMassMatrix = std::make_shared<Matrix<double>>(48);
    LocalCoordinateStiffnessMatrix = std::make_shared<Matrix<double>>(48);
    LocalCoordinateDampingMatrix = std::make_shared<Matrix<double>>(48);
    RotationMatrix = std::make_shared<Matrix<double>>(48);
    GlobalCoordinateStiffnessMatrix = std::make_shared<Matrix<double>>(48);
    GlobalCoordinateMassMatrix = std::make_shared<Matrix<double>>(48);
    GlobalCoordinateDampingMatrix = std::make_shared<Matrix<double>>(48);

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
    ijNode->ConnectedElements.push_back(elmIndex);
    jkNode->ConnectedElements.push_back(elmIndex);
    klNode->ConnectedElements.push_back(elmIndex);
    liNode->ConnectedElements.push_back(elmIndex);
}

SerendipityShell::SerendipityShell()
{
}

SerendipityShell::~SerendipityShell()
{
}

unsigned int SerendipityShell::GetElementIndex()
{
    return ElementIndex;
}

unsigned int SerendipityShell::GetNumberOfDoF()
{
    return 48;
}

ElmType::ElementType SerendipityShell::GetElementTdsaype()
{
    return ElmType::ElementType::SerendipityShell;
}

std::shared_ptr<Matrix<double>> SerendipityShell::GetLocalCoordinateMassMatrix()
{
    return LocalCoordinateMassMatrix;
}

std::shared_ptr<Matrix<double>> SerendipityShell::GetLocalCoordinateStiffnessMatrix()
{
    return LocalCoordinateStiffnessMatrix;
}

std::shared_ptr<Matrix<double>> SerendipityShell::GetLocalCoordinateDampingMatrix()
{
    return LocalCoordinateDampingMatrix;
}

std::shared_ptr<Matrix<double>> SerendipityShell::GetGlobalCoordinateMassMatrix()
{
    return GlobalCoordinateMassMatrix;
}

std::shared_ptr<Matrix<double>> SerendipityShell::GetGlobalCoordinateStiffnessMatrix()
{
    return GlobalCoordinateStiffnessMatrix;
}

std::shared_ptr<Matrix<double>> SerendipityShell::GetGlobalCoordinateDampingMatrix()
{
    return GlobalCoordinateDampingMatrix;
}

std::shared_ptr<Matrix<double>> SerendipityShell::GetRotationMatrix()
{
    return RotationMatrix;
}

std::shared_ptr<Matrix<double>> SerendipityShell::GetElementLoads()
{
    return std::shared_ptr<Matrix<double>>();
}

std::vector<std::shared_ptr<Node>> SerendipityShell::GelElementNodes()
{
    return Nodes;
}

void SerendipityShell::AssembleElementLocalMassMatrix()
{
    auto asd = std::make_shared<Matrix<double>>(48, 48);
    this->LocalCoordinateMassMatrix = asd;
}

void SerendipityShell::AssembleElementLocalStiffnessMatrix()
{
    // Stiffness matrix is calculated at four Gauss points using Gauss Quadrature
    // It is assumed that nodes are oriented in counter-clock wise direction

    // Map coordinates of flat plate to 2-D surface
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
    auto x5 = (x1 + x2) / 2; auto y5 = (y1 + y2) / 2;
    auto x6 = (x2 + x3) / 2; auto y6 = (y2 + y3) / 2;
    auto x7 = (x3 + x4) / 2; auto y7 = (y3 + y4) / 2;
    auto x8 = (x4 + x1) / 2; auto y8 = (y4 + y1) / 2;

    Matrix<double> mappedCoords(8, 2);
    mappedCoords(0, 0) = x1; mappedCoords(0, 1) = y1;
    mappedCoords(1, 0) = x2; mappedCoords(1, 1) = y2;
    mappedCoords(2, 0) = x3; mappedCoords(2, 1) = y3;
    mappedCoords(3, 0) = x4; mappedCoords(3, 1) = y4;
    mappedCoords(4, 0) = x5; mappedCoords(4, 1) = y5;
    mappedCoords(5, 0) = x6; mappedCoords(5, 1) = y6;
    mappedCoords(6, 0) = x7; mappedCoords(6, 1) = y7;
    mappedCoords(7, 0) = x8; mappedCoords(7, 1) = y8;

    // Stiffness matrices
    Matrix<double> kElm(48, 48);
    Matrix<double> kPl(24, 24);
    Matrix<double> kBending(24, 24);
    Matrix<double> kShear(24, 24);

    // Constitutive matrices
    Matrix<double> flexuralRigidity(3, 3);
    auto elas = this->SerendipityMaterial->E;
    auto pois = this->SerendipityMaterial->PoissonsRatio;
    auto thick = this->Thickness;
    auto fRMult = elas * (thick * thick * thick) / (12 * (1 - (pois * pois)));
    flexuralRigidity(0, 0) = fRMult * 1; flexuralRigidity(0, 1) = fRMult * pois;
    flexuralRigidity(1, 0) = fRMult * pois; flexuralRigidity(1, 1) = fRMult * 1;
    flexuralRigidity(2, 2) = fRMult * (1 - pois) / 2;
    Matrix<double> shearRigidity(2, 2);
    double sR = (5.0 / 6.0) * this->SerendipityMaterial->G * thick;
    shearRigidity(0, 0) = sR; shearRigidity(1, 1) = sR;
    
    // To reduce number of mechanism, calculate kBending by using 3x3 Gauss points and
    // calculate kShear by using 2x2 Gauss points. (Gauss points are ordered in CCW orientation)
    auto gp3 = sqrt(3.0 / 5.0);
    double gp3Points[9][2] = { {-gp3, -gp3}, {0, -gp3}, {gp3, -gp3}, {gp3, 0}, {gp3, gp3}, {0, gp3}, {-gp3, gp3}, {-gp3, 0}, {0, 0} };
    
    auto gp2 = 1 / sqrt(3.0);
    double gp2Points[4][2] = { {-gp2, -gp2}, {gp2, -gp2}, {gp2, gp2}, {-gp2, gp2} };
    auto gp2Weight = 1.0;

    // Calculate kBending for each Gauss Points
    for (size_t i = 0; i < 9; i++)
    {
        auto ksi = gp3Points[i][0]; auto eta = gp3Points[i][1];
        
        // Shape functions
        auto N1 = -0.25 * (-1 + ksi) * (-1 + eta) * (ksi + eta + 1);
        auto N2 = 0.25 * (1 + ksi) * (-1 + eta) * (-ksi + eta + 1);
        auto N3 = 0.25 * (1 + ksi) * (1 + eta) * (ksi + eta - 1);
        auto N4 = -0.25 * (-1 + ksi) * (1 + eta) * (-ksi + eta - 1);
        auto N5 = 0.50 * (-1 + (ksi * ksi)) * (-1 + eta);
        auto N6 = -0.50 * (1 + ksi) * (-1 + (eta * eta));
        auto N7 = -0.50 * (-1 + (ksi * ksi)) * (1 + eta);
        auto N8 = 0.50 * (-1 + ksi) * (-1 + (eta * eta));

        // Derivatives of shape functions with respect to natural coordinates
        auto dN1Ksi = 0.25 * ((2 * ksi) - (2 * eta * ksi) - (eta * eta) + (eta));
        auto dN2Ksi = 0.25 * ((2 * ksi) - (2 * eta * ksi) + (eta * eta) - (eta));
        auto dN3Ksi = 0.25 * ((2 * ksi) + (2 * eta * ksi) + (eta * eta) + (eta));
        auto dN4Ksi = 0.25 * ((2 * ksi) + (2 * eta * ksi) - (eta * eta) - (eta));
        auto dN5Ksi = 0.25 * (4 * (ksi) * (-1 + eta));
        auto dN6Ksi = 0.25 * (2 - (2 * eta * eta));
        auto dN7Ksi = 0.25 * (-4 * (ksi) * (1 + eta));
        auto dN8Ksi = 0.25 * (-2 + (2 * eta * eta));

        auto dN1Eta = 0.25 * ((2 * eta) - (ksi * ksi) - (2 * eta * ksi) + (ksi));
        auto dN2Eta = 0.25 * ((2 * eta) - (ksi * ksi) + (2 * eta * ksi) - (ksi));
        auto dN3Eta = 0.25 * ((2 * eta) + (ksi * ksi) + (2 * eta * ksi) + (ksi));
        auto dN4Eta = 0.25 * ((2 * eta) + (ksi * ksi) - (2 * eta * ksi) - (ksi));
        auto dN5Eta = 0.25 * (-2 + (2 * ksi * ksi));
        auto dN6Eta = 0.25 * (-4 * (1 + ksi) * (eta));
        auto dN7Eta = 0.25 * (2 - (2 * ksi * ksi));
        auto dN8Eta = 0.25 * (4 * (-1 + ksi) * eta);

        // Calculate jacobian and inverse jacobian matrix for coordinate transformation
        Matrix<double> j1(2, 8);
        j1(0, 0) = dN1Ksi; j1(0, 1) = dN2Ksi; j1(0, 2) = dN3Ksi; j1(0, 3) = dN4Ksi; j1(0, 4) = dN5Ksi; j1(0, 5) = dN6Ksi; j1(0, 6) = dN7Ksi; j1(0, 7) = dN8Ksi;
        j1(1, 0) = dN1Eta; j1(1, 1) = dN2Eta; j1(1, 2) = dN3Eta; j1(1, 3) = dN4Eta; j1(1, 4) = dN5Eta; j1(1, 5) = dN6Eta; j1(1, 6) = dN7Eta; j1(1, 7) = dN8Eta;
        auto jacobi = j1 * mappedCoords;

        Matrix<double> inversejacobi(2, 2);
        auto detjacobi = (jacobi(0, 0) * jacobi(1, 1)) - (jacobi(0, 1) * jacobi(1, 0));
        inversejacobi(0, 0) = jacobi(1, 1) / detjacobi; inversejacobi(0, 1) = -1 * jacobi(0, 1) / detjacobi;
        inversejacobi(1, 0) = -1 * jacobi(1, 0) / detjacobi; inversejacobi(1, 1) = jacobi(0, 0) / detjacobi;

        // Derivatives of shape functions with respect to natural coordinates
        auto dN1X = ((inversejacobi(0, 0) * dN1Ksi) + (inversejacobi(0, 1) * dN1Eta));
        auto dN2X = ((inversejacobi(0, 0) * dN2Ksi) + (inversejacobi(0, 1) * dN2Eta));
        auto dN3X = ((inversejacobi(0, 0) * dN3Ksi) + (inversejacobi(0, 1) * dN3Eta));
        auto dN4X = ((inversejacobi(0, 0) * dN4Ksi) + (inversejacobi(0, 1) * dN4Eta));
        auto dN5X = ((inversejacobi(0, 0) * dN5Ksi) + (inversejacobi(0, 1) * dN5Eta));
        auto dN6X = ((inversejacobi(0, 0) * dN6Ksi) + (inversejacobi(0, 1) * dN6Eta));
        auto dN7X = ((inversejacobi(0, 0) * dN7Ksi) + (inversejacobi(0, 1) * dN7Eta));
        auto dN8X = ((inversejacobi(0, 0) * dN8Ksi) + (inversejacobi(0, 1) * dN8Eta));

        auto dN1Y = ((inversejacobi(1, 0) * dN1Ksi) + (inversejacobi(1, 1) * dN1Eta));
        auto dN2Y = ((inversejacobi(1, 0) * dN2Ksi) + (inversejacobi(1, 1) * dN2Eta));
        auto dN3Y = ((inversejacobi(1, 0) * dN3Ksi) + (inversejacobi(1, 1) * dN3Eta));
        auto dN4Y = ((inversejacobi(1, 0) * dN4Ksi) + (inversejacobi(1, 1) * dN4Eta));
        auto dN5Y = ((inversejacobi(1, 0) * dN5Ksi) + (inversejacobi(1, 1) * dN5Eta));
        auto dN6Y = ((inversejacobi(1, 0) * dN6Ksi) + (inversejacobi(1, 1) * dN6Eta));
        auto dN7Y = ((inversejacobi(1, 0) * dN7Ksi) + (inversejacobi(1, 1) * dN7Eta));
        auto dN8Y = ((inversejacobi(1, 0) * dN8Ksi) + (inversejacobi(1, 1) * dN8Eta));

        // Create curvature-displacement relation
        Matrix<double> bB(3, 24);
        bB(0, 2) = dN1X; bB(0, 5) = dN2X; bB(0, 8) = dN3X; bB(0, 11) = dN4X; bB(0, 14) = dN5X; bB(0, 17) = dN6X; bB(0, 20) = dN7X; bB(0, 23) = dN8X;
        bB(1, 1) = -dN1Y; bB(1, 4) = -dN2Y; bB(1, 7) = -dN3Y; bB(1, 10) = -dN4Y; bB(1, 13) = -dN5Y; bB(1, 16) = -dN6Y; bB(1, 19) = -dN7Y; bB(1, 22) = -dN8Y;
        bB(2, 1) = -dN1X; bB(2, 4) = -dN2X; bB(2, 7) = -dN3X; bB(2, 10) = -dN4X; bB(2, 13) = -dN5X; bB(2, 16) = -dN6X; bB(2, 19) = -dN7X; bB(2, 22) = -dN8X;
        bB(2, 2) = dN1Y; bB(2, 5) = dN2Y; bB(2, 8) = dN3Y; bB(2, 11) = dN4Y; bB(2, 14) = dN5Y; bB(2, 17) = dN6Y; bB(2, 20) = dN7Y; bB(2, 23) = dN8Y;
        bB *= -1.0;

        // Calculate bending stiffness matrix for Gauss Point
        auto ksiWeight = Utils::AreEqual(abs(ksi), sqrt(0.6)) ? (5.0 / 9.0) : (8.0 / 9.0);
        auto etaWeight = Utils::AreEqual(abs(eta), sqrt(0.6)) ? (5.0 / 9.0) : (8.0 / 9.0);
        auto gaussPointBendingStiffness = bB.transpose() * flexuralRigidity * bB * detjacobi * ksiWeight * etaWeight;

        // Update bending stiffness
        kBending += gaussPointBendingStiffness;
    }

    // Calculate kShear for each Gauss Points
    for (size_t i = 0; i < 4; i++)
    {
        auto ksi = gp2Points[i][0]; auto eta = gp2Points[i][1];

        // Shape functions
        auto N1 = -0.25 * (-1 + ksi) * (-1 + eta) * (ksi + eta + 1);
        auto N2 = 0.25 * (1 + ksi) * (-1 + eta) * (-ksi + eta + 1);
        auto N3 = 0.25 * (1 + ksi) * (1 + eta) * (ksi + eta - 1);
        auto N4 = -0.25 * (-1 + ksi) * (1 + eta) * (-ksi + eta - 1);
        auto N5 = 0.50 * (-1 + (ksi * ksi)) * (-1 + eta);
        auto N6 = -0.50 * (1 + ksi) * (-1 + (eta * eta));
        auto N7 = -0.50 * (-1 + (ksi * ksi)) * (1 + eta);
        auto N8 = 0.50 * (-1 + ksi) * (-1 + (eta * eta));

        // Derivatives of shape functions with respect to natural coordinates
        auto dN1Ksi = 0.25 * ((2 * ksi) - (2 * eta * ksi) - (eta * eta) + (eta));
        auto dN2Ksi = 0.25 * ((2 * ksi) - (2 * eta * ksi) + (eta * eta) - (eta));
        auto dN3Ksi = 0.25 * ((2 * ksi) + (2 * eta * ksi) + (eta * eta) + (eta));
        auto dN4Ksi = 0.25 * ((2 * ksi) + (2 * eta * ksi) - (eta * eta) - (eta));
        auto dN5Ksi = 0.25 * (4 * (ksi) * (-1 + eta));
        auto dN6Ksi = 0.25 * (2 - (2 * eta * eta));
        auto dN7Ksi = 0.25 * (-4 * (ksi) * (1 + eta));
        auto dN8Ksi = 0.25 * (-2 + (2 * eta * eta));

        auto dN1Eta = 0.25 * ((2 * eta) - (ksi * ksi) - (2 * eta * ksi) + (ksi));
        auto dN2Eta = 0.25 * ((2 * eta) - (ksi * ksi) + (2 * eta * ksi) - (ksi));
        auto dN3Eta = 0.25 * ((2 * eta) + (ksi * ksi) + (2 * eta * ksi) + (ksi));
        auto dN4Eta = 0.25 * ((2 * eta) + (ksi * ksi) - (2 * eta * ksi) - (ksi));
        auto dN5Eta = 0.25 * (-2 + (2 * ksi * ksi));
        auto dN6Eta = 0.25 * (-4 * (1 + ksi) * (eta));
        auto dN7Eta = 0.25 * (2 - (2 * ksi * ksi));
        auto dN8Eta = 0.25 * (4 * (-1 + ksi) * eta);

        // Calculate jacobian and inverse jacobian matrix for coordinate transformation
        Matrix<double> j1(2, 8);
        j1(0, 0) = dN1Ksi; j1(0, 1) = dN2Ksi; j1(0, 2) = dN3Ksi; j1(0, 3) = dN4Ksi; j1(0, 4) = dN5Ksi; j1(0, 5) = dN6Ksi; j1(0, 6) = dN7Ksi; j1(0, 7) = dN8Ksi;
        j1(1, 0) = dN1Eta; j1(1, 1) = dN2Eta; j1(1, 2) = dN3Eta; j1(1, 3) = dN4Eta; j1(1, 4) = dN5Eta; j1(1, 5) = dN6Eta; j1(1, 6) = dN7Eta; j1(1, 7) = dN8Eta;
        auto jacobi = j1 * mappedCoords;

        Matrix<double> inversejacobi(2, 2);
        auto detjacobi = (jacobi(0, 0) * jacobi(1, 1)) - (jacobi(0, 1) * jacobi(1, 0));
        inversejacobi(0, 0) = jacobi(1, 1) / detjacobi; inversejacobi(0, 1) = -1 * jacobi(0, 1) / detjacobi;
        inversejacobi(1, 0) = -1 * jacobi(1, 0) / detjacobi; inversejacobi(1, 1) = jacobi(0, 0) / detjacobi;

        // Derivatives of shape functions with respect to natural coordinates
        auto dN1X = ((inversejacobi(0, 0) * dN1Ksi) + (inversejacobi(0, 1) * dN1Eta));
        auto dN2X = ((inversejacobi(0, 0) * dN2Ksi) + (inversejacobi(0, 1) * dN2Eta));
        auto dN3X = ((inversejacobi(0, 0) * dN3Ksi) + (inversejacobi(0, 1) * dN3Eta));
        auto dN4X = ((inversejacobi(0, 0) * dN4Ksi) + (inversejacobi(0, 1) * dN4Eta));
        auto dN5X = ((inversejacobi(0, 0) * dN5Ksi) + (inversejacobi(0, 1) * dN5Eta));
        auto dN6X = ((inversejacobi(0, 0) * dN6Ksi) + (inversejacobi(0, 1) * dN6Eta));
        auto dN7X = ((inversejacobi(0, 0) * dN7Ksi) + (inversejacobi(0, 1) * dN7Eta));
        auto dN8X = ((inversejacobi(0, 0) * dN8Ksi) + (inversejacobi(0, 1) * dN8Eta));

        auto dN1Y = ((inversejacobi(1, 0) * dN1Ksi) + (inversejacobi(1, 1) * dN1Eta));
        auto dN2Y = ((inversejacobi(1, 0) * dN2Ksi) + (inversejacobi(1, 1) * dN2Eta));
        auto dN3Y = ((inversejacobi(1, 0) * dN3Ksi) + (inversejacobi(1, 1) * dN3Eta));
        auto dN4Y = ((inversejacobi(1, 0) * dN4Ksi) + (inversejacobi(1, 1) * dN4Eta));
        auto dN5Y = ((inversejacobi(1, 0) * dN5Ksi) + (inversejacobi(1, 1) * dN5Eta));
        auto dN6Y = ((inversejacobi(1, 0) * dN6Ksi) + (inversejacobi(1, 1) * dN6Eta));
        auto dN7Y = ((inversejacobi(1, 0) * dN7Ksi) + (inversejacobi(1, 1) * dN7Eta));
        auto dN8Y = ((inversejacobi(1, 0) * dN8Ksi) + (inversejacobi(1, 1) * dN8Eta));

        // Create shaer strain-displacement relation
        Matrix<double> bS(2, 24);
        bS(0, 0) = dN1X; bS(0, 2) = N1; bS(0, 3) = dN2X; bS(0, 5) = N2; bS(0, 6) = dN3X; bS(0, 8) = N3; bS(0, 9) = dN4X; bS(0, 11) = N4; bS(0, 12) = dN5X; bS(0, 14) = N5; bS(0, 15) = dN6X; bS(0, 17) = N6; bS(0, 18) = dN7X; bS(0, 20) = N7; bS(0, 21) = dN8X; bS(0, 23) = N8;
        bS(1, 0) = -dN1Y;  bS(1, 1) = N1; bS(1, 3) = -dN2Y;  bS(1, 4) = N2; bS(1, 6) = -dN3Y;  bS(1, 7) = N3; bS(1, 9) = -dN4Y;  bS(1, 10) = N4; bS(1, 12) = -dN5Y;  bS(1, 13) = N5; bS(1, 15) = -dN6Y;  bS(1, 16) = N6; bS(1, 18) = -dN7Y;  bS(1, 19) = N7; bS(1, 21) = -dN8Y;  bS(1, 22) = N8;

        // Calculate shear stiffness matrix for Gauss Point
        auto gaussPointShearStiffness = bS.transpose() * shearRigidity * bS * detjacobi * gp2Weight;

        // Update shear stiffness
        kShear += gaussPointShearStiffness;
    }

    // Sum two stiffness matrices to obtain plate stiffness
    kPl = kBending + kShear;

    // Map plate stifness to element stiffness (For a plate at XY-plane, it resists translation-Z, rotation-X and rotation-Y)
    int mapper1 = 2;
    for (size_t i = 0; i < 24; i++)
    {
        int mapper2 = 2;
        for (size_t j = 0; j < 24; j++)
        {
            kElm(i + mapper1, j + mapper2) = kPl(i, j);
            if (((j + 1) % 3) == 0)
                mapper2 += 3;
        }
        if (((i + 1) % 3) == 0)
            mapper1 += 3;
    }

    this->LocalCoordinateStiffnessMatrix = std::make_shared<Matrix<double>>(kElm);
}

void SerendipityShell::AssembleElementLocalDampingMatrix(double mult1, double mult2)
{
    auto mat1 = *LocalCoordinateMassMatrix * mult1;
    auto mat2 = *LocalCoordinateStiffnessMatrix * mult2;
    this->LocalCoordinateDampingMatrix =
        std::make_shared<Matrix<double>>(mat1 + mat2);
}

void SerendipityShell::AssembleElementRotationMatrix()
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

    for (unsigned int i = 24; i < 27; i++)
        for (unsigned int j = 24; j < 27; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 24, j - 24);

    for (unsigned int i = 27; i < 30; i++)
        for (unsigned int j = 27; j < 30; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 27, j - 27);

    for (unsigned int i = 30; i < 33; i++)
        for (unsigned int j = 30; j < 33; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 30, j - 30);

    for (unsigned int i = 33; i < 36; i++)
        for (unsigned int j = 33; j < 36; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 33, j - 33);

    for (unsigned int i = 36; i < 39; i++)
        for (unsigned int j = 36; j < 39; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 36, j - 36);

    for (unsigned int i = 39; i < 42; i++)
        for (unsigned int j = 39; j < 42; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 39, j - 39);

    for (unsigned int i = 42; i < 45; i++)
        for (unsigned int j = 42; j < 45; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 42, j - 42);

    for (unsigned int i = 45; i < 48; i++)
        for (unsigned int j = 45; j < 48; j++)
            (*this->RotationMatrix)(i, j) = minorRotMat(i - 45, j - 45);
}

void SerendipityShell::AssembleElementGlobalMassMatrix()
{
    auto rotTrans = (*this->RotationMatrix).transpose();
    *this->GlobalCoordinateMassMatrix =
        (rotTrans * (*this->LocalCoordinateMassMatrix)) * (*this->RotationMatrix);
}

void SerendipityShell::AssembleElementGlobalStiffnessMatrix()
{
    auto rotTrans = (*this->RotationMatrix).transpose();
    *this->GlobalCoordinateStiffnessMatrix =
        (rotTrans * (*this->LocalCoordinateStiffnessMatrix)) * (*this->RotationMatrix);
}

void SerendipityShell::AssembleElementGlobalDampingMatrix(double mult1, double mult2)
{
    auto mat1 = *GlobalCoordinateMassMatrix * mult1;
    auto mat2 = *GlobalCoordinateStiffnessMatrix * mult2;
    this->GlobalCoordinateDampingMatrix =
        std::make_shared<Matrix<double>>(mat1 + mat2);
}
