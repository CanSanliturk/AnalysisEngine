#include <armadillo>
#include <iostream>
#include <memory>
#include <Eigen>
#include <vector>
#include <map>
#include "Eigen/Eigenvalues"
#include "StructureSolver.h"
#include "GeometryHelper.h"
#include "Matrix.h"

constexpr double pi = 3.141592653589793;

// Solves K * u = F for unknown u vector (Considering support displacements)
Matrix<double> StructureSolver::GetDisplacementForStaticCase(const Structure& str, SolverChoice solverChoice)
{
    // Define variables
    auto nDofUnrestrained = str.nUnrestrainedDOF;
    auto nDof = str.nDOF;
    auto nDofRestrained = nDof - nDofUnrestrained;

    // Get necessary parts for computations
    auto kUU = str.StiffnessMatrix->getSubmatrix(0, str.nUnrestrainedDOF - 1, 0, str.nUnrestrainedDOF - 1);
    auto kUK = str.StiffnessMatrix->getSubmatrix(0, str.nUnrestrainedDOF - 1, str.nUnrestrainedDOF, str.nDOF - 1);
    auto kKU = kUK.transpose();
    auto kKK = str.StiffnessMatrix->getSubmatrix(str.nUnrestrainedDOF, str.nDOF - 1, str.nUnrestrainedDOF, str.nDOF - 1);
    auto fK = str.ForceVector->getSubmatrix(0, str.nUnrestrainedDOF - 1, 0, 0);

    // Subtract support settlements from force vector
    Matrix<double> uK(str.nDOF - str.nUnrestrainedDOF, 1);

    auto& restraints = *str.Restraints;
    std::map<unsigned int, std::shared_ptr<Restraint>>::iterator iter = restraints.begin();

    // Iterate over the map using Iterator till end.
    while (iter != restraints.end())
    {
        auto restraint = iter->second;
        auto restrainedNode = restraint->RestrainedNode;

        if (restraint->IsRestraintTranslationX)
            uK(restrainedNode->DofIndexTX - nDofUnrestrained - 1, 0) = restraint->TranslationX;
        if (restraint->IsRestraintTranslationY)
            uK(restrainedNode->DofIndexTY - nDofUnrestrained - 1, 0) = restraint->TranslationY;
        if (restraint->IsRestraintTranslationZ)
            uK(restrainedNode->DofIndexTZ - nDofUnrestrained - 1, 0) = restraint->TranslationZ;
        if (restraint->IsRestraintRotationX)
            uK(restrainedNode->DofIndexRX - nDofUnrestrained - 1, 0) = restraint->RotationX;
        if (restraint->IsRestraintRotationY)
            uK(restrainedNode->DofIndexRY - nDofUnrestrained - 1, 0) = restraint->RotationY;
        if (restraint->IsRestraintRotationZ)
            uK(restrainedNode->DofIndexRZ - nDofUnrestrained - 1, 0) = restraint->RotationZ;
        iter++;
    }

    auto subtractVal = kUK * uK;
    auto condFVec = fK - subtractVal;

    // Solve system
    auto resData = LinearEquationSolver(kUU, condFVec, solverChoice);

    // Define a return value. First nUnrestrainedDOF elements will come from resData. Rest is 
    // restrained displacements
    Matrix<double> retVal(nDof, 1);

    for (size_t i = 0; i < nDofUnrestrained; i++)
        retVal(i, 0) = resData(i, 0);

    std::map<unsigned int, std::shared_ptr<Restraint>>::iterator it = restraints.begin();

    // Iterate over the map using Iterator till end.
    while (it != restraints.end())
    {
        auto restraint = it->second;
        auto restrainedNode = restraint->RestrainedNode;

        if (restraint->IsRestraintTranslationX)
            retVal(restrainedNode->DofIndexTX - 1, 0) = restraint->TranslationX;
        if (restraint->IsRestraintTranslationY)
            retVal(restrainedNode->DofIndexTY - 1, 0) = restraint->TranslationY;
        if (restraint->IsRestraintTranslationZ)
            retVal(restrainedNode->DofIndexTZ - 1, 0) = restraint->TranslationZ;
        if (restraint->IsRestraintRotationX)
            retVal(restrainedNode->DofIndexRX - 1, 0) = restraint->RotationX;
        if (restraint->IsRestraintRotationY)
            retVal(restrainedNode->DofIndexRY - 1, 0) = restraint->RotationY;
        if (restraint->IsRestraintRotationZ)
            retVal(restrainedNode->DofIndexRZ - 1, 0) = restraint->RotationZ;

        it++;
    }

    return retVal;
}

// Returns member forces at local coordinates
Matrix<double> StructureSolver::GetMemberEndForcesForLocalCoordinates(Element& elm, Matrix<double>& displacements)
{
    // To get elements end forces in local coordinates, displacement vector of element(which is in global coordinates)
    // should be converted to local coordinates by multiplying it by rotation matrix of the given element.
    auto& rotMat = (*elm.GetRotationMatrix());
    auto nDofElm = elm.GetNumberOfDoF();

    // Retrieve displacements of element end nodes
    auto elmNodes = elm.GelElementNodes();
    Matrix<double> disps(nDofElm, 1);

    unsigned short counter = 0;

    for (auto node : elmNodes)
    {
        disps((counter * 6) + 0, 0) = displacements(node->DofIndexTX - 1, 0);
        disps((counter * 6) + 1, 0) = displacements(node->DofIndexTY - 1, 0);
        disps((counter * 6) + 2, 0) = displacements(node->DofIndexTZ - 1, 0);
        disps((counter * 6) + 3, 0) = displacements(node->DofIndexRX - 1, 0);
        disps((counter * 6) + 4, 0) = displacements(node->DofIndexRY - 1, 0);
        disps((counter * 6) + 5, 0) = displacements(node->DofIndexRZ - 1, 0);
        counter++;
    }

    // Get displacements at local coordinates
    auto localizedDisps = rotMat * disps;

    return (*elm.GetLocalCoordinateStiffnessMatrix()) * localizedDisps;
}

// Returns member forces at global coordinates
Matrix<double> StructureSolver::GetMemberEndForcesForGlobalCoordinates(Element& elm, Matrix<double>& displacements)
{
    // To get elements end forces in local coordinates, displacement vector of element(which is in global coordinates)
    // should be converted to local coordinates by multiplying it by rotation matrix of the given element.
    auto& rotMat = (*elm.GetRotationMatrix());
    auto nDofElm = elm.GetNumberOfDoF();

    // Retrieve displacements of element end nodes
    auto elmNodes = elm.GelElementNodes();
    Matrix<double> disps(nDofElm, 1);

    unsigned short counter = 0;

    for (auto node : elmNodes)
    {
        disps((counter * 6) + 0, 0) = displacements(node->DofIndexTX - 1, 0);
        disps((counter * 6) + 1, 0) = displacements(node->DofIndexTY - 1, 0);
        disps((counter * 6) + 2, 0) = displacements(node->DofIndexTZ - 1, 0);
        disps((counter * 6) + 3, 0) = displacements(node->DofIndexRX - 1, 0);
        disps((counter * 6) + 4, 0) = displacements(node->DofIndexRY - 1, 0);
        disps((counter * 6) + 5, 0) = displacements(node->DofIndexRZ - 1, 0);
        counter++;
    }

    // Get displacements at local coordinates
    auto localizedDisps = rotMat * disps;

    return (*elm.GetGlobalCoordinateStiffnessMatrix()) * localizedDisps;
}

// Returns displacements for given node
Matrix<double> StructureSolver::GetNodalDisplacements(Node& node, Matrix<double>& displacements)
{
    Matrix<double> nodalDips(6, 1);

    nodalDips(0, 0) = displacements(node.DofIndexTX - 1, 0);
    nodalDips(1, 0) = displacements(node.DofIndexTY - 1, 0);
    nodalDips(2, 0) = displacements(node.DofIndexTZ - 1, 0);
    nodalDips(3, 0) = displacements(node.DofIndexRX - 1, 0);
    nodalDips(4, 0) = displacements(node.DofIndexRY - 1, 0);
    nodalDips(5, 0) = displacements(node.DofIndexRZ - 1, 0);

    return nodalDips;
}

// Returns support reactions
Matrix<double> StructureSolver::GetSupportReactions(const Structure& str, Matrix<double>& disps, const Restraint& res, SolverChoice solverChoice)
{
    auto& fullStiffnessMatrix = *str.StiffnessMatrix;
    Matrix<double> retVal(str.nDOF, 1);

    // Perform multiplication using libraries for performance issues
    if (solverChoice == SolverChoice::Eigen)
    {
        // Convert matrices to eigen matrices
        Eigen::MatrixXd aEigenMat(fullStiffnessMatrix.RowCount, fullStiffnessMatrix.ColCount);
        Eigen::VectorXd bEigenVec(str.nDOF);
        for (unsigned int i = 0; i < fullStiffnessMatrix.RowCount; i++)
        {
            for (unsigned int j = 0; j < fullStiffnessMatrix.ColCount; j++)
            {
                aEigenMat(i, j) = fullStiffnessMatrix(i, j);
            }
            bEigenVec(i) = disps(i, 0);
        }

        // Multiply
        Eigen::VectorXd xEigenVec = aEigenMat * bEigenVec;

        // Fill in return value
        for (size_t i = 0; i < fullStiffnessMatrix.RowCount; i++)
            retVal(i, 0) = xEigenVec(i);
    }
    else if (solverChoice == SolverChoice::Armadillo)
    {
        // Convert matrices to armadillo matrices
        arma::mat aArmaMat(fullStiffnessMatrix.RowCount, fullStiffnessMatrix.ColCount);
        arma::vec bArmaVec(fullStiffnessMatrix.RowCount);
        for (unsigned int i = 0; i < fullStiffnessMatrix.RowCount; i++)
        {
            for (unsigned int j = 0; j < fullStiffnessMatrix.ColCount; j++)
            {
                aArmaMat(i, j) = fullStiffnessMatrix(i, j);
            }
            bArmaVec(i) = disps(i, 0);
        }

        // Solve
        arma::vec xArmaVec = aArmaMat * bArmaVec;

        // Fill in return value
        for (size_t i = 0; i < fullStiffnessMatrix.RowCount; i++)
            retVal(i, 0) = xArmaVec(i);
    }
    else
    {
    }

    return retVal;
}

// Return modal periods in decresing order
Matrix<double> StructureSolver::GetModalPeriods(const Structure& str, SolverChoice solverChoice)
{
    // Return value
    Matrix<double> periods(str.nUnrestrainedDOF, 1);

    // Temporary return value (to be able to sort data)
    std::vector<double> t;

    if (solverChoice == SolverChoice::Eigen)
    {

        // Unrestrained parts of stiffness and mass matrices
        Eigen::MatrixXd M(str.nUnrestrainedDOF, str.nUnrestrainedDOF);
        Eigen::MatrixXd K(str.nUnrestrainedDOF, str.nUnrestrainedDOF);

        for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
        {
            for (size_t j = 0; j < str.nUnrestrainedDOF; j++)
            {
                K(i, j) = (*str.StiffnessMatrix)(i, j);
                M(i, j) = (*str.MassMatrix)(i, j);
            }
        }

        Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
        ges.compute(K, M);
        auto eigs = ges.eigenvalues().real();

        for (int i = eigs.size() - 1; -1 < i; i--)
            t.push_back(2 * pi / sqrt(eigs(i)));
    }
    else if (solverChoice == SolverChoice::Armadillo)
    {
        // Unrestrained parts of stiffness and mass matrices
        arma::mat M(str.nUnrestrainedDOF, str.nUnrestrainedDOF);
        arma::mat K(str.nUnrestrainedDOF, str.nUnrestrainedDOF);

        for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
        {
            for (size_t j = 0; j < str.nUnrestrainedDOF; j++)
            {
                K(i, j) = (*str.StiffnessMatrix)(i, j);
                M(i, j) = (*str.MassMatrix)(i, j);
            }
        }

        arma::cx_vec eigVal;
        arma::cx_mat eigVec;
        arma::eig_pair(eigVal, eigVec, K, M);

        for (int i = eigVal.size() - 1; -1 < i; i--)
            t.push_back(2 * pi / sqrt(*(eigVal.at(i)._Val)));
    }
    else {}

    std::sort(t.begin(), t.end(), [&](double a, double b) {return b < a; });

    for (int i = 0; i < t.size(); i++)
        periods(i, 0) = t[i];

    return periods;
}

Matrix<double> StructureSolver::CalculateMembraneNodalStresses(const ShellMember& elm, Matrix<double>& disps, int nodeIndex)
{
    // For now, it is assumed that membrane is at XY-plane. If displacement vector is transformed into plane
    // where shell element lies, it will work for any case.

    // Return value
    Matrix<double> stresses(3, 1);

    // Get nodal displacements
    auto iNodeDisps = GetNodalDisplacements(*elm.Nodes[0], disps);
    auto jNodeDisps = GetNodalDisplacements(*elm.Nodes[1], disps);
    auto kNodeDisps = GetNodalDisplacements(*elm.Nodes[2], disps);
    auto lNodeDisps = GetNodalDisplacements(*elm.Nodes[3], disps);
    Matrix<double> dispVector(8, 1);
    dispVector(0, 0) = iNodeDisps(0, 0);
    dispVector(1, 0) = iNodeDisps(1, 0);
    dispVector(2, 0) = jNodeDisps(0, 0);
    dispVector(3, 0) = jNodeDisps(1, 0);
    dispVector(4, 0) = kNodeDisps(0, 0);
    dispVector(5, 0) = kNodeDisps(1, 0);
    dispVector(6, 0) = lNodeDisps(0, 0);
    dispVector(7, 0) = lNodeDisps(1, 0);

    // Elasticity matrix
    auto e = elm.ShellMaterial->E;
    auto v = elm.ShellMaterial->PoissonsRatio;
    auto eMult = e / (1 - (v * v));
    Matrix<double> eMat(3, 3);
    eMat(0, 0) = eMult * 1; eMat(0, 1) = eMult * v;
    eMat(1, 0) = eMult * v; eMat(1, 1) = eMult * 1;
    eMat(2, 2) = eMult * (1 - v) / 2;

    // Thickness
    auto thickness = elm.Thickness;

    // Map coordinates of flat plane to 2-D surface
    auto d1 = elm.Nodes[0]->Coordinate.DistanceTo(elm.Nodes[1]->Coordinate);
    auto d2 = elm.Nodes[1]->Coordinate.DistanceTo(elm.Nodes[2]->Coordinate);
    auto d3 = elm.Nodes[2]->Coordinate.DistanceTo(elm.Nodes[3]->Coordinate);

    Vector p1V(elm.Nodes[0]->Coordinate);
    Vector p2V(elm.Nodes[1]->Coordinate);
    Vector p3V(elm.Nodes[2]->Coordinate);
    Vector p4V(elm.Nodes[3]->Coordinate);

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

    Matrix<double> bMatrix(3, 8);

    auto rowCounter = 0;
    // Get Gauss point
    double ksi = 0;
    double eta = 0;

    // Select mid-point for bilinear element to avoid shear locking
    if (elm.membraneType != MembraneType::Bilinear)
    {
        if (nodeIndex == 1)
        {
            ksi = -1;
            eta = -1;
        }
        else if (nodeIndex == 2)
        {
            ksi = 1;
            eta = -1;
        }
        else if (nodeIndex == 3)
        {
            ksi = 1;
            eta = 1;
        }
        else if (nodeIndex == 4)
        {
            ksi = -1;
            eta = 1;
        }
    }

    if (elm.membraneType == MembraneType::Bilinear)
    {
        // Calculate jacobi
        Matrix<double> j1(2, 4);
        j1(0, 0) = eta - 1.0; j1(0, 1) = 1.0 - eta; j1(0, 2) = eta + 1.0; j1(0, 3) = -eta - 1.0;
        j1(1, 0) = ksi - 1.0; j1(1, 1) = -ksi - 1.0; j1(1, 2) = ksi + 1.0; j1(1, 3) = 1.0 - ksi;
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
        mat3(0, 0) = eta - 1.0; mat3(0, 2) = 1.0 - eta; mat3(0, 4) = eta + 1.0; mat3(0, 6) = -eta - 1.0;
        mat3(1, 0) = ksi - 1.0; mat3(1, 2) = -ksi - 1.0; mat3(1, 4) = ksi + 1.0; mat3(1, 6) = 1.0 - ksi;
        mat3(2, 1) = eta - 1.0; mat3(2, 3) = 1.0 - eta; mat3(2, 5) = eta + 1.0; mat3(2, 7) = -eta - 1.0;
        mat3(3, 1) = ksi - 1.0; mat3(3, 3) = -ksi - 1.0; mat3(3, 5) = ksi + 1.0; mat3(3, 7) = 1.0 - ksi;
        mat3 *= 0.25;

        auto b = mat1 * mat2 * mat3;
        stresses = eMat * b * dispVector;
    }
    else if (elm.membraneType == MembraneType::Incompatible)
    {
        // Calculate jacobi
        Matrix<double> j1(2, 4);
        j1(0, 0) = eta - 1.0; j1(0, 1) = 1.0 - eta; j1(0, 2) = eta + 1.0; j1(0, 3) = -eta - 1.0;
        j1(1, 0) = ksi - 1.0; j1(1, 1) = -ksi - 1.0; j1(1, 2) = ksi + 1.0; j1(1, 3) = 1.0 - ksi;
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
        mat3(0, 0) = eta - 1.0; mat3(0, 2) = 1.0 - eta; mat3(0, 4) = eta + 1.0; mat3(0, 6) = -eta - 1.0;
        mat3(1, 0) = ksi - 1.0; mat3(1, 2) = -ksi - 1.0; mat3(1, 4) = ksi + 1.0; mat3(1, 6) = 1.0 - ksi;
        mat3(2, 1) = eta - 1.0; mat3(2, 3) = 1.0 - eta; mat3(2, 5) = eta + 1.0; mat3(2, 7) = -eta - 1.0;
        mat3(3, 1) = ksi - 1.0; mat3(3, 3) = -ksi - 1.0; mat3(3, 5) = ksi + 1.0; mat3(3, 7) = 1.0 - ksi;
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
        
        auto invert = [](Matrix<double> m)
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
        };

        // Calculate generalized displacements
        auto& elmK = elm.K1212;;
        auto kRR = elmK->getSubmatrix(0, 7, 0, 7);
        auto kRC = elmK->getSubmatrix(0, 7, 8, 11);
        auto kCR = kRC.transpose();
        auto kCC = elmK->getSubmatrix(8, 11, 8, 11);
        auto invKCC = invert(kCC);
        auto minusInvKCC = invKCC * -1;
        
        auto dc = minusInvKCC * kCR * dispVector;
        Matrix<double> updateDispVector(12, 1);
        
        for (size_t i = 0; i < 8; i++)
            updateDispVector(i, 0) = dispVector(i, 0);
        for (size_t i = 8; i < 12; i++)
            updateDispVector(i, 0) = dc(i - 8, 0);

        stresses = eMat * b * updateDispVector;
    }



    return stresses;
}

double StructureSolver::CondenseStiffnessMatrixForSpecificDOF(const Structure& str, unsigned int dofIndex, SolverChoice solverChoice)
{
    dofIndex--;
    auto rearrangedStiffnessMatrix = (*str.StiffnessMatrix).getSubmatrix(0, str.nUnrestrainedDOF - 1, 0, str.nUnrestrainedDOF - 1).
        sendToCornerForSquareMatrix(dofIndex, dofIndex, false);
    auto kRR = rearrangedStiffnessMatrix.getSubmatrix(0, 0, 0, 0);
    auto kRC = rearrangedStiffnessMatrix.getSubmatrix(0, 0, 1, rearrangedStiffnessMatrix.ColCount - 1);
    auto kCR = kRC.transpose();
    auto kCC = rearrangedStiffnessMatrix.getSubmatrix(1, rearrangedStiffnessMatrix.RowCount - 1, 1, rearrangedStiffnessMatrix.ColCount - 1);

    auto kccInv = GetInverse(kCC, solverChoice);
    auto step1 = kRC * kccInv;
    auto right = step1 * kCR;
    auto kC = kRR - right;
    double retVal = kC(0, 0);

    return retVal;
}

double StructureSolver::CondenseForceVectorForSpecificDOF(const Structure& str, unsigned int dofIndex, SolverChoice solverChoice)
{
    dofIndex--;
    auto rearrangedForceVector = (*str.ForceVector).getSubmatrix(0, str.nUnrestrainedDOF - 1, 0, 0).
        sendItemToBoundVector(dofIndex, false);

    auto rearrangedStiffnessMatrix = (*str.StiffnessMatrix).getSubmatrix(0, str.nUnrestrainedDOF - 1, 0, str.nUnrestrainedDOF - 1).
        sendToCornerForSquareMatrix(dofIndex, dofIndex, false);
    auto kRR = rearrangedStiffnessMatrix.getSubmatrix(0, 0, 0, 0);
    auto kRC = rearrangedStiffnessMatrix.getSubmatrix(0, 0, 1, rearrangedStiffnessMatrix.ColCount - 1);
    auto kCR = kRC.transpose();
    auto kCC = rearrangedStiffnessMatrix.getSubmatrix(1, rearrangedStiffnessMatrix.RowCount - 1, 1, rearrangedStiffnessMatrix.ColCount - 1);

    auto fR = rearrangedForceVector.getSubmatrix(0, 0, 0, 0);
    auto fC = rearrangedForceVector.getSubmatrix(1, str.nUnrestrainedDOF - 1, 0, 0);

    auto kccInv = GetInverse(kCC, solverChoice);
    auto step1 = kRC * kccInv;
    auto right = step1 * fC;
    auto FC = fR - right;
    double retVal = FC(0, 0);

    return retVal;
}

// Returns x as in A * x = b (Assumption is A is square matrix and b is vector)
Matrix<double> StructureSolver::LinearEquationSolver(Matrix<double>& A, Matrix<double>& b, SolverChoice solverChoice)
{
    // Check if matrix sizes do not match or assumptions are not satisfied
    if ((A.RowCount != A.ColCount) || (A.ColCount != b.RowCount) || (b.ColCount != 1))
        throw std::runtime_error("Solver Error: Matrix sizes do not match\n");

    // Create return value
    Matrix<double> retVal(A.RowCount, 1);

    if (solverChoice == SolverChoice::Eigen)
    {
        // Convert matrices to eigen matrices
        Eigen::MatrixXd aEigenMat(A.RowCount, A.ColCount);
        Eigen::VectorXd bEigenVec(b.RowCount);
        for (unsigned int i = 0; i < A.RowCount; i++)
        {
            for (unsigned int j = 0; j < A.ColCount; j++)
            {
                aEigenMat(i, j) = A(i, j);
            }
            bEigenVec(i) = b(i, 0);
        }

        // Solve
        Eigen::VectorXd xEigenVec = aEigenMat.colPivHouseholderQr().solve(bEigenVec);

        // Fill in return value
        for (size_t i = 0; i < A.RowCount; i++)
            retVal(i, 0) = xEigenVec(i);
        return retVal;
    }
    else if (solverChoice == SolverChoice::Armadillo)
    {
        // Convert matrices to armadillo matrices
        arma::mat aArmaMat(A.RowCount, A.ColCount);
        arma::vec bArmaVec(b.RowCount);
        for (unsigned int i = 0; i < A.RowCount; i++)
        {
            for (unsigned int j = 0; j < A.ColCount; j++)
            {
                aArmaMat(i, j) = A(i, j);
            }
            bArmaVec(i) = b(i, 0);
        }

        // Solve
        arma::vec xArmaVec = arma::solve(aArmaMat, bArmaVec);

        // Fill in return value
        for (size_t i = 0; i < A.RowCount; i++)
            retVal(i, 0) = xArmaVec(i);
        return retVal;
    }
    else
    {
        throw std::runtime_error("Solver Error: Unknown solver selection\n");
        return retVal;
    }
    return retVal;
}

// Returns A^-1 provided that A is square matrix
Matrix<double> StructureSolver::GetInverse(Matrix<double>& A, SolverChoice solverChoice)
{
    // Check if input is square matrix or not
    if (A.RowCount != A.ColCount)
        throw std::runtime_error("Solver Error: Provided matrix is not square\n");

    Matrix<double> invA(A.RowCount);

    if (solverChoice == SolverChoice::Eigen)
    {
        // Convert matrices to eigen matrices
        Eigen::MatrixXd aEigenMat(A.RowCount, A.ColCount);
        for (unsigned int i = 0; i < A.RowCount; i++)
            for (unsigned int j = 0; j < A.ColCount; j++)
                aEigenMat(i, j) = A(i, j);

        // Solve
        auto xEigenVec = aEigenMat.inverse();;

        // Fill in return value
        for (size_t i = 0; i < A.RowCount; i++)
            for (size_t j = 0; j < A.ColCount; j++)
                invA(i, j) = xEigenVec(i, j);

        return invA;
    }
    else if (solverChoice == SolverChoice::Armadillo)
    {
        // Convert matrices to armadillo matrices
        arma::mat aArmaMat(A.RowCount, A.ColCount);
        for (unsigned int i = 0; i < A.RowCount; i++)
            for (unsigned int j = 0; j < A.ColCount; j++)
                aArmaMat(i, j) = A(i, j);

        // Solve
        arma::vec xArmaVec = arma::inv(aArmaMat);

        // Fill in return value
        for (size_t i = 0; i < A.RowCount; i++)
            for (size_t j = 0; j < A.ColCount; j++)
                invA(i, j) = xArmaVec(i, j);

        return invA;
    }
    else
    {
        throw std::runtime_error("Solver Error: Unknown solver selection\n");
        return invA;
    }

    return invA;
}
