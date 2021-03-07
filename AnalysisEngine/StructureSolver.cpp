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
