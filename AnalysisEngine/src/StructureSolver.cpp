#include <armadillo>
#include <iostream>
#include <memory>
#include <Eigen>
#include <vector>
#include <map>
#include "../Eigen/Eigenvalues"
#include "StructureSolver.h"
#include "GeometryHelper.h"
#include "Matrix.h"
#include "UtilMethods.h"

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

Matrix<double> StructureSolver::GetDisplacementsForStaticCaseWithForceVector(const Structure& str, Matrix<double>& fVec, SolverChoice solverChoice)
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
    auto fK = fVec.getSubmatrix(0, str.nUnrestrainedDOF - 1, 0, 0);

    // Subtract support settlements from force vector
    Matrix<double> uK(str.nDOF - str.nUnrestrainedDOF, 1);

    auto& restraints = *str.Restraints;
    std::map<unsigned int, std::shared_ptr<Restraint>>::iterator iter = restraints.begin();

    // Iterate over the map using Iterator till end.
    while (iter != restraints.end())
    {
        auto& restraint = iter->second;
        auto& restrainedNode = restraint->RestrainedNode;

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
        auto& restraint = it->second;
        auto& restrainedNode = restraint->RestrainedNode;

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

Matrix<double> StructureSolver::CalculateDisplacements(Matrix<double>& kMat, Matrix<double>& fVec, int nDof, int nUnrestainedDof, SolverChoice solverChoice)
{
    auto kNdof = kMat.getSubmatrix(0, nUnrestainedDof - 1, 0, nUnrestainedDof - 1);
    auto fNdof = fVec.getSubmatrix(0, nUnrestainedDof - 1, 0, 0);
    auto smallerResult = StructureSolver::LinearEquationSolver(kNdof, fNdof, solverChoice);
    Matrix<double> retVal(nDof, 1);
    for (int i = 0; i < nUnrestainedDof; i++)
        retVal(i, 0) = smallerResult(i, 0);
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

std::tuple<std::vector<Matrix<double>>, std::vector<Matrix<double>>, std::vector<Matrix<double>>>
StructureSolver::ImplicitNewmark(const Structure& str, std::vector<Matrix<double>> f, double tMax, double dT,
    double a0, double a1, SolverChoice solverChoice)
{
    // Create fields
    std::vector<Matrix<double>> displacements(f.size());
    std::vector<Matrix<double>> velocities(f.size());
    std::vector<Matrix<double>> accelerations(f.size());

    // Impose boundary conditions
    Matrix<double> initDisps(str.nDOF);
    Matrix<double> initVelos(str.nDOF);
    Matrix<double> initAccs(str.nDOF);

    displacements[0] = initDisps;
    velocities[0] = initVelos;
    accelerations[0] = initAccs;

    // Start solving    
    // Solve for unrestrained dofs. Add restrained dofs after that.
    auto nUDOF = str.nUnrestrainedDOF;
    auto kOr = str.StiffnessMatrix->getSubmatrix(0, nUDOF - 1, 0, nUDOF - 1);
    auto mOr = str.MassMatrix->getSubmatrix(0, nUDOF - 1, 0, nUDOF - 1);

    arma::mat k(nUDOF, nUDOF);
    arma::mat m(nUDOF, nUDOF);

    for (size_t i = 0; i < nUDOF; i++)
    {
        for (size_t j = 0; j < nUDOF; j++)
        {
            k(i, j) = kOr(i, j);
            m(i, j) = mOr(i, j);
        }
    }

    arma::mat c = (a0 * m) + (a1 * k);
    arma::mat mm = arma::inv(m + (c * dT * 0.5) + (k * dT * dT * 0.25));
    arma::mat f1 = (c * dT * 0.5) + (k * dT * dT * 0.25);
    arma::mat f2 = c + (k * dT);

    auto nStep = (int)(tMax / dT);

    for (size_t i = 1; i <= nStep; i++)
    {
        // Get current force vector for given time step
        auto&& fOrCurr = f[i - 1];

        // Get displacements, velocities and accelerations of previous step
        auto&& uOrPrev = displacements[i - 1];
        auto&& vOrPrev = velocities[i - 1];
        auto&& aOrPrev = accelerations[i - 1];

        arma::vec fArmaCurr(nUDOF);
        arma::vec uArmaPrev(nUDOF);
        arma::vec vArmaPrev(nUDOF);
        arma::vec aArmaPrev(nUDOF);

        for (size_t j = 0; j < nUDOF; j++)
        {
            fArmaCurr(j) = fOrCurr(j, 0);
            uArmaPrev(j) = uOrPrev(j, 0);
            vArmaPrev(j) = vOrPrev(j, 0);
            aArmaPrev(j) = aOrPrev(j, 0);
        }

        arma::vec fm = fArmaCurr - (f1 * aArmaPrev) - (f2 * vArmaPrev) - (k * uArmaPrev);

        arma::vec aArmaCurr = mm * fm;

        arma::vec vArmaCurr = vArmaPrev + (0.5 * dT * (aArmaPrev + aArmaCurr));

        arma::vec uArmaCurr = uArmaPrev + (vArmaPrev * dT) + (0.25 * dT * dT * (aArmaPrev + aArmaCurr));

        Matrix<double> uOrForward(str.nDOF, 1);
        Matrix<double> vOrForward(str.nDOF, 1);
        Matrix<double> aOrForward(str.nDOF, 1);

        for (size_t j = 0; j < nUDOF; j++)
        {
            uOrForward(j, 0) = uArmaCurr(j);
            vOrForward(j, 0) = vArmaCurr(j);
            aOrForward(j, 0) = aArmaCurr(j);
        }

        displacements[i] = uOrForward;
        velocities[i] = vOrForward;
        accelerations[i] = aOrForward;
    }

    return std::make_tuple(displacements, velocities, accelerations);
}

/// <summary>
/// THIS IS ONLY VALID FOR LATTICE STRUCTURES
/// </summary>
Matrix<double> StructureSolver::PerformPlasticPushoverForLatticeModel(Structure& str, Node& dispControlNode, double controlDisp, unsigned int controlDofIndex,
    Node& reactionControlNode, double dispIncrement, std::vector<bool> universalRestraintCondition, SolverChoice solverSelection, double modifier)
{
    // Create a new structure
    std::map<unsigned int, std::shared_ptr<Node>> newNodes;
    std::map<unsigned int, std::shared_ptr<Element>> newElements;
    std::map<unsigned int, std::shared_ptr<Restraint>> newRestraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> newNodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> newDistLoads;

    for (size_t nodeIndex = 1; nodeIndex <= str.Nodes->size(); nodeIndex++)
        newNodes[nodeIndex] = std::make_shared<Node>(nodeIndex, str.Nodes->at(nodeIndex)->Coordinate);

    for (size_t elIndex = 1; elIndex <= str.Elements->size(); elIndex++)
        newElements[elIndex] = std::make_shared<TrussMember>(
            elIndex, 
            newNodes.at(str.Elements->at(elIndex)->GelElementNodes().at(0)->NodeIndex),
            newNodes.at(str.Elements->at(elIndex)->GelElementNodes().at(1)->NodeIndex),
            (dynamic_cast<TrussMember*>(&*str.Elements->at(elIndex)))->TrussSection,
            (dynamic_cast<TrussMember*>(&*str.Elements->at(elIndex)))->TrussMaterial
            );

    auto restIterator = str.Restraints->begin();
    unsigned int controlRestraintIndex = 0;
    while (restIterator != str.Restraints->end())
    {
        auto index = restIterator->first;
        auto& originalRestraint = restIterator->second;
        auto& originalRestraintCondition = originalRestraint->IsRestrainedVector;
        auto& originalRestraintAmount = originalRestraint->RestrainedCondition;

        std::vector<bool> restraintCondition;
        std::vector<double> howMuchRestrained;

        for (auto&& cond : originalRestraintCondition)
            restraintCondition.push_back(cond);

        for (size_t i = 0; i < 6; i++)
            if (universalRestraintCondition[i])
                restraintCondition[i] = true;

        for (auto&& amount : originalRestraintAmount)
            howMuchRestrained.push_back(amount);

        if (originalRestraint->RestrainedNode->NodeIndex == dispControlNode.NodeIndex)
        {
            restraintCondition.at(controlDofIndex - 1) = true;
            controlRestraintIndex = index;
        }

        // auto res1 = std::make_shared<Restraint>(node1, isRest, rest); restraints[1] = res1;
        newRestraints[index] = std::make_shared<Restraint>(newNodes.at(originalRestraint->RestrainedNode->NodeIndex), restraintCondition, howMuchRestrained);


        restIterator++;
    }

    auto newStructure = std::make_shared<Structure>(&newNodes, &newElements, &newRestraints, &newNodalLoads, &newDistLoads);
    auto actualControlDofIndex = newStructure->Nodes->at(dispControlNode.NodeIndex)->DofIndexTX - 1 + controlDofIndex - 1;

    // Find number of increments
    auto numOfIncrements = static_cast<int>(controlDisp / dispIncrement);

    // Initialize return value including a data point for (0, 0)
    Matrix<double> retVal(numOfIncrements + 1, 2);

    auto updateRestraint = [&](double settlement) {

        auto& newRestraint = newStructure->Restraints->at(controlRestraintIndex);
        newRestraint->RestrainedCondition.at(controlDofIndex - 1) = true;
        switch (controlDofIndex)
        {
        case 1:
            newRestraint->TranslationX = settlement;
            break;
        case 2:
            newRestraint->TranslationY = settlement;
            break;
        case 3:
            newRestraint->TranslationZ = settlement;
            break;
        case 4:
            newRestraint->RotationX = settlement;
            break;
        case 5:
            newRestraint->RotationY = settlement;
            break;
        case 6:
            newRestraint->RotationZ = settlement;
            break;
        default:
            break;
        }
    };

    // Create a force vector which always be equal to 0
    Matrix<double> fVec(newStructure->nDOF, 1);

    for (size_t i = 1; i < newStructure->Elements->size(); i++)
    {
        auto t = dynamic_cast<TrussMember*>(&*newStructure->Elements->at(i));
        t->updateStiffness(modifier);
    }

    // At each increment, apply displacement and find reaction forces.
    for (size_t i = 1; i <= numOfIncrements; i++)
    {
        // Calculate control nodes displacement for current increment
        auto incrementalDisplacement = i * dispIncrement;
        updateRestraint(incrementalDisplacement);

        // Solve the system
        auto&& currDisplacements = StructureSolver::GetDisplacementsForStaticCaseWithForceVector(*newStructure, fVec, solverSelection);

        // Update structures stiffness matrix for nonlinear elements
        // Update the trusses stiffnesses values using the secant stiffness obtained from material model
        for (size_t elmIndex = 1; elmIndex <= newStructure->Elements->size(); elmIndex++)
        {
            // Get truss member
            auto trussMem = dynamic_cast<TrussMember*>(&*newStructure->Elements->at(elmIndex));

            if (trussMem)
            {
                // Calculate the strain of the member
                auto trussDeformation = trussMem->getTrussDeformation(currDisplacements);
                auto trussStrain = trussDeformation / trussMem->Length;

                // Get current secant modulus at the calculated strain
                auto currSecantModulus = trussMem->TrussMaterial->getSecantModulusAt(trussStrain);

                // Get ratio of current secant modulus of the truss to the last secant modulus of the truss in order to
                // update stiffness.
                auto stiffnessUpdateMultiplier = currSecantModulus / trussMem->elasticityModulusFromMaterialModel;

                // Update the trusses stiffness matrix
                trussMem->updateStiffness(stiffnessUpdateMultiplier);

                // Update elastic modulus for truss
                trussMem->elasticityModulusFromMaterialModel = currSecantModulus;
            }
        }

        // Read data for control nodes and save
        retVal(i, 0) = currDisplacements(newStructure->Nodes->at(dispControlNode.NodeIndex)->DofIndexTY - 1, 0);

        auto&& asd = (*newStructure->StiffnessMatrix) * currDisplacements;
        //auto&& asd = GetSupportReactions(*newStructure, currDisplacements, *newStructure->Restraints->at(controlRestraintIndex), solverSelection);

        auto react = 0.0;
        auto restraintIterator = newStructure->Restraints->begin();
        while (restraintIterator != newStructure->Restraints->end())
        {
            react += asd(restraintIterator->second->RestrainedNode->DofIndexTY - 1, 0);
            restraintIterator++;
        }

        retVal(i, 1) = react;

        // Update stiffness matrix of the structure
        newStructure->updateStiffnessMatrix();
    }

    return retVal;
}

// This will be moved to shell class
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

// This will be moved to shell class
Matrix<double> StructureSolver::CalculatePlateForces(const ShellMember& elm, Matrix<double>& disps)
{
    // Thickness
    auto thickness = elm.Thickness;

    // Map coordinates of flat plane to 2-D surface
    auto d1 = elm.Nodes[0]->Coordinate.DistanceTo(elm.Nodes[1]->Coordinate);
    auto d2 = elm.Nodes[1]->Coordinate.DistanceTo(elm.Nodes[2]->Coordinate);
    auto d3 = elm.Nodes[2]->Coordinate.DistanceTo(elm.Nodes[3]->Coordinate);
    auto d4 = elm.Nodes[3]->Coordinate.DistanceTo(elm.Nodes[0]->Coordinate);

    Vector p1V(elm.Nodes[0]->Coordinate);
    Vector p2V(elm.Nodes[1]->Coordinate);
    Vector p3V(elm.Nodes[2]->Coordinate);
    Vector p4V(elm.Nodes[3]->Coordinate);

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


    // Get displacement vector for element (assume element is at XY-plane)
    auto iNodeDisps = GetNodalDisplacements(*elm.Nodes[0], disps);
    auto jNodeDisps = GetNodalDisplacements(*elm.Nodes[1], disps);
    auto kNodeDisps = GetNodalDisplacements(*elm.Nodes[2], disps);
    auto lNodeDisps = GetNodalDisplacements(*elm.Nodes[3], disps);

    Matrix<double> dispVector(12, 1);
    dispVector(0, 0) = iNodeDisps(2, 0);
    dispVector(1, 0) = iNodeDisps(3, 0);
    dispVector(2, 0) = iNodeDisps(4, 0);
    dispVector(3, 0) = jNodeDisps(2, 0);
    dispVector(4, 0) = jNodeDisps(3, 0);
    dispVector(5, 0) = jNodeDisps(4, 0);
    dispVector(6, 0) = kNodeDisps(2, 0);
    dispVector(7, 0) = kNodeDisps(3, 0);
    dispVector(8, 0) = kNodeDisps(4, 0);
    dispVector(9, 0) = lNodeDisps(2, 0);
    dispVector(10, 0) = lNodeDisps(3, 0);
    dispVector(11, 0) = lNodeDisps(4, 0);

    // Bending forces
    auto gpCoeff2 = 1 / sqrt(3);
    double gaussPoints2[4][2] = { {-gpCoeff2, -gpCoeff2}, {gpCoeff2, -gpCoeff2}, {gpCoeff2, gpCoeff2}, {-gpCoeff2, gpCoeff2} };

    Matrix<double> flexuralRigidity(3, 3);
    auto elas = elm.ShellMaterial->E;
    auto pois = elm.ShellMaterial->PoissonsRatio;
    auto thick = elm.Thickness;
    auto fRMult = elas * (thick * thick * thick) / (12 * (1 - (pois * pois)));
    flexuralRigidity(0, 0) = fRMult * 1; flexuralRigidity(0, 1) = fRMult * pois;
    flexuralRigidity(1, 0) = fRMult * pois; flexuralRigidity(1, 1) = fRMult * 1;
    flexuralRigidity(2, 2) = fRMult * (1 - pois) / 2;
    Matrix<double> bendingMoments(3, 4);

    for (int colCounter = 0; colCounter < 4; colCounter++)
    {
        auto ksi = gaussPoints2[colCounter][0]; auto eta = gaussPoints2[colCounter][1];

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

        // Bilinear shape functions
        auto n1 = 0.25 * (1 - ksi) * (1 - eta);
        auto n2 = 0.25 * (1 + ksi) * (1 - eta);
        auto n3 = 0.25 * (1 + ksi) * (1 + eta);
        auto n4 = 0.25 * (1 - ksi) * (1 + eta);

        // Derivative of shape functions with respect to ksi
        auto dN1Ksi = -0.25 * (1 - eta);
        auto dN2Ksi = 0.25 * (1 - eta);
        auto dN3Ksi = 0.25 * (1 + eta);
        auto dN4Ksi = -0.25 * (1 + eta);

        // Derivative of shape functions with respect to eta
        auto dN1Eta = -0.25 * (1 - ksi);
        auto dN2Eta = -0.25 * (1 + ksi);
        auto dN3Eta = 0.25 * (1 + ksi);
        auto dN4Eta = 0.25 * (1 - ksi);

        // Derivative of shape functions with respect to x
        auto dN1X = (inversejacobi(0, 0) * dN1Ksi) + (inversejacobi(0, 1) * dN1Eta);
        auto dN2X = (inversejacobi(0, 0) * dN2Ksi) + (inversejacobi(0, 1) * dN2Eta);
        auto dN3X = (inversejacobi(0, 0) * dN3Ksi) + (inversejacobi(0, 1) * dN3Eta);
        auto dN4X = (inversejacobi(0, 0) * dN4Ksi) + (inversejacobi(0, 1) * dN4Eta);

        // Derivative of shape functions with respect to y
        auto dN1Y = (inversejacobi(1, 0) * dN1Ksi) + (inversejacobi(1, 1) * dN1Eta);
        auto dN2Y = (inversejacobi(1, 0) * dN2Ksi) + (inversejacobi(1, 1) * dN2Eta);
        auto dN3Y = (inversejacobi(1, 0) * dN3Ksi) + (inversejacobi(1, 1) * dN3Eta);
        auto dN4Y = (inversejacobi(1, 0) * dN4Ksi) + (inversejacobi(1, 1) * dN4Eta);

        Matrix<double> bB(3, 12);
        bB(0, 2) = dN1X; bB(0, 5) = dN2X; bB(0, 8) = dN3X; bB(0, 11) = dN4X;
        bB(1, 1) = -dN1Y; bB(1, 4) = -dN2Y; bB(1, 7) = -dN3Y; bB(1, 10) = -dN4Y;
        bB(2, 1) = -dN1X; bB(2, 4) = -dN2X; bB(2, 7) = -dN3X; bB(2, 10) = -dN4X;
        bB(2, 2) = dN1Y; bB(2, 5) = dN2Y; bB(2, 8) = dN3Y; bB(2, 11) = dN4Y;
        bB *= -1;

        auto momsForPt = flexuralRigidity * bB * dispVector;
        bendingMoments(0, colCounter) = momsForPt(0, 0);
        bendingMoments(1, colCounter) = momsForPt(1, 0);
        bendingMoments(2, colCounter) = momsForPt(2, 0);
    }

    // Extrapolate bending moments
    auto be4Coeff = sqrt(3);
    double be4[4][2] = { {-be4Coeff, -be4Coeff}, {be4Coeff, -be4Coeff}, {be4Coeff, be4Coeff}, {-be4Coeff, be4Coeff} };

    Matrix<double> modifiedBendingMoments(3, 4);
    for (size_t i = 0; i < 4; i++)
    {
        auto ksi = be4[i][0]; auto eta = be4[i][1];

        // Bilinear shape functions
        Matrix<double> N(4, 1);
        N(0, 0) = 0.25 * (1 - ksi) * (1 - eta);
        N(1, 0) = 0.25 * (1 + ksi) * (1 - eta);
        N(2, 0) = 0.25 * (1 + ksi) * (1 + eta);
        N(3, 0) = 0.25 * (1 - ksi) * (1 + eta);

        for (size_t i2 = 0; i2 < 3; i2++)
        {
            auto s = 0.0;

            for (size_t i3 = 0; i3 < 4; i3++)
                s += bendingMoments(i2, i3) * N(i3, 0);

            modifiedBendingMoments(i2, i) = s;
        }

    }

    bendingMoments = modifiedBendingMoments;

    Matrix<double> shearRigidity(2, 2);
    double sR = (5.0 / 6.0) * elm.ShellMaterial->G * elm.Thickness;
    shearRigidity(0, 0) = sR;
    shearRigidity(1, 1) = sR;
    Matrix<double> shearReactions(2, 1);
    for (size_t j = 0; j < 1; j++)
    {
        auto ksi = 0.0; auto eta = 0.0;

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

        // Bilinear shape functions
        auto n1 = 0.25 * (1 - ksi) * (1 - eta);
        auto n2 = 0.25 * (1 + ksi) * (1 - eta);
        auto n3 = 0.25 * (1 + ksi) * (1 + eta);
        auto n4 = 0.25 * (1 - ksi) * (1 + eta);

        // Derivative of shape functions with respect to ksi
        auto dN1Ksi = -0.25 * (1 - eta);
        auto dN2Ksi = 0.25 * (1 - eta);
        auto dN3Ksi = 0.25 * (1 + eta);
        auto dN4Ksi = -0.25 * (1 + eta);

        // Derivative of shape functions with respect to eta
        auto dN1Eta = -0.25 * (1 - ksi);
        auto dN2Eta = -0.25 * (1 + ksi);
        auto dN3Eta = 0.25 * (1 + ksi);
        auto dN4Eta = 0.25 * (1 - ksi);

        // Derivative of shape functions with respect to x
        auto dN1X = (inversejacobi(0, 0) * dN1Ksi) + (inversejacobi(0, 1) * dN1Eta);
        auto dN2X = (inversejacobi(0, 0) * dN2Ksi) + (inversejacobi(0, 1) * dN2Eta);
        auto dN3X = (inversejacobi(0, 0) * dN3Ksi) + (inversejacobi(0, 1) * dN3Eta);
        auto dN4X = (inversejacobi(0, 0) * dN4Ksi) + (inversejacobi(0, 1) * dN4Eta);

        // Derivative of shape functions with respect to y
        auto dN1Y = (inversejacobi(1, 0) * dN1Ksi) + (inversejacobi(1, 1) * dN1Eta);
        auto dN2Y = (inversejacobi(1, 0) * dN2Ksi) + (inversejacobi(1, 1) * dN2Eta);
        auto dN3Y = (inversejacobi(1, 0) * dN3Ksi) + (inversejacobi(1, 1) * dN3Eta);
        auto dN4Y = (inversejacobi(1, 0) * dN4Ksi) + (inversejacobi(1, 1) * dN4Eta);

        Matrix<double> bS(2, 12);
        bS(0, 0) = dN1X; bS(0, 2) = n1;
        bS(0, 3) = dN2X; bS(0, 5) = n2;
        bS(0, 6) = dN3X; bS(0, 8) = n3;
        bS(0, 9) = dN4X; bS(0, 11) = n4;

        bS(1, 0) = -dN1Y; bS(1, 1) = n1;
        bS(1, 3) = -dN2Y; bS(1, 4) = n2;
        bS(1, 6) = -dN3Y; bS(1, 7) = n3;
        bS(1, 9) = -dN4Y; bS(1, 10) = n4;

        shearReactions = shearRigidity * bS * dispVector;
    }

    Matrix<double> reactions(5, 4);

    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 4; j++)
            reactions(i, j) = bendingMoments(i, j);

    reactions(3, 0) = shearReactions(0, 0);
    reactions(4, 0) = shearReactions(1, 0);
    reactions(3, 1) = shearReactions(0, 0);
    reactions(4, 1) = shearReactions(1, 0);
    reactions(3, 2) = shearReactions(0, 0);
    reactions(4, 3) = shearReactions(1, 0);
    reactions(3, 3) = shearReactions(0, 0);
    reactions(4, 3) = shearReactions(1, 0);

    return reactions;
}

// This will be moved to serendipity shell class
Matrix<double> StructureSolver::CalculateSerendipityPlateForces(const SerendipityShell& elm, Matrix<double>& disps)
{
    Matrix<double> reactions(5, 4);
    Matrix<double> bendingMoments(3, 4);
    Matrix<double> shearReactions(2, 1);

    // Map coordinates of flat plate to 2-D surface
    auto d1 = elm.Nodes[0]->Coordinate.DistanceTo(elm.Nodes[1]->Coordinate);
    auto d2 = elm.Nodes[1]->Coordinate.DistanceTo(elm.Nodes[2]->Coordinate);
    auto d3 = elm.Nodes[2]->Coordinate.DistanceTo(elm.Nodes[3]->Coordinate);
    auto d4 = elm.Nodes[3]->Coordinate.DistanceTo(elm.Nodes[0]->Coordinate);

    Vector p1V(elm.Nodes[0]->Coordinate);
    Vector p2V(elm.Nodes[1]->Coordinate);
    Vector p3V(elm.Nodes[2]->Coordinate);
    Vector p4V(elm.Nodes[3]->Coordinate);

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

    // Get displacement vector for element (assume element is at XY-plane)
    auto iNodeDisps = GetNodalDisplacements(*elm.Nodes[0], disps);
    auto jNodeDisps = GetNodalDisplacements(*elm.Nodes[1], disps);
    auto kNodeDisps = GetNodalDisplacements(*elm.Nodes[2], disps);
    auto lNodeDisps = GetNodalDisplacements(*elm.Nodes[3], disps);
    auto ijNodeDisps = GetNodalDisplacements(*elm.Nodes[4], disps);
    auto jkNodeDisps = GetNodalDisplacements(*elm.Nodes[5], disps);
    auto klNodeDisps = GetNodalDisplacements(*elm.Nodes[6], disps);
    auto liNodeDisps = GetNodalDisplacements(*elm.Nodes[7], disps);

    Matrix<double> dispVector(24, 1);
    dispVector(0, 0) = iNodeDisps(2, 0);
    dispVector(1, 0) = iNodeDisps(3, 0);
    dispVector(2, 0) = iNodeDisps(4, 0);

    dispVector(3, 0) = jNodeDisps(2, 0);
    dispVector(4, 0) = jNodeDisps(3, 0);
    dispVector(5, 0) = jNodeDisps(4, 0);

    dispVector(6, 0) = kNodeDisps(2, 0);
    dispVector(7, 0) = kNodeDisps(3, 0);
    dispVector(8, 0) = kNodeDisps(4, 0);

    dispVector(9, 0) = lNodeDisps(2, 0);
    dispVector(10, 0) = lNodeDisps(3, 0);
    dispVector(11, 0) = lNodeDisps(4, 0);

    dispVector(12, 0) = ijNodeDisps(2, 0);
    dispVector(13, 0) = ijNodeDisps(3, 0);
    dispVector(14, 0) = ijNodeDisps(4, 0);

    dispVector(15, 0) = jkNodeDisps(2, 0);
    dispVector(16, 0) = jkNodeDisps(3, 0);
    dispVector(17, 0) = jkNodeDisps(4, 0);

    dispVector(18, 0) = klNodeDisps(2, 0);
    dispVector(19, 0) = klNodeDisps(3, 0);
    dispVector(20, 0) = klNodeDisps(4, 0);

    dispVector(21, 0) = liNodeDisps(2, 0);
    dispVector(22, 0) = liNodeDisps(3, 0);
    dispVector(23, 0) = liNodeDisps(4, 0);

    // Constitutive matrices
    Matrix<double> flexuralRigidity(3, 3);
    auto elas = elm.SerendipityMaterial->E;
    auto pois = elm.SerendipityMaterial->PoissonsRatio;
    auto thick = elm.Thickness;
    auto fRMult = elas * (thick * thick * thick) / (12 * (1 - (pois * pois)));
    flexuralRigidity(0, 0) = fRMult * 1; flexuralRigidity(0, 1) = fRMult * pois;
    flexuralRigidity(1, 0) = fRMult * pois; flexuralRigidity(1, 1) = fRMult * 1;
    flexuralRigidity(2, 2) = fRMult * (1 - pois) / 2;
    Matrix<double> shearRigidity(2, 2);
    double sR = (5.0 / 6.0) * elm.SerendipityMaterial->G * thick;
    shearRigidity(0, 0) = sR; shearRigidity(1, 1) = sR;

    // To reduce number of mechanism, calculate kBending by using 3x3 Gauss points and
    // calculate kShear by using 2x2 Gauss points. (Gauss points are ordered in CCW orientation)
    auto gp2 = 1.0;
    double gp2Points[4][2] = { {-gp2, -gp2}, {gp2, -gp2}, {gp2, gp2}, {-gp2, gp2} };

    // Calculate kBending for each Gauss Points
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

        // Create curvature-displacement relation
        Matrix<double> bB(3, 24);
        bB(0, 2) = dN1X; bB(0, 5) = dN2X; bB(0, 8) = dN3X; bB(0, 11) = dN4X; bB(0, 14) = dN5X; bB(0, 17) = dN6X; bB(0, 20) = dN7X; bB(0, 23) = dN8X;
        bB(1, 1) = -dN1Y; bB(1, 4) = -dN2Y; bB(1, 7) = -dN3Y; bB(1, 10) = -dN4Y; bB(1, 13) = -dN5Y; bB(1, 16) = -dN6Y; bB(1, 19) = -dN7Y; bB(1, 22) = -dN8Y;
        bB(2, 1) = -dN1X; bB(2, 4) = -dN2X; bB(2, 7) = -dN3X; bB(2, 10) = -dN4X; bB(2, 13) = -dN5X; bB(2, 16) = -dN6X; bB(2, 19) = -dN7X; bB(2, 22) = -dN8X;
        bB(2, 2) = dN1Y; bB(2, 5) = dN2Y; bB(2, 8) = dN3Y; bB(2, 11) = dN4Y; bB(2, 14) = dN5Y; bB(2, 17) = dN6Y; bB(2, 20) = dN7Y; bB(2, 23) = dN8Y;
        bB *= -1;

        auto momsForPt = flexuralRigidity * bB * dispVector;
        bendingMoments(0, i) = momsForPt(0, 0);
        bendingMoments(1, i) = momsForPt(1, 0);
        bendingMoments(2, i) = momsForPt(2, 0);
    }

    auto gpShear = 1 / sqrt(3);
    double gpShearPoints[4][2] = { {-gpShear, -gpShear}, {gpShear, -gpShear}, {gpShear, gpShear}, {-gpShear, gpShear} };
    // Calculate kShear for each Gauss Points
    for (size_t i = 0; i < 4; i++)
    {
        auto ksi = gpShearPoints[i][0]; auto eta = gpShearPoints[i][1];

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

        shearReactions = shearRigidity * bS * dispVector;
    }

    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 4; j++)
            reactions(i, j) = bendingMoments(i, j);

    reactions(3, 0) = shearReactions(0, 0);
    reactions(4, 0) = shearReactions(1, 0);
    reactions(3, 1) = shearReactions(0, 0);
    reactions(4, 1) = shearReactions(1, 0);
    reactions(3, 2) = shearReactions(0, 0);
    reactions(4, 3) = shearReactions(1, 0);
    reactions(3, 3) = shearReactions(0, 0);
    reactions(4, 3) = shearReactions(1, 0);

    return reactions;
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

        // auto format = [](double d) {return abs(d) < 1 ? 0 : d; };
        // 
        // auto asd = arma::eig_gen(aArmaMat);
        // for (size_t i = 0; i < asd.size(); i++)
        //     std::cout << " Eigenvalue #" << i + 1 << ": " << format(asd(i).real()) << "\n";

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

double StructureSolver::CalculateInternalEnergy(const Structure& str, Matrix<double>& disps, SolverChoice)
{
    auto strK = str.StiffnessMatrix->getSubmatrix(0, str.nUnrestrainedDOF - 1, 0, str.nUnrestrainedDOF - 1);
    auto disp = disps.getSubmatrix(0, str.nUnrestrainedDOF - 1, 0, 0);
    auto tdisp = disp.transpose();
    auto fir = tdisp * strK;
    auto intEner = fir * disp;
    return intEner(0, 0);
}
