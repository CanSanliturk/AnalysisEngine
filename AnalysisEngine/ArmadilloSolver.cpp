#include <armadillo>
#include <iostream>
#include <memory>
#include <map>
#include "ArmadilloSolver.h"
#include "GeometryHelper.h"
#include "Matrix.h"

constexpr double pi = 3.141592653589793;

std::vector<double> ArmadilloSolver::GetDisplacementForStaticCase(const Structure& str)
{
    // Create armadillo matrices
    arma::mat k(str.nUnrestrainedDOF, str.nUnrestrainedDOF);
    arma::vec f(str.nUnrestrainedDOF);

    auto nDofUnrestrained = str.nUnrestrainedDOF;
    auto nDof = str.nDOF;
    auto nDofRestrained = nDof - nDofUnrestrained;

    k.zeros(nDofUnrestrained, nDofUnrestrained);
    f.zeros(nDofUnrestrained);

    // Fill zeros armadillo matrices
    for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
    {
        f[i] = 0.0;
        for (size_t j = 0; j < str.nUnrestrainedDOF; j++)
        {
            k[(i * str.nUnrestrainedDOF) + j] = 0.0;
        }
    }

    // Fill armadillo matrices
    for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
    {
        f(i) = (*str.ForceVector)(i, 0);
        for (size_t j = 0; j < str.nUnrestrainedDOF; j++)
        {
            k(i, j) = (*str.StiffnessMatrix)(i, j);
        }
    }

    // Subtract support settlements from force vector
    arma::mat kUk(nDofUnrestrained, nDofRestrained);
    arma::vec uK(nDofRestrained);

    kUk.zeros(nDofUnrestrained, nDofRestrained);
    uK.zeros(nDofRestrained);

    for (size_t i = 0; i < nDofUnrestrained; i++)
        for (size_t j = nDofUnrestrained; j < nDof; j++)
            kUk(i, j - nDofUnrestrained) = (*str.StiffnessMatrix)(i, j);

    auto restraints = *str.Restraints;

    std::map<unsigned int, std::shared_ptr<Restraint>>::iterator iter = restraints.begin();

    // Iterate over the map using Iterator till end.
    while (iter != restraints.end())
    {
        auto restraint = iter->second;
        auto restrainedNode = restraint->RestrainedNode;

        if (restraint->IsRestraintTranslationX)
            uK(restrainedNode->DofIndexTX - nDofUnrestrained - 1) = restraint->TranslationX;
        if (restraint->IsRestraintTranslationY)
            uK(restrainedNode->DofIndexTY - nDofUnrestrained - 1) = restraint->TranslationY;
        if (restraint->IsRestraintTranslationZ)
            uK(restrainedNode->DofIndexTZ - nDofUnrestrained - 1) = restraint->TranslationZ;
        if (restraint->IsRestraintRotationX)
            uK(restrainedNode->DofIndexRX - nDofUnrestrained - 1) = restraint->RotationX;
        if (restraint->IsRestraintRotationY)
            uK(restrainedNode->DofIndexRY - nDofUnrestrained - 1) = restraint->RotationY;
        if (restraint->IsRestraintRotationZ)
            uK(restrainedNode->DofIndexRZ - nDofUnrestrained - 1) = restraint->RotationZ;
        iter++;
    }


    auto subtractVal = kUk * uK;
    auto condFVec = f - subtractVal;

    // Solve system
    arma::vec resData = arma::solve(k, condFVec);

    // Store data to return value
    std::vector<double> retVal(str.nDOF);

    for (size_t i = 0; i < str.nDOF; i++)
        retVal.at(i) = 0;

    for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
        retVal.at(i) = resData(i);

    std::map<unsigned int, std::shared_ptr<Restraint>>::iterator it = restraints.begin();

    // Iterate over the map using Iterator till end.
    while (it != restraints.end())
    {
        auto restraint = it->second;
        auto restrainedNode = restraint->RestrainedNode;

        if (restraint->IsRestraintTranslationX)
            retVal.at(restrainedNode->DofIndexTX - 1) = restraint->TranslationX;
        if (restraint->IsRestraintTranslationY)
            retVal.at(restrainedNode->DofIndexTY - 1) = restraint->TranslationY;
        if (restraint->IsRestraintTranslationZ)
            retVal.at(restrainedNode->DofIndexTZ - 1) = restraint->TranslationZ;
        if (restraint->IsRestraintRotationX)
            retVal.at(restrainedNode->DofIndexRX - 1) = restraint->RotationX;
        if (restraint->IsRestraintRotationY)
            retVal.at(restrainedNode->DofIndexRY - 1) = restraint->RotationY;
        if (restraint->IsRestraintRotationZ)
            retVal.at(restrainedNode->DofIndexRZ - 1) = restraint->RotationZ;

        it++;
    }

    return retVal;
}

std::vector<double> ArmadilloSolver::GetMemberEndForcesForLocalCoordinates(Element& elm, std::vector<double> displacements)
{
    // To get elements end forces in local coordinates, displacement vector of element (which is in global coordinates) 
    // should be converted to local coordinates by multiplying it by rotation matrix of the given element.

    // Create arma matrix for rotation matrix
    auto rotMatrix = elm.GetRotationMatrix();
    auto nDofElm = elm.GetNumberOfDoF();

    arma::mat rotMat(nDofElm, nDofElm);
    for (size_t i = 0; i < nDofElm; i++)
        for (size_t j = 0; j < nDofElm; j++)
            rotMat(i, j) = (*rotMatrix)(i, j);

    // Retrieve displacements of element end nodes
    auto elmNodes = elm.GelElementNodes();
    arma::vec disps(nDofElm);

    unsigned short counter = 0;

    for (auto node : elmNodes)
    {
        disps((counter * 6) + 0) = displacements.at(node->DofIndexTX - 1);
        disps((counter * 6) + 1) = displacements.at(node->DofIndexTY - 1);
        disps((counter * 6) + 2) = displacements.at(node->DofIndexTZ - 1);
        disps((counter * 6) + 3) = displacements.at(node->DofIndexRX - 1);
        disps((counter * 6) + 4) = displacements.at(node->DofIndexRY - 1);
        disps((counter * 6) + 5) = displacements.at(node->DofIndexRZ - 1);
        counter++;
    }

    // Get localized displacements
    arma::vec localizedDisps = rotMat * disps;

    // Get local stiffness matrix
    auto localKMat = elm.GetLocalCoordinateStiffnessMatrix();
    arma::mat kMat(nDofElm, nDofElm);

    for (size_t i = 0; i < nDofElm; i++)
        for (size_t j = 0; j < nDofElm; j++)
            kMat(i, j) = (*localKMat)(i, j);

    // Multiply local stiffness matrix and localized displacement vector to obtain element end forces
    arma::vec localForces = kMat * localizedDisps;

    // Convert arma::vec to std::vector
    std::vector<double> retVal(nDofElm);

    for (size_t i = 0; i < nDofElm; i++)
        retVal.at(i) = localForces(i);

    return retVal;
}

std::vector<double> ArmadilloSolver::GetMemberEndForcesForGlobalCoordinates(Element& elm, std::vector<double> displacements)
{
    // Multiplying elements stiffness matrix at global coordinates and elements nodes displacements will returns 
    // member end forces in global coordinates

    // Retrieve displacements of element end nodes
    auto nDofElm = elm.GetNumberOfDoF();
    auto elmNodes = elm.GelElementNodes();
    arma::vec disps(nDofElm);
    unsigned short counter = 0;

    for (auto nodePtr : elmNodes)
    {
        disps((counter * 6) + 0) = displacements.at(nodePtr->DofIndexTX - 1);
        disps((counter * 6) + 1) = displacements.at(nodePtr->DofIndexTY - 1);
        disps((counter * 6) + 2) = displacements.at(nodePtr->DofIndexTZ - 1);
        disps((counter * 6) + 3) = displacements.at(nodePtr->DofIndexRX - 1);
        disps((counter * 6) + 4) = displacements.at(nodePtr->DofIndexRY - 1);
        disps((counter * 6) + 5) = displacements.at(nodePtr->DofIndexRZ - 1);
        counter++;
    }

    // Get global stiffness matrix
    auto globalKMat = elm.GetLocalCoordinateStiffnessMatrix();
    arma::mat kMat(nDofElm, nDofElm);

    for (size_t i = 0; i < nDofElm; i++)
        for (size_t j = 0; j < nDofElm; j++)
            kMat(i, j) = (*globalKMat)(i, j);

    // Multiply global stiffness matrix and displacement vector to obtain element end forces
    arma::vec globalForces = kMat * disps;

    // Convert arma::vec to std::vector
    std::vector<double> retVal(nDofElm);

    for (size_t i = 0; i < nDofElm; i++)
        retVal.at(i) = globalForces(i);

    return retVal;
}

std::vector<double> ArmadilloSolver::GetNodalDisplacements(Node& node, std::vector<double>& displacements)
{
    std::vector<double> nodalDips(6);

    for (size_t i = 0; i < 6; i++)
        nodalDips[i] = 0;

    nodalDips[0] = displacements[node.DofIndexTX - 1];
    nodalDips[1] = displacements[node.DofIndexTY - 1];
    nodalDips[2] = displacements[node.DofIndexTZ - 1];
    nodalDips[3] = displacements[node.DofIndexRX - 1];
    nodalDips[4] = displacements[node.DofIndexRY - 1];
    nodalDips[5] = displacements[node.DofIndexRZ - 1];

    return nodalDips;
}

std::vector<double> ArmadilloSolver::GetSupportReactions(const Structure& str, const std::vector<double> disps, const Restraint& res)
{
    // Multiply full stiffness matrix with full displacement vector. Elements of resulting vector with indices starting from 
    // number of unrestrained dofs are all support reactions in global coordinates

    // Convert stiffness matrix to arma matrix
    // Create armadillo matrices
    arma::mat k(str.nDOF, str.nDOF);
    arma::vec d(str.nDOF);

    k.zeros(str.nDOF, str.nDOF);
    d.zeros(str.nDOF);

    // Fill zeros armadillo matrices
    for (size_t i = 0; i < str.nDOF; i++)
    {
        d[i] = 0.0;
        for (size_t j = 0; j < str.nDOF; j++)
        {
            k[(i * str.nDOF) + j] = 0.0;
        }
    }

    // Fill armadillo matrices
    for (size_t i = 0; i < str.nDOF; i++)
    {
        d(i) = disps[i];
        for (size_t j = 0; j < str.nDOF; j++)
        {
            k(i, j) = (*str.StiffnessMatrix)(i, j);
        }
    }

    // Multiply
    arma::vec forces(str.nDOF);
    forces.zeros();
    forces = k * d;

    // Reactions are obtained above. Now, get the given restraints reactions
    std::vector<double> retVal;
    auto restNode = res.RestrainedNode;

    if (res.IsRestraintTranslationX) retVal.push_back(forces(restNode->DofIndexTX - 1));
    if (res.IsRestraintTranslationY) retVal.push_back(forces(restNode->DofIndexTY - 1));
    if (res.IsRestraintTranslationZ) retVal.push_back(forces(restNode->DofIndexTZ - 1));
    if (res.IsRestraintRotationX) retVal.push_back(forces(restNode->DofIndexRX - 1));
    if (res.IsRestraintRotationY) retVal.push_back(forces(restNode->DofIndexRY - 1));
    if (res.IsRestraintRotationZ) retVal.push_back(forces(restNode->DofIndexRZ - 1));

    return retVal;
}

std::vector<double> ArmadilloSolver::GetModalPeriods(const Structure& str)
{
    // Return value
    std::vector<double> periods;

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

    auto& firstMode = eigVec[0];

    for (int i = eigVal.size() - 1; -1 < i; i--)
        periods.push_back(2 * pi / sqrt(*(eigVal.at(i)._Val)));

    std::sort(periods.begin(), periods.end(), [&](double a, double b) {return b < a; });

    return periods;
}