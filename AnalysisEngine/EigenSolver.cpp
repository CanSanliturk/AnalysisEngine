#include <iostream>
#include <memory>
#include <Eigen>
#include <map>
#include "EigenSolver.h"
#include "MatrixHelper.h"

constexpr auto pi = 3.141592653589793;

std::vector<double> EigenSolver::GetDisplacementForStaticCase(const Structure& str)
{
    std::vector<double> retVal;

    return retVal;
}

std::vector<double> EigenSolver::GetMemberEndForcesForLocalCoordinates(Element& elm, std::vector<double> displacements)
{
    std::vector<double> retVal;

    return retVal;
}

std::vector<double> EigenSolver::GetMemberEndForcesForGlobalCoordinates(Element& elm, std::vector<double> displacements)
{
    std::vector<double> retVal;

    return retVal;
}

std::vector<double> EigenSolver::GetNodalDisplacements(Node& node, std::vector<double>& displacements)
{
    std::vector<double> retVal;

    return retVal;
}

std::vector<double> EigenSolver::GetSupportReactions(const Structure& str, const std::vector<double> disps, const Restraint& res)
{
    std::vector<double> retVal;

    return retVal;
}

std::vector<double> EigenSolver::GetModalPeriods(const Structure& str)
{
    std::vector<double> retVal;

    return retVal;
}