#pragma once
#include <memory>
#include "Structure.h"
#include "Matrix.h"

enum class SolverChoice
{
    None = 0,
    Eigen = 1,
    Armadillo = 2
};

namespace StructureSolver
{
    Matrix<double> GetDisplacementForStaticCase(const Structure& str, SolverChoice solverChoice);
    Matrix<double> GetMemberEndForcesForLocalCoordinates(Element& elm, Matrix<double>& displacements);
    Matrix<double> GetMemberEndForcesForGlobalCoordinates(Element& elm, Matrix<double>& displacements);
    Matrix<double> GetNodalDisplacements(Node& node, Matrix<double>& displacements);
    Matrix<double> GetSupportReactions(const Structure& str, Matrix<double>& disps, const Restraint& res, SolverChoice solverChoice);
    Matrix<double> GetModalPeriods(const Structure& str, SolverChoice solverChoice);
    Matrix<double> CalculateMembraneNodalStresses(const ShellMember& elm, Matrix<double>& disps, int nodeIndex);
    Matrix<double> CalculatePlateForces(const ShellMember& elm, Matrix<double>& disps);
    double CondenseStiffnessMatrixForSpecificDOF(const Structure& str, unsigned int dofIndex, SolverChoice solverChoice);
    double CondenseForceVectorForSpecificDOF(const Structure& str, unsigned int dofIndex, SolverChoice solverChoice);
    Matrix<double> LinearEquationSolver(Matrix<double>& A, Matrix<double>& b, SolverChoice solverChoice);
    Matrix<double> GetInverse(Matrix<double>& A, SolverChoice solverChoice);
}
