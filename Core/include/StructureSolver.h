#pragma once
#include <memory>
#include <tuple>
#include <vector>
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
    Matrix<double> GetDisplacementForStaticCase(Structure& str, SolverChoice solverChoice);
    Matrix<double> GetDisplacementsForStaticCaseWithForceVector(Structure& str, Matrix<double>& fVec, SolverChoice solverChoice);
    Matrix<double> GetDisplacementsWithStaticNodalLoads(Structure& str, std::vector<NodalLoad> nodalLoads, SolverChoice solverChoice);
    Matrix<double> CalculateDisplacements(Matrix<double>& kMat, Matrix<double>& fVec, int nDof, int nUnrestainedDof, SolverChoice solverChoice);
    Matrix<double> GetMemberEndForcesForLocalCoordinates(Element& elm, Matrix<double>& displacements);
    Matrix<double> GetMemberEndForcesForGlobalCoordinates(Element& elm, Matrix<double>& displacements);
    Matrix<double> GetNodalDisplacements(Node& node, Matrix<double>& displacements);
    Matrix<double> GetSupportReactions(Structure& str, Matrix<double>& disps, Restraint& res, SolverChoice solverChoice);
    Matrix<double> GetModalPeriods(Structure& str, SolverChoice solverChoice);

    /// <summary>
    /// Executes impilict newmark integration according to given force vector
    /// </summary>
    /// <param name="str">Structure</param>
    /// <param name="f">Force vectors for ALL time steps</param>
    /// <param name="tMax">Maximum time</param>
    /// <param name="dT">Time step</param>
    /// <returns>
    /// Returns a tuple for displacements, velocities and accelerations
    /// for all time steps
    /// </returns>
    std::tuple<std::vector<Matrix<double>>, std::vector<Matrix<double>>, std::vector<Matrix<double>>>
        ImplicitNewmark(Structure& str, std::vector<Matrix<double>> f, double tMax, double dT, 
            double a0, double a1, SolverChoice solverChoice);



    Matrix<double> PerformPlasticPushoverForLatticeModel(Structure& str, Node& dispControlNode, double controlDisp, unsigned int controlDofIndex,
        Node& reactionControlNode, double dispIncrements, std::vector<bool> universalRestraintCondition, SolverChoice solverSelection, double modifier);


    Matrix<double> CalculateMembraneNodalStresses(ShellMember& elm, Matrix<double>& disps, int nodeIndex);
    Matrix<double> CalculatePlateForces(ShellMember& elm, Matrix<double>& disps);
    Matrix<double> CalculateSerendipityPlateForces(SerendipityShell& elm, Matrix<double>& disps);
    double CondenseStiffnessMatrixForSpecificDOF(Structure& str, unsigned int dofIndex, SolverChoice solverChoice);
    double CondenseForceVectorForSpecificDOF(Structure& str, unsigned int dofIndex, SolverChoice solverChoice);
    Matrix<double> LinearEquationSolver(Matrix<double>& A, Matrix<double>& b, SolverChoice solverChoice);
    Matrix<double> GetInverse(Matrix<double>& A, SolverChoice solverChoice);

    double CalculateInternalEnergy(Structure& str, Matrix<double>& disps, SolverChoice);
}
