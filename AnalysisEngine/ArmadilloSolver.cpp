#include "ArmadilloSolver.h"
#include <armadillo>
#include <iostream>
#include <map>

std::vector<double> ArmadilloSolver::GetDisplacementForStaticCase(const Structure& str)
{
	// Create armadillo matrices
	arma::mat k(str.nUnrestrainedDOF, str.nUnrestrainedDOF);
	arma::vec f(str.nUnrestrainedDOF);

	k.fill(0.0);
	f.fill(0.0);

	// Fill armadillo matrices
	for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
	{
		f(i) = str.ForceVector.at(i);
		for (size_t j = 0; j < str.nUnrestrainedDOF; j++)
		{
			k(i, j) = str.StiffnessMatrix.at(i).at(j);
		}
	}

	// Solve system
	arma::vec resData = arma::solve(k, f);

	// Store data to return value
	std::vector<double> retVal(str.nDOF);
	
	for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
		retVal.at(i) = resData(i);
	for (size_t i = str.nUnrestrainedDOF; i < str.nDOF; i++)
		retVal.at(i) = 0;

	
	auto restraints = *str.Restraints;

	std::map<unsigned int, Restraint*>::iterator it = restraints.begin();

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

std::vector<double> ArmadilloSolver::GetMemberEndForcesForLocalCoordinates(Element elm, std::vector<double> displacements)
{
	// To get elements end forces  


	return std::vector<double>();
}

std::vector<double> ArmadilloSolver::GetMemberEndForcesForGlobalCoordinates(Element elm, std::vector<double> displacements)
{
	return std::vector<double>();
}
