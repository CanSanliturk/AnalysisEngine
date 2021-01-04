#pragma once
#include <vector>
#include "Structure.h"
#include "Element.h"

namespace ArmadilloSolver
{
	std::vector<double> GetDisplacementForStaticCase(const Structure& str);
	std::vector<double> GetMemberEndForcesForLocalCoordinates(Element elm, std::vector<double> displacements);
	std::vector<double> GetMemberEndForcesForGlobalCoordinates(Element elm, std::vector<double> displacements);
}
