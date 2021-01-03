#include "ArmadilloSolver.h"
#include <armadillo>
#include <iostream>

std::vector<double> ArmadilloSolver::GetResult(Structure str)
{
	arma::mat k(str.nUnrestrainedDOF, str.nUnrestrainedDOF);
	arma::vec f(str.nUnrestrainedDOF);

	k.fill(0.0);
	f.fill(0.0);

	for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
	{
		f(i) = str.ForceVector.at(i);
		for (size_t j = 0; j < str.nUnrestrainedDOF; j++)
		{
			k(i, j) = str.StiffnessMatrix.at(i).at(j);
		}
	}


	std::cout << "Stiffness matrix for active DOFs \n";

	for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
	{
		for (size_t j = 0; j < str.nUnrestrainedDOF; j++)
		{
			std::cout << str.StiffnessMatrix.at(i).at(j) << " ";
		}

		std::cout << "\n";
	}

	std::cout << "\n Mass matrix for active DOFs \n";

	for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
	{
		for (size_t j = 0; j < str.nUnrestrainedDOF; j++)
		{
			std::cout << str.MassMatrix.at(i).at(j) << " ";
		}
		std::cout << "\n";
	}

	std::cout << "\n Force vector for active DOFs \n";

	for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
	{
		std::cout << str.ForceVector.at(i) << "\n";
	}

	arma::vec resData = arma::solve(k, f);
	std::vector<double> retVal(str.nUnrestrainedDOF);

	for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
		retVal.at(i) = 0;

	for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
		retVal.at(i) = resData(i);

	std::cout << "\n Displacement vector for active DOFs \n";

	for (size_t i = 0; i < str.nUnrestrainedDOF; i++)
	{
		std::cout << retVal.at(i) << "\n";
	}

	return retVal;
}
