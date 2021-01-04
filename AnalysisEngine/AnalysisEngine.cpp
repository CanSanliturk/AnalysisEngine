#include <iostream>
#include <vector>
#include "XYZPoint.h"
#include "Restraint.h"
#include "Node.h"
#include "Section.h"
#include "FrameMember.h"
#include "Element.h"
#include "MatrixHelper.h"
#include "Structure.h"
#include "ArmadilloSolver.h"

int main()
{
	// Define two nodes
	XYZPoint pt1(0.0, 0.0, 0.0);
	XYZPoint pt2(1.0, 0.0, 0.0);

	Node node1(1, pt1);
	Node node2(2, pt2);
	std::map<unsigned int, Node*> nodes;
	nodes[node1.NodeIndex] = &node1;
	nodes[node2.NodeIndex] = &node2;

	// Assign boundary condition to the first node
	// Fixed for all directions

	std::vector<bool> fix;
	std::vector<double> fixity;

	for (unsigned int i = 0; i < 6; i++)
	{
		fix.push_back(true);
		fixity.push_back(0.0);
	}

	Restraint res(&node1, fix, fixity);
	std::map<unsigned int, Restraint*> restraints;
	restraints[1] = &res;

	// Create section
	auto area = 1.0;
	auto inertia11 = 1.0 / 12.0;
	auto inertia22 = inertia11;
	auto inertia12 = 0.14083333333;

	Section sect(area, inertia11, inertia22, inertia12);

	// Create material
	auto e = 1.0;
	auto v = 0.3;
	auto rho = 1.0;
	Material mat(e, v, rho);

	// Create member
	FrameMember member(1, &node1, &node2, &sect, &mat);
	std::map<unsigned int, Element*> members;
	members[member.ElementIndex] = &member;

	// Add tip load
	double loads[6];
	loads[0] = 0;
	loads[1] = 0;
	loads[2] = -1;
	loads[3] = 0;
	loads[4] = 0;
	loads[5] = 0;

	NodalLoad nL(&node2, loads);
	std::map<unsigned int, NodalLoad*> nload;
	nload[1] = &nL;

	std::map<unsigned int, DistributedLoad*> dist;

	Structure str(&nodes, &members, &restraints, &nload, &dist);

	std::cout << "Model is created successfully" << std::endl;
	std::cout << "Solver starts" << std::endl;

	auto disps = ArmadilloSolver::GetDisplacementForStaticCase(str);
	ArmadilloSolver::GetMemberEndForcesForLocalCoordinates(member, disps);
	ArmadilloSolver::GetMemberEndForcesForGlobalCoordinates(member, disps);

	std::cin.get();

	return 0;
}