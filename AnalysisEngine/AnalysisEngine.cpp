#include <iostream>
#include <vector>
#include <Windows.h>
#include "XYZPoint.h"
#include "Restraint.h"
#include "Node.h"
#include "Section.h"
#include "FrameMember.h"
#include "Element.h"
#include "MatrixHelper.h"
#include "Structure.h"
#include "ArmadilloSolver.h"

#pragma comment(lib, "user32")

#define LOG(x) std::cout << x << "\n"

int main()
{
	// Full screen
	//::SendMessage(::GetConsoleWindow(), WM_SYSKEYDOWN, VK_RETURN, 0x20000000);
	
	// Create points for nodes
	XYZPoint coord1(0.0, 0.0, 0.0);
	XYZPoint coord2(0.0, 0.0, 10.0);
	XYZPoint coord3(10.0, 0.0, 10.0);

	// Define nodes
	Node node1(1, coord1); 
	Node node2(2, coord2);
	Node node3(3, coord3);

	// Add nodes to 
	std::map<unsigned int, Node*> nodes;
	nodes[node1.NodeIndex] = &node1;
	nodes[node2.NodeIndex] = &node2;
	nodes[node3.NodeIndex] = &node3;

	// Assign boundary condition to the first node"
	// and fourth node
	// Fixed for all directions
	std::vector<bool> fix;
	std::vector<double> fixity;

	for (size_t i = 0; i < 6; i++)
	{
		fix.push_back(true);
		fixity.push_back(0.0);
	}

	Restraint restrainNode1(&node1, fix, fixity);

	std::map<unsigned int, Restraint*> restraints;
	restraints[1] = &restrainNode1;

	// Create section
	auto area = 0.16;
	auto inertia11 = 2.133 * 0.001;
	auto inertia22 = 2.133 * 0.001;
	auto inertia12 = 0.0036;
	Section frameSection(area, inertia11, inertia22, inertia12);

	// Create material
	auto e = 200e9;
	auto v = 0.3;
	auto rho = 1.0;
	Material mat(e, v, rho);

	// Create members
	FrameMember column(1, &node1, &node2, &frameSection, &mat);
	FrameMember beam(2, &node2, &node3, &frameSection, &mat);

	std::map<unsigned int, Element*> members;
	members[column.ElementIndex] = &column;
	members[beam.ElementIndex] = &beam;
	
	// Add tip load
	double loads[6]{ 0.0, -5000.0, 0.0, 0.0, 0.0, 0.0 };

	NodalLoad nLFront(&node3, loads);
	std::map<unsigned int, NodalLoad*> nload;
	nload[1] = &nLFront;

	std::map<unsigned int, DistributedLoad*> dist;

	auto str = new Structure(&nodes, &members, &restraints, &nload, &dist);

	LOG("__________________________________________________");
	LOG("|                                                |");
	LOG("|    3-Dimensional Structural Analysis Engine    |");
	LOG("|       Created by Mustafa Can Sanliturk         |");
	LOG("|           All rights reserved ©                |");
	LOG("|________________________________________________|");
	LOG("");

	auto disps = ArmadilloSolver::GetDisplacementForStaticCase(*str);

	for (auto& nodePair : nodes)
	{
		auto node = nodePair.second;
		
		LOG("");
		LOG(" Node Index: ");
		std::cout << " " << node->NodeIndex << "\n";

		auto nodalDisps = ArmadilloSolver::GetNodalDisplacements(*node, disps);

		for (size_t i = 0; i < 6; i++)
			std::cout << " DOF Index: " << i  + 1 << ", Displacement = " << nodalDisps[i] << "\n";
	}

	std::cin.get();

	delete str;

	return 0;
}