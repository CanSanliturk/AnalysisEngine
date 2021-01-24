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
void TableDisplacements();
void CantileverDisplacements();

int main()
{
	// Full screen
	//::SendMessage(::GetConsoleWindow(), WM_SYSKEYDOWN, VK_RETURN, 0x20000000);
	LOG("__________________________________________________");
	LOG("|                                                |");
	LOG("|    3-Dimensional Structural Analysis Engine    |");
	LOG("|       Created by Mustafa Can Sanliturk         |");
	LOG("|           All rights reserved ©                |");
	LOG("|________________________________________________|");
	LOG("");

	CantileverDisplacements();

	std::cin.get();
	return 0;
}

void CantileverDisplacements()
{
	// Coordinates
	XYZPoint pt1(0, 0, 0); // Origin
	XYZPoint pt2(5, 5, 5);

	// Nodes
	std::map<unsigned int, Node*> nodes;
	Node node1(1, pt1); nodes[node1.NodeIndex] = &node1;
	Node node2(2, pt2); nodes[node2.NodeIndex] = &node2;

	// Section
	auto area = 0.16;
	auto inertia11 = 2.133e-3;
	auto inertia22 = 2.133e-3;
	auto inertia12 = 0.0036;
	Section sect(area, inertia11, inertia22, inertia12);

	// Material
	Material mat(200e9, 0.3, 0);

	// Members
	std::map<unsigned int, Element*> elements;
	FrameMember beam(1, &node1, &node2, &sect, &mat); elements[beam.ElementIndex] = &beam;
	
	// Restraints
	std::vector<bool> isRest;
	std::vector<double> rest;

	for (int i = 0; i < 6; i++)
	{
		isRest.push_back(true);
		rest.push_back(0.0);
	}
	std::map<unsigned int, Restraint*> restraints;
	Restraint res1(&node1, isRest, rest); restraints[1] = &res1;

	// Nodal loads
	std::map<unsigned int, NodalLoad*> nodalLoads;
	double nodalLoad[6] = { 0, 0, -5000e3, 0, 0, 0 };
	NodalLoad nl1(&node2, nodalLoad); nodalLoads[1] = &nl1;

	// Distributed loads
	std::map<unsigned int, DistributedLoad*> distLoads;

	// Create structure
	Structure str(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

	// Solve displacement
	auto disps = ArmadilloSolver::GetDisplacementForStaticCase(str);
	for (auto& nodePair : nodes)
	{
		auto node = nodePair.second;

		LOG("");
		LOG(" Node Index: ");
		std::cout << " " << node->NodeIndex << "\n";

		auto nodalDisps = ArmadilloSolver::GetNodalDisplacements(*node, disps);

		for (size_t i = 0; i < 6; i++)
			std::cout << " DOF Index: " << i + 1 << ", Displacement = " << nodalDisps[i] << "\n";
	}
	return;
}

void TableDisplacements()
{
	// Coordinates
	XYZPoint bottomPt1(0, 0, 0); // Origin
	XYZPoint bottomPt2(5, 0, 0);
	XYZPoint bottomPt3(5, 5, 0);
	XYZPoint bottomPt4(0, 5, 0);
	XYZPoint topPt1(0, 0, 5); // Top of the origin
	XYZPoint topPt2(5, 0, 5);
	XYZPoint topPt3(5, 5, 5);
	XYZPoint topPt4(0, 5, 5);

	// Nodes
	std::map<unsigned int, Node*> nodes;
	Node bottomNode1(1, bottomPt1); nodes[1] = &bottomNode1;
	Node bottomNode2(2, bottomPt2); nodes[2] = &bottomNode2;
	Node bottomNode3(3, bottomPt3);	nodes[3] = &bottomNode3;
	Node bottomNode4(4, bottomPt4);	nodes[4] = &bottomNode4;
	Node topNode1(5, topPt1); nodes[5] = &topNode1;
	Node topNode2(6, topPt2); nodes[6] = &topNode2;
	Node topNode3(7, topPt3); nodes[7] = &topNode3;
	Node topNode4(8, topPt4); nodes[8] = &topNode4;

	// Section
	auto area = 0.16;
	auto inertia11 = 2.133 * 0.001;
	auto inertia22 = 2.133 * 0.001;
	auto inertia12 = 0.0036;
	Section sect(area, inertia11, inertia22, inertia12);

	// Material
	Material mat(200e9, 0.3, 0);

	// Members
	std::map<unsigned int, Element*> elements;
	FrameMember col1(1, &bottomNode1, &topNode1, &sect, &mat); elements[1] = &col1;
	FrameMember col2(2, &bottomNode2, &topNode2, &sect, &mat); elements[2] = &col2;
	FrameMember col3(3, &bottomNode3, &topNode3, &sect, &mat); elements[3] = &col3;
	FrameMember col4(4, &bottomNode4, &topNode4, &sect, &mat); elements[4] = &col4;
	FrameMember beam1(5, &topNode1, &topNode2, &sect, &mat); elements[5] = &beam1;
	FrameMember beam2(6, &topNode2, &topNode3, &sect, &mat); elements[6] = &beam2;
	FrameMember beam3(7, &topNode3, &topNode4, &sect, &mat); elements[7] = &beam3;
	FrameMember beam4(7, &topNode4, &topNode1, &sect, &mat); elements[8] = &beam4;

	// Restraints
	std::vector<bool> isRest;
	std::vector<double> rest;

	for (int i = 0; i < 6; i++)
	{
		isRest.push_back(true);
		rest.push_back(0.0);
	}

	std::map<unsigned int, Restraint*> restraints;
	Restraint res1(&bottomNode1, isRest, rest); restraints[1] = &res1;
	Restraint res2(&bottomNode2, isRest, rest); restraints[2] = &res2;
	Restraint res3(&bottomNode3, isRest, rest);	restraints[3] = &res3;
	Restraint res4(&bottomNode4, isRest, rest);	restraints[4] = &res4;

	// Nodal loads
	std::map<unsigned int, NodalLoad*> nodalLoads;
	double nodalLoad[6] = { 5000, 0, 0, 0, 0, 0 };
	double nodalLoad2[6] = { 5000, 0, 0, 0, 0, 0 };
	NodalLoad nl1(&topNode1, nodalLoad); nodalLoads[1] = &nl1;
	NodalLoad nl2(&topNode4, nodalLoad2); nodalLoads[2] = &nl2;

	// Distributed loads
	std::map<unsigned int, DistributedLoad*> distLoads;

	// Create structure
	Structure str(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

	// Solve displacement
	auto disps = ArmadilloSolver::GetDisplacementForStaticCase(str);
	for (auto& nodePair : nodes)
	{
		auto node = nodePair.second;

		LOG("");
		LOG(" Node Index: ");
		std::cout << " " << node->NodeIndex << "\n";

		auto nodalDisps = ArmadilloSolver::GetNodalDisplacements(*node, disps);

		for (size_t i = 0; i < 6; i++)
			std::cout << " DOF Index: " << i + 1 << ", Displacement = " << nodalDisps[i] << "\n";
	}
	return;
}