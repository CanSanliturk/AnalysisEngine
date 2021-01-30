#include <iostream>
#include <vector>
#include <memory>
#include <Windows.h>
#include "XYZPoint.h"
#include "Restraint.h"
#include "Node.h"
#include "Section.h"
#include "Element.h"
#include "FrameMember.h"
#include "TrussMember.h"
#include "MatrixHelper.h"
#include "Structure.h"
#include "ArmadilloSolver.h"

#pragma comment(lib, "user32")

#define LOG(x) std::cout << x << "\n"
void TableDisplacements();
void CantileverDisplacements();
void LDisplacements();
void TriangleTruss();

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

	// Call test function (Later on, these guys will be moved to a unit test project)
	TriangleTruss();
	
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

	// Hinge
	std::map<unsigned int, Hinge*> hinges;

	// Nodal loads
	std::map<unsigned int, NodalLoad*> nodalLoads;
	double nodalLoad[6] = { 0, 0, -5000e3, 0, 0, 0 };
	NodalLoad nl1(&node2, nodalLoad); nodalLoads[1] = &nl1;

	// Distributed loads
	std::map<unsigned int, DistributedLoad*> distLoads;

	// Create structure
	Structure str(&nodes, &elements, &restraints, &hinges, &nodalLoads, &distLoads);

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
	FrameMember beam4(8, &topNode4, &topNode1, &sect, &mat); elements[8] = &beam4;

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

	// Hinge
	std::map<unsigned int, Hinge*> hinges;

	// Nodal loads
	std::map<unsigned int, NodalLoad*> nodalLoads;
	double nodalLoad[6] = { 5000e3, 0, 0, 0, 0, 0 };
	double nodalLoad2[6] = { 5000e3, 0, 0, 0, 0, 0 };
	NodalLoad nl1(&topNode1, nodalLoad); nodalLoads[1] = &nl1;
	NodalLoad nl2(&topNode4, nodalLoad2); nodalLoads[2] = &nl2;

	// Distributed loads
	std::map<unsigned int, DistributedLoad*> distLoads;

	// Create structure
	auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &hinges, &nodalLoads, &distLoads);

	// Solve displacement
	auto disps = ArmadilloSolver::GetDisplacementForStaticCase(*str);
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

void LDisplacements()
{
	// Coordinates
	XYZPoint botLeft(-5, 0, 0);
	XYZPoint botMid(0, 0, 0);
	XYZPoint botRight(0, 5, 0);
	XYZPoint topLeft(-5, 0, 5);
	XYZPoint topMid(0, 0, 5);
	XYZPoint topRight(0, 5, 5);

	// Nodes
	std::map<unsigned int, Node*> nodes;
	Node botLeftNode(1, botLeft); nodes[1] = &botLeftNode;
	Node botMidNode(2, botMid); nodes[2] = &botMidNode;
	Node botRightNode(3, botRight); nodes[3] = &botRightNode;
	Node topLeftNode(4, topLeft); nodes[4] = &topLeftNode;
	Node topMidNode(5, topMid); nodes[5] = &topMidNode;
	Node topRightNode(6, topRight); nodes[6] = &topRightNode;

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
	FrameMember colLeft(1, &botLeftNode, &topLeftNode, &sect, &mat); elements[1] = &colLeft;
	FrameMember colMid(2, &botMidNode, &topMidNode, &sect, &mat); elements[2] = &colMid;
	FrameMember colRight(3, &botRightNode, &topRightNode, &sect, &mat); elements[3] = &colRight;
	FrameMember beamLeft(4, &topLeftNode, &topMidNode, &sect, &mat); elements[4] = &beamLeft;
	FrameMember beamRight(5, &topMidNode, &topRightNode, &sect, &mat); elements[5] = &beamRight;

	// Restraints
	std::vector<bool> isRest;
	std::vector<double> rest;

	for (int i = 0; i < 6; i++)
	{
		isRest.push_back(true);
		rest.push_back(0.0);
	}

	// Hinge
	std::map<unsigned int, Hinge*> hinges;

	std::map<unsigned int, Restraint*> restraints;
	Restraint res1(&botLeftNode, isRest, rest); restraints[1] = &res1;
	Restraint res2(&botMidNode, isRest, rest); restraints[2] = &res2;
	Restraint res3(&botRightNode, isRest, rest); restraints[3] = &res3;

	// Nodal loads
	std::map<unsigned int, NodalLoad*> nodalLoads;
	double nodalLoad[6] = { 0, 0, -5000e3, 0, 0, 0 };
	NodalLoad nl1(&topMidNode, nodalLoad); nodalLoads[1] = &nl1;

	// Distributed loads
	std::map<unsigned int, DistributedLoad*> distLoads;

	// Create structure
	Structure str(&nodes, &elements, &restraints, &hinges, &nodalLoads, &distLoads);

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

void TriangleTruss()
{
	// Coordinates
	XYZPoint leftPt(0, 0, 0); // Origin
	XYZPoint rightPt(10, 0, 0);
	XYZPoint topPt(5, 5, 0);
	
	// Nodes
	std::map<unsigned int, Node*> nodes;
	Node leftNode(1, leftPt); nodes[1] = &leftNode;
	Node rightNode(2, rightPt); nodes[2] = &rightNode;
	Node topNode(3, topPt);	nodes[3] = &topNode;
	
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
	FrameMember bottom(1, &leftNode, &rightNode, &sect, &mat); elements[1] = &bottom;
	FrameMember rightMember(2, &rightNode, &topNode, &sect, &mat); elements[2] = &rightMember;
	FrameMember leftMember(3, &topNode, &leftNode, &sect, &mat); elements[3] = &leftMember;
	
	// Restraints
	std::vector<bool> isRest;
	std::vector<double> rest;

	for (int i = 0; i < 6; i++)
	{
		if (i < 3) isRest.push_back(true);
		else isRest.push_back(false);

		rest.push_back(0.0);
	}

	std::map<unsigned int, Restraint*> restraints;
	Restraint res1(&leftNode, isRest, rest); restraints[1] = &res1;
	Restraint res2(&rightNode, isRest, rest); restraints[2] = &res2;

	// Hinge
	std::vector<bool> isForceSet;;
	std::vector<double> force;

	for (size_t i = 0; i < 6; i++)
	{
		if (i < 3) isForceSet.push_back(false); // Moment release
		else isForceSet.push_back(true);

		force.push_back(0.0);
	}

	std::map<unsigned int, Hinge*> hinges;
	Hinge h1(&leftNode, isForceSet, force); hinges[1] = &h1;
	Hinge h2(&rightNode, isForceSet, force); hinges[2] = &h2;
	Hinge h3(&topNode, isForceSet, force); hinges[3] = &h3;

	// Nodal loads
	std::map<unsigned int, NodalLoad*> nodalLoads;
	double nodalLoad[6] = { 5000e3, 0, 0, 0, 0, 0 };
	NodalLoad nl1(&topNode, nodalLoad); nodalLoads[1] = &nl1;

	// Distributed loads
	std::map<unsigned int, DistributedLoad*> distLoads;

	// Create structure
	Structure str(&nodes, &elements, &restraints, &hinges, &nodalLoads, &distLoads);

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

	auto rightMemForces = ArmadilloSolver::GetMemberEndForcesForLocalCoordinates(rightMember, disps);
	LOG("");
	for (auto& f : rightMemForces) LOG(f);

	return;
}