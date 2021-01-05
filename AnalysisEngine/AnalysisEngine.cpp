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
	::SendMessage(::GetConsoleWindow(), WM_SYSKEYDOWN, VK_RETURN, 0x20000000);

	// Create points for nodes
	XYZPoint bottomLeft(0.0, 0.0, 0.0);
	XYZPoint topLeft(0.0, 0.0, 1.0);
	XYZPoint topRight(1.0, 0.0, 1.0);
	XYZPoint bottomRight(1.0, 0.0, 0);
	XYZPoint bottomLeftBack(0.0, 1.0, 0.0);
	XYZPoint topLeftBack(0.0, 1.0, 1.0);
	XYZPoint topRightBack(1.0, 1.0, 1.0);
	XYZPoint bottomRightBack(1.0, 1.0, 0);

	// Define nodes
	Node node1(1, bottomLeft);
	Node node2(2, topLeft);
	Node node3(3, topRight);
	Node node4(4, bottomRight);
	Node node5(5, bottomLeftBack);
	Node node6(6, topLeftBack);
	Node node7(7, topRightBack);
	Node node8(8, bottomRightBack);

	// Add nodes to 
	std::map<unsigned int, Node*> nodes;
	nodes[node1.NodeIndex] = &node1;
	nodes[node2.NodeIndex] = &node2;
	nodes[node3.NodeIndex] = &node3;
	nodes[node4.NodeIndex] = &node4;
	nodes[node5.NodeIndex] = &node5;
	nodes[node6.NodeIndex] = &node6;
	nodes[node7.NodeIndex] = &node7;
	nodes[node8.NodeIndex] = &node8;

	// Assign boundary condition to the first node
	// and fourth node
	// Fixed for all directions
	std::vector<bool> fix;
	std::vector<double> fixity;

	for (unsigned int i = 0; i < 6; i++)
	{
		fix.push_back(true);
		fixity.push_back(0.0);
	}

	Restraint leftRestraint(&node1, fix, fixity);
	Restraint rightRestraint(&node4, fix, fixity);
	Restraint leftRestraintBack(&node5, fix, fixity);
	Restraint rightRestraintBack(&node8, fix, fixity);

	std::map<unsigned int, Restraint*> restraints;
	restraints[1] = &leftRestraint;
	restraints[2] = &rightRestraint;
	restraints[3] = &leftRestraintBack;
	restraints[4] = &rightRestraintBack;

	// Create section
	auto area = 1.0;
	auto inertia11 = 1.0 / 12.0;
	auto inertia22 = inertia11;
	auto inertia12 = 0.14083333333;

	Section columnSection(area, inertia11, inertia22, inertia12);

	// To have axially rigid beams
	Section beamSection(area * 1000000, inertia11, inertia22, inertia12);

	// Create material
	auto e = 1.0;
	auto v = 0.3;
	auto rho = 1.0;
	Material mat(e, v, rho);

	// Create members
	FrameMember leftColumn(1, &node1, &node2, &columnSection, &mat);
	FrameMember rightColumn(2, &node3, &node4, &columnSection, &mat);
	FrameMember beam(3, &node2, &node3, &beamSection, &mat);
	FrameMember leftColumnBack(4, &node5, &node6, &columnSection, &mat);
	FrameMember rightColumnBack(5, &node7, &node8, &columnSection, &mat);
	FrameMember beamBack(6, &node6, &node7, &beamSection, &mat);
	FrameMember leftConnectionBeam(7, &node2, &node6, &beamSection, &mat);
	FrameMember rightConnectionBeam(8, &node3, &node7, &beamSection, &mat);

	std::map<unsigned int, Element*> members;
	members[leftColumn.ElementIndex] = &leftColumn;
	members[rightColumn.ElementIndex] = &rightColumn;
	members[beam.ElementIndex] = &beam;
	members[leftColumnBack.ElementIndex] = &leftColumnBack;
	members[rightColumnBack.ElementIndex] = &rightColumnBack;
	members[beamBack.ElementIndex] = &beamBack;
	members[leftConnectionBeam.ElementIndex] = &leftConnectionBeam;
	members[rightConnectionBeam.ElementIndex] = &rightConnectionBeam;

	// Add tip load
	double loads[6]{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double loads2[6]{ -1.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	NodalLoad nLFront(&node2, loads);
	NodalLoad nLBack(&node7, loads2);
	std::map<unsigned int, NodalLoad*> nload;
	nload[1] = &nLFront;
	nload[2] = &nLBack;

	std::map<unsigned int, DistributedLoad*> dist;

	//Structure str(&nodes, &members, &restraints, &nload, &dist);
	auto str = new Structure(&nodes, &members, &restraints, &nload, &dist);

	LOG("__________________________________________________");
	LOG("|                                                |");
	LOG("|    3-Dimensional Structural Analysis Engine    |");
	LOG("|       Created by Mustafa Can Sanliturk         |");
	LOG("|           All rights reserved ©                |");
	LOG("|________________________________________________|");
	LOG("");

	auto disps = ArmadilloSolver::GetDisplacementForStaticCase(*str);
	ArmadilloSolver::GetMemberEndForcesForLocalCoordinates(leftColumn, disps);
	ArmadilloSolver::GetMemberEndForcesForGlobalCoordinates(leftColumn, disps);

	for (auto& nodePair : nodes)
	{
		auto node = nodePair.second;
		
		LOG("");
		LOG(" Node Index: ");
		std::cout << " " << node->NodeIndex << "\n";

		auto nodalDisps = ArmadilloSolver::GetNodalDisplacements(*node, disps);

		for (size_t i = 0; i < 6; i++)
			std::cout << " DOF Index: " << i << ", Displacement = " << nodalDisps[i] << "\n";
	}
	
	std::cin.get();

	delete str;

	return 0;
}