#include <iostream>
#include <vector>
#include <memory>
#include <Windows.h>
#include "XYZPoint.h"
#include "Vector.h"
#include "Restraint.h"
#include "Node.h"
#include "Section.h"
#include "Element.h"
#include "FrameMember.h"
#include "TrussMember.h"
#include "MatrixHelper.h"
#include "Structure.h"
#include "ArmadilloSolver.h"
#include "TimoshenkoBeam.h"
#include <iostream>
#include <Eigen>

#pragma comment(lib, "user32")

#define LOG(x) std::cout << x << "\n"
#define MAKESHARED(x) std::make_shared<x>

void TableDisplacements();
void CantileverDisplacements();
void TriangleTruss();
void CantileverColumn();

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
    CantileverColumn();
    LOG("\n Analysis completed without errors....");
    std::cin.get();
    return 0;
}

void CantileverColumn()
{
    // Units are in N & m

    // Coordinates
    XYZPoint pt1(0.0, 0.0, 0); // Origin
    XYZPoint pt2(0.0, 10.0, 0.0);

    // Nodes
    std::map<unsigned int, Node*> nodes;
    Node node1(1, pt1); nodes[node1.NodeIndex] = &node1;
    Node node2(2, pt2); nodes[node2.NodeIndex] = &node2;

    // Section
    auto area = 0.16;
    auto inertia11 = 2.133e-3;
    auto inertia22 = 2.133e-3;
    auto inertia12 = 0.003605;
    Section sect(area, inertia11, inertia22, inertia12);

    // Material
    Material mat(200e9, 0.3, 78500);

    // Releases
    Hinge h1;
    Hinge h2;

    // Members
    std::map<unsigned int, Element*> elements;
    FrameMember beam(1, &node1, &node2, &sect, &mat, &h1, &h2); elements[beam.ElementIndex] = &beam;

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
    double nodalLoad[6] = { 100000, 0, 0, 0, 0, 0 };
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

    // Modal periods
    // auto modalPeriods = ArmadilloSolver::GetModalPeriods(str);
    // LOG("");
    // for (size_t i = 0; i < modalPeriods.size(); i++)
    //     std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods[i] << "\n";

    return;
}

void CantileverDisplacements()
{
    // Units are in N & m

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
    auto inertia12 = 0.003605;
    Section sect(area, inertia11, inertia22, inertia12);

    // Material
    Material mat(200e9, 0.3, 78500);

    // Releases
    Hinge h1;
    Hinge h2;

    // Members
    std::map<unsigned int, Element*> elements;
    FrameMember beam(1, &node1, &node2, &sect, &mat, &h1, &h2); elements[beam.ElementIndex] = &beam;

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

    // Modal periods
    // auto modalPeriods = ArmadilloSolver::GetModalPeriods(str);
    // LOG("");
    // for (size_t i = 0; i < modalPeriods.size(); i++)
    //     std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods[i] << "\n";

    return;
}

void TableDisplacements()
{
    // Coordinates
    XYZPoint bottomPt1(0, 0, 0); // Origin
    XYZPoint bottomPt2(10, 0, 0);
    XYZPoint bottomPt3(10, 10, 0);
    XYZPoint bottomPt4(0, 10, 0);
    XYZPoint topPt1(0, 0, 10); // Top of the origin
    XYZPoint topPt2(10, 0, 10);
    XYZPoint topPt3(10, 10, 10);
    XYZPoint topPt4(0, 10, 10);

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
    auto inertia12 = 0.003605;
    Section sect(area, inertia11, inertia22, inertia12);

    // Material
    Material mat(200e9, 0.3, 78500);

    // Releases
    Hinge h1i, h1j;
    Hinge h2i, h2j;
    Hinge h3i, h3j;
    Hinge h4i, h4j;
    Hinge h5i, h5j;
    Hinge h6i, h6j;
    Hinge h7i, h7j;
    Hinge h8i, h8j;

    // Members
    std::map<unsigned int, Element*> elements;
    FrameMember col1(1, &bottomNode1, &topNode1, &sect, &mat, &h1i, &h1j); elements[1] = &col1;
    FrameMember col2(2, &bottomNode2, &topNode2, &sect, &mat, &h2i, &h2j); elements[2] = &col2;
    FrameMember col3(3, &bottomNode3, &topNode3, &sect, &mat, &h3i, &h3j); elements[3] = &col3;
    FrameMember col4(4, &bottomNode4, &topNode4, &sect, &mat, &h4i, &h4j); elements[4] = &col4;
    FrameMember beam1(5, &topNode1, &topNode2, &sect, &mat, &h5i, &h5j); elements[5] = &beam1;
    FrameMember beam2(6, &topNode2, &topNode3, &sect, &mat, &h6i, &h6j); elements[6] = &beam2;
    FrameMember beam3(7, &topNode3, &topNode4, &sect, &mat, &h7i, &h7j); elements[7] = &beam3;
    FrameMember beam4(8, &topNode4, &topNode1, &sect, &mat, &h8i, &h8j); elements[8] = &beam4;

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
    double nodalLoad[6] = { 5000e3, 0, 0, 0, 0, 0 };
    double nodalLoad2[6] = { 5000e3, 0, 0, 0, 0, 0 };
    NodalLoad nl1(&topNode1, nodalLoad); nodalLoads[1] = &nl1;
    NodalLoad nl2(&topNode4, nodalLoad2); nodalLoads[2] = &nl2;

    // Distributed loads
    std::map<unsigned int, DistributedLoad*> distLoads;

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    /*auto disps = ArmadilloSolver::GetDisplacementForStaticCase(*str);
    for (auto& nodePair : nodes)
    {
        auto node = nodePair.second;

        LOG("");
        LOG(" Node Index: ");
        std::cout << " " << node->NodeIndex << "\n";

        auto nodalDisps = ArmadilloSolver::GetNodalDisplacements(*node, disps);

        for (size_t i = 0; i < 6; i++)
            std::cout << " DOF Index: " << i + 1 << ", Displacement = " << nodalDisps[i] << "\n";
    }*/

    // Modal periods
    auto modalPeriods = ArmadilloSolver::GetModalPeriods(*str);
    LOG("");
    for (size_t i = 0; i < modalPeriods.size(); i++)
        std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods[i] << "\n";

    return;
}

void TriangleTruss()
{
    // Coordinates
    XYZPoint leftPt(0, 0, 0); // Origin
    XYZPoint rightPt(10, 0, 0);
    XYZPoint topPt(5, 0, 5);

    // Nodes
    std::map<unsigned int, Node*> nodes;
    Node leftNode(1, leftPt); nodes[1] = &leftNode;
    Node rightNode(2, rightPt); nodes[2] = &rightNode;
    Node topNode(3, topPt);	nodes[3] = &topNode;

    // Section
    auto area = 0.16;
    auto inertia11 = (1 / 12) * 0.4 * 0.4 * 0.4 * 0.4;
    auto inertia22 = inertia11;
    auto inertia12 = 0.003605333333;
    Section sect(area, inertia11, inertia22, inertia12);

    // Material
    Material mat(200e9, 0.3, 0);

    // Hinge
    std::vector<bool> isReleased{ false, false, false, false, true, true };
    std::vector<double> release{ 0.0,0.0 ,0.0 ,0.0 ,0.0 ,0.0 };

    Hinge h1i(isReleased, release);
    Hinge h1j(isReleased, release);

    Hinge h2i(isReleased, release);
    Hinge h2j(isReleased, release);

    Hinge h3i(isReleased, release);
    Hinge h3j(isReleased, release);

    // Members
    std::map<unsigned int, Element*> elements;
    FrameMember bottomMember(1, &leftNode, &rightNode, &sect, &mat, &h1i, &h1j); elements[1] = &bottomMember;
    FrameMember rightMember(2, &rightNode, &topNode, &sect, &mat, &h2i, &h2j); elements[2] = &rightMember;
    FrameMember leftMember(3, &topNode, &leftNode, &sect, &mat, &h3i, &h3j); elements[3] = &leftMember;

    // Restraints
    std::vector<bool> isRestPin{ true, true, true, true, true, true };
    std::vector<double> restPin{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    std::vector<bool> isRestRoller{ false, true, true, true, true, true };
    std::vector<double> restRoller{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    std::map<unsigned int, Restraint*> restraints;
    Restraint res1(&leftNode, isRestPin, restPin); restraints[1] = &res1;
    Restraint res2(&rightNode, isRestRoller, restRoller); restraints[2] = &res2;

    // Nodal loads
    std::map<unsigned int, NodalLoad*> nodalLoads;
    double nodalLoad[6] = { 0, 0, -5000e3, 0, 0, 0 };
    NodalLoad nl1(&topNode, nodalLoad); nodalLoads[1] = &nl1;

    // Distributed loads
    std::map<unsigned int, DistributedLoad*> distLoads;

    // Create structure
    Structure str(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Modal periods
   /* auto modalPeriods = ArmadilloSolver::GetModalPeriods(str);
    LOG("");
    for (size_t i = 0; i < modalPeriods.size(); i++)
        std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods[i] << "\n";*/

    LOG("");
    LOG("");

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

    auto leftMemForces = ArmadilloSolver::GetMemberEndForcesForLocalCoordinates(leftMember, disps);
    LOG("---------------------------------------");
    LOG("Left member forces");
    for (auto& f : leftMemForces) LOG(f);

    // Check support reactions
    auto leftSupportReactions = ArmadilloSolver::GetSupportReactions(str, disps, res1);
    LOG("---------------------------------------");
    LOG("Left support Reaction");
    for (auto& r : leftSupportReactions) LOG(r);
}