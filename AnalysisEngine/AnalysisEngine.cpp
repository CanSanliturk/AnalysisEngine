#include <iostream>
#include <vector>
#include <memory>
#include <chrono>
#include <Eigen>
#include "StructureSolver.h"
#include "GeometryHelper.h"
#include "FrameMember.h"
#include "Restraint.h"
#include "Structure.h"
#include "XYZPoint.h"
#include "Section.h"
#include "Element.h"
#include "Vector.h"
#include "Matrix.h"
#include "Node.h"

#define LOG(x) std::cout << x << "\n"

void CantileverDisplacements();
void TableDisplacements();

int main()
{
    LOG("__________________________________________________");
    LOG("|                                                |");
    LOG("|    3-Dimensional Structural Analysis Engine    |");
    LOG("|       Created by Mustafa Can Sanliturk         |");
    LOG("|           All rights reserved ©                |");
    LOG("|________________________________________________|");
    LOG("");

    // Start timer
    auto timenow =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    // Call test function (Later on, these guys will be moved to a unit test project)
    TableDisplacements();
    LOG("\n Analysis completed without errors....");

    // Log duration
    auto timenow2 =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    LOG(" Elapsed Time = " << timenow2 - timenow << " seconds\n");

    std::cin.get();
    return 0;
}

void CantileverDisplacements()
{
    // Units are in N & m

    // Coordinates
    XYZPoint pt1(0, 0, 0); // Origin
    XYZPoint pt2(5, 5, 5);

    // Nodes
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    auto node1 = std::make_shared<Node>(1, pt1); nodes[node1->NodeIndex] = node1;
    auto node2 = std::make_shared<Node>(2, pt2); nodes[node2->NodeIndex] = node2;

    // Section
    auto area = 0.16;
    auto inertia11 = 2.133e-3;
    auto inertia22 = 2.133e-3;
    auto inertia12 = 3.605e-3;
    auto sect = std::make_shared<Section>(area, inertia11, inertia22, inertia12);

    // Material
    auto mat = std::make_shared<Material>(200e9, 0.3, 78500);

    // Releases
    auto h1 = std::make_shared<Hinge>();
    auto h2 = std::make_shared<Hinge>();

    // Members
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    auto beam = std::make_shared<FrameMember>(1, node1, node2, sect, mat, h1, h2); elements[beam->ElementIndex] = beam;

    // Restraints
    std::vector<bool> isRest;
    std::vector<double> rest;

    for (int i = 0; i < 6; i++)
    {
        isRest.push_back(true);
        rest.push_back(0.0);
    }
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    auto res1 = std::make_shared<Restraint>(node1, isRest, rest); restraints[1] = res1;

    // Nodal loads
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    double nodalLoad[6] = { 0, 0, -5000e3, 0, 0, 0 };
    auto nl1 = std::make_shared<NodalLoad>(node2, nodalLoad); nodalLoads[1] = nl1;

    // Distributed loads
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Armadillo);
    for (auto& nodePair : nodes)
    {
        auto node = nodePair.second;

        LOG("");
        LOG(" Node Index: ");
        std::cout << " " << node->NodeIndex << "\n";

        auto nodalDisps = StructureSolver::GetNodalDisplacements(*node, disps);

        for (size_t i = 0; i < 6; i++)
            std::cout << " DOF Index: " << i + 1 << ", Displacement = " << nodalDisps(i, 0) << "\n";
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
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    auto  bottomNode1 = std::make_shared<Node>(1, bottomPt1); nodes[1] = bottomNode1;
    auto  bottomNode2 = std::make_shared<Node>(2, bottomPt2); nodes[2] = bottomNode2;
    auto  bottomNode3 = std::make_shared<Node>(3, bottomPt3);	nodes[3] = bottomNode3;
    auto  bottomNode4 = std::make_shared<Node>(4, bottomPt4);	nodes[4] = bottomNode4;
    auto  topNode1 = std::make_shared<Node>(5, topPt1); nodes[5] = topNode1;
    auto  topNode2 = std::make_shared<Node>(6, topPt2); nodes[6] = topNode2;
    auto  topNode3 = std::make_shared<Node>(7, topPt3); nodes[7] = topNode3;
    auto  topNode4 = std::make_shared<Node>(8, topPt4); nodes[8] = topNode4;

    // Section
    auto area = 0.16;
    auto inertia11 = 2.133 * 0.001;
    auto inertia22 = 2.133 * 0.001;
    auto inertia12 = 0.003605;
    auto sect = std::make_shared<Section>(area, inertia11, inertia22, inertia12);

    // Material
    auto mat = std::make_shared<Material>(200e9, 0.3, 78500);

    // Releases
    std::shared_ptr<Hinge> h1i, h1j;
    std::shared_ptr<Hinge> h2i, h2j;
    std::shared_ptr<Hinge> h3i, h3j;
    std::shared_ptr<Hinge> h4i, h4j;
    std::shared_ptr<Hinge> h5i, h5j;
    std::shared_ptr<Hinge> h6i, h6j;
    std::shared_ptr<Hinge> h7i, h7j;
    std::shared_ptr<Hinge> h8i, h8j;

    // Members
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    auto col1 = std::make_shared<FrameMember >(1, bottomNode1, topNode1, sect, mat, h1i, h1j); elements[1] = col1;
    auto col2 = std::make_shared<FrameMember >(2, bottomNode2, topNode2, sect, mat, h2i, h2j); elements[2] = col2;
    auto col3 = std::make_shared<FrameMember >(3, bottomNode3, topNode3, sect, mat, h3i, h3j); elements[3] = col3;
    auto col4 = std::make_shared<FrameMember >(4, bottomNode4, topNode4, sect, mat, h4i, h4j); elements[4] = col4;
    auto beam1 = std::make_shared<FrameMember >(5, topNode1, topNode2, sect, mat, h5i, h5j); elements[5] = beam1;
    auto beam2 = std::make_shared<FrameMember >(6, topNode2, topNode3, sect, mat, h6i, h6j); elements[6] = beam2;
    auto beam3 = std::make_shared<FrameMember >(7, topNode3, topNode4, sect, mat, h7i, h7j); elements[7] = beam3;
    auto beam4 = std::make_shared<FrameMember >(8, topNode4, topNode1, sect, mat, h8i, h8j); elements[8] = beam4;

    // Restraints
    std::vector<bool> isRest;
    std::vector<double> rest;

    for (int i = 0; i < 6; i++)
    {
        isRest.push_back(true);
        rest.push_back(0.0);
    }

    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    auto res1 = std::make_shared<Restraint>(bottomNode1, isRest, rest); restraints[1] = res1;
    auto res2 = std::make_shared<Restraint>(bottomNode2, isRest, rest); restraints[2] = res2;
    auto res3 = std::make_shared<Restraint>(bottomNode3, isRest, rest);	restraints[3] = res3;
    auto res4 = std::make_shared<Restraint>(bottomNode4, isRest, rest);	restraints[4] = res4;

    // Nodal loads
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    double nodalLoad[6] = { 5000e3, 0, 0, 0, 0, 0 };
    double nodalLoad2[6] = { 5000e3, 0, 0, 0, 0, 0 };
    auto nl1 = std::make_shared<NodalLoad>(topNode1, nodalLoad); nodalLoads[1] = nl1;
    auto nl2 = std::make_shared<NodalLoad>(topNode4, nodalLoad2); nodalLoads[2] = nl2;

    // Distributed loads
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Eigen);
    LOG(" Displacement values:");
    for (auto& nodePair : nodes)
    {
        auto node = nodePair.second;

        LOG("");
        LOG(" Node Index: ");
        std::cout << " " << node->NodeIndex << "\n";

        auto nodalDisps = StructureSolver::GetNodalDisplacements(*node, disps);

        for (size_t i = 0; i < 6; i++)
            std::cout << " DOF Index: " << i + 1 << ", Displacement = " << nodalDisps(i, 0) << "\n";
    }

    // Modal periods
    auto modalPeriods = StructureSolver::GetModalPeriods(*str, SolverChoice::Eigen);
    LOG("\n Modal periods:");
    for (size_t i = 0; i < modalPeriods.RowCount; i++)
        std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods(i, 0) << " s\n";

    return;
}
