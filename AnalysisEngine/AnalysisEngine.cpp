#include <iostream>
#include <vector>
#include <memory>
#include <chrono>
#include "StructureSolver.h"
#include "Vector.h"
#include "UtilMethods.h"


#define LOG(x) std::cout << x << "\n"
#define MATRIX Matrix<double>
constexpr double pi = 3.141592653589793;

void CantileverDisplacements3D();
void CantileverDisplacements2D();
void TableDisplacements();
void CE583Sample();
void CE583Assignment1_3();
void CE583Assignment5();
void CE583Assignment6_2();
void CE583Assignment6_3();
void TrussExample();

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
    try
    {
        CE583Assignment6_3();
    }
    catch (const std::runtime_error& e)
    {
        LOG(e.what());
    }
    LOG("\n Analysis completed without errors....");

    // Log duration
    auto timenow2 =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    LOG(" Elapsed Time = " << timenow2 - timenow << " seconds\n");

    std::cin.get();
    return 0;
}

void CantileverDisplacements3D()
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

    // Members
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    bool isLumpedMassMatrix = true;
    auto beam = std::make_shared<FrameMember>(1, node1, node2, sect, mat, isLumpedMassMatrix); elements[beam->ElementIndex] = beam;

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

void CantileverDisplacements2D()
{
    // Units are in N & m

    // Coordinates
    XYZPoint pt1(0, 0, 0); // Origin
    XYZPoint pt2(10, 0, 0);

    // Nodes
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    auto node1 = std::make_shared<Node>(1, pt1); nodes[node1->NodeIndex] = node1;
    auto node2 = std::make_shared<Node>(2, pt2); nodes[node2->NodeIndex] = node2;

    // Section
    auto area = 10;
    auto inertia11 = 100;
    auto inertia22 = 100;
    auto inertia12 = 100;
    auto sect = std::make_shared<Section>(area, inertia11, inertia22, inertia12);

    // Material
    auto mat = std::make_shared<Material>(1000, 0.3, 78500);

    // Members
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    bool isLumpedMassMatrix = true;
    auto beam = std::make_shared<FrameMember>(1, node1, node2, sect, mat, isLumpedMassMatrix); elements[beam->ElementIndex] = beam;

    // Restraints
    std::vector<bool> isRest = { true, true, true, true, true, true };
    std::vector<double> rest = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    std::vector<bool> secondIsRest = { false, false, true, true, true, false };

    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    auto res1 = std::make_shared<Restraint>(node1, isRest, rest); restraints[1] = res1;
    auto res2 = std::make_shared<Restraint>(node2, secondIsRest, rest); restraints[2] = res2;

    // Nodal loads
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    double nodalLoad[6] = { 0, 100, 0, 0, 0, 0 };
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
    XYZPoint bottomPt3(10, 0, 10);
    XYZPoint bottomPt4(0, 0, 10);
    XYZPoint topPt1(0, 10, 0); // Top of the origin
    XYZPoint topPt2(10, 10, 0);
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

    // Members
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    bool isLumpedMassMatrix = true;
    auto col1 = std::make_shared<FrameMember >(1, bottomNode1, topNode1, sect, mat, isLumpedMassMatrix); elements[1] = col1;
    auto col2 = std::make_shared<FrameMember >(2, bottomNode2, topNode2, sect, mat, isLumpedMassMatrix); elements[2] = col2;
    auto col3 = std::make_shared<FrameMember >(3, bottomNode3, topNode3, sect, mat, isLumpedMassMatrix); elements[3] = col3;
    auto col4 = std::make_shared<FrameMember >(4, bottomNode4, topNode4, sect, mat, isLumpedMassMatrix); elements[4] = col4;
    auto beam1 = std::make_shared<FrameMember >(5, topNode1, topNode2, sect, mat, isLumpedMassMatrix); elements[5] = beam1;
    auto beam2 = std::make_shared<FrameMember >(6, topNode2, topNode3, sect, mat, isLumpedMassMatrix); elements[6] = beam2;
    auto beam3 = std::make_shared<FrameMember >(7, topNode3, topNode4, sect, mat, isLumpedMassMatrix); elements[7] = beam3;
    auto beam4 = std::make_shared<FrameMember >(8, topNode4, topNode1, sect, mat, isLumpedMassMatrix); elements[8] = beam4;

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
        if (modalPeriods(i, 0))
            std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods(i, 0) << " s\n";

    return;
}

void CE583Sample()
{
    // Units are in N & m

    // Coordinates
    XYZPoint pt1(0, 0, 0); // Origin
    XYZPoint pt2(0, 3, 0);
    XYZPoint pt3(4, 3, 0);
    XYZPoint pt4(4, 0, 0);

    // Nodes
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    auto node1 = std::make_shared<Node>(1, pt1); nodes[node1->NodeIndex] = node1;
    auto node2 = std::make_shared<Node>(2, pt2); nodes[node2->NodeIndex] = node2;
    auto node3 = std::make_shared<Node>(3, pt3); nodes[node3->NodeIndex] = node3;
    auto node4 = std::make_shared<Node>(4, pt4); nodes[node4->NodeIndex] = node4;

    // Section
    auto sect1 = std::make_shared<Section>(0.02, 0.08, 0.08, 0.08);
    auto sect2 = std::make_shared<Section>(0.01, 0.01, 0.01, 0.01);

    // Material
    auto mat = std::make_shared<Material>(2e5, 0.3, 78500);

    // Members
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    bool isLumpedMassMatrix = true;
    auto elm1 = std::make_shared<FrameMember>(1, node1, node2, sect1, mat, isLumpedMassMatrix); elements[elm1->ElementIndex] = elm1;
    auto elm2 = std::make_shared<FrameMember>(2, node2, node3, sect1, mat, isLumpedMassMatrix); elements[elm2->ElementIndex] = elm2;
    auto elm3 = std::make_shared<FrameMember>(3, node4, node3, sect1, mat, isLumpedMassMatrix); elements[elm3->ElementIndex] = elm3;
    auto elm4 = std::make_shared<FrameMember>(4, node1, node3, sect2, mat, isLumpedMassMatrix); elements[elm4->ElementIndex] = elm4;

    // Restraints
    std::vector<double> zeros = { 0, 0, 0, 0, 0, 0 };
    std::vector<bool> node1RestConds = { true, true, true, true, true, false };
    std::vector<bool> node2RestConds = { false, false, true, true, true, false };
    std::vector<bool> node3RestConds = { false, false, true, true, true, false };
    std::vector<bool> node4RestConds = { false, true, true, true, true, false };

    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    auto rest1 = std::make_shared<Restraint>(node1, node1RestConds, zeros); restraints[1] = rest1;
    auto rest2 = std::make_shared<Restraint>(node2, node2RestConds, zeros); restraints[2] = rest2;
    auto rest3 = std::make_shared<Restraint>(node3, node3RestConds, zeros); restraints[3] = rest3;
    auto rest4 = std::make_shared<Restraint>(node4, node4RestConds, zeros); restraints[4] = rest4;

    // Nodal loads
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    double nodalLoad1[6] = { 10, -10, 0, 0, 0, 0 };
    double nodalLoad2[6] = { 10, -10, 0, 0, 0, 0 };
    auto nl1 = std::make_shared<NodalLoad>(node2, nodalLoad1); nodalLoads[1] = nl1;
    auto nl2 = std::make_shared<NodalLoad>(node3, nodalLoad2); nodalLoads[2] = nl2;

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
    auto modalPeriods = StructureSolver::GetModalPeriods(*str, SolverChoice::Armadillo);
    LOG("");
    for (size_t i = 0; i < modalPeriods.ColCount; i++)
        std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods(i, 0) << "\n";

    return;
}

void CE583Assignment1_3()
{
    // Units are N & m
    XYZPoint leftBottom(0.0, 0.0, 0.0);
    XYZPoint leftTop(0.0, 6.0, 0.0);
    XYZPoint rightTop(6.0, 6.0, 0.0);
    XYZPoint rightBottom(6.0, 0.0, 0.0);

    // Nodes
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    auto node1 = std::make_shared<Node>(1, leftBottom); nodes[node1->NodeIndex] = node1;
    auto node2 = std::make_shared<Node>(2, leftTop); nodes[node2->NodeIndex] = node2;
    auto node3 = std::make_shared<Node>(3, rightTop); nodes[node3->NodeIndex] = node3;
    auto node4 = std::make_shared<Node>(4, rightBottom); nodes[node4->NodeIndex] = node4;

    // Section
    auto area = 1;
    auto inertia11 = 1;
    auto inertia22 = 1;
    auto inertia12 = 1;
    auto sect1 = std::make_shared<Section>(area, inertia11, inertia22, inertia12);

    // Material
    auto mat = std::make_shared<Material>(1, 0, 78500);

    // Members
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    bool isLumpedMassMatrix = true;
    auto elm1 = std::make_shared<FrameMember>(1, node1, node2, sect1, mat, isLumpedMassMatrix); elements[elm1->ElementIndex] = elm1;
    auto elm2 = std::make_shared<FrameMember>(2, node2, node3, sect1, mat, isLumpedMassMatrix); elements[elm2->ElementIndex] = elm2;
    auto elm3 = std::make_shared<FrameMember>(3, node3, node4, sect1, mat, isLumpedMassMatrix); elements[elm3->ElementIndex] = elm3;

    // Restraints
    std::vector<double> zeros = { 0, 0, 0, 0, 0, 0 };
    std::vector<bool> node1RestConds = { true, true, true, true, true, true };
    std::vector<bool> node2RestConds = { false, false, true, true, true, false };
    std::vector<bool> node3RestConds = { false, false, true, true, true, false };
    std::vector<bool> node4RestConds = { false, true, true, true, true, false };

    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    auto rest1 = std::make_shared<Restraint>(node1, node1RestConds, zeros); restraints[1] = rest1;
    auto rest2 = std::make_shared<Restraint>(node2, node2RestConds, zeros); restraints[2] = rest2;
    auto rest3 = std::make_shared<Restraint>(node3, node3RestConds, zeros); restraints[3] = rest3;
    auto rest4 = std::make_shared<Restraint>(node4, node4RestConds, zeros); restraints[4] = rest4;

    // Nodal loads
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    double nodalLoad1[6] = { 0, -15000, 0, 0, 0, -15000 };
    double nodalLoad2[6] = { 0, -15000, 0, 0, 0, 15000 };
    auto nl1 = std::make_shared<NodalLoad>(node2, nodalLoad1); nodalLoads[1] = nl1;
    auto nl2 = std::make_shared<NodalLoad>(node3, nodalLoad2); nodalLoads[2] = nl2;

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
        LOG(" Node Index: " << node->NodeIndex << "\n");

        auto nodalDisps = StructureSolver::GetNodalDisplacements(*node, disps);

        for (size_t i = 0; i < 6; i++)
            std::cout << " DOF Index: " << i + 1 << ", Displacement = " << nodalDisps(i, 0) << "\n";
    }

    // Modal periods
    auto modalPeriods = StructureSolver::GetModalPeriods(*str, SolverChoice::Eigen);
    LOG("\n Modal periods:");
    for (size_t i = 0; i < modalPeriods.RowCount; i++)
        if (modalPeriods(i, 0) != 0)
            std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods(i, 0) << " s\n";

    // Find equivalent stiffness matrix for translation in x-direction for third node.
    auto equivalentStiffness = StructureSolver::CondenseStiffnessMatrixForSpecificDOF(*str, node3->DofIndexTX, SolverChoice::Eigen);
    auto equivalentForceVector = StructureSolver::CondenseForceVectorForSpecificDOF(*str, node3->DofIndexTX, SolverChoice::Eigen);
    auto dispFromEquivalentVals = equivalentForceVector / equivalentStiffness;

    LOG("");
    LOG("Equivalent translational stiffness: " << equivalentStiffness);
    LOG("Equivalent horizontal force: " << equivalentForceVector);
    LOG("Displacement calculated using equivalent value: " << dispFromEquivalentVals);
}

void CE583Assignment5()
{
    // Input Card (Units are in N & m)
    int nElmX = 8;
    int nElmY = 1;
    double lX = 4;
    double lY = 0.6;
    double thickness = 0.4;
    double e = 25e9;
    double v = 0;
    XYZPoint restraintStart(0.0, 0.0, 0.0);
    XYZPoint restraintEnd(0.0, 0.6, 0.0);

    // Solve
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Populate nodes
    auto hX = lX / nElmX;
    auto hY = lY / nElmY;
    auto nNodeX = nElmX + 1;
    auto nNodeY = nElmY + 1;

    // Create a vector representing support face
    Vector supportVector(restraintStart, restraintEnd);
    auto supportUnitVector = supportVector.getUnitVector();

    std::vector<bool> isRest = { true, true, true, true, true, true };
    std::vector<bool> universal = { false, false, true, true, true, true };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };
    auto addRestaint = [&](std::shared_ptr<Node> nR)
    {
        // Create a vector starting from the start of restraint to 
        // node. If unit vectors of new vector and support vector matches
        // and if length of new unit vector is smaller than support unit vector,
        // then the node is restrained. If not, assign universal restraint. Also, 
        // is given node is at the same location with start of restraint, it should
        // be restrained as well.

        Vector supportToNodeVector(restraintStart, nR->Coordinate);
        if ((supportUnitVector == supportToNodeVector.getUnitVector())
            && (supportToNodeVector.Length < (supportVector.Length + sqrt((hX * hX) + (hY * hY))))
            || (Utils::AreEqual(nR->Coordinate.X, restraintStart.X) && Utils::AreEqual(nR->Coordinate.Y, restraintStart.Y)))
        {
            auto restraint = std::make_shared<Restraint>(nR, isRest, rest);
            restraints[nR->NodeIndex] = restraint;
        }
        else
        {
            auto restraint = std::make_shared<Restraint>(nR, universal, rest);
            restraints[nR->NodeIndex] = restraint;
        }
    };

    int maxNodeIdx = 0;
    // Create nodes
    for (size_t i = 0; i < nNodeY; i++)
    {
        for (size_t j = 0; j < nNodeX; j++)
        {
            auto idx = (i * nNodeX) + j + 1;
            XYZPoint pt(j * hX, i * hY, 0);
            auto n = std::make_shared<Node>(idx, pt);
            nodes[n->NodeIndex] = n;
            addRestaint(n);
            maxNodeIdx = idx;
        }
    }

    // Create elements
    auto mt = std::make_shared<Material>(e, v, 0);
    auto isMembrane = true;
    std::vector<ShellMember> membersAtSupport;
    for (size_t i = 0; i < nElmY; i++)
    {
        for (size_t j = 0; j < nElmX; j++)
        {
            auto idx = (i * nElmX) + j + 1;
            auto& iNode = nodes[(i * nNodeX) + j + 1];
            auto& jNode = nodes[iNode->NodeIndex + 1];
            auto& kNode = nodes[jNode->NodeIndex + nNodeX];
            auto& lNode = nodes[kNode->NodeIndex - 1];
            auto sm = std::make_shared<ShellMember>(idx, iNode, jNode, kNode, lNode, mt, thickness, MembraneType::Bilinear, PlateType::NONE); elements[sm->ElementIndex] = sm;
            if (j == 0)
                membersAtSupport.push_back(*sm);
        }
    }

    // Nodal loads
    // Tip load is -20000 kN. Divide to tip nodes
    double nodalForce = static_cast<double>(-20000) / nNodeY;
    double nodalLoad[6] = { 0, nodalForce, 0, 0, 0, 0 };
    for (int i = nNodeX; i <= (nNodeX * nNodeY); i += nNodeX)
    {
        auto nl1 = std::make_shared<NodalLoad>(nodes[i], nodalLoad);
        nodalLoads[i] = nl1;
    }

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Armadillo);

    //for (auto& mem : membersAtSupport)
    //{
    //    for (size_t i = 0; i < 4; i += 3)
    //    {
    //        auto stresses = StructureSolver::CalculateMembraneNodalStresses(mem, disps, i + 1);
    //        //LOG(" Node Index:" << mem.Nodes[i]->NodeIndex << "\n");
    //        //LOG(" Node Index:" << mem.Nodes[i]->NodeIndex << "\n");
    //        //LOG(" Node Coordinate: " << mem.Nodes[i]->Coordinate.X << ", " << mem.Nodes[i]->Coordinate.Y);
    //        //LOG(" SigmaXX: " << stresses(0, 0));
    //        //LOG(" SigmaYY: " << stresses(1, 0));
    //        //LOG(" SigmaXY: " << stresses(2, 0));
    //        //LOG("");
    //        LOG(mem.Nodes[i]->Coordinate.Y << " " << stresses(2, 0) / 1000);
    //    }
    //}

    //for (auto& n : nodes)
    //{
    //    // Get top surface element
    //    if (Utils::AreEqual(n.second->Coordinate.Y, 0.6))
    //    {
    //        auto nodalDisp = StructureSolver::GetNodalDisplacements(*n.second, disps);
    //        LOG(n.second->Coordinate.X << " " << nodalDisp(1, 0));
    //    }
    //}

    auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[nodes.size()], disps);
    LOG(" Node Index: " << nodes[nodes.size()]->NodeIndex);
    LOG(" Node Location: " << nodes[nodes.size()]->Coordinate.X << " m, " << nodes[nodes.size()]->Coordinate.Y << " m");
    LOG(" Vertical Displacement: " << nodalDisp(1, 0) << " m");

    // auto& node = nodes[nNodeX * nNodeY];
    // LOG("");
    // LOG(" Node Index: ");
    // std::cout << " " << node->NodeIndex << "\n";
    // 
    // auto nodalDisps = StructureSolver::GetNodalDisplacements(*node, disps);
    // 
    // for (size_t i = 0; i < 6; i++)
    //     std::cout << " DOF Index: " << i + 1 << ", Displacement = " << nodalDisps(i, 0) << "\n";
}

void CE583Assignment6_2()
{
    // Input Card (Units are in N & m)
    int nElmX = 20;
    int nElmY = 5;
    bool isAxial = false;
    double stressLoc = 0;
    double lX = 4;
    double lY = 0.6;
    double thickness = 0.4;
    double e = 25e9;
    double v = 0;
    XYZPoint restraintStart(0.0, 0.0, 0.0);
    XYZPoint restraintEnd(0.0, 0.6, 0.0);

    // Solve
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Populate nodes
    auto hX = lX / nElmX;
    auto hY = lY / nElmY;
    auto nNodeX = nElmX + 1;
    auto nNodeY = nElmY + 1;

    // Create a vector representing support face
    Vector supportVector(restraintStart, restraintEnd);
    auto supportUnitVector = supportVector.getUnitVector();

    std::vector<bool> isRest = { true, true, true, true, true, true };
    std::vector<bool> universal = { false, false, true, true, true, true };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };
    auto addRestaint = [&](std::shared_ptr<Node> nR)
    {
        // Create a vector starting from the start of restraint to 
        // node. If unit vectors of new vector and support vector matches
        // and if length of new unit vector is smaller than support unit vector,
        // then the node is restrained. If not, assign universal restraint. Also, 
        // is given node is at the same location with start of restraint, it should
        // be restrained as well.

        Vector supportToNodeVector(restraintStart, nR->Coordinate);
        if ((supportUnitVector == supportToNodeVector.getUnitVector())
            && (supportToNodeVector.Length < (supportVector.Length + sqrt((hX * hX) + (hY * hY))))
            || (Utils::AreEqual(nR->Coordinate.X, restraintStart.X) && Utils::AreEqual(nR->Coordinate.Y, restraintStart.Y)))
        {
            auto restraint = std::make_shared<Restraint>(nR, isRest, rest);
            restraints[nR->NodeIndex] = restraint;
        }
        else
        {
            auto restraint = std::make_shared<Restraint>(nR, universal, rest);
            restraints[nR->NodeIndex] = restraint;
        }
    };

    int maxNodeIdx = 0;
    // Create nodes
    for (size_t i = 0; i < nNodeY; i++)
    {
        for (size_t j = 0; j < nNodeX; j++)
        {
            auto idx = (i * nNodeX) + j + 1;
            XYZPoint pt(j * hX, i * hY, 0);
            auto n = std::make_shared<Node>(idx, pt);
            nodes[n->NodeIndex] = n;
            addRestaint(n);
            maxNodeIdx = idx;
        }
    }

    // Create elements
    auto mt = std::make_shared<Material>(e, v, 0);
    auto isMembrane = true;
    std::vector<ShellMember> membersAtSupport;
    for (size_t i = 0; i < nElmY; i++)
    {
        for (size_t j = 0; j < nElmX; j++)
        {
            auto idx = (i * nElmX) + j + 1;
            auto& iNode = nodes[(i * nNodeX) + j + 1];
            auto& jNode = nodes[iNode->NodeIndex + 1];
            auto& kNode = nodes[jNode->NodeIndex + nNodeX];
            auto& lNode = nodes[kNode->NodeIndex - 1];
            auto sm = std::make_shared<ShellMember>(idx, iNode, jNode, kNode, lNode, mt, thickness, MembraneType::Incompatible, PlateType::NONE); elements[sm->ElementIndex] = sm;
            if (Utils::AreEqual(jNode->Coordinate.X, stressLoc) || Utils::AreEqual(iNode->Coordinate.X, stressLoc))
                membersAtSupport.push_back(*sm);
        }
    }

    // Nodal loads
    // Tip load is -20000 kN. Divide to tip nodes
    double nodalForce = static_cast<double>(-20000) / nNodeY;
    double nodalLoad[6] = { 0, nodalForce, 0, 0, 0, 0 };
    for (int i = nNodeX; i <= (nNodeX * nNodeY); i += nNodeX)
    {
        auto nl1 = std::make_shared<NodalLoad>(nodes[i], nodalLoad);
        nodalLoads[i] = nl1;
    }

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Armadillo);

    /*for (auto& mem : membersAtSupport)
    {
        auto stresses = StructureSolver::CalculateMembraneNodalStresses(mem, disps, 1);
        int printIndex = isAxial ? 0 : 2;
        LOG(mem.Nodes[0]->Coordinate.Y - 0.24 << " " << stresses(printIndex, 0) / 1000);
    }*/

    //for (auto& n : nodes)
    //{
    //    // Get top surface element
    //    if (Utils::AreEqual(n.second->Coordinate.Y, 0.6))
    //    {
    //        auto nodalDisp = StructureSolver::GetNodalDisplacements(*n.second, disps);
    //        LOG(n.second->Coordinate.X << " " << nodalDisp(1, 0));
    //    }
    //}

    auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[nodes.size()], disps);
    LOG(" Node Index: " << nodes[nodes.size()]->NodeIndex);
    LOG(" Node Location: " << nodes[nodes.size()]->Coordinate.X << " m, " << nodes[nodes.size()]->Coordinate.Y << " m");
    LOG(" Vertical Displacement: " << nodalDisp(1, 0) << " m");

    // auto& node = nodes[nNodeX * nNodeY];
    // LOG("");
    // LOG(" Node Index: ");
    // std::cout << " " << node->NodeIndex << "\n";
    // 
    // auto nodalDisps = StructureSolver::GetNodalDisplacements(*node, disps);
    // 
    // for (size_t i = 0; i < 6; i++)
    //     std::cout << " DOF Index: " << i + 1 << ", Displacement = " << nodalDisps(i, 0) << "\n";
}

void CE583Assignment6_3()
{
    bool isAxial = true;

    // Units are in N & m
    // Solve
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;
    double e = 25e9;
    double v = 0.0;
    double t = 0.4;

    // Coordinates and nodes
    auto addNode = [&](int nodeIndex, double x, double y, double z) {
        XYZPoint pt(x, y, z);
        nodes[nodeIndex] = std::make_shared<Node>(nodeIndex, pt);
    };

    addNode(1, 0, 0, 0);
    addNode(2, 0.4, 0, 0);
    addNode(3, 1.1, 0, 0);
    addNode(4, 1.4, 0, 0);
    addNode(5, 2.1, 0, 0);
    addNode(6, 2.4, 0, 0);
    addNode(7, 3.1, 0, 0);
    addNode(8, 3.4, 0, 0);
    addNode(9, 4, 0, 0);
    addNode(10, 0, 0.3, 0);
    addNode(11, 0.5, 0.3, 0);
    addNode(12, 1.0, 0.3, 0);
    addNode(13, 1.5, 0.3, 0);
    addNode(14, 2.0, 0.3, 0);
    addNode(15, 2.5, 0.3, 0);
    addNode(16, 3.0, 0.3, 0);
    addNode(17, 3.5, 0.3, 0);
    addNode(18, 4.0, 0.3, 0);
    addNode(19, 0.0, 0.6, 0);
    addNode(20, 0.6, 0.6, 0);
    addNode(21, 0.9, 0.6, 0);
    addNode(22, 1.6, 0.6, 0);
    addNode(23, 1.9, 0.6, 0);
    addNode(24, 2.6, 0.6, 0);
    addNode(25, 2.9, 0.6, 0);
    addNode(26, 3.6, 0.6, 0);
    addNode(27, 4.0, 0.6, 0);

    // Material
    auto mat = std::make_shared<Material>(e, v, 0);

    std::vector<ShellMember> membersAtSupport;
    std::vector<ShellMember> specMembers;
    // Elements
    auto addElem = [&](int elmIndex, int iNodeIndex, int jNodeIndex, int kNodeIndex, int lNodeIndex) {
        auto elm = std::make_shared<ShellMember>(elmIndex, nodes[iNodeIndex], nodes[jNodeIndex], nodes[kNodeIndex], nodes[lNodeIndex], mat, t, MembraneType::Incompatible, PlateType::NONE);
        elements[elmIndex] = elm;
        if (iNodeIndex == 1 || iNodeIndex == 10)
            membersAtSupport.push_back(*elm);
        if ((elmIndex == 1) || (elmIndex == 2) || (elmIndex == 9))
            specMembers.push_back(*elm);
    };

    addElem(1, 1, 2, 11, 10);
    addElem(2, 2, 3, 12, 11);
    addElem(3, 3, 4, 13, 12);
    addElem(4, 4, 5, 14, 13);
    addElem(5, 5, 6, 15, 14);
    addElem(6, 6, 7, 16, 15);
    addElem(7, 7, 8, 17, 16);
    addElem(8, 8, 9, 18, 17);
    addElem(9, 10, 11, 20, 19);
    addElem(10, 11, 12, 21, 20);
    addElem(11, 12, 13, 22, 21);
    addElem(12, 13, 14, 23, 22);
    addElem(13, 14, 15, 24, 23);
    addElem(14, 15, 16, 25, 24);
    addElem(15, 16, 17, 26, 25);
    addElem(16, 17, 18, 27, 26);

    // Restraints and nodal loads
    std::vector<bool> isRest = { true, true, true, true, true, true };
    std::vector<bool> universal = { false, false, true, true, true, true };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };
    auto f = -20000.0 / 3.0;
    double nodalLoad[6] = { 0, f, 0, 0, 0, 0 };
    for (auto& nPair : nodes)
    {
        if (Utils::AreEqual(nPair.second->Coordinate.X, 0))
            restraints[nPair.first] = std::make_shared<Restraint>(nPair.second, isRest, rest);
        else
            restraints[nPair.first] = std::make_shared<Restraint>(nPair.second, universal, rest);

        if (Utils::AreEqual(nPair.second->Coordinate.X, 4))
            nodalLoads[nPair.first] = std::make_shared<NodalLoad>(nPair.second, nodalLoad);
    }

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Armadillo);

    // Print tip displacement
    auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[nodes.size()], disps);
    LOG(" Tip Displacement");
    LOG(" Node Index: " << nodes[nodes.size()]->NodeIndex);
    LOG(" Node Location: " << nodes[nodes.size()]->Coordinate.X << " m, " << nodes[nodes.size()]->Coordinate.Y << " m");
    LOG(" Vertical Displacement: " << nodalDisp(1, 0) << " m");
    LOG("");

    // Axial
    LOG(" Axial Stresses due to Bending at Support");
    for (auto& mem : membersAtSupport)
    {
        int printIndex = 0; // : 2;
        auto stresses1 = StructureSolver::CalculateMembraneNodalStresses(mem, disps, 1);
        auto stresses2 = StructureSolver::CalculateMembraneNodalStresses(mem, disps, 4);
        LOG(" " << mem.Nodes[0]->Coordinate.Y << " " << stresses1(printIndex, 0) / 1000 << " kPa");
        LOG(" " << mem.Nodes[3]->Coordinate.Y << " " << stresses2(printIndex, 0) / 1000 << " kPa");
    }
    LOG("");

    // Shear
    LOG(" Shear Stresses at Support");
    for (auto& mem : membersAtSupport)
    {
        int printIndex = 2;
        auto stresses1 = StructureSolver::CalculateMembraneNodalStresses(mem, disps, 1);
        auto stresses2 = StructureSolver::CalculateMembraneNodalStresses(mem, disps, 4);
        LOG(" " << mem.Nodes[0]->Coordinate.Y << " " << stresses1(printIndex, 0) / 1000 << " kPa");
        LOG(" " << mem.Nodes[3]->Coordinate.Y << " " << stresses2(printIndex, 0) / 1000 << " kPa");
    }

    // First edge
    // Calculate for the first element
    bool isX = false;
    int elemeIndex = 9;
    double eta = -1;
    for (double ksi = -1; ksi <= 1; ksi += 0.1)
    {
        ShellMember elm;
        for (auto& elem : specMembers)
        {
            if (elem.ElementIndex == elemeIndex)
            {
                elm = elem;
                break;
            }
        }

        auto& k1212 = *elm.K1212;
        auto kCR = k1212.getSubmatrix(8, 11, 0, 7);
        auto kCC = k1212.getSubmatrix(8, 11, 8, 11);
        auto invKCC = elm.InvertMatrix4(kCC);
        auto minvKCC = invKCC * -1;

        auto iNodeDisp = StructureSolver::GetNodalDisplacements(*elm.Nodes[0], disps);
        auto jNodeDisp = StructureSolver::GetNodalDisplacements(*elm.Nodes[1], disps);
        auto kNodeDisp = StructureSolver::GetNodalDisplacements(*elm.Nodes[2], disps);
        auto lNodeDisp = StructureSolver::GetNodalDisplacements(*elm.Nodes[3], disps);

        Matrix<double> nodalDispVector(8, 1);
        nodalDispVector(0, 0) = iNodeDisp(0, 0);
        nodalDispVector(1, 0) = iNodeDisp(1, 0);
        nodalDispVector(2, 0) = jNodeDisp(0, 0);
        nodalDispVector(3, 0) = jNodeDisp(1, 0);
        nodalDispVector(4, 0) = kNodeDisp(0, 0);
        nodalDispVector(5, 0) = kNodeDisp(1, 0);
        nodalDispVector(6, 0) = lNodeDisp(0, 0);
        nodalDispVector(7, 0) = lNodeDisp(1, 0);

        auto generalizedDisps = minvKCC * kCR * nodalDispVector;

        Matrix<double> dispVector(12, 1);
        for (size_t i = 0; i < 8; i++)
            dispVector(i, 0) = nodalDispVector(i, 0);
        for (size_t i = 8; i < 11; i++)
            dispVector(i, 0) = generalizedDisps(i - 8, 0);

        auto n1 = 0.25 * (ksi - 1) * (eta - 1);
        auto n2 = -0.25 * (ksi + 1) * (eta - 1);
        auto n3 = 0.25 * (ksi + 1) * (eta + 1);
        auto n4 = -0.25 * (ksi - 1) * (eta + 1);
        auto mult1 = 1 - (ksi * ksi);
        auto mult2 = 1 - (eta * eta);

        double disp = 0;
        if (isX)
            disp = (n1 * dispVector(0, 0)) + (n2 * dispVector(2, 0)) + (n3 * dispVector(4, 0)) + (n4 * dispVector(6, 0)) + (mult1 * generalizedDisps(0, 0)) + (mult2 * generalizedDisps(1, 0));
        else
            disp = (n1 * dispVector(1, 0)) + (n2 * dispVector(3, 0)) + (n3 * dispVector(5, 0)) + (n4 * dispVector(7, 0)) + (mult1 * generalizedDisps(2, 0)) + (mult2 * generalizedDisps(3, 0));

        LOG(ksi << " " << disp * 1000);
    }



}

void TrussExample()
{
    // Units are N & m
    XYZPoint left(0.0, 0.0, 0.0);
    XYZPoint right(10.0, 0.0, 0.0);
    XYZPoint top(5.0, 5.0, 0.0);

    // Nodes
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    auto node1 = std::make_shared<Node>(1, left); nodes[node1->NodeIndex] = node1;
    auto node2 = std::make_shared<Node>(2, right); nodes[node2->NodeIndex] = node2;
    auto node3 = std::make_shared<Node>(3, top); nodes[node3->NodeIndex] = node3;

    // Section
    auto area = 0.16;
    auto inertia11 = 2.13333333333333E-03;
    auto inertia22 = 2.13333333333333E-03;
    auto inertia12 = 0.0036096;
    auto sect1 = std::make_shared<Section>(area, inertia11, inertia22, inertia12);

    // Material
    auto mat = std::make_shared<Material>(200000000000, 0.3, 78500);

    // Members
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    bool isLumpedMassMatrix = true;
    auto elm1 = std::make_shared<TrussMember>(1, node1, node2, sect1, mat, isLumpedMassMatrix); elements[elm1->ElementIndex] = elm1;
    auto elm2 = std::make_shared<TrussMember>(2, node2, node3, sect1, mat, isLumpedMassMatrix); elements[elm2->ElementIndex] = elm2;
    auto elm3 = std::make_shared<TrussMember>(3, node3, node1, sect1, mat, isLumpedMassMatrix); elements[elm3->ElementIndex] = elm3;

    // Restraints
    std::vector<double> zeros = { 0, 0, 0, 0, 0, 0 };
    std::vector<bool> node1RestConds = { true, true, true, false, false, false };
    std::vector<bool> node2RestConds = { false, true, true, false, false, false };

    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    auto rest1 = std::make_shared<Restraint>(node1, node1RestConds, zeros); restraints[1] = rest1;
    auto rest2 = std::make_shared<Restraint>(node2, node2RestConds, zeros); restraints[2] = rest2;

    // Nodal loads
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    double nodalLoad1[6] = { 0, -5e6, 0, 0, 0, 0 };
    auto nl1 = std::make_shared<NodalLoad>(node3, nodalLoad1); nodalLoads[1] = nl1;

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
        LOG(" Node Index: " << node->NodeIndex << "\n");

        auto nodalDisps = StructureSolver::GetNodalDisplacements(*node, disps);

        for (size_t i = 0; i < 6; i++)
            std::cout << " DOF Index: " << i + 1 << ", Displacement = " << nodalDisps(i, 0) << "\n";
    }
}
