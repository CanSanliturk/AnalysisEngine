#include <chrono>
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include "Shape.h"
#include "Vector.h"
#include "UtilMethods.h"
#include "StructureSolver.h"

#define LOG(x) std::cout << x << "\n"
#define LOGW(x) std::wcout << x << "\n"
#define LOGIF(x, cond) if (cond) LOG(x)

void CantileverBeam_PlateModel();

int main()
{
    LOG(" __________________________________________________");
    LOG(" |                                                |");
    LOG(" |  3-Dimensional Finite Element Analysis Engine  |");
    LOG(" |________________________________________________|");
    LOG("");

    // Start timer
    auto timenow =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    // Call test function (Later on, these guys will be moved to a unit test project)
    try
    {
        CantileverBeam_PlateModel();
        LOG(" Analysis completed without errors...");
    }
    catch (const std::runtime_error& e)
    {
        LOG(" Analysis cannot be completed...");
        LOG(" " << e.what());
    }

    // Log duration
    auto timenow2 =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    LOG(" Elapsed Time = " << timenow2 - timenow << " seconds\n");

    return 0;
}

void CantileverBeam_PlateModel()
{
    // Units are in N & m
    int nElmX = 10;
    int nElmY = 40;
    double lX = 0.4;
    double lY = 4.0;
    double e = 25e9;
    double v = 0.2;
    double t = 0.6;
    double verticalForce = -20000.0;

    // Solve
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Material
    auto mat = std::make_shared<Material>(e, v, 0);

    // Helper functions
    // Coordinates, nodes, restraints and loads
    std::vector<bool> isRest = { true, true, true, true, true, true };
    std::vector<bool> universal = { false, false, false, false, false, false };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };
    double f = verticalForce / ((double)nElmX + 1);
    double nodalLoad[6] = { 0, 0, f, 0, 0, 0 };
    int loadIdx = 1;
    std::vector<int> tipNodeIndices;
    int tipMidNodeIndex = -1;
    auto addNode = [&](int nodeIndex, double x, double y, double z) {
        XYZPoint pt(x, y, z);
        nodes[nodeIndex] = std::make_shared<Node>(nodeIndex, pt);

        // Restraint
        if (Utils::AreEqual(y, 0))
            restraints[nodeIndex] = std::make_shared<Restraint>(nodes[nodeIndex], isRest, rest);
        else
            restraints[nodeIndex] = std::make_shared<Restraint>(nodes[nodeIndex], universal, rest);

        // Loads
        if (Utils::AreEqual(y, lY))
        {
            if (Utils::AreEqual(x, lX / 2.0))
                tipMidNodeIndex = nodeIndex;
            nodalLoads[loadIdx++] = std::make_shared<NodalLoad>(nodes[nodeIndex], nodalLoad);
            tipNodeIndices.push_back(nodeIndex);
        }
    };

    // Elements
    std::vector<std::shared_ptr<ShellMember>> membersAtSupport;
    auto addElem = [&](int elmIndex, int iNodeIndex, int jNodeIndex, int kNodeIndex, int lNodeIndex) {
        auto elm = std::make_shared<ShellMember>(elmIndex, nodes[iNodeIndex], nodes[jNodeIndex], nodes[kNodeIndex], nodes[lNodeIndex], mat, t, MembraneType::NONE, PlateType::MindlinFourNode);
        elements[elmIndex] = elm;

        if (elm->Nodes[0]->Coordinate.Y == 0)
            membersAtSupport.push_back(elm);
    };

    // Add nodes
    auto hX = lX / nElmX;
    auto hY = lY / nElmY;
    auto nNodeX = nElmX + 1;
    auto nNodeY = nElmY + 1;

    for (size_t i = 0; i < nNodeY; i++)
    {
        for (size_t j = 0; j < nNodeX; j++)
        {
            auto idx = (i * nNodeX) + j + 1;
            addNode(idx, j * hX, i * hY, 0);
        }
    }

    // Add elements. Nodes should be in CCW-orientation
    for (size_t i = 0; i < nElmY; i++)
    {
        for (size_t j = 0; j < nElmX; j++)
        {
            auto idx = (i * nElmX) + j + 1;
            auto iNodeIdx = (i * nNodeX) + j + 1;
            auto jNodeIdx = iNodeIdx + 1;
            auto kNodeIdx = jNodeIdx + nNodeX;
            auto lNodeIdx = kNodeIdx - 1;
            addElem(idx, iNodeIdx, jNodeIdx, kNodeIdx, lNodeIdx);
            idx++;
        }
    }

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Armadillo);

    // Print tip displacement
    auto checkTipMid = tipMidNodeIndex != -1 ? tipMidNodeIndex : tipNodeIndices[0];
    auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[checkTipMid], disps);


    LOG(" Tip Displacement");
    LOG(" Node Index: " << nodes[checkTipMid]->NodeIndex);
    LOG(" Node Location: " << nodes[checkTipMid]->Coordinate.X << " m, " << nodes[checkTipMid]->Coordinate.Y << " m");
    LOG(" Vertical Displacement: " << nodalDisp(2, 0) << " m");
    LOG("");
}
