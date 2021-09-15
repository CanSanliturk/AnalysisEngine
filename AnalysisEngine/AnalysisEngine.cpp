#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <chrono>
#include "StructureSolver.h"
#include "Vector.h"
#include "UtilMethods.h"
#include "Shape.h"

#define LOG(x) std::cout << x << "\n"
#define LOGW(x) std::wcout << x << "\n"
#define LOGIF(x, cond) if (cond) LOG(x)
constexpr double pi = 3.141592653589793;

void StrutTieDesign();
void StrutTieDesignWithLattice();
void StrutTieDesignWithLatticeAndEnergyConservation();
void CantileverDisplacements3D();
void CantileverDisplacements2D();
void TableModalAnalysis();
void TrussExample();
void CE583Sample();
void CE583Assignment1_3();
void CE583Assignment5();
void CE583Assignment6_2();
void CE583Assignment6_3();
void CE583Assignment7_2();
void CE583Assignment7_3();
void CE583Assignment8_2_MembraneAction();
void CE583Assignment8_2_PlateAction();
void CE583Assignment8_3();
void CE583Assignment9_2_1();
void CE583Assignment9_2_3();
void CE583Assignment9_3();

struct NodeData {
    int nodeIndex;
    XYZPoint coordinate;
};

struct TrussData {
    int trussIndex;
    NodeData* iNode;
    NodeData* jNode;
};

int main()
{
    LOG(" __________________________________________________");
    LOG(" |                                                |");
    LOG(" |  3-Dimensional Finite Element Analysis Engine  |");
    LOG(" |          Created by Can Sanliturk              |");
    LOG(" |________________________________________________|");
    LOG("");

    // Start timer
    auto timenow =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    // Call test function (Later on, these guys will be moved to a unit test project)
    try
    {
        StrutTieDesignWithLatticeAndEnergyConservation();
        LOG("\n Analysis completed without errors...");
    }
    catch (const std::runtime_error& e)
    {
        LOG("\n Analysis cannot be completed...");
        LOG(" " << e.what());
    }

    // Log duration
    auto timenow2 =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    LOG(" Elapsed Time = " << timenow2 - timenow << " seconds\n");

    std::cin.get();
    return 0;
}

void StrutTieDesign()
{
    // INPUT CARD
    // Length: m, Force: N
    // Surface Outer Dimensions
    auto lx = 1.0; // length in x
    auto ly = 3.0; // length in y
    auto thickness = 0.5;
    auto eMod = 30e9; // youngs modulus (xe9 means x gigapascals)
    auto v = 0.3; // poissons ratio
    auto fc = 30e6; // concrete compressive strength
    auto fy = 420e6; // steel yield strength
    auto meshSize = 0.1;
    auto merger = 5;
    double performanceRatioCriteria = 0.05; // minimum performance ratio of elements to be considered
    auto maxIterationCount = 10000; // max iterations for topology optimization
    auto iterationStoppingCount = 30;
    auto solverSelection = SolverChoice::Eigen; // library to be used at linear algebraic equation solvings
    auto isPrintMeshInfo = false;
    auto isPrintForceInfo = true;
    auto isPrintElements = true;
    auto isHorizontalLoad = true;

    // Necessary fields
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Helper function to add restraints
    std::vector<bool> pin = { true, true, true, true, true, true };
    std::vector<bool> roller = { false, true, true, true, true, true };
    std::vector<bool> universal = { false, false, true, true, true, true };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };
    std::vector<int> realRestrainedNodes;
    auto addRestraint = [&](std::shared_ptr<Node> nR)
    {
        if (Utils::AreEqual(nR->Coordinate.X, 0.0, meshSize * 0.1) && Utils::AreEqual(nR->Coordinate.Y, 0.0, meshSize * 0.1))
        {
            restraints[nR->NodeIndex] = std::make_shared<Restraint>(nR, pin, rest);
            realRestrainedNodes.push_back(nR->NodeIndex);
        }
        else if (Utils::AreEqual(nR->Coordinate.X, lx, meshSize * 0.1) && Utils::AreEqual(nR->Coordinate.Y, 0.0, meshSize * 0.1))
        {
            restraints[nR->NodeIndex] = std::make_shared<Restraint>(nR, pin, rest);
            realRestrainedNodes.push_back(nR->NodeIndex);
        }
        else
            restraints[nR->NodeIndex] = std::make_shared<Restraint>(nR, universal, rest);
    };

    // Populate nodes
    int nodeIndex = 1;
    int nNodeX = ((int)(lx / meshSize)) + 1;
    int nNodeY = ((int)(ly / meshSize)) + 1;
    LOGIF(" NODES (Node Index, X-Coordinate, Y-Coordinate)", isPrintMeshInfo);
    for (auto yCoord = 0.0; yCoord < ly + meshSize * 0.1; yCoord += meshSize)
    {
        for (auto xCoord = 0.0; xCoord < lx + meshSize * 0.1; xCoord += meshSize)
        {
            LOGIF(" " << nodeIndex << ", " << xCoord << ", " << yCoord, isPrintMeshInfo);
            XYZPoint pt(xCoord, yCoord, 0);
            auto n = std::make_shared<Node>(nodeIndex++, pt);
            nodes[n->NodeIndex] = n;
            addRestraint(n);
        }
    }

    // Populate elements
    // Initially, populate horizontal members
    int memberIndex = 1;
    auto&& mat = std::make_shared<Material>(eMod, v, 0);
    auto&& sect = std::make_shared<Section>(0.01, 0, 0, 0);
    LOGIF("\n ELEMENT-NODE CONNECTIVITY (Member Index, i-Node Index, j-Node Index)", isPrintMeshInfo);
    for (int i = 1; i <= nNodeY; i++)
    {
        for (int j = 1; j < nNodeX; j++)
        {
            auto index = j + (nNodeX * (i - 1));
            auto& elmINode = nodes[index];
            auto& elmJNode = nodes[index + 1];
            LOGIF(" Horizontal Member, Member Index: " << memberIndex << ", " << elmINode->NodeIndex << ", " << elmJNode->NodeIndex, isPrintMeshInfo);
            elements[memberIndex] = std::make_shared<TrussMember>(memberIndex, elmINode, elmJNode, sect, mat);
            memberIndex++;
        }
    }

    // Then, populate vertical members
    for (size_t i = 1; i < nNodeY; i++)
    {
        for (size_t j = 1; j <= nNodeX; j++)
        {
            auto iIdx = j + (nNodeX * (i - 1));
            auto jIdx = iIdx + nNodeX;
            LOGIF(" Vertical Member, Memmber Index: " << memberIndex << ", " << iIdx << ", " << jIdx, isPrintMeshInfo);
            elements[memberIndex] = std::make_shared<TrussMember>(memberIndex, nodes[iIdx], nodes[jIdx], sect, mat);
            memberIndex++;
        }
    }

    // Finally, populate diagonal members
    for (size_t i = 1; i < nNodeY; i++)
    {
        for (size_t j = 1; j < nNodeX; j++)
        {
            // Two cross-truss will be added at this stage
            auto lowerLeftNodeIndex = j + (nNodeX * (i - 1));
            auto lowerRightNodeIndex = lowerLeftNodeIndex + 1;
            auto upperLeftNodeIndex = lowerLeftNodeIndex + nNodeX;
            auto upperRightNodeIndex = upperLeftNodeIndex + 1;

            LOGIF(" Diagonal Member, Member Index: " << memberIndex << ", " << lowerLeftNodeIndex << ", " << upperRightNodeIndex, isPrintMeshInfo);
            elements[memberIndex] = std::make_shared<TrussMember>(memberIndex, nodes[lowerLeftNodeIndex], nodes[upperRightNodeIndex], sect, mat);
            memberIndex++;

            LOGIF(" Diagonal Member, Member Index: " << memberIndex << ", " << lowerRightNodeIndex << ", " << upperLeftNodeIndex, isPrintMeshInfo);
            elements[memberIndex] = std::make_shared<TrussMember>(memberIndex, nodes[lowerRightNodeIndex], nodes[upperLeftNodeIndex], sect, mat);
            memberIndex++;
        }
    }

    LOGIF(" Total member count: " << elements.size(), isPrintMeshInfo);

    // Force vector
    for (auto& n : nodes)
    {
        if ((Utils::AreEqual(n.second->Coordinate.X, lx / 2.0, meshSize * 0.1)) && (Utils::AreEqual(n.second->Coordinate.Y, ly, meshSize * 0.1)))
        {
            if (!isHorizontalLoad)
            {
                double load[] = { 0, -5000000000000, 0, 0, 0, 0 };
                nodalLoads[1] = std::make_shared<NodalLoad>(n.second, load);
            }
            else
            {
                double load[] = { 5000000000000, 0, 0, 0, 0, 0 };
                nodalLoads[1] = std::make_shared<NodalLoad>(n.second, load);
            }
            break;
        }
    }

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);
    LOGIF(" Total DOF: " << str->nUnrestrainedDOF, false);
    auto fVec = str->getForceVector(&nodalLoads);

    if (isPrintElements)
    {
        std::ofstream elementOutputStream;

        elementOutputStream.open("plot\\elements.dat");
        elementOutputStream << "BEGIN SCENE\n";
        for (auto&& e : elements)
        {
            auto&& firstNode = e.second->GelElementNodes().at(0);
            auto&& secondNode = e.second->GelElementNodes().at(1);
            elementOutputStream << firstNode->Coordinate.X << " " << firstNode->Coordinate.Y << "\n";
            elementOutputStream << secondNode->Coordinate.X << " " << secondNode->Coordinate.Y << "\n";
            elementOutputStream << "7\n";
        }

        elementOutputStream << "END SCENE\n";
        elementOutputStream.close();
    }

    // Start iterations
    auto comp = 0.0, tens = 0.0;
    auto prevNumOfElmStsfyCrt = 0;
    auto unchangedIterCount = 0;

    std::ofstream animationStream;
    animationStream.open("plot\\animation.dat");

    for (int i = 0; i < maxIterationCount; i++)
    {
        auto currNumOfElmStsfyCrt = 0;
        LOGIF(" Iteration #" << i + 1, true);
        // Fields to be used for storing extreme forces on members
        auto maxCompression = 0.0, maxTension = 0.0;
        str->updateStiffnessMatrix();

        // Perform initial analysis
        //auto disps = StructureSolver::CalculateDisplacements(*(str->StiffnessMatrix), fVec, str->nDOF, str->nUnrestrainedDOF, solverSelection);
        auto disps = StructureSolver::GetDisplacementForStaticCase(*str, solverSelection);

        auto leftTopNode = str->getNodeAt(0, ly, 0);
        auto rightTopNode = str->getNodeAt(lx, ly, 0);

        auto leftTopDofTX = leftTopNode->DofIndexTX;
        auto leftTopDofTY = leftTopNode->DofIndexTY;
        auto rightTopDofTX = rightTopNode->DofIndexTX;
        auto rightTopDofTY = rightTopNode->DofIndexTY;

        auto leftTopTX = disps(leftTopDofTX - 1, 0);
        auto leftTopTY = disps(leftTopDofTY - 1, 0);
        auto rightTopTX = disps(rightTopDofTX - 1, 0);
        auto rightTopTY = disps(rightTopDofTY - 1, 0);

        // Find internal forces on elements
        for (auto& elm : elements)
        {
            auto trMem = dynamic_cast<TrussMember*>(&*elm.second);
            auto axialForce = trMem->getAxialForce(disps);
            if (axialForce < maxCompression) maxCompression = axialForce;
            if (maxTension < axialForce) maxTension = axialForce;
            LOGIF(" Comp: " << maxCompression << ", Tens:" << maxTension, false);
        }

        animationStream << "BEGIN SCENE\n";
        // Modify elements
        for (auto& elm : elements)
        {
            auto trMem = dynamic_cast<TrussMember*>(&*elm.second);
            auto axialForce = trMem->getAxialForce(disps);
            auto isCompression = axialForce < 0;
            auto isTension = 0 < axialForce;
            auto ratio = axialForce / (axialForce < 0 ? maxCompression : maxTension);
            if (performanceRatioCriteria < ratio) currNumOfElmStsfyCrt++;
            ratio = pow(ratio, 1.0 / (i + 1));
            elm.second->updateStiffness(ratio);

            int colour = 0;
            if (isCompression)
            {
                if (ratio < 0.2) colour = 2;
                else if (ratio < 0.4) colour = 3;
                else if (ratio < 0.6) colour = 4;
                else if (ratio < 0.8) colour = 5;
                else if (ratio <= 1.0) colour = 6;
            }
            else if (isTension)
            {
                if (ratio < 0.2) colour = 7;
                else if (ratio < 0.4) colour = 8;
                else if (ratio < 0.6) colour = 9;
                else if (ratio < 0.8) colour = 10;
                else if (ratio <= 1.0) colour = 11;
            }
            else
            {
                colour = 0;
            }

            animationStream << trMem->Nodes[0]->Coordinate.X << " " << trMem->Nodes[0]->Coordinate.Y << "\n";
            animationStream << trMem->Nodes[1]->Coordinate.X << " " << trMem->Nodes[1]->Coordinate.Y << "\n";
            animationStream << colour << "\n";
        }

        animationStream << "END SCENE\n";

        comp = maxCompression;
        tens = maxTension;

        if (currNumOfElmStsfyCrt == prevNumOfElmStsfyCrt)
            unchangedIterCount++;

        prevNumOfElmStsfyCrt = currNumOfElmStsfyCrt;

        if (unchangedIterCount == iterationStoppingCount)
            break;
    }
    animationStream.close();

    LOG(" Iterations for optimization finished.");
    auto disps = StructureSolver::CalculateDisplacements(*(str->StiffnessMatrix), fVec, str->nDOF, str->nUnrestrainedDOF, solverSelection);

    // Identify nodes
    std::vector<int> nodeIndices;
    for (auto& nodePair : nodes)
    {
        // Get node
        auto& node = nodePair.second;

        // If there is a boundary condition at a node (either natural or essential), identify node 
        // as strut&tie system node and continue
        bool isStop = false;
        for (auto& res : realRestrainedNodes)
        {
            if (res == node->NodeIndex)
            {
                nodeIndices.push_back(node->NodeIndex);
                isStop = true;
                break;
            }
        }

        if (isStop) continue;

        for (auto& nLoad : nodalLoads)
        {
            if (nLoad.second->ActingNode->NodeIndex == node->NodeIndex)
            {
                nodeIndices.push_back(node->NodeIndex);
                isStop = true;
                break;
            }
        }

        if (isStop) continue;

        // Get list of indices of connected elements
        auto& listOfConnectedElements = node->ConnectedElements;

        // Calculate axial forces on connected elements. If their ratio are higher than performance criteria,
        // store them
        std::vector<int> elementsSatisfyPerformanceCriteria;
        for (auto& elmIndex : listOfConnectedElements)
        {
            auto trussElm = dynamic_cast<TrussMember*>(&*elements[elmIndex]);
            auto axialForce = trussElm->getAxialForce(disps);
            auto ratio = axialForce / (axialForce < 0 ? comp : tens);
            if (performanceRatioCriteria <= ratio)
                elementsSatisfyPerformanceCriteria.push_back(elmIndex);
        }

        // If there are at least three load carrying elements connecting to node, identify node as 
        // strut&tie system cantidate node and continue
        if (elementsSatisfyPerformanceCriteria.size() > 2)
        {
            nodeIndices.push_back(node->NodeIndex);
            continue;
        }
        // Else if there are two elements connecting to the node and if they are not on the same line,
        // identify node as strut&tie candidate node and continue
        else if (elementsSatisfyPerformanceCriteria.size() == 2)
        {
            auto firstElm = dynamic_cast<TrussMember*>(&*elements[elementsSatisfyPerformanceCriteria[0]]);
            auto secondElm = dynamic_cast<TrussMember*>(&*elements[elementsSatisfyPerformanceCriteria[1]]);
            Vector firstElmVector(firstElm->Nodes[0]->Coordinate, firstElm->Nodes[1]->Coordinate);
            Vector secondElmVector(secondElm->Nodes[0]->Coordinate, secondElm->Nodes[1]->Coordinate);

            auto firstElmUnitVec = firstElmVector.getUnitVector();
            auto secondElmUnitVec = secondElmVector.getUnitVector();
            auto asdasd = abs(firstElmUnitVec.dotProduct(secondElmUnitVec));
            // If dot-product of element vectors results in -1, they are on the same line.
            if (!Utils::AreEqual(abs(firstElmUnitVec.dotProduct(secondElmUnitVec)), 1))
                nodeIndices.push_back(node->NodeIndex);
        }
    }

    // After obtaining nodes, start merging nodes that are too close
    // Create node data structs
    std::vector<NodeData> nodeDataVector;
    for (size_t i = 0; i < nodeIndices.size(); i++)
    {
        NodeData nd;
        nd.nodeIndex = i;
        nd.coordinate = nodes[nodeIndices[i]]->Coordinate;
        nodeDataVector.push_back(nd);
    }

    LOG("\n STRUT AND TIE SYSTEM INITIAL NODES");
    for (auto& nd : nodeDataVector)
        LOG(" Node Index: " << (nd.nodeIndex + 1) << ", Coordinates(x,y): (" << nd.coordinate.X << ", " << nd.coordinate.Y << ")");

    // Merge the nodes that are close to each other
    std::vector<NodeData> mergedNodes;
    std::vector<NodeData> nodesToBeSkipped;
    auto nodeIndexer = 1;
    for (auto& nd : nodeDataVector)
    {
        // Skip the already encountered nodes
        bool isSkipOutside = false;
        for (auto& nskp : nodesToBeSkipped)
        {
            if (nskp.nodeIndex == nd.nodeIndex)
            {
                isSkipOutside = true;
                break;
            }
        }

        if (isSkipOutside)
            continue;

        std::vector<NodeData> nodesWillBeMerged;
        bool isFirstEncountering = true;
        for (auto& nnd : nodeDataVector)
        {
            // Skip the same nodes
            if (nd.nodeIndex == nnd.nodeIndex)
                continue;

            // Skip the already encountered nodes
            bool isSkipInside = false;
            for (auto& nskp : nodesToBeSkipped)
            {
                if (nskp.nodeIndex == nnd.nodeIndex)
                {
                    isSkipInside = true;
                    break;
                }
            }

            if (isSkipInside)
                continue;

            auto dist = nd.coordinate.DistanceTo(nnd.coordinate);
            if (dist < merger * meshSize)
            {
                if (isFirstEncountering)
                {
                    nodesWillBeMerged.push_back(nd);
                    isFirstEncountering = false;
                }

                nodesWillBeMerged.push_back(nnd);
            }
        }

        for (auto& nwbm : nodesWillBeMerged)
        {
            bool isInsert = true;

            for (auto& ntbs : nodesToBeSkipped)
            {
                if (nwbm.nodeIndex == ntbs.nodeIndex)
                {
                    isInsert = false;
                    break;
                }
            }

            if (isInsert)
                nodesToBeSkipped.push_back(nwbm);
        }

        if (1 < nodesWillBeMerged.size())
        {
            auto sumXCoord = 0.0, sumYCoord = 0.0;
            for (auto& nid : nodesWillBeMerged)
            {
                sumXCoord += nid.coordinate.X;
                sumYCoord += nid.coordinate.Y;
            }

            auto xCd = sumXCoord / (double(nodesWillBeMerged.size()));
            auto yCd = sumYCoord / (double(nodesWillBeMerged.size()));
            XYZPoint pt(xCd, yCd, 0.0);
            NodeData newNode;
            newNode.nodeIndex = nodeIndexer;
            newNode.coordinate = pt;
            mergedNodes.push_back(newNode);
            nodeIndexer++;
        }
        else
        {
            NodeData newNode;
            newNode.nodeIndex = nodeIndexer;
            newNode.coordinate = nd.coordinate;
            mergedNodes.push_back(newNode);
            nodeIndexer++;
        }
    }

    LOG("\n STRUT AND TIE SYSTEM MERGED NODES");
    for (auto& nd : mergedNodes)
        LOG(" Node Index: " << nd.nodeIndex << ", Coordinates(x,y): (" << nd.coordinate.X << ", " << nd.coordinate.Y << ")");

    // Since nodes of the system are found, create truss analogy by linking all the nodes to each other
    std::vector<TrussData> trussDataVector;
    int trussIndexer = 1;
    LOG("\n STRUT AND TIE SYSTEM TRUSSES");
    std::ofstream trussDataPrinter;

    if (isPrintElements)
    {
        trussDataPrinter.open("plot\\strutTieSystem.dat");
        trussDataPrinter << "BEGIN SCENE\n";
    }

    for (size_t i = 0; i < mergedNodes.size(); i++)
    {
        auto iNode = mergedNodes[i];
        for (size_t j = i + 1; j < mergedNodes.size(); j++)
        {
            auto jNode = mergedNodes[j];
            TrussData t;
            t.trussIndex = trussIndexer;
            t.iNode = &iNode;
            t.jNode = &jNode;
            trussDataVector.push_back(t);
            trussIndexer++;
            LOG(" Truss Index: " << t.trussIndex << ", i-End Node Index: " << t.iNode->nodeIndex << ", j-End Node Index: " << t.jNode->nodeIndex);
            if (isPrintElements)
            {
                trussDataPrinter << iNode.coordinate.X << " " << iNode.coordinate.Y << "\n";
                trussDataPrinter << jNode.coordinate.X << " " << jNode.coordinate.Y << "\n";
                trussDataPrinter << "7\n";
            }
        }

    }

    if (isPrintElements)
    {
        trussDataPrinter << "END SCENE";
        trussDataPrinter.close();
    }
}

void StrutTieDesignWithLattice()
{
    // INPUT CARD
    // Length: m, Force: N
    // Surface Outer Dimensions
    auto lx = 1.0; // length in x
    auto ly = 1.0; // length in y
    auto thickness = 0.5;
    auto eMod = 30e9; // youngs modulus (xe9 means x gigapascals)
    auto v = 0.3; // poissons ratio
    auto fc = 30e6; // concrete compressive strength
    auto fy = 420e6; // steel yield strength
    auto meshSize = 0.05;
    auto horizon = 3.5 * meshSize;
    auto merger = 3.5;
    double performanceRatioCriteria = 0.05; // minimum performance ratio of elements to be considered
    auto maxIterationCount = 10000; // max iterations for topology optimization
    auto iterationStoppingCount = 350;
    auto solverSelection = SolverChoice::Eigen; // library to be used at linear algebraic equation solvings
    auto isPrintMeshInfo = false;
    auto isPrintForceInfo = true;
    auto isPrintElements = true;
    auto isHorizontalLoad = false;

    // Necessary fields
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Helper function to add restraints
    std::vector<bool> pin = { true, true, true, true, true, true };
    std::vector<bool> roller = { false, true, true, true, true, true };
    std::vector<bool> universal = { false, false, true, true, true, true };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };
    std::vector<int> realRestrainedNodes;
    auto addRestraint = [&](std::shared_ptr<Node> nR)
    {
        if (Utils::AreEqual(nR->Coordinate.X, 0.0, meshSize * 0.1) && Utils::AreEqual(nR->Coordinate.Y, 0.0, meshSize * 0.1))
        {
            restraints[nR->NodeIndex] = std::make_shared<Restraint>(nR, pin, rest);
            realRestrainedNodes.push_back(nR->NodeIndex);
        }
        else if (Utils::AreEqual(nR->Coordinate.X, lx, meshSize * 0.1) && Utils::AreEqual(nR->Coordinate.Y, 0.0, meshSize * 0.1))
        {
            restraints[nR->NodeIndex] = std::make_shared<Restraint>(nR, pin, rest);
            realRestrainedNodes.push_back(nR->NodeIndex);
        }
        else
            restraints[nR->NodeIndex] = std::make_shared<Restraint>(nR, universal, rest);
    };

    // Populate nodes
    int nodeIndex = 1;

    XYPoint firstPt(0, 0);
    XYPoint secondPt(lx, 0);
    XYPoint thirdPt(lx, ly);
    XYPoint fourthPt(0, ly);
    std::vector<XYPoint> vertices;
    vertices.push_back(firstPt);
    vertices.push_back(secondPt);
    vertices.push_back(thirdPt);
    vertices.push_back(fourthPt);

    XYPoint hole1(0.4, 0.4);
    XYPoint hole2(0.6, 0.4);
    XYPoint hole3(0.6, 0.6);
    XYPoint hole4(0.4, 0.6);
    std::vector<XYPoint> holeVertices;
    holeVertices.push_back(hole1);
    holeVertices.push_back(hole2);
    holeVertices.push_back(hole3);
    holeVertices.push_back(hole4);
    Shape hole(holeVertices);
    std::vector<Shape> holes;
    holes.push_back(hole);

    Shape sh(vertices, holes);
    auto meshPoints = sh.getPoints(meshSize);

    int nNodeX = ((int)(lx / meshSize)) + 1;
    int nNodeY = ((int)(ly / meshSize)) + 1;
    for (auto&& pt : meshPoints)
    {
        XYZPoint pt(pt.X, pt.Y, 0);
        auto n = std::make_shared<Node>(nodeIndex++, pt);
        nodes[n->NodeIndex] = n;
        addRestraint(n);
    }

    // Populate elements
    // Initially, populate horizontal members
    int memberIndex = 1;
    auto&& mat = std::make_shared<Material>(eMod, v, 0);
    auto&& sect = std::make_shared<Section>(0.01, 0, 0, 0);
    for (auto&& thisNodePair : nodes)
    {
        for (auto&& thatNodePair : nodes)
        {
            if ((thisNodePair.first != thatNodePair.first) && (thisNodePair.first < thatNodePair.first))
            {
                if (thisNodePair.second->Coordinate.DistanceTo(thatNodePair.second->Coordinate) < horizon)
                {
                    auto midX = (thisNodePair.second->Coordinate.X + thatNodePair.second->Coordinate.X) / 2.0;
                    auto midY = (thisNodePair.second->Coordinate.Y + thatNodePair.second->Coordinate.Y) / 2.0;
                    XYPoint midPt(midX, midY);
                    if (!hole.isInside(midPt))
                    {
                        elements[memberIndex] = std::make_shared<TrussMember>(memberIndex, thisNodePair.second, thatNodePair.second, sect, mat);
                        memberIndex++;
                    }
                }
            }
        }
    }

    // Force vector
    for (auto& n : nodes)
    {
        if ((Utils::AreEqual(n.second->Coordinate.X, lx / 2.0, meshSize * 0.1)) && (Utils::AreEqual(n.second->Coordinate.Y, ly, meshSize * 0.1)))
        {
            if (!isHorizontalLoad)
            {
                double load[] = { 0, -5000000000000, 0, 0, 0, 0 };
                nodalLoads[1] = std::make_shared<NodalLoad>(n.second, load);
            }
            else
            {
                double load[] = { 5000000000000, 0, 0, 0, 0, 0 };
                nodalLoads[1] = std::make_shared<NodalLoad>(n.second, load);
            }
            break;
        }
    }

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);
    LOGIF(" Total DOF: " << str->nUnrestrainedDOF, false);
    auto fVec = str->getForceVector(&nodalLoads);

    if (isPrintElements)
    {
        std::ofstream elementOutputStream;

        elementOutputStream.open("plot\\elements.dat");
        elementOutputStream << "BEGIN SCENE\n";
        for (auto&& e : elements)
        {
            auto&& firstNode = e.second->GelElementNodes().at(0);
            auto&& secondNode = e.second->GelElementNodes().at(1);
            elementOutputStream << firstNode->Coordinate.X << " " << firstNode->Coordinate.Y << "\n";
            elementOutputStream << secondNode->Coordinate.X << " " << secondNode->Coordinate.Y << "\n";
            elementOutputStream << "7\n";
        }

        elementOutputStream << "END SCENE\n";
        elementOutputStream.close();
    }

    // Start iterations
    auto comp = 0.0, tens = 0.0;
    auto prevNumOfElmStsfyCrt = 0;
    auto unchangedIterCount = 0;

    std::ofstream animationStream;
    animationStream.open("plot\\animation.dat");

    for (int i = 0; i < maxIterationCount; i++)
    {
        auto currNumOfElmStsfyCrt = 0;
        LOGIF(" Iteration #" << i + 1, true);
        // Fields to be used for storing extreme forces on members
        auto maxCompression = 0.0, maxTension = 0.0;
        str->updateStiffnessMatrix();

        // Perform initial analysis
        //auto disps = StructureSolver::CalculateDisplacements(*(str->StiffnessMatrix), fVec, str->nDOF, str->nUnrestrainedDOF, solverSelection);
        auto disps = StructureSolver::GetDisplacementForStaticCase(*str, solverSelection);

        auto leftTopNode = str->getNodeAt(0, ly, 0);
        auto rightTopNode = str->getNodeAt(lx, ly, 0);

        auto leftTopDofTX = leftTopNode->DofIndexTX;
        auto leftTopDofTY = leftTopNode->DofIndexTY;
        auto rightTopDofTX = rightTopNode->DofIndexTX;
        auto rightTopDofTY = rightTopNode->DofIndexTY;

        auto leftTopTX = disps(leftTopDofTX - 1, 0);
        auto leftTopTY = disps(leftTopDofTY - 1, 0);
        auto rightTopTX = disps(rightTopDofTX - 1, 0);
        auto rightTopTY = disps(rightTopDofTY - 1, 0);

        // Find internal forces on elements
        for (auto& elm : elements)
        {
            auto trMem = dynamic_cast<TrussMember*>(&*elm.second);
            auto axialForce = trMem->getAxialForce(disps);
            if (axialForce < maxCompression) maxCompression = axialForce;
            if (maxTension < axialForce) maxTension = axialForce;
            LOGIF(" Comp: " << maxCompression << ", Tens:" << maxTension, false);
        }

        animationStream << "BEGIN SCENE\n";
        // Modify elements
        for (auto& elm : elements)
        {
            auto trMem = dynamic_cast<TrussMember*>(&*elm.second);
            auto axialForce = trMem->getAxialForce(disps);
            auto isCompression = axialForce < 0;
            auto isTension = 0 < axialForce;
            auto ratio = axialForce / (axialForce < 0 ? maxCompression : maxTension);
            if (performanceRatioCriteria < ratio) currNumOfElmStsfyCrt++;
            ratio = pow(ratio, 1.0 / (i + 1));
            elm.second->updateStiffness(ratio);

            int colour = 0;
            if (isCompression)
            {
                if (ratio < 0.2) colour = 2;
                else if (ratio < 0.4) colour = 3;
                else if (ratio < 0.6) colour = 4;
                else if (ratio < 0.8) colour = 5;
                else if (ratio <= 1.0) colour = 6;
            }
            else if (isTension)
            {
                if (ratio < 0.2) colour = 7;
                else if (ratio < 0.4) colour = 8;
                else if (ratio < 0.6) colour = 9;
                else if (ratio < 0.8) colour = 10;
                else if (ratio <= 1.0) colour = 11;
            }
            else
            {
                colour = 0;
            }

            animationStream << trMem->Nodes[0]->Coordinate.X << " " << trMem->Nodes[0]->Coordinate.Y << "\n";
            animationStream << trMem->Nodes[1]->Coordinate.X << " " << trMem->Nodes[1]->Coordinate.Y << "\n";
            animationStream << colour << "\n";
        }

        animationStream << "END SCENE\n";

        comp = maxCompression;
        tens = maxTension;

        if (currNumOfElmStsfyCrt == prevNumOfElmStsfyCrt)
            unchangedIterCount++;

        prevNumOfElmStsfyCrt = currNumOfElmStsfyCrt;

        if (unchangedIterCount == iterationStoppingCount)
            break;
    }
    animationStream.close();

    LOG(" Iterations for optimization finished.");
    auto disps = StructureSolver::CalculateDisplacements(*(str->StiffnessMatrix), fVec, str->nDOF, str->nUnrestrainedDOF, solverSelection);

    // Identify nodes
    std::vector<int> nodeIndices;
    for (auto& nodePair : nodes)
    {
        // Get node
        auto& node = nodePair.second;

        // If there is a boundary condition at a node (either natural or essential), identify node 
        // as strut&tie system node and continue
        bool isStop = false;
        for (auto& res : realRestrainedNodes)
        {
            if (res == node->NodeIndex)
            {
                nodeIndices.push_back(node->NodeIndex);
                isStop = true;
                break;
            }
        }

        if (isStop) continue;

        for (auto& nLoad : nodalLoads)
        {
            if (nLoad.second->ActingNode->NodeIndex == node->NodeIndex)
            {
                nodeIndices.push_back(node->NodeIndex);
                isStop = true;
                break;
            }
        }

        if (isStop) continue;

        // Get list of indices of connected elements
        auto& listOfConnectedElements = node->ConnectedElements;

        // Calculate axial forces on connected elements. If their ratio are higher than performance criteria,
        // store them
        std::vector<int> elementsSatisfyPerformanceCriteria;
        for (auto& elmIndex : listOfConnectedElements)
        {
            auto trussElm = dynamic_cast<TrussMember*>(&*elements[elmIndex]);
            auto axialForce = trussElm->getAxialForce(disps);
            auto ratio = axialForce / (axialForce < 0 ? comp : tens);
            if (performanceRatioCriteria <= ratio)
                elementsSatisfyPerformanceCriteria.push_back(elmIndex);
        }

        // If there are at least three load carrying elements connecting to node, identify node as 
        // strut&tie system cantidate node and continue
        if (elementsSatisfyPerformanceCriteria.size() > 2)
        {
            nodeIndices.push_back(node->NodeIndex);
            continue;
        }
        // Else if there are two elements connecting to the node and if they are not on the same line,
        // identify node as strut&tie candidate node and continue
        else if (elementsSatisfyPerformanceCriteria.size() == 2)
        {
            auto firstElm = dynamic_cast<TrussMember*>(&*elements[elementsSatisfyPerformanceCriteria[0]]);
            auto secondElm = dynamic_cast<TrussMember*>(&*elements[elementsSatisfyPerformanceCriteria[1]]);
            Vector firstElmVector(firstElm->Nodes[0]->Coordinate, firstElm->Nodes[1]->Coordinate);
            Vector secondElmVector(secondElm->Nodes[0]->Coordinate, secondElm->Nodes[1]->Coordinate);

            auto firstElmUnitVec = firstElmVector.getUnitVector();
            auto secondElmUnitVec = secondElmVector.getUnitVector();
            auto asdasd = abs(firstElmUnitVec.dotProduct(secondElmUnitVec));
            // If dot-product of element vectors results in -1, they are on the same line.
            if (!Utils::AreEqual(abs(firstElmUnitVec.dotProduct(secondElmUnitVec)), 1))
                nodeIndices.push_back(node->NodeIndex);
        }
    }

    // After obtaining nodes, start merging nodes that are too close
    // Create node data structs
    std::vector<NodeData> nodeDataVector;
    for (size_t i = 0; i < nodeIndices.size(); i++)
    {
        NodeData nd;
        nd.nodeIndex = i;
        nd.coordinate = nodes[nodeIndices[i]]->Coordinate;
        nodeDataVector.push_back(nd);
    }

    LOG("\n STRUT AND TIE SYSTEM INITIAL NODES");
    for (auto& nd : nodeDataVector)
        LOG(" Node Index: " << (nd.nodeIndex + 1) << ", Coordinates(x,y): (" << nd.coordinate.X << ", " << nd.coordinate.Y << ")");

    // Merge the nodes that are close to each other
    std::vector<NodeData> mergedNodes;
    std::vector<NodeData> nodesToBeSkipped;
    auto nodeIndexer = 1;
    for (auto& nd : nodeDataVector)
    {
        // Skip the already encountered nodes
        bool isSkipOutside = false;
        for (auto& nskp : nodesToBeSkipped)
        {
            if (nskp.nodeIndex == nd.nodeIndex)
            {
                isSkipOutside = true;
                break;
            }
        }

        if (isSkipOutside)
            continue;

        std::vector<NodeData> nodesWillBeMerged;
        bool isFirstEncountering = true;
        for (auto& nnd : nodeDataVector)
        {
            // Skip the same nodes
            if (nd.nodeIndex == nnd.nodeIndex)
                continue;

            // Skip the already encountered nodes
            bool isSkipInside = false;
            for (auto& nskp : nodesToBeSkipped)
            {
                if (nskp.nodeIndex == nnd.nodeIndex)
                {
                    isSkipInside = true;
                    break;
                }
            }

            if (isSkipInside)
                continue;

            auto dist = nd.coordinate.DistanceTo(nnd.coordinate);
            if (dist < merger * meshSize)
            {
                if (isFirstEncountering)
                {
                    nodesWillBeMerged.push_back(nd);
                    isFirstEncountering = false;
                }

                nodesWillBeMerged.push_back(nnd);
            }
        }

        for (auto& nwbm : nodesWillBeMerged)
        {
            bool isInsert = true;

            for (auto& ntbs : nodesToBeSkipped)
            {
                if (nwbm.nodeIndex == ntbs.nodeIndex)
                {
                    isInsert = false;
                    break;
                }
            }

            if (isInsert)
                nodesToBeSkipped.push_back(nwbm);
        }

        if (1 < nodesWillBeMerged.size())
        {
            auto sumXCoord = 0.0, sumYCoord = 0.0;
            for (auto& nid : nodesWillBeMerged)
            {
                sumXCoord += nid.coordinate.X;
                sumYCoord += nid.coordinate.Y;
            }

            auto xCd = sumXCoord / (double(nodesWillBeMerged.size()));
            auto yCd = sumYCoord / (double(nodesWillBeMerged.size()));
            XYZPoint pt(xCd, yCd, 0.0);
            NodeData newNode;
            newNode.nodeIndex = nodeIndexer;
            newNode.coordinate = pt;
            mergedNodes.push_back(newNode);
            nodeIndexer++;
        }
        else
        {
            NodeData newNode;
            newNode.nodeIndex = nodeIndexer;
            newNode.coordinate = nd.coordinate;
            mergedNodes.push_back(newNode);
            nodeIndexer++;
        }
    }

    LOG("\n STRUT AND TIE SYSTEM MERGED NODES");
    for (auto& nd : mergedNodes)
        LOG(" Node Index: " << nd.nodeIndex << ", Coordinates(x,y): (" << nd.coordinate.X << ", " << nd.coordinate.Y << ")");

    // Since nodes of the system are found, create truss analogy by linking all the nodes to each other
    std::vector<TrussData> trussDataVector;
    int trussIndexer = 1;
    LOG("\n STRUT AND TIE SYSTEM TRUSSES");
    std::ofstream trussDataPrinter;

    if (isPrintElements)
    {
        trussDataPrinter.open("plot\\strutTieSystem.dat");
        trussDataPrinter << "BEGIN SCENE\n";
    }

    for (size_t i = 0; i < mergedNodes.size(); i++)
    {
        auto iNode = mergedNodes[i];
        for (size_t j = i + 1; j < mergedNodes.size(); j++)
        {
            auto jNode = mergedNodes[j];
            TrussData t;
            t.trussIndex = trussIndexer;
            t.iNode = &iNode;
            t.jNode = &jNode;
            trussDataVector.push_back(t);
            trussIndexer++;
            LOG(" Truss Index: " << t.trussIndex << ", i-End Node Index: " << t.iNode->nodeIndex << ", j-End Node Index: " << t.jNode->nodeIndex);
            if (isPrintElements)
            {
                trussDataPrinter << iNode.coordinate.X << " " << iNode.coordinate.Y << "\n";
                trussDataPrinter << jNode.coordinate.X << " " << jNode.coordinate.Y << "\n";
                trussDataPrinter << "7\n";
            }
        }

    }

    if (isPrintElements)
    {
        trussDataPrinter << "END SCENE";
        trussDataPrinter.close();
    }
}

void StrutTieDesignWithLatticeAndEnergyConservation()
{
    // INPUT CARD
    // Length: m, Force: N
    // Surface Outer Dimensions
    auto lx = 1.0; // length in x
    auto ly = 1.0; // length in y
    auto thickness = 0.25;
    auto eMod = 30e9; // youngs modulus (xe9 means x gigapascals)
    auto v = 0.3; // poissons ratio
    auto fc = 30e6; // concrete compressive strength
    auto fy = 420e6; // steel yield strength
    auto meshSize = 0.05;
    auto horizon = 3.5 * meshSize;
    auto merger = 3.5;
    double performanceRatioCriteria = 0.30; // minimum performance ratio of elements to be considered
    auto maxIterationCount = 10000; // max iterations for topology optimization
    auto iterationStoppingCount = 250;
    auto solverSelection = SolverChoice::Eigen; // library to be used at linear algebraic equation solvings
    auto isPrintMeshInfo = false;
    auto isPrintForceInfo = true;
    auto isPrintElements = true;
    auto isHorizontalLoad = false;
    auto isConserveEnergy = true;

    // Necessary fields
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Helper function to add restraints
    std::vector<bool> pin = { true, true, true, true, true, true };
    std::vector<bool> roller = { false, true, true, true, true, true };
    std::vector<bool> universal = { false, false, true, true, true, true };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };
    std::vector<int> realRestrainedNodes;
    auto addRestraint = [&](std::shared_ptr<Node> nR)
    {
        if (Utils::AreEqual(nR->Coordinate.X, 0.0, meshSize * 0.1) && Utils::AreEqual(nR->Coordinate.Y, 0.0, meshSize * 0.1))
        {
            restraints[nR->NodeIndex] = std::make_shared<Restraint>(nR, pin, rest);
            realRestrainedNodes.push_back(nR->NodeIndex);
        }
        else if (Utils::AreEqual(nR->Coordinate.X, lx, meshSize * 0.1) && Utils::AreEqual(nR->Coordinate.Y, 0.0, meshSize * 0.1))
        {
            restraints[nR->NodeIndex] = std::make_shared<Restraint>(nR, pin, rest);
            realRestrainedNodes.push_back(nR->NodeIndex);
        }
        else
            restraints[nR->NodeIndex] = std::make_shared<Restraint>(nR, universal, rest);
    };

    // Populate nodes
    int nodeIndex = 1;

    XYPoint firstPt(0, 0);
    XYPoint secondPt(lx, 0);
    XYPoint thirdPt(lx, ly);
    XYPoint fourthPt(0, ly);
    std::vector<XYPoint> vertices;
    vertices.push_back(firstPt);
    vertices.push_back(secondPt);
    vertices.push_back(thirdPt);
    vertices.push_back(fourthPt);

    XYPoint hole1(0.2, 0.2);
    XYPoint hole2(0.8, 0.2);
    XYPoint hole3(0.8, 0.8);
    XYPoint hole4(0.2, 0.8);
    std::vector<XYPoint> holeVertices;
    holeVertices.push_back(hole1);
    holeVertices.push_back(hole2);
    holeVertices.push_back(hole3);
    holeVertices.push_back(hole4);
    Shape hole(holeVertices);
    std::vector<Shape> holes;
    holes.push_back(hole);

    Shape sh(vertices, holes);
    auto meshPoints = sh.getPoints(meshSize);

    int nNodeX = ((int)(lx / meshSize)) + 1;
    int nNodeY = ((int)(ly / meshSize)) + 1;
    for (auto&& pt : meshPoints)
    {
        XYZPoint nodePt(pt.X, pt.Y, 0);
        auto n = std::make_shared<Node>(nodeIndex++, nodePt);
        nodes[n->NodeIndex] = n;
        addRestraint(n);
    }

    // Populate elements
    // Initially, populate horizontal members
    int memberIndex = 1;
    auto&& mat = std::make_shared<Material>(eMod, v, 0);
    auto&& sect = std::make_shared<Section>(0.01, 0, 0, 0);
    for (auto&& thisNodePair : nodes)
    {
        for (auto&& thatNodePair : nodes)
        {
            //if ((thisNodePair.first != thatNodePair.first) && (thisNodePair.first < thatNodePair.first))
            if (thisNodePair.first != thatNodePair.first)
            {
                if (thisNodePair.second->Coordinate.DistanceTo(thatNodePair.second->Coordinate) < horizon)
                {
                    auto midX = (thisNodePair.second->Coordinate.X + thatNodePair.second->Coordinate.X) / 2.0;
                    auto midY = (thisNodePair.second->Coordinate.Y + thatNodePair.second->Coordinate.Y) / 2.0;
                    XYPoint midPt(midX, midY);
                    if (!hole.isInside(midPt))
                    {
                        elements[memberIndex] = std::make_shared<TrussMember>(memberIndex, thisNodePair.second, thatNodePair.second, sect, mat);
                        memberIndex++;
                    }
                }
            }
        }
    }

    // Force vector
    for (auto& n : nodes)
    {
        if ((Utils::AreEqual(n.second->Coordinate.X, lx / 2.0, meshSize * 0.1)) && (Utils::AreEqual(n.second->Coordinate.Y, ly, meshSize * 0.1)))
        {
            if (!isHorizontalLoad)
            {
                double load[] = { 0, -5000000000000, 0, 0, 0, 0 };
                nodalLoads[1] = std::make_shared<NodalLoad>(n.second, load);
            }
            else
            {
                double load[] = { 5000000000000, 0, 0, 0, 0, 0 };
                nodalLoads[1] = std::make_shared<NodalLoad>(n.second, load);
            }
            break;
        }
    }

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);
    LOGIF(" Total DOF: " << str->nUnrestrainedDOF, false);
    auto fVec = str->getForceVector(&nodalLoads);

    if (isPrintElements)
    {
        std::ofstream elementOutputStream;

        elementOutputStream.open("plot\\elements.dat");
        elementOutputStream << "BEGIN SCENE\n";
        for (size_t i = 1; i <= str->Elements->size(); i++)
        {
            auto e = dynamic_cast<TrussMember*>(&*str->Elements->at(i));
            auto& firstNode = e->Nodes[0];
            auto& secondNode = e->Nodes[1];
            elementOutputStream << firstNode->Coordinate.X << " " << firstNode->Coordinate.Y << "\n";
            elementOutputStream << secondNode->Coordinate.X << " " << secondNode->Coordinate.Y << "\n";
            elementOutputStream << "7\n";
        }

        elementOutputStream << "END SCENE\n";
        elementOutputStream.close();
    }

    // Start iterations
    auto comp = 0.0, tens = 0.0;
    auto prevNumOfElmStsfyCrt = 0;
    auto unchangedIterCount = 0;

    std::ofstream animationStream;
    animationStream.open("plot\\animation.dat");

    for (int i = 0; i < maxIterationCount; i++)
    {
        auto currNumOfElmStsfyCrt = 0;
        LOGIF(" Iteration #" << i + 1, true);
        // Fields to be used for storing extreme forces on members
        auto maxCompression = 0.0, maxTension = 0.0;
        str->updateStiffnessMatrix();

        // Perform initial analysis
        auto disps = StructureSolver::CalculateDisplacements(*(str->StiffnessMatrix), fVec, str->nDOF, str->nUnrestrainedDOF, solverSelection);
        //auto disps = StructureSolver::GetDisplacementForStaticCase(*str, solverSelection);

        auto leftTopNode = str->getNodeAt(0, ly, 0);
        auto rightTopNode = str->getNodeAt(lx, ly, 0);

        auto leftTopDofTX = leftTopNode->DofIndexTX;
        auto leftTopDofTY = leftTopNode->DofIndexTY;
        auto rightTopDofTX = rightTopNode->DofIndexTX;
        auto rightTopDofTY = rightTopNode->DofIndexTY;

        auto leftTopTX = disps(leftTopDofTX - 1, 0);
        auto leftTopTY = disps(leftTopDofTY - 1, 0);
        auto rightTopTX = disps(rightTopDofTX - 1, 0);
        auto rightTopTY = disps(rightTopDofTY - 1, 0);

        // Find internal forces on elements
        for (auto& elm : elements)
        {
            auto trMem = dynamic_cast<TrussMember*>(&*elm.second);
            auto axialForce = trMem->getAxialForce(disps);
            if (axialForce < maxCompression) maxCompression = axialForce;
            if (maxTension < axialForce) maxTension = axialForce;
            LOGIF(" Comp: " << maxCompression << ", Tens:" << maxTension, false);
        }

        animationStream << "BEGIN SCENE\n";
        // Modify elements
        for (auto& elm : elements)
        {
            auto trMem = dynamic_cast<TrussMember*>(&*elm.second);
            auto axialForce = trMem->getAxialForce(disps);
            auto isCompression = axialForce < 0;
            auto isTension = 0 < axialForce;
            auto ratio = axialForce / (axialForce < 0 ? maxCompression : maxTension);
            if (performanceRatioCriteria < ratio) currNumOfElmStsfyCrt++;
            ratio = pow(ratio, 1.0 / (i + 1));
            elm.second->updateStiffness(ratio);

            int colour = 0;
            if (isCompression)
            {
                if (ratio < 0.2) colour = 2;
                else if (ratio < 0.4) colour = 3;
                else if (ratio < 0.6) colour = 4;
                else if (ratio < 0.8) colour = 5;
                else if (ratio <= 1.0) colour = 6;
            }
            else if (isTension)
            {
                if (ratio < 0.2) colour = 7;
                else if (ratio < 0.4) colour = 8;
                else if (ratio < 0.6) colour = 9;
                else if (ratio < 0.8) colour = 10;
                else if (ratio <= 1.0) colour = 11;
            }
            else
            {
                colour = 0;
            }

            animationStream << trMem->Nodes[0]->Coordinate.X << " " << trMem->Nodes[0]->Coordinate.Y << "\n";
            animationStream << trMem->Nodes[1]->Coordinate.X << " " << trMem->Nodes[1]->Coordinate.Y << "\n";
            animationStream << colour << "\n";
        }

        animationStream << "END SCENE\n";

        comp = maxCompression;
        tens = maxTension;

        if (currNumOfElmStsfyCrt == prevNumOfElmStsfyCrt)
            unchangedIterCount++;

        prevNumOfElmStsfyCrt = currNumOfElmStsfyCrt;

        if (unchangedIterCount == iterationStoppingCount)
            break;
    }
    animationStream.close();

    LOG(" Iterations for optimization finished.");
    auto disps = StructureSolver::CalculateDisplacements(*(str->StiffnessMatrix), fVec, str->nDOF, str->nUnrestrainedDOF, solverSelection);

    // Identify nodes
    std::vector<int> nodeIndices;
    for (auto& nodePair : nodes)
    {
        // Get node
        auto& node = nodePair.second;

        // If there is a boundary condition at a node (either natural or essential), identify node 
        // as strut&tie system node and continue
        bool isStop = false;
        for (auto& res : realRestrainedNodes)
        {
            if (res == node->NodeIndex)
            {
                nodeIndices.push_back(node->NodeIndex);
                isStop = true;
                break;
            }
        }

        if (isStop) continue;

        for (auto& nLoad : nodalLoads)
        {
            if (nLoad.second->ActingNode->NodeIndex == node->NodeIndex)
            {
                nodeIndices.push_back(node->NodeIndex);
                isStop = true;
                break;
            }
        }

        if (isStop) continue;

        // Get list of indices of connected elements
        auto& listOfConnectedElements = node->ConnectedElements;

        // Calculate axial forces on connected elements. If their ratio are higher than performance criteria,
        // store them
        std::vector<int> elementsSatisfyPerformanceCriteria;
        for (auto& elmIndex : listOfConnectedElements)
        {
            auto trussElm = dynamic_cast<TrussMember*>(&*elements[elmIndex]);
            auto axialForce = trussElm->getAxialForce(disps);
            auto ratio = axialForce / (axialForce < 0 ? comp : tens);
            if (performanceRatioCriteria <= ratio)
                elementsSatisfyPerformanceCriteria.push_back(elmIndex);
        }

        // If there are at least three load carrying elements connecting to node, identify node as 
        // strut&tie system cantidate node and continue
        if (elementsSatisfyPerformanceCriteria.size() > 2)
        {
            nodeIndices.push_back(node->NodeIndex);
            continue;
        }
        // Else if there are two elements connecting to the node and if they are not on the same line,
        // identify node as strut&tie candidate node and continue
        else if (elementsSatisfyPerformanceCriteria.size() == 2)
        {
            auto firstElm = dynamic_cast<TrussMember*>(&*elements[elementsSatisfyPerformanceCriteria[0]]);
            auto secondElm = dynamic_cast<TrussMember*>(&*elements[elementsSatisfyPerformanceCriteria[1]]);
            Vector firstElmVector(firstElm->Nodes[0]->Coordinate, firstElm->Nodes[1]->Coordinate);
            Vector secondElmVector(secondElm->Nodes[0]->Coordinate, secondElm->Nodes[1]->Coordinate);

            auto firstElmUnitVec = firstElmVector.getUnitVector();
            auto secondElmUnitVec = secondElmVector.getUnitVector();
            auto asdasd = abs(firstElmUnitVec.dotProduct(secondElmUnitVec));
            // If dot-product of element vectors results in -1, they are on the same line.
            if (!Utils::AreEqual(abs(firstElmUnitVec.dotProduct(secondElmUnitVec)), 1))
                nodeIndices.push_back(node->NodeIndex);
        }
    }

    // After obtaining nodes, start merging nodes that are too close
    // Create node data structs
    std::vector<NodeData> nodeDataVector;
    for (size_t i = 0; i < nodeIndices.size(); i++)
    {
        NodeData nd;
        nd.nodeIndex = i;
        nd.coordinate = nodes[nodeIndices[i]]->Coordinate;
        nodeDataVector.push_back(nd);
    }

    LOG("\n STRUT AND TIE SYSTEM INITIAL NODES");
    for (auto& nd : nodeDataVector)
        LOG(" Node Index: " << (nd.nodeIndex + 1) << ", Coordinates(x,y): (" << nd.coordinate.X << ", " << nd.coordinate.Y << ")");

    // Merge the nodes that are close to each other
    std::vector<NodeData> mergedNodes;
    std::vector<NodeData> nodesToBeSkipped;
    auto nodeIndexer = 1;
    for (auto& nd : nodeDataVector)
    {
        // Skip the already encountered nodes
        bool isSkipOutside = false;
        for (auto& nskp : nodesToBeSkipped)
        {
            if (nskp.nodeIndex == nd.nodeIndex)
            {
                isSkipOutside = true;
                break;
            }
        }

        if (isSkipOutside)
            continue;

        std::vector<NodeData> nodesWillBeMerged;
        bool isFirstEncountering = true;
        for (auto& nnd : nodeDataVector)
        {
            // Skip the same nodes
            if (nd.nodeIndex == nnd.nodeIndex)
                continue;

            // Skip the already encountered nodes
            bool isSkipInside = false;
            for (auto& nskp : nodesToBeSkipped)
            {
                if (nskp.nodeIndex == nnd.nodeIndex)
                {
                    isSkipInside = true;
                    break;
                }
            }

            if (isSkipInside)
                continue;

            auto dist = nd.coordinate.DistanceTo(nnd.coordinate);
            if (dist < merger * meshSize)
            {
                if (isFirstEncountering)
                {
                    nodesWillBeMerged.push_back(nd);
                    isFirstEncountering = false;
                }

                nodesWillBeMerged.push_back(nnd);
            }
        }

        for (auto& nwbm : nodesWillBeMerged)
        {
            bool isInsert = true;

            for (auto& ntbs : nodesToBeSkipped)
            {
                if (nwbm.nodeIndex == ntbs.nodeIndex)
                {
                    isInsert = false;
                    break;
                }
            }

            if (isInsert)
                nodesToBeSkipped.push_back(nwbm);
        }

        if (1 < nodesWillBeMerged.size())
        {
            auto sumXCoord = 0.0, sumYCoord = 0.0;
            for (auto& nid : nodesWillBeMerged)
            {
                sumXCoord += nid.coordinate.X;
                sumYCoord += nid.coordinate.Y;
            }

            auto xCd = sumXCoord / (double(nodesWillBeMerged.size()));
            auto yCd = sumYCoord / (double(nodesWillBeMerged.size()));
            XYZPoint pt(xCd, yCd, 0.0);
            NodeData newNode;
            newNode.nodeIndex = nodeIndexer;
            newNode.coordinate = pt;
            mergedNodes.push_back(newNode);
            nodeIndexer++;
        }
        else
        {
            NodeData newNode;
            newNode.nodeIndex = nodeIndexer;
            newNode.coordinate = nd.coordinate;
            mergedNodes.push_back(newNode);
            nodeIndexer++;
        }
    }

    LOG("\n STRUT AND TIE SYSTEM MERGED NODES");
    for (auto& nd : mergedNodes)
        LOG(" Node Index: " << nd.nodeIndex << ", Coordinates(x,y): (" << nd.coordinate.X << ", " << nd.coordinate.Y << ")");

    // Since nodes of the system are found, create truss analogy by linking all the nodes to each other
    std::vector<TrussData> trussDataVector;
    int trussIndexer = 1;
    LOG("\n STRUT AND TIE SYSTEM TRUSSES");
    std::ofstream trussDataPrinter;

    if (isPrintElements)
    {
        trussDataPrinter.open("plot\\strutTieSystem.dat");
        trussDataPrinter << "BEGIN SCENE\n";
    }

    for (size_t i = 0; i < mergedNodes.size(); i++)
    {
        auto iNode = mergedNodes[i];
        for (size_t j = i + 1; j < mergedNodes.size(); j++)
        {
            auto jNode = mergedNodes[j];
            TrussData t;
            t.trussIndex = trussIndexer;
            t.iNode = &iNode;
            t.jNode = &jNode;
            trussDataVector.push_back(t);
            trussIndexer++;
            LOG(" Truss Index: " << t.trussIndex << ", i-End Node Index: " << t.iNode->nodeIndex << ", j-End Node Index: " << t.jNode->nodeIndex);
            if (isPrintElements)
            {
                trussDataPrinter << iNode.coordinate.X << " " << iNode.coordinate.Y << "\n";
                trussDataPrinter << jNode.coordinate.X << " " << jNode.coordinate.Y << "\n";
                trussDataPrinter << "7\n";
            }
        }

    }

    if (isPrintElements)
    {
        trussDataPrinter << "END SCENE";
        trussDataPrinter.close();
    }
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
    auto modalPeriods = StructureSolver::GetModalPeriods(*str, SolverChoice::Eigen);
    LOG("");
    for (size_t i = 0; i < modalPeriods.RowCount; i++)
        std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods(i, 0) << "\n";

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
    auto modalPeriods = StructureSolver::GetModalPeriods(*str, SolverChoice::Eigen);
    LOG("");
    for (size_t i = 0; i < modalPeriods.RowCount; i++)
        std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods(i, 0) << "\n";

    return;
}

void TableModalAnalysis()
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
    int nElmX = 1;
    int nElmY = 1;
    auto membraneType = MembraneType::Drilling;
    auto plateType = PlateType::NONE;
    double lX = 4;
    double lY = 0.6;
    double thickness = 0.4;
    double e = 25e9;
    double v = 0.2;
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

    std::vector<bool> isRest = { true, true, true, true, true, false };
    //std::vector<bool> universal = { false, false, true, true, true, false };
    std::vector<bool> universal = { false, false, false, false, false, false };
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
            auto sm = std::make_shared<ShellMember>(idx, iNode, jNode, kNode, lNode, mt, thickness, membraneType, plateType); elements[sm->ElementIndex] = sm;
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
    int nElmX = 8;
    int nElmY = 1;
    auto membraneType = MembraneType::Incompatible;
    auto plateType = PlateType::NONE;
    bool isAxial = false;
    double stressLoc = 0; // At support
    double lX = 4;
    double lY = 0.6;
    double thickness = 0.4;
    double e = 25e9;
    double v = 0.0;
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

    std::vector<bool> isRest = { true, true, true, true, true, false };
    std::vector<bool> universal = { false, false, true, true, true, false };
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
            auto sm = std::make_shared<ShellMember>(idx, iNode, jNode, kNode, lNode, mt, thickness, membraneType, plateType); elements[sm->ElementIndex] = sm;
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

void CE583Assignment7_2()
{
    // Units are in N & m
    int nElmX = 5;
    int nElmY = 20;
    double lX = 0.4;
    double lY = 4.0;
    double e = 25e9;
    double v = 0.0;
    double t = 0.6;

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
    std::vector<bool> universal = { true, true, false, false, false, true };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };
    auto f = -20000.0 / (nElmX + 1.0);
    double nodalLoad[6] = { 0, 0, f, 0, 0, 0 };
    int loadIdx = 1;
    std::vector<int> tipNodeIndices;
    auto addNode = [&](int nodeIndex, double x, double y, double z) {
        XYZPoint pt(x, y, z);
        nodes[nodeIndex] = std::make_shared<Node>(nodeIndex, pt);

        // Restraint
        if (Utils::AreEqual(y, 0))
            restraints[nodeIndex] = std::make_shared<Restraint>(nodes[nodeIndex], isRest, rest);
        else
            restraints[nodeIndex] = std::make_shared<Restraint>(nodes[nodeIndex], universal, rest);

        // Loads
        if (Utils::AreEqual(y, 4))
        {
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

    /*for (auto& n : nodes)
        LOG(" Node Index: " << n.first << ", X,Y = " << n.second->Coordinate.X << ", " << n.second->Coordinate.Y);
    LOG("");
    for (auto& n : nodes)
        LOG(" Node Index: " << n.first << ", DOF Indices = " << n.second->Coordinate.X << ", " << n.second->Coordinate.Y);
    LOG("");
    for (auto& eP : elements)
        LOG(" Element Index: " << eP.first << ", I-Node:" << eP.second->GelElementNodes()[0]->NodeIndex << ", J-Node:" << eP.second->GelElementNodes()[1]->NodeIndex
            << ", K-Node:" << eP.second->GelElementNodes()[2]->NodeIndex << ", L-Node:" << eP.second->GelElementNodes()[3]->NodeIndex);
    LOG("");
    for (auto& r : restraints)
        LOG(" Restrained Node ID: " << r.second->RestrainedNode->NodeIndex);

    LOG("");
    for (auto& nl : nodalLoads)
        LOG(" Loaded Node ID: " << nl.second->ActingNode->NodeIndex << ", Load: " << nl.second->Loads[2] << " kN");*/

        // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Armadillo);

    // Print tip displacement
    auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[tipNodeIndices[0]], disps);

    LOG(" Tip Displacement");
    LOG(" Node Index: " << nodes[tipNodeIndices[0]]->NodeIndex);
    LOG(" Node Location: " << nodes[tipNodeIndices[0]]->Coordinate.X << " m, " << nodes[tipNodeIndices[0]]->Coordinate.Y << " m");
    LOG(" Vertical Displacement: " << nodalDisp(2, 0) << " m");
    LOG("");

    auto format = [](double d) { return abs(d) < 1 ? 0 : d / 1000; };

    auto reactions = StructureSolver::CalculatePlateForces(*membersAtSupport[0], disps);
    LOG(" Reactions");
    LOG("   Moments");
    LOG("     " << format(reactions(1, 0)) << " kN-m/m");
    LOG("     " << format(reactions(1, 1)) << " kN-m/m");
    LOG("     " << format(reactions(1, 2)) << " kN-m/m");
    LOG("     " << format(reactions(1, 3)) << " kN-m/m");
    LOG("   Shear");
    LOG("     " << format(reactions(4, 0)) << " kN/m");

    LOG("");
}

void CE583Assignment7_3()
{
    // Units are in N & m

    // Coordinates
    XYZPoint ptI(0, 0, 0);
    XYZPoint ptJ(0.4, 0, 0);
    XYZPoint ptK(0.4, 4, 0);
    XYZPoint ptL(0, 4, 0);
    XYZPoint ptIJ(0.2, 0, 0);
    XYZPoint ptJK(0.4, 2, 0);
    XYZPoint ptKL(0.2, 4, 0);
    XYZPoint ptLI(0, 2, 0);

    // Nodes
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    auto nodeI = std::make_shared<Node>(1, ptI); nodes[nodeI->NodeIndex] = nodeI;
    auto nodeJ = std::make_shared<Node>(2, ptJ); nodes[nodeJ->NodeIndex] = nodeJ;
    auto nodeK = std::make_shared<Node>(3, ptK); nodes[nodeK->NodeIndex] = nodeK;
    auto nodeL = std::make_shared<Node>(4, ptL); nodes[nodeL->NodeIndex] = nodeL;
    auto nodeIJ = std::make_shared<Node>(5, ptIJ); nodes[nodeIJ->NodeIndex] = nodeIJ;
    auto nodeJK = std::make_shared<Node>(6, ptJK); nodes[nodeJK->NodeIndex] = nodeJK;
    auto nodeKL = std::make_shared<Node>(7, ptKL); nodes[nodeKL->NodeIndex] = nodeKL;
    auto nodeLI = std::make_shared<Node>(8, ptLI); nodes[nodeLI->NodeIndex] = nodeLI;

    // Material
    auto mat = std::make_shared<Material>(25e9, 0.0, 0.0);

    // Members
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    auto elm = std::make_shared<SerendipityShell>(1, nodeI, nodeJ, nodeK, nodeL, nodeIJ, nodeJK, nodeKL, nodeLI, mat, 0.6); elements[elm->ElementIndex] = elm;

    // Restraints
    std::vector<bool> fullyFixed = { true, true, true, true, true, true };
    std::vector<bool> universal = { true, true, false, false, false, true };
    std::vector<double> rest = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    auto fixed1 = std::make_shared<Restraint>(nodeI, fullyFixed, rest); restraints[1] = fixed1;
    auto fixed2 = std::make_shared<Restraint>(nodeIJ, fullyFixed, rest); restraints[2] = fixed2;
    auto fixed3 = std::make_shared<Restraint>(nodeJ, fullyFixed, rest); restraints[3] = fixed3;
    auto univerasllyRestrained1 = std::make_shared<Restraint>(nodeJK, universal, rest); restraints[4] = univerasllyRestrained1;
    auto univerasllyRestrained2 = std::make_shared<Restraint>(nodeK, universal, rest); restraints[5] = univerasllyRestrained2;
    auto univerasllyRestrained3 = std::make_shared<Restraint>(nodeKL, universal, rest); restraints[6] = univerasllyRestrained3;
    auto univerasllyRestrained4 = std::make_shared<Restraint>(nodeL, universal, rest); restraints[7] = univerasllyRestrained4;
    auto univerasllyRestrained5 = std::make_shared<Restraint>(nodeLI, universal, rest); restraints[8] = univerasllyRestrained5;

    // Nodal loads
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    auto f = (-20000.0 / 3.0);
    double nodalLoad[6] = { 0, 0, f, 0, 0, 0 };
    auto nl1 = std::make_shared<NodalLoad>(nodeK, nodalLoad); nodalLoads[1] = nl1;
    auto nl2 = std::make_shared<NodalLoad>(nodeKL, nodalLoad); nodalLoads[2] = nl2;
    auto nl3 = std::make_shared<NodalLoad>(nodeL, nodalLoad); nodalLoads[3] = nl3;

    // Distributed loads
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Armadillo);
    // Print tip displacement
    auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[7], disps);

    LOG(" Tip Displacement");
    LOG(" Node Index: " << nodes[7]->NodeIndex);
    LOG(" Node Location: " << nodes[7]->Coordinate.X << " m, " << nodes[7]->Coordinate.Y << " m");
    LOG(" Vertical Displacement: " << nodalDisp(2, 0) << " m");
    LOG("");

    auto format = [](double d) { return abs(d) < 1 ? 0 : d / 1000; };

    auto reactions = StructureSolver::CalculateSerendipityPlateForces(*elm, disps);
    LOG(" Reactions");
    LOG("   Moments");
    LOG("     " << format(reactions(1, 0)) << " kN-m/m");
    LOG("     " << format(reactions(1, 1)) << " kN-m/m");
    LOG("     " << format(reactions(1, 2)) << " kN-m/m");
    LOG("     " << format(reactions(1, 3)) << " kN-m/m");
    LOG("   Shear");
    LOG("     " << format(reactions(4, 0)) << " kN/m");

    LOG("");
}

void CE583Assignment8_2_MembraneAction()
{
    // Input Card (Units are in N & m)
    auto membraneType = MembraneType::Drilling;
    auto plateType = PlateType::MindlinFourNode;
    double thickness = 0.4;
    double e = 25e9;
    double v = 0.0;

    // Solve
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    XYZPoint pt1(0.0, 0.0, 0.0);
    XYZPoint pt2(4.0, 0.0, 0.0);
    XYZPoint pt3(4.0, 0.6, 0.0);
    XYZPoint pt4(0.0, 0.6, 0.0);

    nodes[1] = std::make_shared<Node>(1, pt1);
    nodes[2] = std::make_shared<Node>(2, pt2);
    nodes[3] = std::make_shared<Node>(3, pt3);
    nodes[4] = std::make_shared<Node>(4, pt4);

    std::vector<bool> fixed = { true, true, true, true, true, false };
    std::vector<bool> universal = { false, false, false, false, false, false };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };

    restraints[1] = std::make_shared<Restraint>(nodes[1], fixed, rest);
    restraints[2] = std::make_shared<Restraint>(nodes[4], fixed, rest);
    restraints[3] = std::make_shared<Restraint>(nodes[2], universal, rest);
    restraints[4] = std::make_shared<Restraint>(nodes[3], universal, rest);

    elements[1] = std::make_shared<ShellMember>(1, nodes[1], nodes[2], nodes[3], nodes[4],
        std::make_shared<Material>(e, v, 0), 0.4, membraneType, plateType);

    // Nodal loads
    // Tip load is -20000 kN. Divide to tip nodes
    double nodalForce = 10000;
    double nodalLoad[6] = { 0, -10000, 0, 0, 0, 0 };
    nodalLoads[1] = std::make_shared<NodalLoad>(nodes[2], nodalLoad);
    nodalLoads[2] = std::make_shared<NodalLoad>(nodes[3], nodalLoad);

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Armadillo);

    auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[2], disps);
    LOG(" Node Index: " << nodes[2]->NodeIndex);
    LOG(" Node Location: " << nodes[2]->Coordinate.X << " m, " << nodes[2]->Coordinate.Y << " m");
    LOG(" Vertical Displacement: " << nodalDisp(1, 0) << " m");
}

void CE583Assignment8_2_PlateAction()
{
    // Input Card (Units are in N & m)
    auto membraneType = MembraneType::Drilling;
    auto plateType = PlateType::MindlinFourNode;
    double thickness = 0.6;
    double e = 25e6;
    double v = 0.0;

    // Solve
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    XYZPoint pt1(0.0, 0.0, 0.2);
    XYZPoint pt2(4.0, 0., 0.2);
    XYZPoint pt3(4.0, 0., -0.2);
    XYZPoint pt4(0.0, 0.0, -0.2);

    nodes[1] = std::make_shared<Node>(1, pt1);
    nodes[2] = std::make_shared<Node>(2, pt2);
    nodes[3] = std::make_shared<Node>(3, pt3);
    nodes[4] = std::make_shared<Node>(4, pt4);

    std::vector<bool> fixed = { true, true, true, true, false,  true };
    std::vector<bool> universal = { false, false, false, false, false, false };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };

    restraints[1] = std::make_shared<Restraint>(nodes[1], fixed, rest);
    restraints[2] = std::make_shared<Restraint>(nodes[4], fixed, rest);
    restraints[3] = std::make_shared<Restraint>(nodes[2], universal, rest);
    restraints[4] = std::make_shared<Restraint>(nodes[3], universal, rest);

    elements[1] = std::make_shared<ShellMember>(1, nodes[1], nodes[2], nodes[3], nodes[4],
        std::make_shared<Material>(e, v, 0), thickness, membraneType, plateType);

    // Nodal loads
    // Tip load is -20000 kN. Divide to tip nodes
    double nodalForce = 10;
    double nodalLoad[6] = { 0, -10, 0, 0, 0, 0 };
    nodalLoads[1] = std::make_shared<NodalLoad>(nodes[2], nodalLoad);
    nodalLoads[2] = std::make_shared<NodalLoad>(nodes[3], nodalLoad);

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Armadillo);

    auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[2], disps);
    LOG(" Node Index: " << nodes[2]->NodeIndex);
    LOG(" Node Location: " << nodes[2]->Coordinate.X << " m, " << nodes[2]->Coordinate.Y << " m");
    LOG(" Vertical Displacement: " << nodalDisp(1, 0) << " m");
}

void CE583Assignment8_3()
{
    // Maps
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Nodes
    XYZPoint pt1(0, 0.195, 0.2);
    XYZPoint pt2(0, 0.195, 0);
    XYZPoint pt3(0, 0.195, -0.2);
    XYZPoint pt4(0.4, 0.195, 0.2);
    XYZPoint pt5(0.4, 0.195, 0);
    XYZPoint pt6(0.4, 0.195, -0.2);
    XYZPoint pt7(0.8, 0.195, 0.2);
    XYZPoint pt8(0.8, 0.195, 0);
    XYZPoint pt9(0.8, 0.195, -0.2);
    XYZPoint pt10(1.2, 0.195, 0.2);
    XYZPoint pt11(1.2, 0.195, 0);
    XYZPoint pt12(1.2, 0.195, -0.2);
    XYZPoint pt13(1.6, 0.195, 0.2);
    XYZPoint pt14(1.6, 0.195, 0);
    XYZPoint pt15(1.6, 0.195, -0.2);
    XYZPoint pt16(2, 0.195, 0.2);
    XYZPoint pt17(2, 0.195, 0);
    XYZPoint pt18(2, 0.195, -0.2);
    XYZPoint pt19(0, -0.195, 0.2);
    XYZPoint pt20(0, -0.195, 0);
    XYZPoint pt21(0, -0.195, -0.2);
    XYZPoint pt22(0.4, -0.195, 0.2);
    XYZPoint pt23(0.4, -0.195, 0);
    XYZPoint pt24(0.4, -0.195, -0.2);
    XYZPoint pt25(0.8, -0.195, 0.2);
    XYZPoint pt26(0.8, -0.195, 0);
    XYZPoint pt27(0.8, -0.195, -0.2);
    XYZPoint pt28(1.2, -0.195, 0.2);
    XYZPoint pt29(1.2, -0.195, 0);
    XYZPoint pt30(1.2, -0.195, -0.2);
    XYZPoint pt31(1.6, -0.195, 0.2);
    XYZPoint pt32(1.6, -0.195, 0);
    XYZPoint pt33(1.6, -0.195, -0.2);
    XYZPoint pt34(2, -0.195, 0.2);
    XYZPoint pt35(2, -0.195, 0);
    XYZPoint pt36(2, -0.195, -0.2);
    XYZPoint pt37(0, 0, 0);
    XYZPoint pt38(0.4, 0, 0);
    XYZPoint pt39(0.8, 0, 0);
    XYZPoint pt40(1.2, 0, 0);
    XYZPoint pt41(1.6, 0, 0);
    XYZPoint pt42(2, 0, 0);

    nodes[1] = std::make_shared<Node>(1, pt1);
    nodes[2] = std::make_shared<Node>(2, pt2);
    nodes[3] = std::make_shared<Node>(3, pt3);
    nodes[4] = std::make_shared<Node>(4, pt4);
    nodes[5] = std::make_shared<Node>(5, pt5);
    nodes[6] = std::make_shared<Node>(6, pt6);
    nodes[7] = std::make_shared<Node>(7, pt7);
    nodes[8] = std::make_shared<Node>(8, pt8);
    nodes[9] = std::make_shared<Node>(9, pt9);
    nodes[10] = std::make_shared<Node>(10, pt10);
    nodes[11] = std::make_shared<Node>(11, pt11);
    nodes[12] = std::make_shared<Node>(12, pt12);
    nodes[13] = std::make_shared<Node>(13, pt13);
    nodes[14] = std::make_shared<Node>(14, pt14);
    nodes[15] = std::make_shared<Node>(15, pt15);
    nodes[16] = std::make_shared<Node>(16, pt16);
    nodes[17] = std::make_shared<Node>(17, pt17);
    nodes[18] = std::make_shared<Node>(18, pt18);
    nodes[19] = std::make_shared<Node>(19, pt19);
    nodes[20] = std::make_shared<Node>(20, pt20);
    nodes[21] = std::make_shared<Node>(21, pt21);
    nodes[22] = std::make_shared<Node>(22, pt22);
    nodes[23] = std::make_shared<Node>(23, pt23);
    nodes[24] = std::make_shared<Node>(24, pt24);
    nodes[25] = std::make_shared<Node>(25, pt25);
    nodes[26] = std::make_shared<Node>(26, pt26);
    nodes[27] = std::make_shared<Node>(27, pt27);
    nodes[28] = std::make_shared<Node>(28, pt28);
    nodes[29] = std::make_shared<Node>(29, pt29);
    nodes[30] = std::make_shared<Node>(30, pt30);
    nodes[31] = std::make_shared<Node>(31, pt31);
    nodes[32] = std::make_shared<Node>(32, pt32);
    nodes[33] = std::make_shared<Node>(33, pt33);
    nodes[34] = std::make_shared<Node>(34, pt34);
    nodes[35] = std::make_shared<Node>(35, pt35);
    nodes[36] = std::make_shared<Node>(36, pt36);
    nodes[37] = std::make_shared<Node>(37, pt37);
    nodes[38] = std::make_shared<Node>(38, pt38);
    nodes[39] = std::make_shared<Node>(39, pt39);
    nodes[40] = std::make_shared<Node>(40, pt40);
    nodes[41] = std::make_shared<Node>(41, pt41);
    nodes[42] = std::make_shared<Node>(42, pt42);

    // Elements
    auto mat = std::make_shared<Material>(70e9, 0.0, 0.0);
    elements[1] = std::make_shared<ShellMember>(1, nodes[1], nodes[4], nodes[5], nodes[2], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[2] = std::make_shared<ShellMember>(2, nodes[2], nodes[5], nodes[6], nodes[3], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[3] = std::make_shared<ShellMember>(3, nodes[4], nodes[7], nodes[8], nodes[5], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[4] = std::make_shared<ShellMember>(4, nodes[5], nodes[8], nodes[9], nodes[6], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[5] = std::make_shared<ShellMember>(5, nodes[7], nodes[10], nodes[11], nodes[8], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[6] = std::make_shared<ShellMember>(6, nodes[8], nodes[11], nodes[12], nodes[9], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[7] = std::make_shared<ShellMember>(7, nodes[10], nodes[13], nodes[14], nodes[11], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[8] = std::make_shared<ShellMember>(8, nodes[11], nodes[14], nodes[15], nodes[12], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[9] = std::make_shared<ShellMember>(9, nodes[13], nodes[16], nodes[17], nodes[14], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[10] = std::make_shared<ShellMember>(10, nodes[14], nodes[17], nodes[18], nodes[15], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[11] = std::make_shared<ShellMember>(11, nodes[19], nodes[22], nodes[23], nodes[20], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[12] = std::make_shared<ShellMember>(12, nodes[20], nodes[23], nodes[24], nodes[21], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[13] = std::make_shared<ShellMember>(13, nodes[22], nodes[25], nodes[26], nodes[23], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[14] = std::make_shared<ShellMember>(14, nodes[23], nodes[26], nodes[27], nodes[24], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[15] = std::make_shared<ShellMember>(15, nodes[25], nodes[28], nodes[29], nodes[26], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[16] = std::make_shared<ShellMember>(16, nodes[26], nodes[29], nodes[30], nodes[27], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[17] = std::make_shared<ShellMember>(17, nodes[28], nodes[31], nodes[32], nodes[29], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[18] = std::make_shared<ShellMember>(18, nodes[29], nodes[32], nodes[33], nodes[30], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[19] = std::make_shared<ShellMember>(19, nodes[31], nodes[34], nodes[35], nodes[32], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[20] = std::make_shared<ShellMember>(20, nodes[32], nodes[35], nodes[36], nodes[33], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[21] = std::make_shared<ShellMember>(21, nodes[20], nodes[23], nodes[38], nodes[37], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[22] = std::make_shared<ShellMember>(22, nodes[37], nodes[38], nodes[5], nodes[2], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[23] = std::make_shared<ShellMember>(23, nodes[23], nodes[26], nodes[39], nodes[38], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[24] = std::make_shared<ShellMember>(24, nodes[38], nodes[39], nodes[8], nodes[5], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[25] = std::make_shared<ShellMember>(25, nodes[26], nodes[29], nodes[40], nodes[39], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[26] = std::make_shared<ShellMember>(26, nodes[39], nodes[40], nodes[11], nodes[8], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[27] = std::make_shared<ShellMember>(27, nodes[29], nodes[32], nodes[41], nodes[40], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[28] = std::make_shared<ShellMember>(28, nodes[40], nodes[41], nodes[14], nodes[11], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[29] = std::make_shared<ShellMember>(29, nodes[32], nodes[35], nodes[42], nodes[41], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);
    elements[30] = std::make_shared<ShellMember>(30, nodes[41], nodes[42], nodes[17], nodes[14], mat, 0.01, MembraneType::Drilling, PlateType::MindlinFourNode);

    // Restraints
    std::vector<bool> fixed = { true, true, true, true, true, true };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };
    restraints[1] = std::make_shared<Restraint>(nodes[1], fixed, rest);
    restraints[2] = std::make_shared<Restraint>(nodes[2], fixed, rest);
    restraints[3] = std::make_shared<Restraint>(nodes[3], fixed, rest);
    restraints[19] = std::make_shared<Restraint>(nodes[19], fixed, rest);
    restraints[20] = std::make_shared<Restraint>(nodes[20], fixed, rest);
    restraints[21] = std::make_shared<Restraint>(nodes[21], fixed, rest);
    restraints[37] = std::make_shared<Restraint>(nodes[37], fixed, rest);

    // Nodal loads
    // double nodalLoad[6] = { 0, -20000.0 / 3.0, 0, 0, 0, 0 };
    // nodalLoads[1] = std::make_shared<NodalLoad>(nodes[17], nodalLoad);
    // nodalLoads[2] = std::make_shared<NodalLoad>(nodes[35], nodalLoad);
    // nodalLoads[3] = std::make_shared<NodalLoad>(nodes[42], nodalLoad);

    double nodalLoadPos[6] = { 0, 0, 40000.0 / 0.39, 0, 0, 0 };
    double nodalLoadNeg[6] = { 0, 0, -40000.0 / 0.39, 0, 0, 0 };
    nodalLoads[1] = std::make_shared<NodalLoad>(nodes[16], nodalLoadPos);
    nodalLoads[2] = std::make_shared<NodalLoad>(nodes[18], nodalLoadPos);
    nodalLoads[3] = std::make_shared<NodalLoad>(nodes[34], nodalLoadNeg);
    nodalLoads[4] = std::make_shared<NodalLoad>(nodes[36], nodalLoadNeg);

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    auto nd = str->nUnrestrainedDOF - 1;

    // Solve displacement
    auto disps = StructureSolver::GetDisplacementForStaticCase(*str, SolverChoice::Eigen);

    auto printDisplacements = [&](int nodeIdx) {
        auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[nodeIdx], disps);
        LOG(" Node Index: " << nodes[nodeIdx]->NodeIndex);
        LOG(" Node Location: " << nodes[nodeIdx]->Coordinate.X << " m, " << nodes[nodeIdx]->Coordinate.Y << " m, " << nodes[nodeIdx]->Coordinate.Z << " m");
        LOG(" Vertical Displacement: " << nodalDisp(1, 0) << " m");
        LOG("");
    };

    auto printAngleOfTwist = [&](int nodeIdx) {
        auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[nodeIdx], disps);
        LOG(" Node Index: " << nodes[nodeIdx]->NodeIndex);
        LOG(" Node Location: " << nodes[nodeIdx]->Coordinate.X << " m, " << nodes[nodeIdx]->Coordinate.Y << " m, " << nodes[nodeIdx]->Coordinate.Z << " m");
        LOG(" Angle of twist: " << nodalDisp(3, 0) << " rad");
        LOG("");
    };

    auto printVerticalDisplacements = [&](int nodeIdx) {
        auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[nodeIdx], disps);
        LOG((nodes[nodeIdx]->Coordinate.Z + nodalDisp(2, 0)) << "   " << (nodes[nodeIdx]->Coordinate.Y + nodalDisp(1, 0)));
    };

    auto printHorizontalDipslacements = [&](int nodeIdx) {
        auto nodalDisp = StructureSolver::GetNodalDisplacements(*nodes[nodeIdx], disps);
        LOG((nodes[nodeIdx]->Coordinate.Z) << "   " << (nodalDisp(0, 0)));
    };

    printHorizontalDipslacements(16);
    printHorizontalDipslacements(17);
    printHorizontalDipslacements(18);
    LOG("---------------");
    printHorizontalDipslacements(34);
    printHorizontalDipslacements(35);
    printHorizontalDipslacements(36);

}

void CE583Assignment9_2_1()
{
    // Fields
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Nodal coordinates
    XYZPoint pt1(0, 0, 0);
    XYZPoint pt2(30.0, 0, 0);

    // Nodes
    nodes[1] = std::make_shared<Node>(1, pt1);
    nodes[2] = std::make_shared<Node>(2, pt2);

    // Section
    auto area = 0.4;
    auto inertia11 = 1.0 / 30.0;
    auto inertia22 = 0.064 / 12;
    auto inertia12 = 0.015936;
    // auto inertia22 = 0.064 / 12.0;
    // auto inertia12 = 0;

    auto sect = std::make_shared<Section>(area, inertia11, inertia22, inertia12);

    // Material
    auto mat = std::make_shared<Material>(25e9, 0.0, 23536.8);

    // Fixities
    std::vector<bool> fix1 = { true, true, true, true, true, false };
    std::vector<bool> fix2 = { false, true, true, true, true, false };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };

    restraints[1] = std::make_shared<Restraint>(nodes[1], fix1, rest);
    restraints[2] = std::make_shared<Restraint>(nodes[2], fix2, rest);

    // Element
    elements[1] = std::make_shared<FrameMember>(1, nodes[1], nodes[2], sect, mat, false);

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    LOG(" Mass Matrix");
    str->MassMatrix->getSubmatrix(0, str->nUnrestrainedDOF - 1, 0, str->nUnrestrainedDOF - 1).printElements();
    LOG("");

    LOG(" Stiffness Matrix");
    str->StiffnessMatrix->getSubmatrix(0, str->nUnrestrainedDOF - 1, 0, str->nUnrestrainedDOF - 1).printElements();
    LOG("");

    // Modal periods
    auto modalPeriods = StructureSolver::GetModalPeriods(*str, SolverChoice::Armadillo);
    LOG(" Modal periods:");
    for (size_t i = 0; i < modalPeriods.RowCount; i++)
        if (modalPeriods(i, 0))
            std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods(i, 0) << " s\n";

    return;
}

void CE583Assignment9_2_3()
{
    // Fields
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Nodal coordinates
    XYZPoint pt1(0, 0, 0);
    XYZPoint pt2(6.0, 0, 0);
    XYZPoint pt3(12.0, 0, 0);
    XYZPoint pt4(18.0, 0, 0);
    XYZPoint pt5(24.0, 0, 0);
    XYZPoint pt6(30.0, 0, 0);

    // Nodes
    nodes[1] = std::make_shared<Node>(1, pt1);
    nodes[2] = std::make_shared<Node>(2, pt2);
    nodes[3] = std::make_shared<Node>(3, pt3);
    nodes[4] = std::make_shared<Node>(4, pt4);
    nodes[5] = std::make_shared<Node>(5, pt5);
    nodes[6] = std::make_shared<Node>(6, pt6);

    // Section
    auto area = 0.4;
    auto inertia11 = 1.0 / 30.0;
    auto inertia22 = 0.064 / 12;
    auto inertia12 = 0.015936;
    // auto inertia22 = 0.064 / 12.0;
    // auto inertia12 = 0;

    auto sect = std::make_shared<Section>(area, inertia11, inertia22, inertia12);

    // Material
    auto mat = std::make_shared<Material>(25e9, 0.0, 23536.8);

    // Fixities
    std::vector<bool> fix1 = { true, true, true, true, true, false };
    std::vector<bool> fix2 = { false, true, true, true, true, false };
    std::vector<bool> fix3 = { false, false, true, true, true, false };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };

    restraints[1] = std::make_shared<Restraint>(nodes[1], fix1, rest);
    restraints[2] = std::make_shared<Restraint>(nodes[2], fix3, rest);
    restraints[3] = std::make_shared<Restraint>(nodes[3], fix3, rest);
    restraints[4] = std::make_shared<Restraint>(nodes[4], fix3, rest);
    restraints[5] = std::make_shared<Restraint>(nodes[5], fix3, rest);
    restraints[6] = std::make_shared<Restraint>(nodes[6], fix2, rest);

    // Element
    elements[1] = std::make_shared<FrameMember>(1, nodes[1], nodes[2], sect, mat, false);
    elements[2] = std::make_shared<FrameMember>(2, nodes[2], nodes[3], sect, mat, false);
    elements[3] = std::make_shared<FrameMember>(3, nodes[3], nodes[4], sect, mat, false);
    elements[4] = std::make_shared<FrameMember>(4, nodes[4], nodes[5], sect, mat, false);
    elements[5] = std::make_shared<FrameMember>(5, nodes[5], nodes[6], sect, mat, false);

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    /*LOG(" Mass Matrix");
    str->MassMatrix->getSubmatrix(0, str->nUnrestrainedDOF - 1, 0, str->nUnrestrainedDOF - 1).printElements();
    LOG("");

    LOG(" Stiffness Matrix");
    str->StiffnessMatrix->getSubmatrix(0, str->nUnrestrainedDOF - 1, 0, str->nUnrestrainedDOF - 1).printElements();
    LOG("");*/

    // Modal periods
    auto modalPeriods = StructureSolver::GetModalPeriods(*str, SolverChoice::Armadillo);
    LOG(" Modal periods:");
    for (size_t i = 0; i < modalPeriods.RowCount; i++)
        if (modalPeriods(i, 0))
            std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods(i, 0) << " s\n";

    return;
}

void CE583Assignment9_3()
{
    // INPUTS
    auto speed = 5; // in km/hr
    auto dt = 0.1; // in seconds

    // Fields
    std::map<unsigned int, std::shared_ptr<Node>> nodes;
    std::map<unsigned int, std::shared_ptr<Element>> elements;
    std::map<unsigned int, std::shared_ptr<Restraint>> restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>> nodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>> distLoads;

    // Generate nodes and put restraints
    int nodeCounter = 0;
    int restraintCounter = 0;

    std::vector<bool> pin = { true, true, true, true, true, true };
    std::vector<bool> roller = { false, true, true, true, true, true };
    std::vector<bool> universal = { false, false, true, true, true, true };
    std::vector<double> rest = { 0, 0, 0, 0, 0, 0 };

    for (size_t i = 0; i < 2; i++)
    {
        double yCoord = i * 1.0;
        for (size_t j = 0; j < 31; j++)
        {
            nodeCounter++;
            double xCoord = j * 1.0;
            XYZPoint coord(xCoord, yCoord, 0);
            nodes[nodeCounter] = std::make_shared<Node>(nodeCounter, coord);

            // Check restraint
            if ((i == 1) && ((j == 0) || (j == 30))) // Top restraints
            {
                restraintCounter++;
                restraints[restraintCounter] = std::make_shared<Restraint>(nodes[nodeCounter], roller, rest);
            }
            else if ((i == 0) && (j == 30)) // Bottom right roller
            {
                restraintCounter++;
                restraints[restraintCounter] = std::make_shared<Restraint>(nodes[nodeCounter], roller, rest);
            }
            else if ((i == 0) && (j == 0)) // Bottom left pin
            {
                restraintCounter++;
                restraints[restraintCounter] = std::make_shared<Restraint>(nodes[nodeCounter], pin, rest);
            }
            else // Universal
            {
                restraintCounter++;
                restraints[restraintCounter] = std::make_shared<Restraint>(nodes[nodeCounter], universal, rest);
            }
        }
    }

    // Create material
    auto mat = std::make_shared<Material>(25e6, 0.0, 23.5368);

    // Place elements
    auto thickness = 0.40;
    auto membraneType = MembraneType::Incompatible;
    auto plateType = PlateType::NONE;

    for (size_t i = 1; i < 31; i++)
    {
        auto iNodeIdx = i;
        auto jNodeIdx = iNodeIdx + 1;
        auto kNodeIdx = jNodeIdx + 31;
        auto lNodeIdx = kNodeIdx - 1;
        elements[i] = std::make_shared<ShellMember>(i, nodes[iNodeIdx], nodes[jNodeIdx], nodes[kNodeIdx], nodes[lNodeIdx],
            mat, thickness, membraneType, plateType);
    }

    // Create structure
    auto str = std::make_shared<Structure>(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

    // Modal periods
    auto modalPeriods = StructureSolver::GetModalPeriods(*str, SolverChoice::Armadillo);
    LOG(" Modal periods:");
    //for (size_t i = 0; i < modalPeriods.RowCount; i++)
    for (size_t i = 0; i < 10; i++)
        std::cout << " Mode Number: " << i + 1 << ", Period = " << modalPeriods(i, 0) << " s\n";

    LOG("\n Modal frequencies:");
    //for (size_t i = 0; i < modalPeriods.RowCount; i++)
    for (size_t i = 0; i < 10; i++)
        std::cout << " Mode Number: " << i + 1 << ", Frequency = " << 2.0 * 3.141593 / modalPeriods(i, 0) << " rad/s\n";

    // Create force vector for all time steps
    // There is a 1 kN moving load with certain velocity. Create force vectors according to location
    auto v = speed / 3.6; // Speed in m/s
    auto tMax = 40.0;
    auto nStep = (int)(tMax / dt) + 1;
    std::vector<Matrix<double>> forceVectors(nStep);
    auto p = 1.0;
    for (size_t i = 0; i < nStep; i++)
    {
        Matrix<double> f(str->nDOF, 1);

        // If force is on a node, take it as point load on that node. If it is in between, 
        // use linear interpolation
        auto forceLocation = ((double)i * dt) * v;

        // Check nodes. There are 62 nodes.
        auto smallerIdx = -1;
        auto largerIdx = -1;
        auto nIdx = -1;

        for (size_t j = 32; j <= 62; j++)
        {
            auto&& n = nodes[j];

            if (Utils::AreEqual(n->Coordinate.X, forceLocation, 0.0001))
            {
                nIdx = n->NodeIndex;
                break;
            }

            if (forceLocation < n->Coordinate.X)
            {
                smallerIdx = n->NodeIndex - 1;
                largerIdx = smallerIdx + 1;
                nIdx = -1;
                break;
            }
        }

        if (nIdx != -1)
        {
            auto&& loadedNode = nodes[nIdx];
            f(loadedNode->DofIndexTY - 1, 0) = -1 * p;
        }
        else if (forceLocation <= 30)
        {
            auto&& leftNode = nodes[smallerIdx];
            auto&& rightNode = nodes[largerIdx];
            f(leftNode->DofIndexTY - 1, 0) = (forceLocation - rightNode->Coordinate.X) * p;
            f(rightNode->DofIndexTY - 1, 0) = (leftNode->Coordinate.X - forceLocation) * p;
        }

        forceVectors[i] = f;
    }

    // Calculate a0 and a1

    // Call solver
    auto wi = 2 * 3.14159262 / modalPeriods(0, 0);
    auto wj = 2 * 3.14159262 / modalPeriods(2, 0);

    // auto a0 = 0.0251;
    // auto a1 = 0.00013736;

    auto a0 = 0.01 * 2 * wi * wj / (wi + wj);
    auto a1 = 0.01 * 2 / (wi + wj);

    LOG("\n Dynamic Parameters");
    LOG(" a0 = " << a0);
    LOG(" a1 = " << a1);

    auto res = StructureSolver::ImplicitNewmark(*str, forceVectors, tMax, dt, a0, a1, SolverChoice::Armadillo);

    auto& disps = std::get<0>(res);

    auto midNodeVerticalDispIndex = nodes[16]->DofIndexTY - 1;
    auto absMaxDisp = 0;
    auto maxDisp = 0;
    auto t = 0.0;

    LOG("");
    for (size_t i = 0; i < disps.size(); i++)
    {
        auto& disp = disps[i];
        auto midSpanDisp = disp(midNodeVerticalDispIndex, 0);

        if (absMaxDisp < abs(midSpanDisp))
        {
            absMaxDisp = abs(midSpanDisp);
            maxDisp = midSpanDisp;
        }
        LOG(" " << t << "," << midSpanDisp);
        t += dt;
    }

    //LOG(" " << absMaxDisp);


    return;
}
