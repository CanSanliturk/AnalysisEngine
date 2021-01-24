#include "..\AnalysisEngine\ArmadilloSolver.h"
#include "..\AnalysisEngine\DistributedLoad.h"
#include "..\AnalysisEngine\EigenSolver.h"
#include "..\AnalysisEngine\Element.h"
#include "..\AnalysisEngine\FrameMember.h"
#include "..\AnalysisEngine\Material.h"
#include "..\AnalysisEngine\MatrixHelper.h"
#include "..\AnalysisEngine\NodalLoad.h"
#include "..\AnalysisEngine\Node.h"
#include "..\AnalysisEngine\Piece.h"
#include "..\AnalysisEngine\Restraint.h"
#include "..\AnalysisEngine\Section.h"
#include "..\AnalysisEngine\Structure.h"
#include "..\AnalysisEngine\TrussMember.h"
#include "..\AnalysisEngine\Vector.h"
#include "..\AnalysisEngine\XYPoint.h"
#include "..\AnalysisEngine\XYZPoint.h"
#include "pch.h"
#include "CppUnitTest.h"
#include <iostream>

#define LOG(x) std::cout << x << "\n"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace EngineTests
{
	TEST_CLASS(EntitiyTests)
	{
	public:

		TEST_METHOD(TableStructure)
		{
			// Structure used in example is a table. 
			// All elements are 5 meter long and have 
			// cross-section of 0.5m*0.5m with E = 200 GPa.
			// (Units are in N & m)

			// Coordinates
			XYZPoint bottomPt1(0, 0, 0); // Origin
			XYZPoint bottomPt2(5, 0, 0);
			XYZPoint bottomPt3(5, 5, 0);
			XYZPoint bottomPt4(5, 0, 0);
			XYZPoint topPt1(0, 0, 5); // Top of the origin
			XYZPoint topPt2(5, 0, 5);
			XYZPoint topPt3(5, 5, 5);
			XYZPoint topPt4(5, 0, 5);

			// Nodes
			std::map<unsigned int, Node*> nodes;
			Node bottomNode1(1, bottomPt1); nodes[bottomNode1.NodeIndex] = &bottomNode1;
			Node bottomNode2(2, bottomPt2); nodes[bottomNode2.NodeIndex] = &bottomNode2;
			Node bottomNode3(3, bottomPt3);	nodes[bottomNode3.NodeIndex] = &bottomNode3;
			Node bottomNode4(4, bottomPt4);	nodes[bottomNode4.NodeIndex] = &bottomNode4;
			Node topNode1(5, topPt1); nodes[topNode1.NodeIndex] = &topNode1;
			Node topNode2(6, topPt2); nodes[topNode2.NodeIndex] = &topNode2;
			Node topNode3(7, topPt3); nodes[topNode3.NodeIndex] = &topNode3;
			Node topNode4(8, topPt4); nodes[topNode4.NodeIndex] = &topNode4;

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
			FrameMember col1(1, &bottomNode1, &topNode1, &sect, &mat); elements[col1.ElementIndex] = &col1;
			FrameMember col2(2, &bottomNode2, &topNode2, &sect, &mat); elements[col2.ElementIndex] = &col2;
			FrameMember col3(3, &bottomNode3, &topNode3, &sect, &mat); elements[col3.ElementIndex] = &col3;
			FrameMember col4(4, &bottomNode4, &topNode4, &sect, &mat); elements[col4.ElementIndex] = &col4;
			FrameMember beam1(5, &topNode1, &topNode2, &sect, &mat); elements[beam1.ElementIndex] = &beam1;
			FrameMember beam2(6, &topNode2, &topNode3, &sect, &mat); elements[beam2.ElementIndex] = &beam2;
			FrameMember beam3(7, &topNode3, &topNode4, &sect, &mat); elements[beam3.ElementIndex] = &beam3;
			FrameMember beam4(7, &topNode4, &topNode1, &sect, &mat); elements[beam4.ElementIndex] = &beam4;

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
			NodalLoad nl1(&topNode1, nodalLoad); nodalLoads[1] = &nl1;
			NodalLoad nl2(&topNode4, nodalLoad); nodalLoads[2] = &nl2;

			// Distributed loads
			std::map<unsigned int, DistributedLoad*> distLoads;

			// Create structure
			Structure str(&nodes, &elements, &restraints, &nodalLoads, &distLoads);

			// Solve displacement
			auto displacements = ArmadilloSolver::GetDisplacementForStaticCase(str);
			for (auto& nodePair : nodes)
			{
				auto node = nodePair.second;

				LOG("");
				LOG(" Node Index: ");
				std::cout << " " << node->NodeIndex << "\n";

				auto nodalDisps = ArmadilloSolver::GetNodalDisplacements(*node, displacements);

				for (size_t i = 0; i < 6; i++)
					std::cout << " DOF Index: " << i + 1 << ", Displacement = " << nodalDisps[i] << "\n";
			}

			Assert::IsTrue(true);
		}
	};
}
