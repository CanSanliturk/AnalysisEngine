#pragma once
#include <vector>
#include <memory>
#include <Map>
#include "Node.h"
#include "SerendipityShell.h"
#include "DistributedLoad.h"
#include "FrameMember.h"
#include "TrussMember.h"
#include "Restraint.h"
#include "NodalLoad.h"
#include "Element.h"
#include "Matrix.h"
#include "Hinge.h"
#include "Shell.h"

class Structure
{
public:
    std::map<unsigned int, std::shared_ptr<Node>>* Nodes;
    std::map<unsigned int, std::shared_ptr<Element>>* Elements;
    std::map<unsigned int, std::shared_ptr<Restraint>>* Restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>>* NodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>>* DistributedLoads;

    std::shared_ptr<Matrix<double>> MassMatrix;
    std::shared_ptr<Matrix<double>> StiffnessMatrix;
    std::shared_ptr<Matrix<double>> ForceVector;

    unsigned int nDOF;
    unsigned int nUnrestrainedDOF;

    Structure(std::map<unsigned int, std::shared_ptr<Node>>* nodeMap, std::map<unsigned int, std::shared_ptr<Element>>* elementMap, std::map<unsigned int,
        std::shared_ptr<Restraint>>*restraintMap, std::map<unsigned int, std::shared_ptr<NodalLoad>>* nodalLoadMap, std::map<unsigned int, 
        std::shared_ptr<DistributedLoad>>* distLoadMap);
    Structure();
    ~Structure();

    void updateStiffnessMatrix();
    void updateDofIndicesAndMatrices();
    void updateMassMatrix();
    Matrix<double> getForceVector(std::map<unsigned int, std::shared_ptr<NodalLoad>>* nodalLoadMap);    
    std::shared_ptr<Node> getNodeAt(XYZPoint coord);
    std::shared_ptr<Node> getNodeAt(double x, double y, double z);

private:
    void AssignDegreesOfFreedom(unsigned int& unrestDofCount, unsigned int& totalDofCount);
    void AssembleMassMatrix(unsigned int totalDofCount);
    void AssembleStiffnessMatrix(unsigned int totalDofCount);
    void AssembleForceVector(unsigned int totalDofCount);
};