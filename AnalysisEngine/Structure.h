#pragma once
#include <Map>
#include <vector>
#include <memory>
#include "Node.h"
#include "Restraint.h"
#include "FrameMember.h"
#include "NodalLoad.h"
#include "DistributedLoad.h"
#include "Element.h"
#include "Hinge.h"

class Structure
{
public:
    std::map<unsigned int, std::shared_ptr<Node>>* Nodes;
    std::map<unsigned int, std::shared_ptr<Element>>* Elements;
    std::map<unsigned int, std::shared_ptr<Restraint>>* Restraints;
    std::map<unsigned int, std::shared_ptr<NodalLoad>>* NodalLoads;
    std::map<unsigned int, std::shared_ptr<DistributedLoad>>* DistributedLoads;

    std::vector<std::vector<double>> MassMatrix;
    std::vector<std::vector<double>> DampingMatrix;
    std::vector<std::vector<double>> StiffnessMatrix;
    std::vector<double> ForceVector;

    unsigned int nDOF;
    unsigned int nUnrestrainedDOF;

    Structure(std::map<unsigned int, std::shared_ptr<Node>>* nodeMap, std::map<unsigned int, std::shared_ptr<Element>>* elementMap, std::map<unsigned int,
        std::shared_ptr<Restraint>>*restraintMap, std::map<unsigned int, std::shared_ptr<NodalLoad>>* nodalLoadMap, std::map<unsigned int, 
        std::shared_ptr<DistributedLoad>>* distLoadMap);
    Structure();
    ~Structure();

private:
    void AssignDegreesOfFreedom(unsigned int& unrestDofCount, unsigned int& totalDofCount);
    void AssembleStiffnessMatrix(unsigned int totalDofCount);
    void AssembleMassMatrix(unsigned int totalDofCount);
    void AssembleForceVector(unsigned int totalDofCount);
};