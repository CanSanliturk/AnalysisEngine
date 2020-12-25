#pragma once
#include <Map>
#include <vector>
#include "Node.h"
#include "Restraint.h"
#include "FrameMember.h"
#include "NodalLoad.h"
#include "DistributedLoad.h"
#include "Element.h"

class Structure
{
public:	
 	std::map<unsigned int, Node*>* Nodes;
	std::map<unsigned int, Element*>* Elements;
	std::map<unsigned int, Restraint*>* Restraints;
	std::map<unsigned int, NodalLoad*>* NodalLoads;
	std::map<unsigned int, DistributedLoad*>* DistributedLoads;

	std::vector<std::vector<double>> StiffnessMatrix;
	std::vector<std::vector<double>> MassMatrix;
	std::vector<double> ForceVector;

	Structure(std::map<unsigned int, Node*>* nodeMap, std::map<unsigned int, Element*>* elementMap, std::map<unsigned int, Restraint*>* restraintMap, 
		std::map<unsigned int, NodalLoad*>* nodalLoadMap, std::map<unsigned int, DistributedLoad*>* distLoadMap);
	Structure();
	~Structure();
	

private:
	void AssignDegreesOfFreedom(unsigned int &unrestDofCount, unsigned int &totalDofCount);
	void AssembleStiffnessMatrix(unsigned int totalDofCount);
	void AssembleMassMatrix(unsigned int totalDofCount);
	void AssembleForceVector(unsigned int totalDofCount);
};

