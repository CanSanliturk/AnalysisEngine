#include "Structure.h"
#include <Map>
#include <vector>
#include <iostream>

Structure::Structure(std::map<unsigned int, Node*>* nodeMap, std::map<unsigned int, Element*>* elementMap, std::map<unsigned int, Restraint*>* restraintMap,
	std::map<unsigned int, Hinge*>* hingeMap, std::map<unsigned int, NodalLoad*>* nodalLoadMap, std::map<unsigned int, DistributedLoad*>* distLoadMap)
	: nDOF(0), nUnrestrainedDOF(0), Nodes(nodeMap), Elements(elementMap), Hinges(hingeMap), Restraints(restraintMap), NodalLoads(nodalLoadMap), 
	DistributedLoads(distLoadMap)
{
	unsigned int totalDofCount = 0;
	unsigned int unrestDofCount = 0;
	this->AssignDegreesOfFreedom(unrestDofCount, totalDofCount);
	this->AssembleStiffnessMatrix(totalDofCount);
	this->AssembleMassMatrix(totalDofCount);
	this->AssembleForceVector(totalDofCount);
}

Structure::Structure()
{
}

Structure::~Structure()
{
}

void Structure::AssignDegreesOfFreedom(unsigned int& unrestDofCount, unsigned int& totalDofCount)
{
	unsigned int dofIdx = 0;

	for (auto& nodePair : *this->Nodes)
	{
		auto node = nodePair.second;
		
		bool isTransXRest = false;
		bool isTransYRest = false;
		bool isTransZRest = false;
		bool isRotXRest = false;
		bool isRotYRest = false;
		bool isRotZRest = false;

		for (auto& restPair : *this->Restraints)
		{
			auto rest = restPair.second;
			if (rest->RestrainedNode->NodeIndex != nodePair.first)
				continue; // Look for other restraint

			isTransXRest = isTransXRest || (rest->IsRestraintTranslationX);
			isTransYRest = isTransYRest || (rest->IsRestraintTranslationY);
			isTransZRest = isTransZRest || (rest->IsRestraintTranslationZ);
			isRotXRest = isRotXRest || (rest->IsRestraintRotationX);
			isRotYRest = isRotYRest || (rest->IsRestraintRotationY);
			isRotZRest = isRotZRest || (rest->IsRestraintRotationZ);
		}

		if (!isTransXRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexTX = dofIdx;
		}

		if (!isTransYRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexTY = dofIdx;
		}

		if (!isTransZRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexTZ = dofIdx;
		}

		if (!isRotXRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexRX = dofIdx;
		}

		if (!isRotYRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexRY = dofIdx;
		}

		if (!isRotZRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexRZ = dofIdx;
		}
	}

	unrestDofCount = dofIdx;
	this->nUnrestrainedDOF = unrestDofCount;

	for (auto nodePair : *this->Nodes)
	{
		auto node = nodePair.second;

		bool isTransXRest = false;
		bool isTransYRest = false;
		bool isTransZRest = false;
		bool isRotXRest = false;
		bool isRotYRest = false;
		bool isRotZRest = false;

		for (auto& restPair : *this->Restraints)
		{
			auto rest = restPair.second;
			if (rest->RestrainedNode->NodeIndex != nodePair.first)
				continue; // Look for other restraint

			isTransXRest = isTransXRest || (rest->IsRestraintTranslationX);
			isTransYRest = isTransYRest || (rest->IsRestraintTranslationY);
			isTransZRest = isTransZRest || (rest->IsRestraintTranslationZ);
			isRotXRest = isRotXRest || (rest->IsRestraintRotationX);
			isRotYRest = isRotYRest || (rest->IsRestraintRotationY);
			isRotZRest = isRotZRest || (rest->IsRestraintRotationZ);
		}

		if (isTransXRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexTX = dofIdx;
		}

		if (isTransYRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexTY = dofIdx;
		}

		if (isTransZRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexTZ = dofIdx;
		}

		if (isRotXRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexRX = dofIdx;
		}

		if (isRotYRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexRY = dofIdx;
		}

		if (isRotZRest)
		{
			dofIdx++;
			this->Nodes->at(nodePair.first)->DofIndexRZ = dofIdx;
		}
	}

	totalDofCount = dofIdx;
	this->nDOF = totalDofCount;
}

void Structure::AssembleStiffnessMatrix(unsigned int totalDofCount)
{
	std::vector<double> zeroVector(totalDofCount, 0.0);
	for (unsigned int i = 0; i < totalDofCount; i++) this->StiffnessMatrix.push_back(zeroVector);

	for (auto& elmPair : *this->Elements)
	{
		auto elm = elmPair.second;
		auto nodes = elm->GelElementNodes();
		std::vector<unsigned int> steerVector;

		for (auto& n : nodes)
		{
			steerVector.push_back(n->DofIndexTX);
			steerVector.push_back(n->DofIndexTY);
			steerVector.push_back(n->DofIndexTZ);
			steerVector.push_back(n->DofIndexRX);
			steerVector.push_back(n->DofIndexRY);
			steerVector.push_back(n->DofIndexRZ);
		}

		auto nDof = elm->GetNumberOfDoF();
		auto elmStiffnessMat = (double*)elm->GetGlobalCoordinateStiffnessMatrix();

		for (size_t i = 0; i < nDof; i++)
		{
			for (size_t j = 0; j < nDof; j++)
			{
				auto idx1 = steerVector.at(i);
				auto idx2 = steerVector.at(j);
				this->StiffnessMatrix.at(idx1 - 1).at(idx2 - 1) += *elmStiffnessMat;
				elmStiffnessMat++;
			}
		}
	}
}

void Structure::AssembleMassMatrix(unsigned int totalDofCount)
{
	std::vector<double> zeroVector(totalDofCount, 0.0);
	for (unsigned int i = 0; i < totalDofCount; i++) this->MassMatrix.push_back(zeroVector);

	for (auto& elmPair : *this->Elements)
	{
		auto elm = elmPair.second;
		auto nodes = elm->GelElementNodes();
		std::vector<unsigned int> steerVector;

		for (auto n : nodes)
		{
			steerVector.push_back(n->DofIndexTX);
			steerVector.push_back(n->DofIndexTY);
			steerVector.push_back(n->DofIndexTZ);
			steerVector.push_back(n->DofIndexRX);
			steerVector.push_back(n->DofIndexRY);
			steerVector.push_back(n->DofIndexRZ);
		}

		auto nDof = elm->GetNumberOfDoF();
		auto elmMassMat = (double*)elm->GetGlobalCoordinateMassMatrix();

		for (size_t i = 0; i < nDof; i++)
		{
			for (size_t j = 0; j < nDof; j++)
			{
				auto idx1 = steerVector.at(i);
				auto idx2 = steerVector.at(j);
				this->MassMatrix.at(idx1 - 1).at(idx2 - 1) += *elmMassMat;
				elmMassMat++;
			}
		}
	}
}

void Structure::AssembleForceVector(unsigned int totalDofCount)
{
	for (size_t i = 0; i < totalDofCount; i++) this->ForceVector.push_back(0.0);
	
	for (auto& nodalLoadPair : *this->NodalLoads)
	{
		auto load = nodalLoadPair.second;
		auto node = load->ActingNode;

		this->ForceVector[node->DofIndexTX - 1] = load->Loads[0];
		this->ForceVector[node->DofIndexTY - 1] = load->Loads[1];
		this->ForceVector[node->DofIndexTZ - 1] = load->Loads[2];
		this->ForceVector[node->DofIndexRX - 1] = load->Loads[3];
		this->ForceVector[node->DofIndexRY - 1] = load->Loads[4];
		this->ForceVector[node->DofIndexRZ - 1] = load->Loads[5];
	}

}
