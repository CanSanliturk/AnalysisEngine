#pragma once

#include <vector>
#include <Map>
#include "Node.h"

namespace ElmType {
	enum ElementType
	{
		None = 0,
		Frame = 1
	};
}

class Element {
public:
	virtual std::vector<Node*> GelElementNodes() 
	{
		Node n;
		std::vector<Node*> asd;
		asd.push_back(&n);
		return asd; 
	};
	virtual unsigned int GetElementIndex() { return 0; };
	virtual ElmType::ElementType GetElementType() { return ElmType::ElementType::None; };
	virtual unsigned int GetNumberOfDoF() { return 0; };
	virtual void* GetLocalCoordinateStiffnessMatrix() { void* vPtr = 0x00000000;  return vPtr; };
	virtual void* GetLocalCoordinateMassMatrix() { void* vPtr = 0x00000000;  return vPtr; };
	virtual void* GetGlobalCoordinateStiffnessMatrix() { void* vPtr = 0x00000000;  return vPtr; };
	virtual void* GetGlobalCoordinateMassMatrix() { void* vPtr = 0x00000000;  return vPtr; };
	virtual void* GetRotationMatrix() { void* vPtr = 0x00000000;  return vPtr; };

	virtual void GetElm() { };

	ElmType::ElementType Type;
};
