#pragma once
#include "XYZPoint.h"
#include <vector>

class Node
{
public:
	unsigned int NodeIndex;

	unsigned int DofIndexTX = 0;
	unsigned int DofIndexTY = 0;
	unsigned int DofIndexTZ = 0;
	unsigned int DofIndexRX = 0;
	unsigned int DofIndexRY = 0;
	unsigned int DofIndexRZ = 0;
	std::vector<unsigned int> ConnectedElements;

	XYZPoint Coordinate;

	Node(unsigned int NodeIndex, XYZPoint pt);
	Node();
	~Node();

private:
};