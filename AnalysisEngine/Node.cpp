#include "Node.h"

Node::Node(unsigned int nodeIdx, XYZPoint pt)
{
	this->NodeIndex = nodeIdx;
	this->Coordinate = pt;
}

Node::~Node()
{

}

Node::Node()
{
	this->NodeIndex = 0;
	XYZPoint pt;
	this->Coordinate = pt;
}