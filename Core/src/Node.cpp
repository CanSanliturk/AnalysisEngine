#include "Node.h"

Node::Node(unsigned int nodeIdx, XYZPoint pt)
    : NodeIndex(nodeIdx), Coordinate(pt)
{ }

Node::Node()
    : NodeIndex(0)
{
    XYZPoint pt;
    this->Coordinate = pt;
}

Node::~Node() 
{

}