#pragma once
#include "Node.h"
#include <vector>

class DistributedLoad
{
public:
    Node* StartNode;
    Node* EndNode;
    double InitialValue[6];
    double FinalValue[6];

    DistributedLoad();
    DistributedLoad(Node* startNode, Node* endNode, double initialValue[6], double finalValue[6]);
    ~DistributedLoad();

private:


};