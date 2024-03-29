#pragma once
#include "Node.h"
#include <vector>
#include <memory>

class NodalLoad
{
public:
    std::shared_ptr<Node> ActingNode;
    double Loads[6] = { 0, 0, 0, 0, 0, 0 }; // Forces and moments in X, Y and Z directions

    NodalLoad(std::shared_ptr<Node> node, double loads[6])
        : ActingNode(node)
    {
        this->Loads[0] = loads[0];
        this->Loads[1] = loads[1];
        this->Loads[2] = loads[2];
        this->Loads[3] = loads[3];
        this->Loads[4] = loads[4];
        this->Loads[5] = loads[5];
    };

    NodalLoad()
        : ActingNode(std::make_shared<Node>())
    { };

    ~NodalLoad()
    { };

private:
};