#pragma once

#include <vector>
#include <memory>
#include "Node.h"

class Restraint
{
public:
    std::shared_ptr<Node> RestrainedNode;

    double TranslationX;
    double TranslationY;
    double TranslationZ;
    double RotationX;
    double RotationY;
    double RotationZ;

    bool IsRestraintTranslationX = false;
    bool IsRestraintTranslationY = false;
    bool IsRestraintTranslationZ = false;
    bool IsRestraintRotationX = false;
    bool IsRestraintRotationY = false;
    bool IsRestraintRotationZ = false;

    std::vector<bool> IsRestrainedVector;
    std::vector<double> RestrainedCondition;


    Restraint(std::shared_ptr<Node> resNode, std::vector<bool> isRest, std::vector<double> rest)
        : RestrainedNode(resNode), 
        IsRestraintTranslationX(isRest.at(0)), IsRestraintTranslationY(isRest.at(1)), IsRestraintTranslationZ(isRest.at(2)),
        IsRestraintRotationX(isRest.at(3)), IsRestraintRotationY(isRest.at(4)), IsRestraintRotationZ(isRest.at(5)),
        TranslationX(rest.at(0)), TranslationY(rest.at(1)), TranslationZ(rest.at(2)),
        RotationX(rest.at(3)), RotationY(rest.at(4)), RotationZ(rest.at(5)), IsRestrainedVector(isRest), RestrainedCondition(rest)
    { }

    Restraint()
        : RestrainedNode(std::make_shared<Node>()),
        IsRestraintTranslationX(false), IsRestraintTranslationY(false), IsRestraintTranslationZ(false),
        IsRestraintRotationX(false), IsRestraintRotationY(false), IsRestraintRotationZ(false),
        TranslationX(0.0), TranslationY(0.0), TranslationZ(0.0),
        RotationX(0.0), RotationY(0.0), RotationZ(0.0)
    { };

    ~Restraint()
    { };

};