#pragma once
#include <memory>
#include <vector>
#include <Map>
#include "Matrix.h"
#include "Node.h"

namespace ElmType {
    enum ElementType
    {
        None = 0,
        Frame = 1,
        Truss = 2,
        Shell = 3,
        SerendipityShell = 4
    };
}

class Element {
public:
    virtual std::vector<std::shared_ptr<Node>> GelElementNodes()
    {
        std::vector<std::shared_ptr<Node>> nodeList;
        return nodeList;
    };
    virtual unsigned int GetElementIndex() { return 0; }
    virtual ElmType::ElementType GetElementType() { return ElmType::ElementType::None; }
    virtual unsigned int GetNumberOfDoF() { return 0; };
    virtual std::shared_ptr<Matrix<double>> GetLocalCoordinateMassMatrix() { auto nullVal = std::make_shared<Matrix<double>>(); return nullVal; }
    virtual std::shared_ptr<Matrix<double>> GetLocalCoordinateStiffnessMatrix() { auto nullVal = std::make_shared<Matrix<double>>(); return nullVal; }
    virtual std::shared_ptr<Matrix<double>> GetLocalCoordinateDampingMatrix() { auto nullVal = std::make_shared<Matrix<double>>(); return nullVal; }
    virtual std::shared_ptr<Matrix<double>> GetRotationMatrix() { auto nullVal = std::make_shared<Matrix<double>>(); return nullVal; }
    virtual std::shared_ptr<Matrix<double>> GetGlobalCoordinateMassMatrix() { auto nullVal = std::make_shared<Matrix<double>>(); return nullVal; }
    virtual std::shared_ptr<Matrix<double>> GetGlobalCoordinateStiffnessMatrix() { auto nullVal = std::make_shared<Matrix<double>>(); return nullVal; }
    virtual std::shared_ptr<Matrix<double>> GetGlobalCoordinateDampingMatrix() { auto nullVal = std::make_shared<Matrix<double>>(); return nullVal; }
    virtual std::shared_ptr<Matrix<double>> GetElementLoads() { auto nullVal = std::make_shared<Matrix<double>>(); return nullVal; }
    virtual void GetElm() { };

    ElmType::ElementType Type;
};
