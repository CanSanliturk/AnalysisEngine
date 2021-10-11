#pragma once
#include "Vector.h"
#include <vector>
#include "XYPoint.h"

class Shape {

public:
    Shape();
    Shape(std::vector<XYPoint> vertices);
    Shape(std::vector<XYPoint> vertices, std::vector<Shape> holes);

    /// <summary>
    /// Returns a collection of coordinates within shape with provided 
    /// interval, approximately
    /// </summary>
    std::vector<XYPoint> getPoints(double interval);

    bool isInside(XYPoint coord);
    bool isOnBoundary(XYPoint coord);
    short getNumberOfHoles();
    std::tuple<XYPoint, XYPoint> getBoundingBoxLowerLeftAndUpperRight();
    std::vector<Shape> getHoles();

private:
    std::vector<XYPoint> vertices;
    std::vector<Shape> holes;
};