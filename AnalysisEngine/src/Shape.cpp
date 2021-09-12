#include "Shape.h"
#include "UtilMethods.h"

Shape::Shape()
{
    std::vector<XYPoint> v;
    this->vertices = v;
    std::vector<Shape> h;
    this->holes = h;
}

Shape::Shape(std::vector<XYPoint> vertices)
{
    this->vertices = vertices;
    std::vector<Shape> h;
    this->holes = h;
}

Shape::Shape(std::vector<XYPoint> vertices, std::vector<Shape> holes)
{
    this->vertices = vertices;
    this->holes = holes;
}

///<summary>
///Returns a collection of coordinates within shape with provided 
///interval, approximately
///</summary>
std::vector<XYPoint> Shape::getPoints(double interval)
{
    std::vector<XYPoint> retVal;

    auto bounds = this->getBoundingBoxLowerleftAndUpperright();

    auto lowerLeftX = std::get<0>(bounds).X;
    auto lowerLeftY = std::get<0>(bounds).Y;

    auto upperRightX = std::get<1>(bounds).X;
    auto upperRightY = std::get<1>(bounds).Y;

    for (auto y = lowerLeftY; !Utils::AreEqual(y, upperRightY + interval, interval * 0.1); y += interval)
    {
        for (auto x = lowerLeftX; !Utils::AreEqual(x, upperRightX + interval, interval * 0.1); x += interval)
        {
            XYPoint pt(x, y);

            if (this->isInside(pt))
            {
                retVal.push_back(pt);
            }
        }
    }


    return retVal;
}

bool Shape::isInside(XYPoint coord)
{
    //bool isInside(std::vector<std::tuple<double, double>> vertices, std::tuple<double, double> pt);

    std::vector<std::tuple<double, double>> points;
    std::tuple<double, double> pt = std::make_tuple(coord.X, coord.Y);

    for (auto&& item : vertices)
        points.push_back(std::make_tuple(item.X, item.Y));

    auto isInsideOuterShape = Utils::isInside(points, pt);
    auto isOnBound = this->isOnBoundary(coord);

    if (!holes.size())
        return isInsideOuterShape || isOnBound;

    for (auto&& hole : holes)
        if (hole.isInside(coord) && !hole.isOnBoundary(coord))
            return false;

    return true;
}

bool Shape::isOnBoundary(XYPoint coord)
{
    auto pt = std::make_tuple(coord.X, coord.Y);

    for (size_t i = 0; i < this->vertices.size() - 1; i++)
    {
        auto iEnd = std::make_tuple(this->vertices[i].X, this->vertices[i].Y);
        auto jEnd = std::make_tuple(this->vertices[i + 1].X, this->vertices[i + 1].Y);

        if (Utils::isPointOnLine(iEnd, jEnd, pt))
            return true;
    }

    auto firstPt = std::make_tuple(this->vertices[0].X, this->vertices[0].Y);
    auto lastPt = std::make_tuple(this->vertices[this->vertices.size() - 1].X, this->vertices[this->vertices.size() - 1].Y);

    return Utils::isPointOnLine(firstPt, lastPt, pt);
}

std::tuple<XYPoint, XYPoint> Shape::getBoundingBoxLowerleftAndUpperright()
{
    auto lowerLeftX = this->vertices[0].X;
    auto lowerLeftY = this->vertices[0].Y;

    auto upperRightX = this->vertices[0].X;
    auto upperRightY = this->vertices[0].Y;

    for (size_t i = 1; i < this->vertices.size(); i++)
    {
        auto x = this->vertices[i].X;
        auto y = this->vertices[i].Y;
        if (x < lowerLeftX) lowerLeftX = x;
        if (y < lowerLeftY) lowerLeftY = y;
        if (upperRightX < x) upperRightX = x;
        if (upperRightY < y) upperRightY = y;
    }

    XYPoint lowerLeft(lowerLeftX, lowerLeftY);
    XYPoint upperRight(upperRightX, upperRightY);

    return std::make_tuple(lowerLeft, upperRight);
}
