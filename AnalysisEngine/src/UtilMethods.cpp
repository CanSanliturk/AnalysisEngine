#include "UtilMethods.h"

bool Utils::AreEqual(const double& firstVal, const double& secondVal, const double& tolerance)
{
    return ((abs(abs(firstVal) - abs(secondVal))) <= (abs(tolerance)));
}

bool Utils::AreEqual(const double& firstVal, const double& secondVal)
{
    auto tol = abs(firstVal) / 500;
    return Utils::AreEqual(firstVal, secondVal, tol);
}

bool Utils::isInside(std::vector<std::tuple<double, double>> vertices, std::tuple<double, double> pt)
{
    auto c = false;
    auto x = std::get<0>(pt);
    auto y = std::get<1>(pt);
   
    for (int m = 0, n = vertices.size() - 1; m < vertices.size(); n = m++)
        if (((std::get<0>(vertices[m]) > x) != (std::get<0>(vertices[n]) > x)) &&
            (y < (std::get<1>(vertices[n]) - std::get<1>(vertices[m])) * (x - std::get<0>(vertices[m])) / (std::get<0>(vertices[n]) - std::get<0>(vertices[m])) + std::get<1>(vertices[m])))
            c = !c;

    return c;
}

bool Utils::isPointOnLine(std::tuple<double, double> iEnd, std::tuple<double, double> jEnd, std::tuple<double, double> pt)
{
    auto iX = std::get<0>(iEnd);
    auto iY = std::get<1>(iEnd);
    auto jX = std::get<0>(jEnd);
    auto jY = std::get<1>(jEnd);
    auto pX = std::get<0>(pt);
    auto pY = std::get<1>(pt);

    if (((Utils::AreEqual(iX, pX)) && (Utils::AreEqual(iY, pY))) ||
        ((Utils::AreEqual(jX, pX)) && (Utils::AreEqual(jY, pY))))
        return true;

    auto originalLineSlope = atan2(jY - iY, jX - iX);
    auto lineSlope = atan2(pY- iY, pX - iX);

    if (!Utils::AreEqual(originalLineSlope, lineSlope))
        return false;

    auto originalDist = sqrt(((jX - iX) * (jX - iX)) + ((jY - iY) * (jY - iY)));
    auto dist = sqrt(((pX - iX) * (pX - iX)) + ((pY - iY) * (pY - iY)));

    return (dist < originalDist) || (Utils::AreEqual(originalDist, dist));
}


