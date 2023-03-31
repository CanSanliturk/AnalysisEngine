#pragma once
#include <cmath>
#include <vector>
#include <tuple>

namespace Utils
{
    /// <summary>
    /// A comparer function that checks equality of two double-precision numbers using provided 
    /// tolarance
    /// </summary>
    bool AreEqual(const double& firstVal, const double& secondVal, const double& tolerance);

    /// <summary>
    /// A comparer function that checks equality of two doubles using a tolerance which is equal to
    /// 1/500 of first value
    /// </summary>
    bool AreEqual(const double& firstVal, const double& secondVal);

    /// <summary>
    /// Method that returns whether provided coordinate is in given vertices
    /// </summary>
    bool isInside(std::vector<std::tuple<double, double>> vertices, std::tuple<double, double> pt);

    bool isPointOnLine(std::tuple<double, double> iEnd, std::tuple<double, double> jEnd, std::tuple<double, double> pt);
}