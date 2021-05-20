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