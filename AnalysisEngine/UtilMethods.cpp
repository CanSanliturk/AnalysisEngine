#include "UtilMethods.h"

bool Utils::AreEqual(double &firstVal, double &secondVal, double &tolerance)
{
	return ((abs(abs(firstVal) - abs(secondVal))) <= (abs(tolerance)));
}