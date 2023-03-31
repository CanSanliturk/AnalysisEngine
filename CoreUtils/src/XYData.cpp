#include <XYData.h>
#include <iostream>

XYData::XYData() { }

XYData::XYData(ExtrapolationMethod extrpltnMethod)
    : extrapolationMethod(extrpltnMethod){ }

XYData::XYData(XYData& that)
    : extrapolationMethod(that.extrapolationMethod),
    xData(that.xData), yData(that.yData) { }

XYData::~XYData() { }

void XYData::addItem(double x, double y)
{
    xData.push_back(x);
    yData.push_back(y);
}

double XYData::getMaxNegative() {
    auto maxNeg = 0.0;
    auto i = 0;
    while ((i < xData.size()) && (xData[i] < 0))  {
        maxNeg = xData[i];
        ++i;
    }
    return maxNeg;
}

double XYData::getMinPositive(){
    auto minPos = xData[xData.size() - 1];
    auto i = xData.size() - 1;
    while ((0 <= i) && (xData[i] > 0))  {
        minPos = xData[i];
        --i;
    }
    return minPos;
}

double XYData::integrate(double xStart, double xEnd)
{
    auto area = 0.0;

    /*/// <summary>
    /// Perform calculation if integration interval starts before data
    /// </summary>
    if (xStart < this->xData[0])
    {
        switch (extrapolationMethod)
        {
        case ExtrapolationMethod::Zero:
            break;
        case ExtrapolationMethod::Horizontal:
            area += (this->xData[0] - xStart) * this->yData[0];
            break;
        default:
            break;
        }
    }

    /// <summary>
    /// Perform calculation if integration interval ends after data
    /// </summary>
    if (this->xData[this->getCount() - 1] < xEnd)
    {
        switch (extrapolationMethod)
        {
        case ExtrapolationMethod::Zero:
            break;
        case ExtrapolationMethod::Horizontal:
            area += (xEnd - this->xData[this->getCount() - 1]) * this->yData[this->getCount() - 1];
            break;
        default:
            break;
        }
    }*/

    return area;
}

double XYData::integrate()
{
    auto area = 0.0;

    for (unsigned int i = 0; i < this->getCount() - 1; i++)
        area += 0.5 * (this->xData[i + 1] - this->xData[i]) * (this->yData[i + 1] + this->yData[i]);

    return area;
}

double XYData::getY(double x)
{
    auto y = 0.0;
    auto count = getCount();
    auto firstX = xData[0];
    auto lastX = xData[count - 1];

    if (x < xData[0])
    {
        // Case where the given x is smaller than the first data point
        switch (extrapolationMethod)
        {
        case ExtrapolationMethod::Zero:
            // Do not update return value
            break;
        case ExtrapolationMethod::Horizontal:
            y = yData[0];
            break;
        default:
            break;
        }
    }
    else if (xData[getCount() - 1] < x)
    {
        // Case where the given x is greater than the last data point
        switch (extrapolationMethod)
        {
        case ExtrapolationMethod::Zero:
            // Do not update the return value
            break;
        case ExtrapolationMethod::Horizontal:
            y = yData[getCount() - 1];
            break;
        default:
            break;
        }
    }
    else
    {
        // Case where the given x is within the data interval
        for (size_t i = 1; i < getCount(); i++)
        {
            if (x < xData[i])
            {
                y = yData[i - 1] + ((yData[i] - yData[i - 1]) * ((x - xData[i - 1]) / (xData[i] - xData[i - 1])));
                break;
            }
        }
    }

    return y;
}

double XYData::getX(double y)
{
    return 0.0;
}

unsigned int XYData::getCount()
{
    return this->xData.size();
}

const std::vector<double> XYData::getXData()
{
    return xData;
}

const std::vector<double> XYData::getYData()
{
    return yData;
}