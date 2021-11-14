#include "..\include\XYData.h"

XYData::XYData() { }

XYData::XYData(ExtrapolationMethod extrpltnMethod)
    : extrapolationMethod(extrapolationMethod) { }

XYData::XYData(XYData& const that)
    : extrapolationMethod(that.extrapolationMethod),
    xData(that.xData), yData(that.yData) { }

XYData::~XYData() { }

void XYData::addItem(double x, double y)
{
    xData.push_back(x);
    yData.push_back(y);
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

    for (int i = 0; i < this->getCount() - 1; i++)
        area += 0.5 * (this->xData[i + 1] - this->xData[i]) * (this->yData[i + 1] + this->yData[i]);

    return area;
}

double XYData::getY(double x)
{
    return 0.0;
}

double XYData::getX(double y)
{
    return 0.0;
}

unsigned int XYData::getCount()
{
    return this->xData.size();
}
