#pragma once
#include <vector>

enum class ExtrapolationMethod {
    /// <summary>
    /// Out of the graph has zero value
    /// </summary>
    Zero = 0,

    /// <summary>
    /// Extends to positive and negative infinity with first and last values.
    /// </summary>
    Horizontal = 1
};


class XYData
{
public:
    /// <summary>
    /// Default constructor
    /// </summary>
    XYData();

    /// <summary>
    /// Constructor
    /// </summary>
    XYData(ExtrapolationMethod extrpltnMethod);

    /// <summary>
    /// Construtor
    /// </summary>
    XYData(XYData& that);

    /// <summary>
    /// Destructor
    /// </summary>
    ~XYData();

    /// <summary>
    /// Adds item
    /// </summary>
    void addItem(double x, double y);

    /// <summary>
    /// Gets integral from given start point to given end point using provided
    /// extrapolation method
    /// </summary>
    double integrate(double xStart, double xEnd);

    /// <summary>
    /// Gets integral from start point of data to end point of data
    /// </summary>
    double integrate();

    /// <summary>
    /// Returns y for given x using linear interpolation if given x is between the range of
    /// data and and given extrapolation method if given x is out of the range of the data
    /// </summary>
    double getY(double x);

    /// <summary>
    /// Returns x for given y using linear interpolation if given y is between the range of
    /// data and and given extrapolation method if given y is out of the range of the data
    /// </summary>
    double getX(double y);

    /// <summary>
    /// Gets the number of items of XY-data
    /// </summary>
    /// <returns></returns>
    unsigned int getCount();

    const std::vector<double> getXData();
    const std::vector<double> getYData();

    double getMaxNegative();
    double getMinPositive();

private:
    std::vector<double> xData;
    std::vector<double> yData;
    ExtrapolationMethod extrapolationMethod = static_cast<ExtrapolationMethod>(0);
};