#include <iostream>
#include <vector>
#include <math.h>
#include "MatrixHelper.h"
#include "UtilMethods.h"

Matrix<double> MatrixHelper::GetRotationMatrix(double thetaX, double thetaY, double thetaZ)
{
    auto rotX = MatrixHelper::GetRotationMatrixAboutX(thetaX);
    auto rotY = MatrixHelper::GetRotationMatrixAboutY(thetaY);
    auto rotZ = MatrixHelper::GetRotationMatrixAboutZ(thetaZ);
    return ((rotZ * rotY) * rotX);
}

Matrix<double> MatrixHelper::GetTranslationalRotationMatrix(Vector elmVector, double rotationAngle)
{
    // Skew angle is ignored. Create section in a skewed way for now.
    //auto alpha = rotationAngle;

    // For vertical members (i.e., both cX and cY are zero, different approach should be followed)
    auto cX = elmVector.X / elmVector.Length;
    auto cY = elmVector.Y / elmVector.Length;
    auto cZ = elmVector.Z / elmVector.Length;
    auto cXY = sqrt((cX * cX) + (cY * cY));

    std::vector<double> firstRow;
    std::vector<double> secondRow;
    std::vector<double> thirdRow;

    double compareVal = 0.0;
    double alpha = 0.0;
    double tol = 1e-8;

    Matrix<double> retVal(3);

    if (Utils::AreEqual(cXY, compareVal, tol)) // Vertical members
    {
        retVal(0, 0) = 0.0;
        retVal(0, 1) = cZ;
        retVal(0, 2) = 0.0;

        retVal(1, 0) = -1 * cZ;
        retVal(1, 1) = 0;
        retVal(1, 2) = 0;

        retVal(2, 0) = 0.0;
        retVal(2, 1) = 0.0;
        retVal(2, 2) = 1.0;
    }
    else // General case
    {
        retVal(0, 0) = cX;
        retVal(0, 1) = cZ;
        retVal(0, 2) = cY;

        retVal(1, 0) = -1 * ((cX * cZ * cos(alpha)) + (cY * sin(alpha))) / cXY;
        retVal(1, 1) = cXY * cos(alpha);
        retVal(1, 2) = ((-cZ * cY * cos(alpha)) + (cX * sin(alpha))) / cXY;

        retVal(2, 0) = ((cX * cZ * sin(alpha)) - (cY * cos(alpha))) / cXY;
        retVal(2, 1) = -1 * cXY * sin(alpha);
        retVal(2, 2) = ((cZ * cY * sin(alpha)) + (cX * cos(alpha))) / cXY;
    }

    return retVal;
}

Matrix<double> MatrixHelper::GetTranspose(Matrix<double>& matrix)
{
    auto rowCount = matrix.RowCount;
    auto colCount = matrix.ColCount;

    Matrix<double> trans(colCount, rowCount);

    for (unsigned int colIdx = 0; colIdx < colCount; colIdx++)
        for (unsigned int rowIdx = 0; rowIdx < rowCount; rowIdx++)
            trans(colIdx, rowIdx) = matrix(rowIdx, colIdx);

    return trans;
}

Matrix<double> MatrixHelper::GetRotationMatrixAboutX(double thetaX)
{
    // Define rotation matrix
    Matrix<double> rotAboutX(3);

    // Define first row
    rotAboutX(0, 0) = 1;
    rotAboutX(0, 1) = 0;
    rotAboutX(0, 2) = 0;

    // Define second row
    rotAboutX(1, 0) = 0;
    rotAboutX(1, 1) = cos(thetaX);
    rotAboutX(1, 2) = -1 * sin(thetaX);

    // Define third row
    rotAboutX(2, 0) = 0;
    rotAboutX(2, 1) = sin(thetaX);
    rotAboutX(2, 2) = cos(thetaX);

    return rotAboutX;
}

Matrix<double> MatrixHelper::GetRotationMatrixAboutY(double thetaY)
{
    // Define rotation matrix
    Matrix<double> rotAboutY(3);

    // Define first row
    rotAboutY(0, 0) = cos(thetaY);
    rotAboutY(0, 1) = 0;
    rotAboutY(0, 2) = sin(thetaY);

    // Define second row
    rotAboutY(1, 0) = 0;
    rotAboutY(1, 1) = 1;
    rotAboutY(1, 2) = 0;

    // Define third row
    rotAboutY(2, 0) = -1 * sin(thetaY);
    rotAboutY(2, 1) = 0;
    rotAboutY(2, 2) = cos(thetaY);

    return rotAboutY;
}

Matrix<double> MatrixHelper::GetRotationMatrixAboutZ(double thetaZ)
{
    // Define rotation matrix
    Matrix<double> rotAboutZ(3);

    // Define first row
    std::vector<double> firstRow;
    rotAboutZ(0, 0) = cos(thetaZ);
    rotAboutZ(0, 1) = -1 * sin(thetaZ);
    rotAboutZ(0, 2) = 0;

    // Define second row
    rotAboutZ(1, 0) = sin(thetaZ);
    rotAboutZ(1, 1) = cos(thetaZ);
    rotAboutZ(1, 2) = 0;

    // Define third row
    rotAboutZ(2, 0) = 0;
    rotAboutZ(2, 1) = 0;
    rotAboutZ(2, 2) = 1;

    return rotAboutZ;
}
