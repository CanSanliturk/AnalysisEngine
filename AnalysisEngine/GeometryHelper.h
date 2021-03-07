#pragma once

#include <vector>
#include <memory>
#include <math.h>
#include "Matrix.h"
#include "vector.h"

namespace GeometryHelper
{
    Matrix<double> GetRotationMatrix(double thetaX, double thetaY, double thetaZs);
    Matrix<double> GetTranslationalRotationMatrix(Vector elmVector, double rotationAngle);
    Matrix<double> GetRotationMatrixAboutX(double thetaX);
    Matrix<double> GetRotationMatrixAboutY(double thetaY);
    Matrix<double> GetRotationMatrixAboutZ(double thetaZ);
};