#pragma once

#include <vector>
#include <memory>
#include <math.h>
#include "Matrix.h"
#include "vector.h"

struct MatrixHelper
{
public:
	static Matrix<double> GetRotationMatrix(double thetaX, double thetaY, double thetaZs);
	static Matrix<double> GetTranslationalRotationMatrix(Vector elmVector, double rotationAngle);
	static Matrix<double> GetRotationMatrixAboutX(double thetaX);
	static Matrix<double> GetRotationMatrixAboutY(double thetaY);
	static Matrix<double> GetRotationMatrixAboutZ(double thetaZ);

private:
	

};