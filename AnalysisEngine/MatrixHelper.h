#pragma once

#include <vector>
#include "vector.h"
#include <math.h>

struct MatrixHelper
{
public:
	static std::vector<std::vector<double>> GetTranspose(std::vector<std::vector<double>> matrix);
	static std::vector<std::vector<double>> MultiplyMatrices(std::vector<std::vector<double>> leftMatrix, std::vector<std::vector<double>> rightMatrix);
	static std::vector<std::vector<double>> GetRotationMatrix(double thetaX, double thetaY, double thetaZs);
	static std::vector<std::vector<double>> GetTranslationalRotationMatrix(Vector elmVector, double rotationAngle);
	static std::vector<std::vector<double>> GetRotationMatrixAboutX(double thetaX);
	static std::vector<std::vector<double>> GetRotationMatrixAboutY(double thetaY);
	static std::vector<std::vector<double>> GetRotationMatrixAboutZ(double thetaZ);

private:
	

};