#include "MatrixHelper.h"
#include <vector>
#include <math.h>
#include <iostream>

std::vector<std::vector<double>> MatrixHelper::GetRotationMatrix(double thetaX, double thetaY, double thetaZ)
{
    auto rotX = MatrixHelper::GetRotationMatrixAboutX(thetaX);
    auto rotY = MatrixHelper::GetRotationMatrixAboutY(thetaY);
    auto rotZ = MatrixHelper::GetRotationMatrixAboutZ(thetaZ);

    return MatrixHelper::MultiplyMatrices((MatrixHelper::MultiplyMatrices(rotZ, rotY)), rotX);
}

std::vector<std::vector<double>> MatrixHelper::GetTranspose(std::vector<std::vector<double>> matrix)
{
    auto rowCount = (unsigned int)matrix.size();
    auto colCount = (unsigned int)matrix.at(0).size();

    std::vector<std::vector<double>> trans;

    for (unsigned int colIdx = 0;  colIdx < colCount;  colIdx++)
    {
        // Read values at columns of matrix. Assign them to row 
        // vector and push it back to transpose matrix as a row
        std::vector<double> row;

        for (unsigned int rowIdx = 0; rowIdx < rowCount; rowIdx++)
            row.push_back(matrix.at(rowIdx).at(colIdx));

        trans.push_back(row);
    }

    return trans;
}

std::vector<std::vector<double>> MatrixHelper::MultiplyMatrices(std::vector<std::vector<double>> leftMatrix, std::vector<std::vector<double>> rightMatrix)
{
    std::vector<std::vector<double>> result;

    // Check if dimensions match
    auto leftMatrixRowCount = (unsigned int)leftMatrix.size();
    auto leftMatrixColCount = (unsigned int)leftMatrix.at(0).size();

    auto rightMatrixRowCount = (unsigned int)rightMatrix.size();
    auto rightMatrixColCount = (unsigned int)rightMatrix.at(0).size();

    if (leftMatrixColCount != rightMatrixRowCount)
    {
        std::cout << "Matrix sizes do not match!" << std::endl;
        return result;
    }

    for (unsigned int i = 0; i < leftMatrixRowCount; i++)
    {
        // Row from left matrix
        auto row = leftMatrix.at(i);
        
        std::vector<double> resultRow;

        for (unsigned int j = 0; j < rightMatrixColCount; j++)
        {
            // Column from right matrix
            std::vector<double> col;
            for (unsigned int k = 0; k < rightMatrixRowCount; k++)
            {
                col.push_back(rightMatrix.at(k).at(j));
            }

            auto val = 0.0;

            for (unsigned int t = 0; t < (unsigned int)row.size(); t++)
            {
                val += row.at(t) * col.at(t);
            }

            resultRow.push_back(val);
        }
        result.push_back(resultRow);
    }

    return result;
}

std::vector<std::vector<double>> MatrixHelper::GetRotationMatrixAboutX(double thetaX)
{
    // Define rotation matrix
    std::vector<std::vector<double>> rotAboutX;

    // Define first row
    std::vector<double> firstRow;
    firstRow.push_back(1);
    firstRow.push_back(0);
    firstRow.push_back(0);

    // Define second row
    std::vector<double> secondRow;
    secondRow.push_back(0);
    secondRow.push_back(cos(thetaX));
    secondRow.push_back(-1 * sin(thetaX));

    // Define third row
    std::vector<double> thirdRow;
    thirdRow.push_back(0);
    thirdRow.push_back(sin(thetaX));
    thirdRow.push_back(cos(thetaX));

    // Push rows into rotation matrix
    rotAboutX.push_back(firstRow);
    rotAboutX.push_back(secondRow);
    rotAboutX.push_back(thirdRow);

    return rotAboutX;
}

std::vector<std::vector<double>> MatrixHelper::GetRotationMatrixAboutY(double thetaY)
{
    // Define rotation matrix
    std::vector<std::vector<double>> rotAboutY;

    // Define first row
    std::vector<double> firstRow;
    firstRow.push_back(cos(thetaY));
    firstRow.push_back(0);
    firstRow.push_back(sin(thetaY));

    // Define second row
    std::vector<double> secondRow;
    secondRow.push_back(0);
    secondRow.push_back(1);
    secondRow.push_back(0);

    // Define third row
    std::vector<double> thirdRow;
    thirdRow.push_back(-1 * sin(thetaY));
    thirdRow.push_back(0);
    thirdRow.push_back(cos(thetaY));

    // Push rows into rotation matrix
    rotAboutY.push_back(firstRow);
    rotAboutY.push_back(secondRow);
    rotAboutY.push_back(thirdRow);

    return rotAboutY;
}

std::vector<std::vector<double>> MatrixHelper::GetRotationMatrixAboutZ(double thetaZ)
{
    // Define rotation matrix
    std::vector<std::vector<double>> rotAboutZ;

    // Define first row
    std::vector<double> firstRow;
    firstRow.push_back(cos(thetaZ));
    firstRow.push_back(-1 * sin(thetaZ));
    firstRow.push_back(0);

    // Define second row
    std::vector<double> secondRow;
    secondRow.push_back(sin(thetaZ));
    secondRow.push_back(cos(thetaZ));
    secondRow.push_back(0);

    // Define third row
    std::vector<double> thirdRow;
    thirdRow.push_back(0);
    thirdRow.push_back(0);
    thirdRow.push_back(1);

    // Push rows into rotation matrix
    rotAboutZ.push_back(firstRow);
    rotAboutZ.push_back(secondRow);
    rotAboutZ.push_back(thirdRow);

    return rotAboutZ;
}
