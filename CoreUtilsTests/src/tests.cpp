#include <iostream>
#include <vector>

#include <Matrix.h>
#include <Shape.h>
#include <XYPoint.h>

#define LOG(x) std::cout << " " << x
#define LOGL(x) LOG(x) << "\n"

int main(){

    LOGL("Core Utilities Tests...");
    LOGL("");

    // create vertices
	XYPoint pt1(0.0, 0.0);
    XYPoint pt2(6.0, 0.0);
    XYPoint pt3(6.0, 1.0);
    XYPoint pt4(0.0, 1.0);
    std::vector<XYPoint> vertices;
    vertices.push_back(pt1);
    vertices.push_back(pt2);
    vertices.push_back(pt3);
    vertices.push_back(pt4);

    // create shape
    Shape shape(vertices);

    // check isInside method of shape
    XYPoint ptInside(3.0, 0.5);
    if (shape.isInside(ptInside)) {
        LOGL("isInside method is OK!..");    
    }
    else {
        LOGL("isInside method failed!..");
        return -1;
    }

    LOGL("");

    // create matrices
    Matrix<double> mat1(2, 2);
    Matrix<double> mat2(2, 2);
    auto num = 1.0;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            mat1(i, j) = num;
            mat2(i, j) = 2.0 * num;
            num++;
        }
    }

    LOGL("Matrix 1 is...");
    mat1.printElements();
    LOGL("");

    LOGL("Matrix 2 is...");
    mat2.printElements();
    LOGL("");

    auto&& sumMat = mat1 + mat2;
    LOGL("Summation of matrices is...");
    sumMat.printElements();
    LOGL("");

    auto&& subsMat = mat2 - mat1;
    LOGL("Difference between matrices is...");
    subsMat.printElements();
    LOGL("");

    auto&& multMat = mat1 * mat2;
    LOGL("Multiplication of matrices is...");
    multMat.printElements();
    LOGL("");

    auto multiplier = 2.0;
    auto&& multWithScalar = mat1 * multiplier;
    LOGL("Multiplication of Matrix 1 by " << multiplier << " is...");
    multWithScalar.printElements();
    LOGL("Success!..");

    multWithScalar *= multiplier;
    LOGL("Matrix 1 is multiplied by multiplier once more...");
    multWithScalar.printElements();
    LOGL("");

    return 0;
}















