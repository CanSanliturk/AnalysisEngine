#include "TrussMember.h"
#include "XYZPoint.h"
#include "Vector.h"
#include "MatrixHelper.h"
#include <math.h>

TrussMember::TrussMember(unsigned int ElmIndex, Node* iNode, Node* jNode, Section* section, Material* material)
	: ElementIndex(ElmIndex), TrussSection(section), TrussMaterial(material), Type(ElmType::ElementType::Truss),
	Length(jNode->Coordinate.DistanceTo(iNode->Coordinate))
{
	this->Nodes[0] = iNode;
	this->Nodes[1] = jNode;
	AssembleElementLocalStiffnessMatrix();
	AssembleElementLocalMassMatrix();
	AssembleElementRotationMatrix();
	AssembleElementGlobalStiffnessMatrix();
	AssembleElementGlobalMassMatrix();
	iNode->ConnectedElements.push_back(ElmIndex);
}

TrussMember::TrussMember()
{
	Node iN;
	Node jN;
	Section s;
	Material m;
	this->ElementIndex = 0;
	this->Nodes[0] = &iN;
	this->Nodes[1] = &jN;
	this->TrussSection = &s;
	this->TrussMaterial = &m;
	this->Length = -123.456789;
}

TrussMember::~TrussMember()
{
}

void TrussMember::AssembleElementLocalStiffnessMatrix()
{
	auto A = this->TrussSection->Area;
	auto E = this->TrussMaterial->E;
	auto L = this->Length;
	double kElm[6][6];

	for (unsigned int i = 0; i < 6; i++)
		for (unsigned int j = 0; j < 6; j++)
			kElm[i][j] = 0;

	kElm[0][0] = E * A / L;
	kElm[3][3] = E * A / L;

	for (unsigned int i = 0; i < 6; i++)
		for (unsigned int j = 0; j < 6; j++)
			this->LocalCoordinateStiffnessMatrix[i][j] = kElm[i][j];
}

void TrussMember::AssembleElementLocalMassMatrix()
{
	auto A = this->TrussSection->Area;
	auto rho = this->TrussMaterial->UnitWeight;
	auto rX2 = this->TrussSection->Inertia11 / A;
	auto a = this->Length / 2.0;
	auto a2 = a * a;
	auto mult = rho * A * a / 210.0;
	double mElm[12][12];

	for (unsigned int i = 0; i < 6; i++)
		for (unsigned int j = 0; j < 6; j++)
			mElm[i][j] = 0;

	for (unsigned int i = 0; i < 6; i++)
		for (unsigned int j = 0; j < 6; j++)
			this->LocalCoordinateMassMatrix[i][j] = mElm[i][j];
}

void TrussMember::AssembleElementRotationMatrix()
{
	double pi = 3.141592653589793;

	Vector elmVector(this->Nodes[0]->Coordinate, this->Nodes[1]->Coordinate);

	auto thetaX = acos(elmVector.X / elmVector.Length);
	auto thetaY = (pi / 2) - acos(elmVector.Y / elmVector.Length);
	auto thetaZ = (pi / 2) - acos(elmVector.Z / elmVector.Length);

	auto minorRotMat = MatrixHelper::GetRotationMatrix(thetaX, thetaY, thetaZ);

	for (unsigned int i = 0; i < 6; i++)
		for (unsigned int j = 0; j < 6; j++)
			this->RotationMatrix[i][j] = 0.0;

	for (unsigned int i = 0; i < 3; i++)
		for (unsigned int j = 0; j < 3; j++)
			this->RotationMatrix[i][j] = minorRotMat.at(i).at(j);

	for (unsigned int i = 3; i < 6; i++)
		for (unsigned int j = 3; j < 6; j++)
			this->RotationMatrix[i][j] = minorRotMat.at(i - 3).at(j - 3);
}

void TrussMember::AssembleElementGlobalStiffnessMatrix()
{
	double kElmGlob[6][6];
	for (unsigned int i = 0; i < 6; i++)
		for (unsigned int j = 0; j < 6; j++)
			kElmGlob[i][j] = 0.0;

	Vector elmVec(this->Nodes[0]->Coordinate, this->Nodes[1]->Coordinate);

	auto cX = (elmVec.X) / elmVec.Length;
	auto cY = (elmVec.Y) / elmVec.Length;
	auto cZ = (elmVec.Z) / elmVec.Length;
	auto mult = this->TrussMaterial->E * this->TrussSection->Area / this->Length;

	kElmGlob[0][0] = mult * cX * cX;
	kElmGlob[0][1] = mult * cX * cY;
	kElmGlob[0][2] = mult * cX * cZ;
	kElmGlob[0][3] = mult * cX * cX * -1;
	kElmGlob[0][4] = mult * cX * cY * -1;
	kElmGlob[0][5] = mult * cX * cZ * -1;

	kElmGlob[1][0] = kElmGlob[0][1];
	kElmGlob[1][1] = mult * cY * cY;
	kElmGlob[1][2] = mult * cY * cZ;
	kElmGlob[1][3] = mult * cX * cY * -1;
	kElmGlob[1][4] = mult * cY * cY * -1;
	kElmGlob[1][5] = mult * cY * cZ * -1;

	kElmGlob[2][0] = kElmGlob[0][2];
	kElmGlob[2][1] = kElmGlob[1][2];
	kElmGlob[2][2] = mult * cZ * cZ;
	kElmGlob[2][3] = mult * cZ * cX * -1;
	kElmGlob[2][4] = mult * cZ * cY * -1;
	kElmGlob[2][5] = mult * cZ * cZ * -1;

	kElmGlob[3][0] = kElmGlob[0][3];
	kElmGlob[3][1] = kElmGlob[1][3];
	kElmGlob[3][2] = kElmGlob[2][3];
	kElmGlob[3][3] = mult * cX * cX;
	kElmGlob[3][4] = mult * cX * cY;
	kElmGlob[3][5] = mult * cX * cZ;

	kElmGlob[4][0] = kElmGlob[0][4];
	kElmGlob[4][1] = kElmGlob[1][4];
	kElmGlob[4][2] = kElmGlob[2][4];
	kElmGlob[4][3] = kElmGlob[3][4];
	kElmGlob[4][4] = mult * cY * cY;
	kElmGlob[4][5] = mult * cY * cZ;

	kElmGlob[5][0] = kElmGlob[0][5];
	kElmGlob[5][1] = kElmGlob[1][5];
	kElmGlob[5][2] = kElmGlob[2][5];
	kElmGlob[5][3] = kElmGlob[3][5];
	kElmGlob[5][3] = kElmGlob[4][5];
	kElmGlob[5][5] = mult * cZ * cZ;
	
	for (unsigned int i = 0; i < 6; i++)
		for (unsigned int j = 0; j < 6; j++)
			this->GlobalCoordinateStiffnessMatrix[i][j] = kElmGlob[i][j];
}

void TrussMember::AssembleElementGlobalMassMatrix()
{
	std::vector<std::vector<double>> mElm;
	std::vector<std::vector<double>> rot;

	for (unsigned int i = 0; i < 6; i++)
	{
		std::vector<double> mElmRow;
		std::vector<double> rotRow;

		for (unsigned int j = 0; j < 6; j++)
		{
			mElmRow.push_back(this->LocalCoordinateMassMatrix[i][j]);
			rotRow.push_back(this->RotationMatrix[i][j]);
		}

		mElm.push_back(mElmRow);
		rot.push_back(rotRow);
	}

	auto rotTrans = MatrixHelper::GetTranspose(rot);

	auto firstStep = MatrixHelper::MultiplyMatrices(rotTrans, mElm);
	auto mElmGlob = MatrixHelper::MultiplyMatrices(firstStep, rot);

	for (unsigned int i = 0; i < 6; i++)
		for (unsigned int j = 0; j < 6; j++)
			this->GlobalCoordinateMassMatrix[i][j] = mElmGlob.at(i).at(j);
}
