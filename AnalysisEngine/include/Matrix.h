#pragma once
#include <iostream>
#include <memory>

// Template matrix class for internal use
template <class T>
class Matrix
{
public:
    // Default constructor for matrix class
    Matrix() {}

    // Initialize a square matrix
    Matrix(unsigned int size)
        : m_rowCount(size), m_colCount(size),
        RowCount(size), ColCount(size)
    {
        this->allocate_memory();
        this->fillZeros();
    }

    // Initialize a matrix of size rowCount * colCount
    Matrix(unsigned int rowCount, unsigned int colCount)
        : m_rowCount(rowCount), m_colCount(colCount),
        RowCount(rowCount), ColCount(colCount)
    {
        this->allocate_memory();
        this->fillZeros();
    }

    // Copy constructor
    Matrix(Matrix& const that)
        : m_rowCount(that.m_rowCount), m_colCount(that.m_colCount),
        RowCount(that.RowCount), ColCount(that.ColCount)
    {
        this->allocate_memory();
        std::memcpy(this->firstElementAdress, that.firstElementAdress, m_rowCount * m_colCount * sizeof(T));
    }

    // Move constructor
    Matrix(Matrix&& const that)
        : m_rowCount(that.m_rowCount), m_colCount(that.m_colCount),
        RowCount(that.RowCount), ColCount(that.ColCount),
        firstElementAdress(that.firstElementAdress)
    {
        that.firstElementAdress = nullptr;
    }

    // Move assignment operator
    Matrix& operator=(Matrix that)
    {
        std::swap(m_rowCount, that.m_rowCount);
        std::swap(m_colCount, that.m_colCount);
        std::swap(RowCount, that.RowCount);
        std::swap(ColCount, that.ColCount);
        std::swap(firstElementAdress, that.firstElementAdress);
        return *this;
    }

    // Destructor
    ~Matrix()
    {
        delete[] firstElementAdress;
    }

    // Index usage
    T& operator()(unsigned int i, unsigned j)
    {
        if ((this->m_rowCount <= i) || (this->m_colCount <= j))
            throw std::runtime_error("Matrix Index Error: Subscript out of range");

        return this->firstElementAdress[(i * m_colCount) + j];
    }

    // Matrix multiplication
    Matrix<T> operator*(Matrix<T>& const that)
    {
        // Check if dimensions match
        if (this->m_colCount != that.m_rowCount)
            throw std::runtime_error("Matrix Multiplication Error: Subscript indices does not match\n");

        Matrix<T> result(this->m_rowCount, that.m_colCount);

        for (size_t i = 0; i < this->m_rowCount; i++)
        {
            for (size_t j = 0; j < that.m_colCount; j++)
            {
                auto sum = 0.0;
                for (size_t k = 0; k < this->m_colCount; k++)
                {
                    sum += (*this)(i, k) * that(k, j);
                }
                result(i, j) = sum;
            }
        }

        return result;
    }

    // Matrix multiplication with scalar
    Matrix<T> operator*(const T& mult)
    {
        Matrix<T> result(*this);
        for (size_t i = 0; i < result.m_rowCount; i++)
            for (size_t j = 0; j < result.m_colCount; j++)
                result(i, j) = result(i, j) * mult;
        return result;
    }

    void operator*=(const T& mult)
    {
        for (size_t i = 0; i < this->m_rowCount; i++)
            for (size_t j = 0; j < this->m_colCount; j++)
                (*this)(i, j) *= mult;
    }

    // Matrix subtraction
    Matrix<T> operator-(Matrix<T>& const that)
    {
        // Check if dimensions match
        if ((this->m_colCount != that.m_colCount) || (this->m_rowCount != that.m_rowCount))
            throw std::runtime_error("Matrix Subtraction Error: Matrix sizes does not match\n");

        Matrix<T> result(this->RowCount, this->ColCount);
        for (size_t i = 0; i < result.m_rowCount; i++)
            for (size_t j = 0; j < result.m_colCount; j++)
                result(i, j) = (*this)(i, j) - that(i, j);
        return result;
    }

    // Matrix summation
    Matrix<T> operator+(Matrix<T>& const that)
    {
        // Check if dimensions match
        if ((this->m_colCount != that.m_colCount) || (this->m_rowCount != that.m_rowCount))
            throw std::runtime_error("Matrix Summation Error: Matrix sizes does not match\n");

        Matrix<T> result(this->RowCount, this->ColCount);
        for (size_t i = 0; i < result.m_rowCount; i++)
            for (size_t j = 0; j < result.m_colCount; j++)
                result(i, j) = (*this)(i, j) + that(i, j);
        return result;
    }

    void operator+=(Matrix<T>& const that)
    {
        // Check if dimensions match
        if ((this->m_colCount != that.m_colCount) || (this->m_rowCount != that.m_rowCount))
            throw std::runtime_error("Matrix Summation Error: Matrix sizes does not match\n");

        for (size_t i = 0; i < this->m_rowCount; i++)
            for (size_t j = 0; j < this->m_colCount; j++)
                (*this)(i, j) += that(i, j);
    }

    // Transpose of matrix
    Matrix<T> transpose()
    {
        Matrix<T> retVal(this->ColCount, this->RowCount);
        for (size_t i = 0; i < this->RowCount; i++)
            for (size_t j = 0; j < this->ColCount; j++)
                retVal(j, i) = (*this)(i, j);
        return retVal;
    }

    // Get submatrix within given 0-based indices
    Matrix<T> getSubmatrix(unsigned int rowStart, unsigned int rowEnd, unsigned int colStart, unsigned int colEnd)
    {
        if ((this->m_rowCount <= rowEnd) || (this->m_colCount <= colEnd) ||
            (rowEnd < rowStart) || (colEnd < colStart))
            throw std::runtime_error("Get Submatrix Error: Check desired submatrix indices\n");

        Matrix<T> retVal(rowEnd - rowStart + 1, colEnd - colStart + 1);

        for (size_t i = rowStart; i <= rowEnd; i++)
        {
            for (size_t j = colStart; j <= colEnd; j++)
                retVal(i - rowStart, j - colStart) = (*this)(i, j);
        }

        return retVal;
    }

    void printElements()
    {
        for (size_t i = 0; i < this->m_rowCount; i++)
        {
            for (size_t j = 0; j < this->m_colCount; j++)
            {
                std::cout << this->firstElementAdress[(i * m_colCount) + j] << " ";
            }
            std::cout << "\n";
        }
    }

    void fill(T filler)
    {
        for (size_t i = 0; i < m_rowCount; i++)
            for (size_t j = 0; j < m_colCount; j++)
                (*this)(i, j) = filler;
    }

    Matrix<T> sendToCornerForSquareMatrix(unsigned int rowIdx, unsigned int colIdx, bool isRightBottom)
    {
        // Check if matrix is square or not
        if (m_rowCount != m_colCount)
            throw std::runtime_error("Send to Corner Error: Matrix is not square\n");

        // Check if indices are within boundaries or not
        if ((m_rowCount <= rowIdx) || (m_colCount <= colIdx))
            throw std::runtime_error("Send to Corner Error: Indices are out of bounds\n");

        // Check if element is diagonal or not
        if (rowIdx != colIdx)
            throw std::runtime_error("Send to Corner Error: Element is not on diagonal\n");

        // Get row and column of element
        auto row = this->getSubmatrix(rowIdx, rowIdx, 0, m_colCount - 1);
        auto col = this->getSubmatrix(0, m_rowCount - 1, colIdx, colIdx);

        // Create a new matrix as return value
        Matrix<T> retVal(m_rowCount, m_colCount);

        if (isRightBottom)
        {
            // Start to store elements of original matrix starting from
            // left top excluding row and column of given element. Store
            // them as last row and column.
            for (size_t i = 0; i < m_rowCount; i++)
            {
                for (size_t j = 0; j < m_colCount; j++)
                {
                    if ((i == rowIdx) || (j == colIdx))
                        continue;

                    unsigned storeRowIdx = i < rowIdx ? i : i - 1;
                    unsigned storeColIdx = j < colIdx ? j : j - 1;
                    retVal(storeRowIdx, storeColIdx) = (*this)(i, j);
                }
            }

            // Fill last column
            for (size_t i = 0; i < m_rowCount; i++)
            {
                if (i == rowIdx)
                    continue;
                unsigned storeRowIdx = i < rowIdx ? i : i - 1;
                retVal(storeRowIdx, m_colCount - 1) = col(i, 0);
            }

            // Fill last row
            for (size_t j = 0; j < m_colCount; j++)
            {
                if (j == colIdx)
                    continue;
                unsigned storeColIdx = j < colIdx ? j : j - 1;
                retVal(m_rowCount - 1, storeColIdx) = row(0, j);
            }

            retVal(m_rowCount - 1, m_colCount - 1) = (*this)(rowIdx, colIdx);
        }
        else
        {

            for (size_t i = 0; i < m_rowCount; i++)
            {
                for (size_t j = 0; j < m_colCount; j++)
                {
                    if ((i == rowIdx) || (j == colIdx))
                        continue;

                    unsigned storeRowIdx = i < rowIdx ? i + 1 : i;
                    unsigned storeColIdx = j < colIdx ? j + 1 : j;
                    retVal(storeRowIdx, storeColIdx) = (*this)(i, j);
                }
            }

            // Fill first column
            for (size_t i = 0; i < m_rowCount; i++)
            {
                if (i == rowIdx)
                    continue;
                unsigned storeRowIdx = i < rowIdx ? i + 1 : i;
                retVal(storeRowIdx, 0) = col(i, 0);
            }

            // Fill last row
            for (size_t j = 0; j < m_colCount; j++)
            {
                if (j == colIdx)
                    continue;
                unsigned storeColIdx = j < colIdx ? j + 1 : j;
                retVal(0, storeColIdx) = row(0, j);
            }

            retVal(0, 0) = (*this)(rowIdx, colIdx);

        }

        return retVal;
    }

    Matrix<T> sendItemToBoundVector(unsigned int rowIdx, bool isBot)
    {
        // Check if given matrix is a vector
        if (m_colCount != 1)
            throw std::runtime_error("Vector Operation Error: Given matrix is not vector\n");

        // Check whether given index is within boundaries
        if (m_rowCount <= rowIdx)
            throw std::runtime_error("Vector Operation Error: Index is out of bound\n");

        Matrix<T> retVal(m_rowCount, 1);

        if (isBot)
        {
            for (size_t i = 0; i < m_rowCount; i++)
            {
                if (i == rowIdx)
                    continue;
                auto storeRowIdx = i < rowIdx ? i : i - 1;
                retVal(storeRowIdx, 0) = (*this)(i, 0);
            }
            retVal(m_rowCount - 1, 0) = (*this)(rowIdx, 0);
        }
        else
        {
            for (size_t i = 0; i < m_rowCount; i++)
            {
                if (i == rowIdx)
                    continue;
                auto storeRowIdx = i < rowIdx ? i + 1 : i;
                retVal(storeRowIdx, 0) = (*this)(i, 0);
            }
            retVal(0, 0) = (*this)(rowIdx, 0);
        }


        return retVal;
    }

    unsigned int RowCount;
    unsigned int ColCount;

private:
    unsigned int m_rowCount = 0;
    unsigned int m_colCount = 0;
    T* firstElementAdress = nullptr;

    void allocate_memory()
    {
        //firstElementAdress = (T*)malloc(m_rowCount * m_colCount * sizeof(T));
        firstElementAdress = new T[m_rowCount * m_colCount];
    };

    void fillZeros()
    {
        for (size_t i = 0; i < this->m_rowCount; i++)
            for (size_t j = 0; j < this->m_colCount; j++)
                this->firstElementAdress[(i * m_colCount) + j] = 0;
    }

};
