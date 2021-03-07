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
    Matrix(const Matrix& that)
        : m_rowCount(that.m_rowCount), m_colCount(that.m_colCount),
        RowCount(that.RowCount), ColCount(that.ColCount)
    {
        this->allocate_memory();
        std::memcpy(this->firstElementAdress, that.firstElementAdress, m_rowCount * m_colCount);
    }

    // Move constructor
    Matrix(Matrix&& that)
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
        free(firstElementAdress);
    }

    // Index usage
    T& operator()(unsigned int i, unsigned j)
    {
        if ((this->m_rowCount <= i) || (this->m_colCount <= j))
            throw std::runtime_error("Matrix Error: Subscript out of range");

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

    // Matrix addition
    Matrix<T> operator+(Matrix<T>& const that)
    {
        // Check if dimensions match
        if ((this->m_colCount != that.m_colCount) || (this->m_rowCount != that.m_rowCount))
            throw std::runtime_error("Matrix Addition Error: Matrix sizes does not match\n");

        Matrix<T> result(*this);
        for (size_t i = 0; i < result.m_rowCount; i++)
            for (size_t j = 0; j < result.m_colCount; j++)
                result(i, j) += that(i, j);
        return result;
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
        if ((rowStart < 0) || (colStart < 0) ||
            (this->m_rowCount <= rowEnd) || (this->m_colCount <= colEnd) ||
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

    unsigned int RowCount;
    unsigned int ColCount;

private:
    unsigned int m_rowCount = 0;
    unsigned int m_colCount = 0;
    T* firstElementAdress = nullptr;

    void allocate_memory()
    {
        firstElementAdress = (T*)malloc(m_rowCount * m_colCount * sizeof(T));
    };

    void fillZeros()
    {
        for (size_t i = 0; i < this->m_rowCount; i++)
            for (size_t j = 0; j < this->m_colCount; j++)
                this->firstElementAdress[(i * m_colCount) + j] = 0;
    }

};
