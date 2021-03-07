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
        RowCount(that.RowCount), ColCount(that.ColCount),
        firstElementAdress(that.firstElementAdress) { }

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
            throw std::runtime_error("Matrix Multiplication Error: Subscript indices does not match!\n");

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

    Matrix<T> transpose()
    {
        Matrix<T> retVal(this->ColCount, this->RowCount);
        for (size_t i = 0; i < this->RowCount; i++)
            for (size_t j = 0; j < this->ColCount; j++)
                retVal(j, i) = (*this)(i, j);
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
