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

    T& operator()(unsigned int i, unsigned j)
    {
        if ((this->m_rowCount <= i) || (this->m_colCount <= j))
            throw std::runtime_error("Matrix Error: Subscript out of range");

        return this->firstElementAdress[(i * m_colCount) + j];
    }

    T& operator()(std::shared_ptr<Matrix> m, unsigned i, unsigned j)
    {
        return (*m)->firstElementAdress[(i * m_colCount) + j];
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
