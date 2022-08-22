#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <omp.h>
//#include "logger.h"


struct MatrixValue{
    size_t y;
    double value;
};

class Matrix
{
    std::vector< std::vector<MatrixValue> > _nonZeroValues;

    inline bool _coordsIsCorrect(size_t x, size_t y) const;
    size_t _size;

public:
    Matrix(size_t size = 0);
    size_t size() const;

    void setValue(size_t x, size_t y, double v, bool checkForCoincidence=true);
    double getValue(size_t x, size_t y) const;

    unsigned long nonZeroElementsSize(size_t i) const;
    MatrixValue getMatrValue(size_t i, size_t num) const;

    Matrix getTranspose();

    /// logging
    void print() const;
};

class Vector
{
    std::vector<double> values;

    inline bool _coordIsCorrect(size_t x) const;
    size_t _size;

public:
    Vector(size_t size = 0);
    size_t size() const;

    void setValue(size_t x, double v);
    void addValue(size_t x, double v);
    double getValue(size_t x) const;

    double normL2() const;
    double normL2sq() const;

    void saveToFile(const std::string filename) const;
    void readFromFile(const std::string filename);
};

// Operations
void multVecToValue(const Vector& a, double alpha, Vector& out); // alpha*a
void summVectors(const Vector& a, const Vector& b, Vector& out); // a + b
void summVectorsAlpha(const Vector& a, const Vector& b, double alpha, Vector& out); // A + aB

double dotProduct(const Vector & a, const Vector & b);
double distL2sq(const Vector & a, const Vector & b);
void copyVector(const Vector& a, Vector& out);

void multMatrixToVec(const Matrix& A, const Vector& x, Vector& out); // Ax
double summOfSquares(const Matrix& A, const Vector& x);

#endif // MATRIX_H
