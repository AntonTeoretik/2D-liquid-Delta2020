#include <matrix.h>
#include <cmath>
#include <iostream>

//Matrix

bool Matrix::_coordsIsCorrect(size_t x, size_t y) const
{
    if(x >= _size || y >= _size )
    {
        return false;
    }
    else
    {
        return true;
    }
}

Matrix Matrix::getTranspose()
{
    Matrix New_Matrix(size());
    for (size_t x = 0; x < size(); x++)
    {
        for (size_t y = 0; y < nonZeroElementsSize(x); y++)
        {
            MatrixValue Element = getMatrValue(x, y);
            //Element - vector that holds the number of column and the coefficient
            New_Matrix.setValue(Element.y, x, Element.value);
        }
    }
    return New_Matrix;
}

void Matrix::setValue(size_t x, size_t y, double v, bool checkForCoincidence)
{
    if (checkForCoincidence)
    {
        for (size_t j = 0; j < nonZeroElementsSize(x); j++)
        {
            //std::cout << nonZeroElementsSize(x) << std::endl;
  //          std::cout << j << std::endl;
            if(_nonZeroValues[x][j].y == y)
            {
                _nonZeroValues[x][j].value = v;
                return;
            }
        }
    }

    MatrixValue New_Cell = {y, v};
    _nonZeroValues[x].push_back(New_Cell);
    //std::cout << "SET VALUE: " << x << " " << y << " " << v << std::endl;
}

double Matrix::getValue(size_t x, size_t y) const
{
    for(size_t i = 0; i < nonZeroElementsSize(x); i++)
    {
       // std::cout << i << std::endl;
        if(_nonZeroValues[x][i].y == y)
        {
 //           std::cout << "test" << std::endl;
            return _nonZeroValues[x][i].value;
        }
    }
    return 0;
}

MatrixValue Matrix::getMatrValue (size_t i, size_t num) const
{
    size_t len = nonZeroElementsSize(i);
    MatrixValue value;
    value.value = 0;

    if(num >= len)
    {
        return(value);
    }
    value = _nonZeroValues[i][num];
    return value;
}

unsigned long Matrix::nonZeroElementsSize(size_t i) const
{
    return _nonZeroValues[i].size();
}

Matrix::Matrix(size_t size)
{
    _nonZeroValues.resize(size, std::vector<MatrixValue>());
    _size = size;
 //   std::cout << "test1" << std::endl;
}

size_t Matrix::size() const
{
    return _size;
}

//Vector

Vector::Vector(size_t size)
{
    values.resize(size, 0);
    _size = size;
}
size_t Vector::size() const
{
    return _size;
}

double Vector::getValue(size_t x) const
{
    return values[x];
}

void Vector::setValue(size_t x, double v)
{
    values[x] = v;
}

void Vector::addValue(size_t x, double v)
{
    values[x] += v;
}

double Vector::normL2sq() const
{
    double sum1 = 0;
    for (size_t i = 0; i < _size; i++)
    {
        sum1 += getValue(i) * getValue(i);
    }
    return sum1;
}

double Vector::normL2() const
{
    return sqrt(normL2sq());
}

inline bool Vector::_coordIsCorrect(size_t x) const
{
    if (x < _size)
        return true;
    return false;
}

//Artem
double summOfSquares(const Matrix& A, const Vector& x)
{
    double res = 0;
    for (size_t i = 0; i < A.size(); i++)
    {
        for (size_t j = 0; j < A.nonZeroElementsSize(i); j++) {
            auto p = A.getMatrValue(i, j);
            double val = p.value * x.getValue(p.y);
            res += val * val;
        }
    }
    return res;
}

double distL2sq(const Vector& a, const Vector& b)
{
    double x1 = a.getValue(1);
    double x2 = b.getValue(1);
    double y1 = a.getValue(0);
    double y2 = b.getValue(0);

    return sqrt(pow((x1-x2),2) + pow((y1-y2),2));
}

//Eugene

void copyVector(const Vector& a, Vector& out)
{
    for(size_t i = 0;i < a.size(); ++i)
    {
        out.setValue(i, a.getValue(i));
    }
}

void summVectors(const Vector& a, const Vector& b, Vector& out)
{
    for(size_t i = 0;i < a.size();++i)
    {
        out.setValue(i, a.getValue(i) + b.getValue(i));
    }
}

void summVectorsAlpha(const Vector& a, const Vector& b, double alpha, Vector& out)
{
    multVecToValue(b, alpha, out);
    summVectors(a, out, out);
}

//Sasha

double dotProduct(const Vector & a, const Vector & b)
{
    double thing = 0;
    size_t sizea = a.size();

    for (size_t i = 0; i < sizea; i++)
    {
 //       double need = a.getValue(i);
        double need = a.getValue(i) * b.getValue(i);
        thing = thing + need;
    }
    return(thing);
}

void multMatrixToVec(const Matrix& A, const Vector& x, Vector& out)
{
    #pragma omp parallel for
    for (size_t i = 0; i < out.size(); i++)
        out.setValue(i, 0.0);

    #pragma omp parallel for
    for (size_t i = 0; i < A.size(); i++) {
        double new_value = 0.0;
        for (size_t j = 0; j < A.nonZeroElementsSize(i); j++) {
            auto p = A.getMatrValue(i, j);
            new_value += p.value * x.getValue(p.y);
        }
        out.setValue(i, new_value);
//        std::cout << new_value << std::endl;
    }
}

void multVecToValue(const Vector& a, double alpha, Vector& out)
{
    size_t sizea = a.size();
    for(size_t i = 0; i < sizea; i++)
    {
        out.setValue(i, a.getValue(i) * alpha);
    }
    return;
}

