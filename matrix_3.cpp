#include "matrix.h"

inline bool Matrix::_coordsIsCorrect(size_t x, size_t y) const
{
    return x < _size and y < _size;
}

Matrix::Matrix(size_t size) : _size(size)
{
    _nonZeroValues.resize(size);
}

size_t Matrix::size() const
{
    return _size;
}

void Matrix::setValue(size_t x, size_t y, double v, bool checkForCoincidence)
{
    // Check for 0
    if(v == 0.0) {
        for(size_t i = 0; i < _nonZeroValues[x].size(); i++)
            if(_nonZeroValues[x][i].y == y) {
                _nonZeroValues[x].erase(_nonZeroValues[x].begin() + long(i));
                return;
            }
        return;
    }

    // Check for coincidence
    if(checkForCoincidence)
        for(size_t i = 0; i < _nonZeroValues[x].size(); i++)
            if(_nonZeroValues[x][i].y == y) {
                _nonZeroValues[x][i].value = v;
                return;
            }

    // Without checking or if was no coincidence
    _nonZeroValues[x].push_back({y, v});
}

double Matrix::getValue(size_t x, size_t y) const
{
    for(size_t i = 0; i < _nonZeroValues[x].size(); i++)
        if(_nonZeroValues[x][i].y == y)
            return _nonZeroValues[x][i].value;

    return 0.0;
}

unsigned long Matrix::nonZeroElementsSize(size_t i) const
{
    return _nonZeroValues[i].size();
}

MatrixValue Matrix::getMatrValue(size_t i, size_t num) const
{
    return _nonZeroValues[i][num];
}

Matrix Matrix::getTranspose()
{
    Matrix At(_size);
    for (size_t i = 0; i < _size; i++) {
        for(size_t j = 0; j < _nonZeroValues[i].size(); j++)
        {
            auto p = _nonZeroValues[i][j];
            At.setValue(p.y, i, p.value, true);
        }
    }
    return  At;
}

void Matrix::print() const
{
    for(size_t i = 0; i < _size; i++) {
        for(size_t j = 0; j < _size; j++)
            std::cout << getValue(i, j) << " ";
        std::cout << std::endl;
    }
}

/// Vector
inline bool Vector::_coordIsCorrect(size_t x) const
{
    return x < _size;
}

Vector::Vector(size_t size): _size(size)
{
    values = std::vector<double>(_size);
}

size_t Vector::size() const
{
    return _size;
}

void Vector::setValue(size_t x, double v)
{
    values[x] = v;
}

void Vector::addValue(size_t x, double v)
{
    values[x] += v;
}

double Vector::getValue(size_t x) const
{
    return values[x];
}

double Vector::normL2() const
{
    return std::sqrt(normL2sq());
}

double Vector::normL2sq() const
{
    double result = 0.0;
    for (size_t i = 0; i < values.size(); ++i) {
        result += values[i] * values[i];
    }
    return result;
}

void Vector::saveToFile(const std::string filename) const
{
    std::ofstream fout(filename);

    fout << size() << std::endl;

    for (size_t i = 0; i < size(); i++)
        fout << i << " " << getValue(i) << std::endl;

    fout.close();
}

void Vector::readFromFile(const std::string filename)
{
    std::ifstream fin(filename);

    size_t s;
    fin >> s;
    if(s != size())
        return;

    for (size_t i = 0; i < size(); i++) {
        size_t pos;
        double val;
        fin >> pos >> val;
        setValue(pos, val);
    }

    fin.close();
}

void summVectors(const Vector &a, const Vector &b, Vector &out)
{
    //#pragma omp parallel for
    for(size_t i = 0; i < a.size(); i++)
        out.setValue(i, a.getValue(i) + b.getValue(i));
}

void summVectorsAlpha(const Vector &a, const Vector &b, double alpha, Vector &out)
{
    //#pragma omp parallel for
    for(size_t i = 0; i < a.size(); i++)
        out.setValue(i, a.getValue(i) + alpha * b.getValue(i));
}

void multMatrixToVec(const Matrix &A, const Vector &x, Vector &out)
{
    //erase out
    #pragma omp parallel for
    for (size_t i = 0; i < out.size(); i++)
        out.setValue(i, 0.0);

    // WARNING
    #pragma omp parallel for
    for (size_t i = 0; i < A.size(); i++) {
        double summ = 0.0;
        for (size_t j = 0; j < A.nonZeroElementsSize(i); j++) {
            auto p = A.getMatrValue(i, j);
            summ += p.value * x.getValue(p.y);
        }
        out.setValue(i, summ);
    }
}

double dotProduct(const Vector &a, const Vector &b)
{
    double result = 0.0;
    for (size_t i = 0; i < a.size(); i++)
        result += a.getValue(i) * b.getValue(i);

    return result;
}

void copyVector(const Vector &a, Vector &out)
{
    for (size_t i = 0; i < a.size(); i++)
        out.setValue(i, a.getValue(i));
}

void multVecToValue(const Vector &a, double alpha, Vector &out)
{
    #pragma omp parallel for
    for(size_t i = 0; i < a.size(); i++)
        out.setValue(i, alpha * a.getValue(i));
}

double distL2sq(const Vector &a, const Vector &b)
{
    double result = 0.0;
    for (size_t i = 0; i < a.size(); i++)
        result += (a.getValue(i) - b.getValue(i)) * (a.getValue(i) - b.getValue(i));

    return result;
}


double summOfSquares(const Matrix &A, const Vector &x)
{
    double res = 0.0;

    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A.nonZeroElementsSize(i); j++) {
            auto p = A.getMatrValue(i, j);
            double val = p.value * x.getValue(p.y);
            res += val * val;
        }
    }
    return res;
}
