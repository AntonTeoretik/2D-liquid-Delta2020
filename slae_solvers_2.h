#ifndef SLAE_SOLVERS_H
#define SLAE_SOLVERS_H

#include "matrix.h"
#include "logger.h"
#include <algorithm>

class AbstractSolver : public Logger
{
protected:
    //Logging
    const Matrix& A;
    const Vector& b;

public:
    AbstractSolver(const Matrix &A, const Vector &b): A(A), b(b) {}
    virtual ~AbstractSolver();

    // Main part
    virtual void makeOneStep() = 0;
    virtual void solve(Vector &out) = 0;
};

class SymmetricConjugateGradients: public AbstractSolver
{
    const Matrix& At;

    double _eps;
    size_t _logIterNumber; // Logging ones per some iterations
    bool _logInsideStep;

    // Data
    Vector x_s, x_sp1;
    Vector r_s, r_sm1;
    Vector q_s, q_sm1;
    Vector p_s, p_sm1;
    Vector temp;

    void _setVectors(size_t size);
    void _makeFirstIteration();

    void _swap();

public:
    SymmetricConjugateGradients(const Matrix &A, const Matrix &At, const Vector &b, double eps = 1e-6); // Need to know At
    virtual ~SymmetricConjugateGradients();

    virtual void makeOneStep();
    virtual void solve(Vector &out);
};

#endif // SOLVERS_H
