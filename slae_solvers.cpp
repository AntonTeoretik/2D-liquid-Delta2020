#include <iostream>
#include "slae_solvers.h"
#include "logger.h"

void SymmetricConjugateGradients::_setVectors(size_t size)
{
    x_s = Vector(size);
    x_sp1 = Vector(size);
    r_s = Vector(size);
    r_sm1 = Vector(size);
    q_s = Vector(size);
    q_sm1 = Vector(size);
    p_s = Vector(size);
    p_sm1 = Vector(size);
    temp = Vector(size);
}

void SymmetricConjugateGradients::_swap()
{
    std::swap(x_s, x_sp1);
    std::swap(r_s, r_sm1);
    std::swap(q_s, q_sm1);
    std::swap(p_s, p_sm1);
}

SymmetricConjugateGradients::SymmetricConjugateGradients(const Matrix &A, const Matrix &At, const Vector &b, double eps) :
    AbstractSolver(A, b), At(At), _eps(eps)
{
    //fin
}

SymmetricConjugateGradients::~SymmetricConjugateGradients() {}

AbstractSolver::~AbstractSolver() {}

void SymmetricConjugateGradients::makeOneStep()
{
    double alpha = -1.0 / dotProduct(q_sm1, p_sm1);
    if(float(alpha) == INFINITY or float(-alpha) == INFINITY) alpha = 0.0;

    summVectorsAlpha(r_sm1, q_sm1, alpha, r_s);

    alpha = 1.0 / dotProduct(r_s, r_s);
    if(float(alpha) == INFINITY or float(-alpha) == INFINITY) alpha = 0.0;

    summVectorsAlpha(p_sm1, r_s, alpha, p_s);

    multMatrixToVec(A, p_s, temp);
    multMatrixToVec(At, temp, q_s);

    alpha = -1.0 / dotProduct(q_s, p_s);
    if(float(alpha) == INFINITY or float(-alpha) == INFINITY) alpha = 0.0;

    summVectorsAlpha(x_s, p_s, alpha, x_sp1);
    _swap();
}

void SymmetricConjugateGradients::_makeFirstIteration()
{
    multMatrixToVec(A, x_s, temp);
    summVectorsAlpha(temp, b, -1, r_sm1); /// ?
    multMatrixToVec(At, r_sm1, r_s);

    double alpha = 1.0 / dotProduct(r_s, r_s);
    if(float(alpha) == INFINITY or float(-alpha) == INFINITY) alpha = 0.0;
    summVectorsAlpha(p_sm1, r_s, alpha, p_s);

    multMatrixToVec(A, p_s, temp);
    multMatrixToVec(At, temp, q_s);

    alpha = -1.0 / dotProduct(q_s, p_s);
    if(float(alpha) == INFINITY or float(-alpha) == INFINITY) alpha = 0.0;
    summVectorsAlpha(x_s, p_s, alpha, x_sp1);
    _swap();
}

void SymmetricConjugateGradients::solve(Vector& out)
{
    _setVectors(out.size());
    _makeFirstIteration();

    double r1 = r_sm1.normL2();

    uint counter = 0;
    while((r_s.normL2() / r1) < 1e-6)
    {
        makeOneStep();

        std::cout << "While 1, i = " + std::to_string(counter) + ", " << r_s.normL2() / r1 << std::endl;
        counter++;
    }

    double z = b.normL2sq() + summOfSquares(A, x_s);

    double epsilon_squared = 1e-20; //    double epsilon = pow(10,-16.2);
    double epsilon_squared_recipocal_times_z = 1.0 / (epsilon_squared * z);
    double one_over_r_s_normLsq = 0.0; //1.0 / r_s.normL2sq();


    counter++;
    while(one_over_r_s_normLsq < epsilon_squared_recipocal_times_z)
    {
        makeOneStep();
        one_over_r_s_normLsq += 1.0 / r_sm1.normL2sq();
        if(counter % 1000 == 0) {
        std::cout << "While 2, i = " + std::to_string(counter) + ", " <<
                     one_over_r_s_normLsq / epsilon_squared_recipocal_times_z << std::endl;
        }

        counter++;

    }

    multMatrixToVec(A, x_s, r_s);
    summVectorsAlpha(r_s, b, -1.0, p_s);
    std::cout << " || r_s || = " << p_s.normL2sq() << std::endl;

    copyVector(x_s, out);
}







