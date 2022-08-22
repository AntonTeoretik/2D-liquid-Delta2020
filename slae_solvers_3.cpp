#include "slae_solvers.h"

void SymmetricConjugateGradients::_setVectors(size_t size)
{
    x_s = Vector(size);
    x_sp1 = Vector(size);

    p_s = Vector(size);
    p_sm1 = Vector(size);

    q_s = Vector(size);
    q_sm1 = Vector(size);

    r_s = Vector(size);
    r_sm1 = Vector(size);

    temp = Vector(size); //Ap_s
}

void SymmetricConjugateGradients::_makeFirstIteration()
{
    // x_1 = 0 (here is a posibility to change x_0)

    _calcR1();
    _calcP1();
    _calcQS();
    _calcXSp1();
    _swap();
}

void SymmetricConjugateGradients::_calcR1()
{
    // r1 = A^t (Ax_1 - b)
    multMatrixToVec(A, x_s, temp); // temp = A * x_1
    summVectorsAlpha(temp, b, -1.0, r_sm1); // r_sm1 = temp - b = A*x_1 - b. r_sm1 uses here as temporary result
    multMatrixToVec(At, r_sm1, r_s); // r_s = A^t * r_sm1 =  A^t (A*x_1 - b)
}

void SymmetricConjugateGradients::_calcP1()
{
    // p1 = r_1 / || r_1 ||^2
    double alpha = 1.0 / r_s.normL2sq();
    if(float(alpha) == INFINITY or float(-alpha) == INFINITY) alpha = 0.0;

    multVecToValue(r_s, alpha, p_s);
}

void SymmetricConjugateGradients::_calcRS()
{
    double alpha = -1.0 / dotProduct(q_sm1, p_sm1);
    //td::cout << alpha << std::endl;
    if(float(alpha) == INFINITY or float(-alpha) == INFINITY) alpha = 0.0;


    summVectorsAlpha(r_sm1, q_sm1, alpha, r_s);
}

void SymmetricConjugateGradients::_calcPS()
{
    double alpha = 1.0 / r_s.normL2sq();
    if(float(alpha) == INFINITY or float(-alpha) == INFINITY) alpha = 0.0;


    summVectorsAlpha(p_sm1, r_s, alpha, p_s);
}

void SymmetricConjugateGradients::_calcQS()
{
    // q_1 = A^t (A p_1)

    multMatrixToVec(A, p_s, temp); // temp = A * p_1
    multMatrixToVec(At, temp, q_s); // q_s = A^t * temp

}

void SymmetricConjugateGradients::_calcXSp1()
{

    // x_2 = x_1 - p_s / <q_s, p_s>
    double alpha = - 1.0 / dotProduct(q_s, p_s);
    if(float(alpha) == INFINITY or float(-alpha) == INFINITY) alpha = 0.0;

    summVectorsAlpha(x_s, p_s, alpha, x_sp1);
}

void SymmetricConjugateGradients::_swap()
{
    std::swap(x_sp1, x_s);
    std::swap(r_s, r_sm1);
    std::swap(p_s, p_sm1);
    std::swap(q_s, q_sm1);
}

SymmetricConjugateGradients::SymmetricConjugateGradients(const Matrix &A, const Matrix &At, const Vector &b, double eps) : AbstractSolver(A, b), At(At)
{
    _eps = eps;
}

SymmetricConjugateGradients::~SymmetricConjugateGradients()
{

}

void SymmetricConjugateGradients::makeOneStep()
{
    _calcRS();
    _calcPS();
    _calcQS();
    _calcXSp1();
    _swap();
}

void SymmetricConjugateGradients::solve(Vector &out)
{
    if(A.size() != b.size())
        log("Wrong sizes... ");

    log("Preparing data... ");

    //preparing data
    size_t size = A.size();
    _setVectors(size);

    log("Makes first iteration... ");

    _makeFirstIteration(); //s = 1

    _logIterNumber = 1000;

    // calculate r1
    double r1normSqEps = r_sm1.normL2sq() * _eps * _eps;
    size_t s = 1;


    // First
    while (r_sm1.normL2sq() > r1normSqEps) {

        std::cout << "  Iteration number: " + std::to_string(s) + ": " << r_sm1.normL2sq() << " " << r1normSqEps << std::endl;
        makeOneStep();
        s++;
    }
    log("End of first part. Starting second...");


    double z = summOfSquares(A, x_s) + b.normL2sq();
    double errorEps = 1e-13; // Magic number
    double factor = 1.0 / (errorEps * errorEps * z);
    double summOfSqrRev = 0.0;

    // Second
    while (summOfSqrRev < factor) {

        makeOneStep(); /// ! Make verbosity here
        s++;

        if(s % _logIterNumber == 0)
            std::cout << "  Iteration number: " + std::to_string(s) + ": " << summOfSqrRev << " " << factor << " " << r_sm1.normL2sq() << std::endl;


        // update sqrRew
        summOfSqrRev += 1.0 / r_sm1.normL2sq();
    }

    log("Copying to out... ");
    //copyVector(x_s, out);
    copyVector(x_sp1, out);
}

AbstractSolver::~AbstractSolver()
{

}
