#ifndef DIFFEQ_SOLVERS_CPP
#define DIFFEQ_SOLVERS_CPP

#include "diffeq_solvers.h"
#include "grids.h"
#include "slae_solvers.h"
#include "diff_operators.h"
#include "logger.h"
#include "matrix.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

PoissonSolver2DRectangleGrid::PoissonSolver2DRectangleGrid(
        const RectangleGrid2D &grid,
        const conditions &cond) :
    PoissonSolver2D(grid, cond), op(grid) {}

PoissonSolver2DRectangleGrid::~PoissonSolver2DRectangleGrid() {}

PoissonSolver2D::PoissonSolver2D(const AbstractGrid2D &grid, const conditions &cond) : _grid(grid), _cond(cond) {}

PoissonSolver2D::~PoissonSolver2D() {}

void PoissonSolver2DRectangleGrid::solve(Vector &out) const
{

    Matrix A(out.size());
    Vector b(out.size());
    //uint counter = 0;

    for (size_t i = 0; i < _cond.size(); i++) {
        b.setValue(_cond[i].index, _cond[i].value);

        dataVec a;

        if(_cond[i].cond == DIRICHLET) {
            a = op.dirichlet(_cond[i].index);
        } else if (_cond[i].cond == NEUMANN) {
            a = op.neumann(_cond[i].index);
            //counter++;
        } else if (_cond[i].cond == LAPLACE) {
            a = op.laplace(_cond[i].index);
        }

        for (size_t u = 0; u < a.size(); ++u) {
            A.setValue(_cond[i].index, a[u].index, a[u].value, false);
        }
    }
    //std::cout << "**********[ " << counter << ": " << _grid.getNumberOfBoundaryCells(true) << " ]********" << std::endl;
    Matrix At = A.getTranspose();
    SymmetricConjugateGradients solver(A, At, b);
    solver.solve(out);
}

void NavierStokesSolver2DRectangleGrid::_loadForces(std::ifstream& fin){

    forces_x = Vector(_grid.getNumberOfCells());
    forces_y = Vector(_grid.getNumberOfCells());

    size_t NumOfForces;
    size_t a, b;
    double f_x, f_y;

    fin >> NumOfForces;

    for (uint i = 0; i < NumOfForces; i++)
    {
        fin >> a >> b >> f_x >> f_y;
        size_t z = _grid.pairToIndex(a, b);
        forces_x.setValue(z, f_x);
        forces_y.setValue(z, f_y);
    }
}

NavierStokesSolver2DRectangleGrid::NavierStokesSolver2DRectangleGrid(
        const RectangleGrid2D &grid,
        double tau,
        double rho,
        double nu) :
    rho(rho),
    nu(nu),
    tau(tau),
    _grid(grid),
    op(grid)
{
    _setVectors(grid.getNumberOfCells());
}

NavierStokesSolver2DRectangleGrid::~NavierStokesSolver2DRectangleGrid()
{

}

void NavierStokesSolver2DRectangleGrid::_loadPressureConditions(std::ifstream& fin)
{
    pressureConditions.resize(0);

    int NumOfConditions;
    fin >> NumOfConditions;
    for (int i = 0; i < NumOfConditions; i++)
    {
        size_t num;
        double value;
        int type;
        fin >> num >> value >> type;
        IndexValueType b = {_grid.getIndexOfBoundaryCell(num, true), value, CondType(type)};
        pressureConditions.push_back(b);
    }
}

void NavierStokesSolver2DRectangleGrid::startSimulation(std::string out_filename_velocity,
                                                        std::string out_filename_pressure,
                                                        std::string forces,
                                                        std::string pressureConditions)
{
    setImportance(ALL);

    std::ifstream fin_f(forces);
    std::ifstream fin_p(pressureConditions);

    uint numberOfSimulations_in_fin_p;
    uint numberOfSimulations_in_fin_f;

    fin_p >> numberOfSimulations_in_fin_p;
    fin_f >> numberOfSimulations_in_fin_f;


    std::ofstream fout_v(out_filename_velocity);
    std::ofstream fout_p(out_filename_pressure);

    uint numberOfSimulations = std::min(numberOfSimulations_in_fin_f, numberOfSimulations_in_fin_p);

    fout_v << numberOfSimulations << std::endl;
    fout_p << numberOfSimulations << std::endl;

    fout_v << _grid.getM() << " " << _grid.getN() << std::endl;
    fout_p << _grid.getM() << " " << _grid.getN() << std::endl;

    for (size_t i = 0; i < numberOfSimulations; i++)
    {
        log("i = ", double(i));
        _loadPressureConditions(fin_p);
        _loadForces(fin_f);

        _makeOneStep();

        _saveVelocityToFile(fout_v);
        _savePressureToFile(fout_p);
    }

    fin_f.close();
    fin_p.close();
    fout_p.close();
    fout_v.close();
}

void NavierStokesSolver2DRectangleGrid::_makeOneStep()
{
    _makeFirstStep();
    _makeSecondStep();
    _makeThirdStep();

    std::swap(vx, vx_next);
    std::swap(vy, vy_next);
}

double NavierStokesSolver2DRectangleGrid::_operate(Vector &vec, dataVec data)
{
    double answer = 0.0;
    for(auto iv : data)
        answer += vec.getValue(iv.index) * iv.value;
    return answer;
}

void NavierStokesSolver2DRectangleGrid::_setVectors(size_t size)
{
    vx = Vector(size);
    vy = Vector(size);

    vx_sub = Vector(size);
    vy_sub = Vector(size);

    vx_next = Vector(size);
    vy_next = Vector(size);

    p = Vector(size);

    forces_x = Vector(size);
    forces_y = Vector(size);
}


void NavierStokesSolver2DRectangleGrid::_saveVelocityToFile(std::ofstream & fout)
{
    for (size_t i = 0; i < vx.size(); i++) {
        Pair coords = _grid.indexToPair(i);
        fout << coords.x << " " << coords.y << " " << vx.getValue(i) << " " << vy.getValue(i) << std::endl;
    }
}

void NavierStokesSolver2DRectangleGrid::_savePressureToFile(std::ofstream & fout)
{
    for (size_t i = 0; i < p.size(); i++) {
        Pair coords  = _grid.indexToPair(i);
        fout << coords.x << " " << coords.y << " " << p.getValue(i) << std::endl;
    }
}

void NavierStokesSolver2DRectangleGrid::_makeFirstStep()
{
    size_t N = _grid.getNumberOfCells();
    for(size_t i = 0; i < N; i++){
        double dvx_dx     = _operate(vx, op.ddx(i));
        double dvx_dy     = _operate(vx, op.ddy(i));
        double d2vx_dxdy  = _operate(vx, op.d2dxdy(i));
        double d2vx_dy2   = _operate(vx, op.d2dy2(i));

        double dvy_dx     = _operate(vy, op.ddx(i));
        double dvy_dy     = _operate(vy, op.ddy(i));
        double d2vy_dxdy  = _operate(vy, op.d2dxdy(i));
        double d2vy_dx2   = _operate(vy, op.d2dx2(i));
        //std::cout << "nu =  " << (d2vy_dxdy - d2vx_dy2) << std::endl;
        vx_sub.setValue(i, vx.getValue(i) +
                        tau * ( (dvx_dx * vx.getValue(i) + dvx_dy * vy.getValue(i)) - nu * (d2vy_dxdy - d2vx_dy2) + forces_x.getValue(i)));
        vy_sub.setValue(i, vy.getValue(i) +
                        tau * ( (dvy_dx * vx.getValue(i) + dvy_dy * vy.getValue(i)) - nu * (d2vx_dxdy - d2vy_dx2) + forces_y.getValue(i)));
    }
}

void NavierStokesSolver2DRectangleGrid::_makeSecondStep()
{
    for(size_t i = 0; i < _grid.getNumberOfCells(); i++)
    {
        if(_grid.getType(i) == CENTRAL) {
            double dvx_dx = _operate(vx, op.ddx(i));
            double dvy_dy = _operate(vy, op.ddy(i));

            double value = (dvx_dx + dvy_dy) / tau;
            //std::cout << value << std::endl;
            IndexValueType cond = {i, value, LAPLACE};
            pressureConditions.push_back(cond);
        }
    }
    PoissonSolver2DRectangleGrid solver(_grid, pressureConditions);
    solver.solve(p);
}

void NavierStokesSolver2DRectangleGrid::_makeThirdStep()
{
    size_t N = _grid.getNumberOfCells();
    for(size_t i = 0; i < N; i++){
        double dp_dx = _operate(p, op.ddx(i));
        double dp_dy = _operate(p, op.ddy(i));

        vx_next.setValue(i, vx_sub.getValue(i) - rho * tau * dp_dx);
        vy_next.setValue(i, vy_sub.getValue(i) - rho * tau * dp_dy);
    }

    for(size_t i = 0; i < _grid.getNumberOfBoundaryCells(true); i++) {
        vx_next.setValue(_grid.getIndexOfBoundaryCell(i), 0);
        vy_next.setValue(_grid.getIndexOfBoundaryCell(i), 0);
    }
}


#endif


