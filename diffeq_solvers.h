#ifndef DIFFEQ_SOLVERS_H
#define DIFFEQ_SOLVERS_H

#include "grids.h"
#include "slae_solvers.h"
#include "diff_operators.h"
#include "logger.h"
#include <vector>
#include <iostream>
#include <fstream>

enum CondType
{
    DIRICHLET,
    NEUMANN,
    LAPLACE // ...
};

struct IndexValueType
{
    size_t index;
    double value;
    CondType cond;
};
typedef std::vector<IndexValueType> conditions;



class PoissonSolver2D : public Logger //For any grid
{
protected:
    // const AbstractDiffOperators2D & op;
    const AbstractGrid2D & _grid;
    const conditions & _cond;

public:
    PoissonSolver2D(const AbstractGrid2D & grid, const conditions& cond);
    virtual ~PoissonSolver2D();
    virtual void solve(Vector& out) const = 0;

};

class PoissonSolver2DRectangleGrid : public PoissonSolver2D  //For any rectangle grid
{
    DiffOperators2DForRectangleGrid op;
public:
    PoissonSolver2DRectangleGrid(const RectangleGrid2D& grid, const conditions& cond);
    virtual ~PoissonSolver2DRectangleGrid ();
    virtual void solve(Vector& out) const;
};


class NavierStokesSolver2DRectangleGrid : public Logger // Only for rectangle grid
{
    // Data
    Vector vx; // velocity_x
    Vector vy; // velocity_y

    Vector p;   // pressure

    double rho; // density
    double nu;  // ...

    double tau; // time step

    Vector forces_x;
    Vector forces_y;

    Vector vx_sub; // velocity_x_s+1/2
    Vector vy_sub; // velocity_y_s+1/2

    Vector vx_next; // velocity_x_s+1
    Vector vy_next; // velocity_y_s+1

    conditions pressureConditions;

    void _setVectors(size_t size); // Initialise vectors

    const RectangleGrid2D & _grid;
    DiffOperators2DForRectangleGrid op;

    void _saveVelocityToFile(std::ofstream & fout);
    void _savePressureToFile(std::ofstream & fout);

    double _operate(Vector &vec, dataVec data);

    void _makeOneStep(); //makes one timestep

    void _makeFirstStep();
    void _makeSecondStep();
    void _makeThirdStep();

    void _loadPressureConditions(std::ifstream& fin);
    void _loadForces(std::ifstream& fin);


public:

    NavierStokesSolver2DRectangleGrid(const RectangleGrid2D & grid,
                                      double tau, double rho = 1.0, double nu = 0.0);

    ~NavierStokesSolver2DRectangleGrid();


    void startSimulation(std::string out_filename_velocity,
                         std::string out_filename_pressure,
                         std::string forces,
                         std::string pressureConditions
                         ); //

};

#endif // DIFFEQ_SOLVERS_H
