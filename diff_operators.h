#ifndef DIFF_OPERATORS_H
#define DIFF_OPERATORS_H

#include "grids.h"
#include <vector>
#include <iostream>

struct IndexAndValue{
    size_t index;
    double value;
};

typedef std::vector<IndexAndValue> dataVec;

class AbstractDiffOperators2D{
protected:
    const AbstractGrid2D & _grid;

public:
    AbstractDiffOperators2D(const AbstractGrid2D & grid);
    virtual ~AbstractDiffOperators2D();

    virtual dataVec dirichlet(size_t index) const = 0;
    virtual dataVec laplace(size_t index) const = 0;
    virtual dataVec neumann(size_t index) const = 0;

    virtual dataVec ddx(size_t index) const = 0;
    virtual dataVec ddy(size_t index) const = 0;

    virtual dataVec d2dxdy(size_t index) const = 0;
    virtual dataVec d2dx2(size_t index) const = 0;
    virtual dataVec d2dy2(size_t index) const = 0;
};

class DiffOperators2DForRectangleGrid : public AbstractDiffOperators2D{

public:
    DiffOperators2DForRectangleGrid(const RectangleGrid2D& grid);
    virtual ~DiffOperators2DForRectangleGrid();

    virtual dataVec dirichlet(size_t index) const;
    virtual dataVec laplace(size_t index) const;
    virtual dataVec neumann(size_t index) const;

    virtual dataVec ddx(size_t index) const;
    virtual dataVec ddy(size_t index) const;

    virtual dataVec d2dxdy(size_t index) const;
    virtual dataVec d2dx2(size_t index) const;
    virtual dataVec d2dy2(size_t index) const;
};

#endif
// DIFF_OPERATORS_H
