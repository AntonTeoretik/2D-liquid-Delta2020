#ifndef DIFF_OPERATORS_CPP
#define DIFF_OPERATORS_CPP

#include "diffeq_solvers.h"
#include "grids.h"
#include <vector>
#include <iostream>

AbstractDiffOperators2D::AbstractDiffOperators2D(const AbstractGrid2D &grid) : _grid(grid)
{

}

AbstractDiffOperators2D::~AbstractDiffOperators2D() {}


DiffOperators2DForRectangleGrid::DiffOperators2DForRectangleGrid(const RectangleGrid2D &grid) :
    AbstractDiffOperators2D(grid)
{

}

DiffOperators2DForRectangleGrid::~DiffOperators2DForRectangleGrid() {}

dataVec DiffOperators2DForRectangleGrid::dirichlet(size_t index) const {
    dataVec data;
    IndexAndValue toData;
    cellType type = _grid.getType(index);
    if(type != ANGLE){return data;}
    toData = {index, 1.0};
    data.push_back(toData);
    return data;
}

dataVec DiffOperators2DForRectangleGrid::neumann(size_t index) const {
    const RectangleGrid2D & gr = dynamic_cast<const RectangleGrid2D&>(_grid);
    double hx = gr.getHx();
    double hy = gr.getHy();


    dataVec data;
    IndexAndValue toData;
    cellType type = _grid.getType(index);
    if(type == LEFT){
        toData = {_grid.pairToIndex(_grid.indexToPair(index).x, _grid.indexToPair(index).y), -1.0 / hx};
        data.push_back(toData);
        toData = {_grid.pairToIndex(_grid.indexToPair(index).x + 1, _grid.indexToPair(index).y), 1.0/hx};
        data.push_back(toData);
        return data;
    }
    else if(type == RIGHT){
        toData = {_grid.pairToIndex(_grid.indexToPair(index).x, _grid.indexToPair(index).y), -1.0/hx};
        data.push_back(toData);
        toData = {_grid.pairToIndex(_grid.indexToPair(index).x - 1, _grid.indexToPair(index).y), 1.0/hx};
        data.push_back(toData);
        return data;
    }
    else if(type == TOP){
        toData = {_grid.pairToIndex(_grid.indexToPair(index).x, _grid.indexToPair(index).y), -1.0/hy};
        data.push_back(toData);
        toData = {_grid.pairToIndex(_grid.indexToPair(index).x, _grid.indexToPair(index).y + 1), 1.0/hy};
        data.push_back(toData);
        return data;
    }
    else if(type == BOTTOM){
        toData = {_grid.pairToIndex(_grid.indexToPair(index).x, _grid.indexToPair(index).y), -1.0/hy};
        data.push_back(toData);
        toData = {_grid.pairToIndex(_grid.indexToPair(index).x, _grid.indexToPair(index).y - 1), 1.0/hy};
        data.push_back(toData);
        return data;
    }
    return data;
}

dataVec DiffOperators2DForRectangleGrid::laplace(size_t index) const
{
    const RectangleGrid2D & gr = dynamic_cast<const RectangleGrid2D&>(_grid);
    double hx = gr.getHx();
    double hy = gr.getHy();

    dataVec data;
    IndexAndValue toData;
    cellType type = _grid.getType(index);
    if(type != CENTRAL){return data;}
    toData = {_grid.pairToIndex(_grid.indexToPair(index).x, _grid.indexToPair(index).y - 1), 1.0/(hy*hy)};
    data.push_back(toData);
    toData = {_grid.pairToIndex(_grid.indexToPair(index).x, _grid.indexToPair(index).y + 1), 1.0/(hy*hy)};
    data.push_back(toData);
    toData = {_grid.pairToIndex(_grid.indexToPair(index).x - 1, _grid.indexToPair(index).y),   1.0/(hx*hx)};
    data.push_back(toData);
    toData = {_grid.pairToIndex(_grid.indexToPair(index).x + 1, _grid.indexToPair(index).y),   1.0/(hx*hx)};
    data.push_back(toData);
    toData = {index, - 2.0/(hy*hy) - 2.0/(hx*hx)};
    data.push_back(toData);
    return data;
}

dataVec DiffOperators2DForRectangleGrid::ddx(size_t index) const
{
    const RectangleGrid2D & gr = dynamic_cast<const RectangleGrid2D&>(_grid);
    double hx = gr.getHx();

    dataVec data;

    cellType type = _grid.getType(index);
    if(type == LEFT or type == RIGHT or type == ANGLE or type == OUT_OF_RANGE){ return data; }

    IndexAndValue toData;
    Pair pair = _grid.indexToPair(index);

    toData = { _grid.pairToIndex(pair.x + 1, pair.y), 0.5 / hx};
    data.push_back(toData);

    toData = { _grid.pairToIndex(pair.x - 1, pair.y), -0.5 / hx};
    data.push_back(toData);

    return data;
}

dataVec DiffOperators2DForRectangleGrid::ddy(size_t index) const{

    const RectangleGrid2D & gr = dynamic_cast<const RectangleGrid2D&>(_grid);
    double hy = gr.getHy();

    dataVec data;
    cellType type = _grid.getType(index);
    if(type == BOTTOM or type == TOP or type == ANGLE or type == OUT_OF_RANGE){return data;}

    IndexAndValue toData;
    Pair pair = _grid.indexToPair(index);

    toData = { _grid.pairToIndex(pair.x, pair.y + 1), 0.5 / hy };
    data.push_back(toData);

    toData = { _grid.pairToIndex(pair.x, pair.y - 1), -0.5 / hy };
    data.push_back(toData);
    return data;
}

dataVec DiffOperators2DForRectangleGrid::d2dx2(size_t index) const{
    const RectangleGrid2D & gr = dynamic_cast<const RectangleGrid2D&>(_grid);
    double hx = gr.getHx();

    dataVec data;
    cellType type = _grid.getType(index);
    if(type == LEFT or type == RIGHT or type == ANGLE or type == OUT_OF_RANGE){return data;}

    IndexAndValue toData;
    toData = {_grid.pairToIndex(_grid.indexToPair(index).x - 1, _grid.indexToPair(index).y), 1.0 /( hx * hx )};
    data.push_back(toData);
    toData = {index, -2.0 / (hx * hx) };
    data.push_back(toData);
    toData = {_grid.pairToIndex(_grid.indexToPair(index).x + 1, _grid.indexToPair(index).y), 1.0 / (hx * hx )};
    data.push_back(toData);
    return data;
}

dataVec DiffOperators2DForRectangleGrid::d2dy2(size_t index) const
{
    const RectangleGrid2D & gr = dynamic_cast<const RectangleGrid2D&>(_grid);
    double hy = gr.getHy();

    dataVec data;
    IndexAndValue toData;
    cellType type = _grid.getType(index);
    if(type == BOTTOM or type == TOP or type == ANGLE or type == OUT_OF_RANGE){return data;}

    toData = { _grid.pairToIndex(_grid.indexToPair(index).x, _grid.indexToPair(index).y - 1), 1.0/(hy*hy)};
    data.push_back(toData);
    toData = {index, -2.0/ ( hy * hy )};
    data.push_back(toData);
    toData = { _grid.pairToIndex(_grid.indexToPair(index).x, _grid.indexToPair(index).y + 1), 1.0/(hy*hy)};
    data.push_back(toData);
    return data;
}

dataVec DiffOperators2DForRectangleGrid::d2dxdy(size_t index) const
{
    const RectangleGrid2D & gr = dynamic_cast<const RectangleGrid2D&>(_grid);
    double hx = gr.getHx();
    double hy = gr.getHy();

    dataVec data;
    cellType type = _grid.getType(index);

    if(type != CENTRAL){return data;}
    auto p = _grid.indexToPair(index);

    data = {
        {_grid.pairToIndex(p.x - 1, p.y - 1),  0.25/(hx*hy)},
        {_grid.pairToIndex(p.x - 1, p.y + 1), -0.25/(hx*hy)},
        {_grid.pairToIndex(p.x + 1, p.y - 1), -0.25/(hx*hy)},
        {_grid.pairToIndex(p.x + 1, p.y + 1),  0.25/(hx*hy)}
    };

    return data;
}


#endif

