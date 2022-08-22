#ifndef GRIDS_CPP
#define GRIDS_CPP

#include "grids.h"
#include <vector>
#include <iostream>



AbstractGrid2D::AbstractGrid2D(){}
AbstractGrid2D::~AbstractGrid2D(){}

size_t AbstractGrid2D::pairToIndex(Pair p) const {return pairToIndex(p.x, p.y);}


//Michael

RectangleGrid2D::RectangleGrid2D(size_t M, size_t N, double hx, double hy)
{
    this->M = M;
    this->N = N;
    this->hx = hx;
    this->hy = hy;
}

RectangleGrid2D::~RectangleGrid2D(){}
bool RectangleGrid2D::_checkPair(size_t x, size_t y) const{

    if(x > M-1){return false;}
    if(y > N-1){return false;}
    return true;
}

bool RectangleGrid2D::_checkIndex(size_t index) const{
    if(index > M * N - 1){return false;}
    return true;
}

bool RectangleGrid2D::_checkNum(size_t num, bool withAngles) const{
    size_t numBoundary = getNumberOfBoundaryCells(withAngles);
    if(num > numBoundary){return false;}
    return true;
}

//Denis
cellType RectangleGrid2D::getType(size_t x, size_t y) const{
    if(not _checkPair(x, y)) {return OUT_OF_RANGE;}

    /* std::vector < std::vector <cellType> > arr =
    {
        {ANGLE, LEFT, ANGLE},
        {TOP, CENTRAL, BOTTOM},
        {ANGLE, RIGHT, ANGLE}
    };

    return arr[ uint(x > 0) + uint(x == M - 1) ][ uint(y > 0) + uint(y == M - 1) ]; */

    if(y == 0) {
        if(0 < x and x < M - 1) {
            return TOP;
        } else {
            return ANGLE;
        }
    }
    if(x == 0){
        if(0 < y and y < N - 1) {
            return LEFT;
        } else {
            return ANGLE;
        }
    }
    if(y == N-1) {
        if(0 < x and x < M - 1){
            return BOTTOM;
        } else {
            return ANGLE;
        }
    }
    if(x == M-1) {
        if(0 < y and y < N - 1) {
            return RIGHT;
        } else {
            return ANGLE;
        }
    }
    return CENTRAL;
}
cellType RectangleGrid2D::getType(size_t index) const{
    if(not _checkIndex(index)){return OUT_OF_RANGE;}
    Pair pair = indexToPair(index);
    cellType type = getType(pair.x, pair.y);
    return type;
}

bool RectangleGrid2D::isBoundary(size_t index, bool withAngles) const
{
    cellType type = getType(index);
    if(withAngles){
        if(type == CENTRAL) {
            return false;
        } else {
            return true;
        }
    } else {
        if(type == CENTRAL or type == ANGLE) {
            return false;
        } else {
            return true;
        }
    }
}
size_t RectangleGrid2D::getIndexOfBoundaryCell(size_t num, bool withAngles) const
{
    if(not _checkNum(num, withAngles)){return getNumberOfCells();}
    if(withAngles) {
        if(num < M){return num;}
        else if(M-2 < num and num < M+N-1){return (num-M + 2) * M - 1;}
        else if(M+N-3 < num and num < 2*M+N-2){return N+M - 3 - num + N*M;}
        else if(2*M+N-4 < num){return (2 * (N+M)- 4 -num) * M;}
    } else {
        if(num < M-2){return num+1;}
        else if(M-3 < num and num < M+N-4){return (num-M + 4) * M - 1;}
        else if(M+N-5 < num and num < 2*M+N-6){return N+M - 4 - num + N*M - 1;}
        else if(2*M+N-5 < num){return (2 * (N+M) - 8 - num) * M;}
    }
    return getNumberOfCells();
}
//Misha
size_t RectangleGrid2D::getNumberOfBoundaryCells(bool withAngles) const{
    if(withAngles){
        return (M - 1) * 2 + (N - 1) * 2;
    }else{
        return (M-2) * 2 + (N-2) * 2;
    }
}
Pair RectangleGrid2D::indexToPair(size_t index) const
{
    return {index / N, index % N};
}
size_t RectangleGrid2D::pairToIndex(size_t x, size_t y) const
{
    return x * N + y;
}
size_t RectangleGrid2D::getNumberOfCells() const{return M * N;}


#endif
