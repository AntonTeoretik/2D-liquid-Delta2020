#ifndef GRIDS_H
#define GRIDS_H

#include <vector>
#include <iostream>

enum cellType{
    CENTRAL,

    // For rectangle
    LEFT,
    RIGHT,
    BOTTOM,
    TOP,
    ANGLE,

    OUT_OF_RANGE

    // for Triangle
};

struct Pair{
    size_t x, y;
};

class AbstractGrid2D{
private:
    virtual bool _checkPair(size_t x, size_t y) const = 0;
    virtual bool _checkIndex(size_t index) const = 0;
    virtual bool _checkNum(size_t num, bool withAngles) const = 0;

public:
    AbstractGrid2D();
    virtual ~AbstractGrid2D();

    virtual Pair indexToPair(size_t index) const = 0;
    virtual size_t pairToIndex(size_t x, size_t y) const = 0;
    virtual size_t pairToIndex(Pair p) const;

    virtual cellType getType(size_t index) const = 0;
    virtual cellType getType(Pair p) const {return getType(p.x, p.y);}
    virtual cellType getType(size_t x, size_t y) const = 0;

    virtual size_t getNumberOfCells() const = 0;

    //boundary cells info
    virtual size_t getNumberOfBoundaryCells(bool withAngles=false) const = 0;
    virtual bool isBoundary(size_t index, bool withAngles=false) const = 0;
    virtual size_t getIndexOfBoundaryCell(size_t num, bool withAngles=false) const = 0;
};

class RectangleGrid2D : public AbstractGrid2D{
    size_t M, N;
    double hx, hy;

private:
    virtual bool _checkPair(size_t x, size_t y) const;
    virtual bool _checkIndex(size_t index) const;
    virtual bool _checkNum(size_t num, bool withAngles) const;

public:
    RectangleGrid2D(size_t M, size_t N, double hx, double hy);
    virtual ~RectangleGrid2D();

    virtual Pair indexToPair(size_t index) const; // check Pair
    virtual size_t pairToIndex(size_t x, size_t y) const; //check Pair

    virtual cellType getType(size_t index) const; // checkpair
    virtual cellType getType(size_t x, size_t y) const;

    virtual size_t getNumberOfCells() const;

    //boundary cells info
    virtual size_t getNumberOfBoundaryCells(bool withAngles=false) const;
    virtual bool isBoundary(size_t index, bool withAngles=false) const;
    virtual size_t getIndexOfBoundaryCell(size_t num, bool withAngles=false) const;

    size_t getM() const {return M;}
    size_t getN() const {return N;}
    double getHx() const {return hx;}
    double getHy() const {return hy;}
};


#endif // GRIDS_H
