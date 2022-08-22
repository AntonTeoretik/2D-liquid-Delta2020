TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS +=  -fopenmp
QMAKE_CXXFLAGS += -O3
QMAKE_CXXFLAGS += -march=native


SOURCES += \
        diff_operators.cpp \
        diffeq_solvers.cpp \
        drawer.cpp \
        grids.cpp \
        logger.cpp \
        main.cpp \
        matrix.cpp \
        scene.cpp \
        slae_solvers.cpp

HEADERS += \
    bitmap_image.hpp \
    diff_operators.h \
    diffeq_solvers.h \
    drawer.h \
    grids.h \
    logger.h \
    matrix.h \
    scene.h \ \
    slae_solvers.h

