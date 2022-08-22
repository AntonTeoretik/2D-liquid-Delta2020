#ifndef DRAWER_H
#define DRAWER_H
#include <bitmap_image.hpp>
#include <vector>
//#include "bitmap_image.hpp"
#include <chrono>
#include <random>

class Drawer
{
    // constants
    bool showPressureFlag = false;
    uint M, N, scale, maxLife, numOfParticles;
    double pressureScale;

    std::vector<std::pair<double, double> > particles;

    std::vector< std::vector <std::pair<double, double> > > vecField;
    std::vector< std::vector <double> > pressField;
    std::vector< uint > lifeTime;

    void initializeRandom();

    void loadVecField(std::ifstream &fout);
    void loadPressField(std::ifstream &fout);

    void drawImg(std::string filename);
    void makeStep(double  timeStep);

    //
    //void _setRandomParticle(size_t indexOfParticle);

public:
    Drawer(uint maxLife, uint numOfParticles, bool showPressure = false, double pressureScale = 1.0);

    void startSimulation(std::string filename_pres,
                         std::string filename_vec,
                         std::string outfilename,
                         uint scale = 5,
                         uint freeSimulations = 30,
                         double timeStep = 0.01);


};

#endif // DRAWER_H
