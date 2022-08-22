///#ifndef DRAWER_H
#include <vector>
#include <bitmap_image.hpp>
#include <bits/stdc++.h>
#include <chrono>
#include <random>
#include <bitset>
#include "drawer.h"
#include "matrix.h"

using namespace std;

void Drawer::makeStep(double  timeStep)
{
    size_t index_y;
    size_t index_x;
    double new_x, new_y;
    //size_t lifetime;

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> distM(scale, (M - 1) * scale);
    std::uniform_real_distribution<double> distN(scale, (N - 1) * scale);
    std::uniform_int_distribution<uint> distT(1, maxLife);

    for(size_t  i = 0; i < particles.size(); i++)
    {
        double x = particles[i].first;
        double y = particles[i].second;

        index_x = size_t(x / scale);
        index_y = size_t(y / scale);
        pair<double, double> a = vecField[index_x][index_y];
        new_x = x + a.first * timeStep;
        new_y = y + a.second * timeStep;

        particles[i].first = new_x;
        particles[i].second = new_y;
        lifeTime[i]--;

        if(lifeTime[i] == 0 or (N - 1)*scale < new_x or (M - 1)*scale < new_y or scale > new_x or scale > new_y)
        {
            particles[i].first = distM(rng);
            particles[i].second = distN(rng);
            lifeTime[i] = distT(rng);
        }
    }
}

Drawer::Drawer(uint maxLife, uint numOfParticles, bool showPressure, double p_scale) :
    showPressureFlag(showPressure), maxLife(maxLife), numOfParticles(numOfParticles), pressureScale(p_scale)
{

}

void Drawer::startSimulation(string filename_pres,
                             string filename_vec,
                             string outfilename,
                             uint _scale,
                             uint freeSimulations,
                             double timeStep)
{
    uint numberOfIteration;
    scale = _scale;

    std::ifstream fin(filename_vec);
    std::ifstream fin_1(filename_pres);

    fin >> numberOfIteration;
    fin_1 >> numberOfIteration;

    fin >> N >> M;
    fin_1 >> N >> M; /// КОСТЫЛЬ

    initializeRandom();
    uint cnt_of_pics = 0;

    for(uint i = 0; i < numberOfIteration; i++)
    {
        std::cout << i << std::endl;
        loadVecField(fin);
        std::cout << i << std::endl;
        loadPressField(fin_1);
        std::cout << i << std::endl;

        for(size_t j = 0; j < freeSimulations; j++)
        {
            makeStep(timeStep);
            string curr_outfilename = outfilename + std::to_string(cnt_of_pics) + ".bmp";
            ++cnt_of_pics;
            drawImg(curr_outfilename);
        }

    }
}

void Drawer::loadVecField(std::ifstream &fin)
{
    vecField.resize(N);

    for (size_t i = 0; i < N; ++i)
        vecField[i].resize(M);

    size_t x, y;
    double xa, ya;
    for(size_t i = 0; i < N * M; ++i) {
        fin >> x >> y;
        fin >> xa >> ya;
        vecField[x][y] = {xa, ya};

   }
}



void Drawer::loadPressField(std::ifstream &fin)
{
    pressField.resize(N);
    for (size_t i = 0; i < N; ++i) {
        pressField[i].resize(M);
    }
    size_t x, y;
    double xa;
    for(size_t i = 0; i < M * N; ++i) {
        fin >> x >> y;
        fin >> xa;
        pressField[x][y] = xa;
    }
}

void Drawer::drawImg(std::string filename)
{
    uint sizeX = N * scale;
    uint sizeY = M * scale; /// M
    double pixelX = 0.0;
    double pixelY = 0.0;
    double presValue = 0.0;

    bitmap_image B(sizeX, sizeY);

    for(size_t i = 0; i < particles.size(); ++i) {
        pixelX = particles[i].first;
        pixelY = particles[i].second;

        uint cnt_x = uint(pixelX / scale),
             cnt_y = uint(pixelY / scale);

        presValue = pressureScale * pressField[cnt_x][cnt_y];

        //u_char presColor = u_char(2.55 * presValue);
        u_char presColor = u_char(255.0 * (1.0 - (1.0 / (1.0 + std::fabs(presValue) ) ) ) );


        if(showPressureFlag)
        {
            if(presValue > 0)
                B.set_pixel(uint(pixelX), uint(pixelY), 255, 255 - presColor, 255 - presColor);
            else
                B.set_pixel(uint(pixelX), uint(pixelY), 255 - presColor, 255 - presColor, 255);
        }
        else
        {
            B.set_pixel(uint(pixelX), uint(pixelY), 255, 255, 255);
        }
    }
    B.save_image(filename);
}

void Drawer::initializeRandom()
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> distM(scale, (M - 1) * scale);
    std::uniform_real_distribution<double> distN(scale, (N - 1) * scale);
    std::uniform_int_distribution<uint> distT(1, maxLife);

    particles.resize(numOfParticles);
    lifeTime.resize(numOfParticles);

    for(size_t i = 0; i < particles.size(); i++){
        particles[i].first = distN(rng);
        particles[i].second = distM(rng);
        lifeTime[i] = distT(rng);
    }
}
