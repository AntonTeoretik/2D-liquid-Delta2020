//#include "grids.h"

#include "diffeq_solvers.h"
#include "drawer.h"

void firstScene(uint N, double rad, uint numOfIterations, std::string forcesFile, std::string pressureFile)
{
    RectangleGrid2D grid(N, N, 1.0 / (N - 1), 1.0 / (N - 1));
    std::ofstream fout_f(forcesFile);
    std::ofstream fout_p(pressureFile);

    fout_f << numOfIterations << std::endl;
    fout_p << numOfIterations << std::endl;

    for(uint sim = 0; sim < numOfIterations; sim++) {

        //forces
        fout_f << N * N << std::endl;

        for (uint i = 0; i < N; i++) {
            for (uint j = 0; j < N; j++) {
                double fx = 1e-1, fy = 2*1e-1;

                double r = (i - N / 2.0)*(i - N / 2.0) + (j - N / 2.0)*(j - N / 2.0);
                double val = std::max(0.0, (rad*rad - r) / (rad*rad));
                fx *= val;
                fy *= val;

                fout_f << i << " " << j << " " << fx << " " << fy << std::endl;

            }
        }

        //pressure
        fout_p << grid.getNumberOfBoundaryCells(true) << std::endl;
        for (uint k = 0; k < grid.getNumberOfBoundaryCells(true); k++) {
            if(k == N/2 or grid.getType( grid.getIndexOfBoundaryCell(k, true) ) == ANGLE) {
                fout_p << k << " " << 0.0 << " " << DIRICHLET << std::endl;
                //std::cout << "D" << std::endl;
            }
            else
                fout_p << k << " " << 0.0 << " " << NEUMANN << std::endl;
        }
        //std::cout << std::endl;
    }
    fout_f.close();
    fout_p.close();

}

void secondScene(uint N, double rad, uint numOfIterations, std::string forcesFile, std::string pressureFile)
{
    RectangleGrid2D grid(N, N, 1.0 / (N - 1), 1.0 / (N - 1));
    std::ofstream fout_f(forcesFile);
    std::ofstream fout_p(pressureFile);

    fout_f << numOfIterations << std::endl;
    fout_p << numOfIterations << std::endl;

    for(uint sim = 0; sim < numOfIterations; sim++) {

        //forces
        fout_f << N * N << std::endl;

        for (uint i = 0; i < N; i++) {
            for (uint j = 0; j < N; j++) {
                double fx = 1e-1, fy = 2*1e-1;

                double r = (i - N / 2.0)*(i - N / 2.0) + (j - N / 2.0)*(j - N / 2.0);
                double val = std::max(0.0, (rad*rad - r) / (rad*rad)) * std::sin(sim * 3.1415 / 4);
                fx *= val;
                fy *= val;

                fout_f << i << " " << j << " " << fx << " " << fy << std::endl;

            }
        }

        //pressure
        fout_p << grid.getNumberOfBoundaryCells(true) << std::endl;
        for (uint k = 0; k < grid.getNumberOfBoundaryCells(true); k++) {
            if(k == N/2 or grid.getType( grid.getIndexOfBoundaryCell(k, true) ) == ANGLE) {
                fout_p << k << " " << 0.0 << " " << DIRICHLET << std::endl;
                //std::cout << "D" << std::endl;
            }
            else
                fout_p << k << " " << 0.0 << " " << NEUMANN << std::endl;
        }
        //std::cout << std::endl;
    }
    fout_f.close();
    fout_p.close();
}

void thirdScene(uint N, double rad, uint numOfIterations, std::string forcesFile, std::string pressureFile)
{
    RectangleGrid2D grid(N, N, 1.0 / (N - 1), 1.0 / (N - 1));
    std::ofstream fout_f(forcesFile);
    std::ofstream fout_p(pressureFile);

    fout_f << numOfIterations << std::endl;
    fout_p << numOfIterations << std::endl;

    for(uint sim = 0; sim < numOfIterations; sim++) {
        std::cout << "Make files : " << sim << " of " << numOfIterations << std::endl;
        //forces
        fout_f << N * N << std::endl;

        for (uint i = 0; i < N; i++) {
            for (uint j = 0; j < N; j++) {
                double fx, fy;

                double center1x = sim - 10.0;
                double center2y = sim - 40.0;

                double r2 = (i - N / 2.0)*(i - N / 2.0) + (j - center1x)*(j - center1x);
                double val = std::max(-1.0, (rad*rad - r2) / (rad * rad));

                if (val > 0) {
                    fy = 1e-1 * val;
                    fx = 0 * val;
                } else {
                    r2 = (i - center2y)*(i - center2y) + (j - N / 2.0)*(j - N / 2.0);
                    val = std::max(-1.0, (rad*rad - r2) / (rad * rad));
                    if(val > 0) {
                        fy = 0 * val;
                        fx = 1e-1 * val;
                    } else {
                        fy = 0.0;
                        fx = 0.0;
                    }
               }
               fout_f << i << " " << j << " " << fx << " " << fy << std::endl;
            }
        }

        //pressure
        fout_p << grid.getNumberOfBoundaryCells(true) << std::endl;
        for (uint k = 0; k < grid.getNumberOfBoundaryCells(true); k++) {
            if(k == N/2 or grid.getType( grid.getIndexOfBoundaryCell(k, true) ) == ANGLE) {
                fout_p << k << " " << 0.0 << " " << DIRICHLET << std::endl;
                //std::cout << "D" << std::endl;
            }
            else
                fout_p << k << " " << 0.0 << " " << NEUMANN << std::endl;
        }
        //std::cout << std::endl;
    }
    fout_f.close();
    fout_p.close();
}

void x4Scene(uint N, double rad, uint numOfIterations, std::string forcesFile, std::string pressureFile)
{
    RectangleGrid2D grid(N, N, 1.0 / (N - 1), 1.0 / (N - 1));
    std::ofstream fout_f(forcesFile);
    std::ofstream fout_p(pressureFile);

    fout_f << numOfIterations << std::endl;
    fout_p << numOfIterations << std::endl;

    for(uint sim = 0; sim < numOfIterations; sim++) {
        std::cout << "Make files : " << sim << " of " << numOfIterations << std::endl;
        //forces
        fout_f << N * N << std::endl;

        for (uint i = 0; i < N; i++) {
            for (uint j = 0; j < N; j++) {
                double fx, fy;

                double center1x = sim - 0.0;
                double center2y = sim - 20.0;

                double r2 = (i - N / 2.0)*(i - N / 2.0) + (j - center1x)*(j - center1x);
                double val = std::max(-1.0, (rad*rad - r2) / (rad * rad));

                if (val > 0) {
                    fy = 10 * val;
                    fx = 0 * val;
                } else {
                    r2 = (i - center2y)*(i - center2y) + (j - N / 2.0)*(j - N / 2.0);
                    val = std::max(-1.0, (rad*rad - r2) / (rad * rad));
                    if(val > 0) {
                        fy = 0 * val;
                        fx = 10 * val;
                    } else {
                        fy = 0.0;
                        fx = 0.0;
                    }
               }
               fout_f << i << " " << j << " " << fx << " " << fy << std::endl;
            }
        }

        //pressure
        fout_p << grid.getNumberOfBoundaryCells(true) << std::endl;
        for (uint k = 0; k < grid.getNumberOfBoundaryCells(true); k++) {
            if(k == N/2 or grid.getType( grid.getIndexOfBoundaryCell(k, true) ) == ANGLE) {
                fout_p << k << " " << 0.0 << " " << DIRICHLET << std::endl;
                //std::cout << "D" << std::endl;
            } else if( grid.getType( grid.getIndexOfBoundaryCell(k, true) ) == TOP and
                       grid.indexToPair( grid.getIndexOfBoundaryCell(k, true) ).x < N/2 ) {
                fout_p << k << " " << 1 << " " << NEUMANN << std::endl;
            } else if( grid.getType( grid.getIndexOfBoundaryCell(k, true) ) == BOTTOM and
                       grid.indexToPair( grid.getIndexOfBoundaryCell(k, true) ).x < N/2) {
                fout_p << k << " " << -1 << " " << NEUMANN << std::endl;
            } else {
                fout_p << k << " " << 0.0 << " " << NEUMANN << std::endl;
            }
        }
        //std::cout << std::endl;
    }
    fout_f.close();
    fout_p.close();
}

void x5Scene(uint N, uint numOfIterations, std::string forcesFile, std::string pressureFile)
{
    RectangleGrid2D grid(N, N, 1.0 / (N - 1), 1.0 / (N - 1));
    std::ofstream fout_f(forcesFile);
    std::ofstream fout_p(pressureFile);

    fout_f << numOfIterations << std::endl;
    fout_p << numOfIterations << std::endl;

    for(uint sim = 0; sim < numOfIterations; sim++) {
        std::cout << "Make files : " << sim << " of " << numOfIterations << std::endl;
        //forces
        fout_f << N * N << std::endl;

        for (uint i = 0; i < N; i++) {
            for (uint j = 0; j < N; j++) {
                double fx = 0.0, fy = 0.0;
                double rad1 = N / 6.0;
                double rad2 = N / 3.0;
                double cx = N / 2.0, cy = N / 2.0;

                double dist = sqrt((i - cx)*(i - cx) + (j - cy)*(j - cy));
                if(dist < rad2 and dist > rad1) {
                    double val = (rad2 - dist)*(dist - rad1);
                    fx = val * (j - cy) / dist / (1.0 + sim / 10);
                    fy = -val * (i - cx) / dist / (1.0 + sim / 10);
                }
                fout_f << i << " " << j << " " << fx << " " << fy << std::endl;
            }
        }
        //pressure
        //pressure
        fout_p << grid.getNumberOfBoundaryCells(true) << std::endl;
        for (uint k = 0; k < grid.getNumberOfBoundaryCells(true); k++) {
            if(k == N/2 or grid.getType( grid.getIndexOfBoundaryCell(k, true) ) == ANGLE) {
                fout_p << k << " " << 0.0 << " " << DIRICHLET << std::endl;
                //std::cout << "D" << std::endl;
            }
            else
                fout_p << k << " " << 0.0 << " " << NEUMANN << std::endl;
        }
        //std::cout << std::endl;
    }
    fout_f.close();
    fout_p.close();
}

void x6Scene(uint N, uint numOfIterations, std::string forcesFile, std::string pressureFile)
{
    RectangleGrid2D grid(N, N, 1.0 / (N - 1), 1.0 / (N - 1));
    std::ofstream fout_f(forcesFile);
    std::ofstream fout_p(pressureFile);

    fout_f << numOfIterations << std::endl;
    fout_p << numOfIterations << std::endl;

    for(uint sim = 0; sim < numOfIterations; sim++) {
        std::cout << "Make files : " << sim << " of " << numOfIterations << std::endl;
        //forces
        fout_f << N * N << std::endl;

        for (uint i = 0; i < N; i++) {
            for (uint j = 0; j < N; j++) {
                double fx = 0.0, fy = 0.0;

                if(i < N/6.0) {
                    fy = 0.001 * i * (N / 6.0 - i) * j * (N - j) / (1.0 + sim / 10);
                } else if( i > 5 * N / 6.0) {
                    fy = 0.001 * (N - i) * (5 * N / 6.0 - i) * j * (N - j) / (1.0 + sim / 10);
                }
                fout_f << i << " " << j << " " << fx << " " << fy << std::endl;
            }
        }
        //pressure
        fout_p << grid.getNumberOfBoundaryCells(true) << std::endl;
        for (uint k = 0; k < grid.getNumberOfBoundaryCells(true); k++) {
            if(k == N/2 or grid.getType( grid.getIndexOfBoundaryCell(k, true) ) == ANGLE) {
                fout_p << k << " " << 0.0 << " " << DIRICHLET << std::endl;
            }
            else
                fout_p << k << " " << 0.0 << " " << NEUMANN << std::endl;
        }
        //std::cout << std::endl;
    }
    fout_f.close();
    fout_p.close();
}

void x7Scene(uint N, uint numOfIterations, std::string forcesFile, std::string pressureFile)
{
    RectangleGrid2D grid(N, N, 1.0 / (N - 1), 1.0 / (N - 1));
    std::ofstream fout_f(forcesFile);
    std::ofstream fout_p(pressureFile);

    fout_f << numOfIterations << std::endl;
    fout_p << numOfIterations << std::endl;

    for(uint sim = 0; sim < numOfIterations; sim++) {
        std::cout << "Make files : " << sim << " of " << numOfIterations << std::endl;
        //forces
        fout_f << N * N << std::endl;

        double rad = N / 7.0;

        for (uint i = 0; i < N; i++) {
            for (uint j = 0; j < N; j++) {
                double fx = 0.0, fy = 0.0;

                double r = sqrt ((i - N / 2.0)*(i - N / 2.0) + (j - N / 2.0)*(j - N / 2.0));
                double val = std::max(0.0, (rad*rad - r*r) / (rad*rad));
                fx = 9 * val * sin( 2 * 3.1415926 * (sim / 40.0) );
                fy = 9 * val * cos( 2 * 3.1415926 * (sim / 40.0) );

                fout_f << i << " " << j << " " << fx << " " << fy << std::endl;

            }
        }
        //pressure
        fout_p << grid.getNumberOfBoundaryCells(true) << std::endl;
        for (uint k = 0; k < grid.getNumberOfBoundaryCells(true); k++) {
            if(k == N/2 or grid.getType( grid.getIndexOfBoundaryCell(k, true) ) == ANGLE) {
                fout_p << k << " " << 0.0 << " " << DIRICHLET << std::endl;
            } else {
                fout_p << k << " " << 1 * sin( 6 * 3.1415926 * double(k) / grid.getNumberOfBoundaryCells(true) ) << " " << NEUMANN << std::endl;
            }
        }
        //std::cout << std::endl;
    }
    fout_f.close();
    fout_p.close();
}


int main()
{
    uint N = 70;
    std::string forcesFile = "test/forces.txt";
    std::string prCondFile = "test/pressure.txt";

    x7Scene(N, 1200, forcesFile, prCondFile);

    RectangleGrid2D grid(N, N, 1.0 / (N - 1), 1.0 / (N - 1));
    NavierStokesSolver2DRectangleGrid solver(grid, 1e-10, 1.0, 1e+6);
    solver.startSimulation("vel.txt", "press.txt", "test/forces.txt", "test/pressure.txt");

    std::cout << "Finished" << std::endl;

    Drawer D(100, 100000, true, 5*1e-1);
    D.startSimulation("press.txt", "vel.txt", "imgs/pic", 10, 6, 24*1e+8);

}
