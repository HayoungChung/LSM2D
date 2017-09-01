#ifndef SLSM_MODULE_H
#define SLSM_MODULE_H

#include <iostream>
#include <vector>
#include <fstream>

#include "src/slsm.h"

class slsm::LevelSet;
class slsm::Boundary;
class slsm::Mesh;

class slsm_Module{
    public:
        slsm_Module();
        double wrapper_optimise(std::vector<slsm::BoundaryPoint>&, const std::vector<double>&,
        std::vector<double>&, std::vector <double>&, double);

    private:
        slsm::Mesh* mesh;
        slsm::LevelSet* levelSet;
        slsm::Boundary* boundary;
        slsm::InputOutput* io;
        slsm::Optimise* optimise;
        

};


#endif