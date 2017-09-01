#ifndef LIN_MODULE_H
#define LIN_MODULE_H

#include <iostream>
#include <vector>
#include <fstream>

#include "feasrc/fea.h"
#include "slsmsrc/slsm.h"


// int main();
class fea::Mesh;
//class fea::CHakElasticity;
class slsm::LevelSet;
class slsm::Boundary;
class slsm::Mesh;

class Lin_Module
{
    public:
        Lin_Module();

        std::vector<double> Field;

        //std::vector<fea::CHakElasticity> Materials;

        fea::Mesh* get_CMesh();
        slsm::Mesh* get_mesh();
        slsm::LevelSet* get_levelset();
        slsm::Boundary* get_boundary();
        fea::ElasticitySensitivities<fea::BilinearQuad>* get_CSensitivities();

        std::vector<double>* get_Field();
        std::vector<double> get_Field2();

        void run(fea::Mesh& CMesh, slsm::Mesh& mesh, slsm::LevelSet& levelSet, 
								slsm::Boundary& boundary, slsm::InputOutput& io, std::vector<double>& Field);
    
    protected:
        fea::Mesh* CMesh;
        fea::ElasticitySensitivities<fea::BilinearQuad>* CSensitivities;

        slsm::Mesh* mesh;
        slsm::LevelSet* levelSet;
        slsm::Boundary* boundary;
        // slsm::BoundaryPoint ;
        
        
        

};

#endif