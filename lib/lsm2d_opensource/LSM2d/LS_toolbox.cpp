#include "LS_toolbox.h"
#include "LinAlg_MDO.cpp"
#include <assert.h>

Levelset2D::Levelset2D(int num_nodes_x, int num_nodes_y, double maxArea)
: mesh(num_nodes_x-1,num_nodes_y-1), levelset(mesh), boundary(levelset)
{
    this->num_nodes_x = num_nodes_x;
    this->num_nodes_y = num_nodes_y;
    this->num_elem_x = num_nodes_x-1;
    this->num_elem_y = num_nodes_y-1;
    this->discretize();
    this->maxArea = maxArea;
    this->meshArea = (double) (num_elem_x*num_elem_y);
}

void Levelset2D::get_area_fractions(double* areafraction){
    // assume that discretization is already called
    for (int ii = 0; ii < (num_elem_x * num_elem_y) ; ii++){
        areafraction[ii] = mesh.elements[ii].area;
    }
}

double Levelset2D::get_num_boundary_coords(){
    return boundary.points.size();
}

void Levelset2D::get_boundary_coords(double* x, double* y){
    // assume that discretization is already called
    for (int jj = 0; jj < boundary.points.size(); jj ++){
            x[jj] = boundary.points[jj].coord.x;
            y[jj] = boundary.points[jj].coord.y;
    }
}

void Levelset2D::get_phi(double* phi){
    for (int ii = 0; ii < num_nodes_x*num_nodes_y; ii++){
        phi[ii] = this->levelset.signedDistance[ii];
    }
}

void Levelset2D::get_delphi(double* delphi){
    double timeStep; // unchanged
    double movelimit = 0.5;
    Optimise optimise(boundary.points,  timeStep, movelimit) ;

    optimise.length_x = mesh.width;
    optimise.length_y = mesh.height;
    optimise.boundary_area = boundary.area; // area of structure
    optimise.mesh_area = meshArea; // area of the entire mesh
    optimise.max_area = maxArea; // maximum area

    // Perform the optimisation.
    optimise.Solve_With_NewtonRaphson() ;
    // LSM::MersenneTwister rng ;
    
    levelset.computeVelocities(boundary.points);//, timeStep, 0, rng) ;
    
    levelset.computeGradients() ;
    for (int ii = 0; ii < num_nodes_x*num_nodes_y; ii++){
        delphi[ii] = levelset.gradient[ii]*levelset.velocity[ii];        
    }
}

void Levelset2D::set_sensitivities(double* FEAsensitivities){
    for (int ii = 0; ii < boundary.points.size(); ii++){
        boundary.points[ii].sensitivities[0] = -FEAsensitivities[ii];
        boundary.points[ii].sensitivities[1] = -1;
    }    
}

void Levelset2D::reinitialize(){
    levelset.reinitialise();
}

void Levelset2D::update(){
    levelset.update(1) ;    
}

void Levelset2D::set_signedDistance(double* phi){
    for (int ii = 0; ii < (num_nodes_x*num_nodes_y) ; ii++){
        this->levelset.signedDistance[ii] = phi[ii];
    }
    this->discretize();
}

void Levelset2D::discretize(){
    this->boundary.discretise(false, 2);
    this->boundary.computeAreaFractions();
}


/*
int main(){
    Levelset2D testclass(160,80,0.3);
    // cout << testclass.mesh.width << endl;
    LSM::InputOutput io;
    io.saveLevelSetVTK('a', testclass.levelset);
    for (int ii = 0;ii < testclass.boundary.points.size(); ii++){
        cout << testclass.boundary.points[ii].coord.x << "\t" << 
         testclass.boundary.points[ii].coord.y << endl;
    }
    
    return 0;
}
*/