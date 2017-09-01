#include <iostream>
#include <vector>
#include <cmath>
#include "./M2DO_LSM.h"

using namespace std;
namespace LSM = M2DO_LSM;

typedef vector<vector<double> > Matrix;
typedef vector<double> Vector;

class Levelset2D {
    public:
    // constructor: get fixed grid mesh & initialize level set to the nodes
    Levelset2D(int num_node_x, int num_nodes_y, int maxArea); // in current form, exy = Lxy (size = 1)
    ~Levelset2D(){};

    void get_area_fractions(double* areafraction);
    void get_suboptim_functions(double* lambdas, double* Aux); 
    void get_velocity(double* lambdas, double* velocities);
    void get_gradphi(double* gradphi);
    // void reinitialize();
    // void update(); // to compare ODE results
    // void subOptimize(); // compare with openMDAO results
    void set_suboptim(Vector lambdas);
    void set_signedDistance(double* phi); // set level set function from Python
    void compute_bndVel(vector<double> lambdas);
        
    class LSM::Mesh mesh; // from mesh class 
    class LSM::LevelSet levelset;
    class LSM::Boundary boundary;
    
    private:
    void discretize();

    Vector areafraction;
    int num_nodes, num_elems;
    int num_nodes_x, num_nodes_y, num_elem_x, num_elem_y;
    double length_x, length_y;
    vector<LSM::BoundaryPoint> boundarypoints;
    vector<double> displacements;
    double maxArea, meshArea;
    Vector Aux;

    // double* compute_scaleFactors();
    
    // // reused
    // LSM::LevelSet& levelset;
    // LSM::Boundary& boundary;
    // Vector nodal_velocities;
    // Vector nodal_gradients;
    // Vector lambdas; 
};
