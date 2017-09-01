// this program is based on M2DO-opensource v1

#include <iostream>
#include <vector>
#include <cmath>
// #include "./M2DO_LSM.h"

using namespace std;
// namespace LSM = M2DO_LSM;

// #include "debug.h"
// #include "min_unit.h"
#include "common.h"

#include "mesh.h"
#include "level_set.h"
#include "boundary.h"
#include "optimise.h"

typedef vector<vector<double> > Matrix;
typedef vector<double> Vector;

class Levelset2D {
    public:
    // constructor: get fixed grid mesh & initialize level set to the nodes
    Levelset2D(int num_node_x, int num_nodes_y, double maxArea); // in current form, exy = Lxy (size = 1)
    ~Levelset2D(){};

    void get_area_fractions(double* areafraction); // toFEA
    double get_num_boundary_coords();
    void get_boundary_coords(double* x, double* y); // toFEA
    
    void get_phi(double* phi);
    void get_delphi(double* delphi); // to ODE

    void set_sensitivities(double* FEAsensitivities); // fromFEA
    void reinitialize();
    void update(); // for a comparison with previous results

    void set_signedDistance(double* phi); // set level set function from Python
    // void compute_velocity();
        
    void get_lambdaLimits(double* negLambdaLim, double* posLambdaLim);

    class Mesh mesh; 
    class LevelSet levelset;
    class Boundary boundary;
    class Optimise* optimise;
    
    void discretize();
    private:

    Vector areafraction;
    int num_nodes, num_elems;
    int num_nodes_x, num_nodes_y, num_elem_x, num_elem_y;
    double length_x, length_y;
    vector<BoundaryPoint> boundarypoints;
    vector<double> displacements;
    double maxArea, meshArea;
};
