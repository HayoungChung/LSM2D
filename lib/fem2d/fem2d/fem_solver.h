#include <iostream>
#include <vector>

// NOTE THAT the structured mesh is assumed herein
// if isNodal = false, then elem-wise desvar

using namespace std;
typedef vector<vector<double> > Matrix;
typedef vector<double> Vector;

class FEMSolver {
public:
  FEMSolver(
    int num_nodes_x, int num_nodes_y, double length_x, double length_y,
    double E, double nu
  );
  ~FEMSolver();
  // compute K_tot with the multipliers as design variables given 
  void get_stiffness_matrix(double* multipliers, double* data, int* rows, int* cols, bool isNodal = true); 
  // compute K_tot with areafraction specified
  void get_stiffness_matrix(double* data, int* rows, int* cols); 
  
  // compute Kij * uj (partial of residuals): Hwang's node-wise desvar
  void get_stiffness_matrix_derivs(double* states, double* data, int* rows, int* cols, bool isNodal = true);
  void set_area_fractions(const double* areafraction);
  void get_sensitivity_LSTO(double* u, double* xpos, double* ypos, double* sens);

private:
  Matrix area_fraction;
  int num_nodes, num_elems;
  int num_nodes_x, num_nodes_y;
  double length_x, length_y;
  double dr_dx, ds_dy;
  vector<vector<vector<double> > > nodes;
  vector<vector<vector<int> > > elems;
  vector<vector<vector<int> > > cnnt;
  Matrix D_voigt;
  Matrix Ke, Ke0, Ke1, Ke2, Ke3; // K matrices (indices 0-3 indicate each gausspoint)
  void compute_nodes();
  void compute_elems();
  void compute_D(double E, double nu);
  void compute_Ke();
};
