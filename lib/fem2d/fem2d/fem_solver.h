#include <iostream>
#include <vector>

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
  void get_stiffness_matrix(double* multipliers, double* data, int* rows, int* cols);
  void get_stiffness_matrix(double* data, int* rows, int* cols);
  void get_stiffness_matrix_derivs(double* states, double* data, int* rows, int* cols);
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
  Matrix Ke, Ke0, Ke1, Ke2, Ke3;
  void compute_nodes();
  void compute_elems();
  void compute_D(double E, double nu);
  void compute_Ke();
};
