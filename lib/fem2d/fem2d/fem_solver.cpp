#include "fem_solver.h"
#include "LinAlg_MDO.cpp"
#include <assert.h>

FEMSolver::FEMSolver(
  int num_nodes_x, int num_nodes_y, double length_x, double length_y,
  double E, double nu){
  this->num_nodes_x = num_nodes_x;
  this->num_nodes_y = num_nodes_y;
  this->length_x = length_x;
  this->length_y = length_y;

  dr_dx = 2.0 / length_x * (num_nodes_x - 1);
  ds_dy = 2.0 / length_y * (num_nodes_y - 1);

  D_voigt.resize(3, vector<double>(3,0.0));

  compute_nodes();
  compute_elems();
  compute_D(E, nu);
  
  compute_Ke();
}

void FEMSolver::compute_D(double E, double nu){
  double sc = E/(1.0  -nu*nu);

  // E / (1 - nu^2) *
  // [ 1, nu, 0       ]
  // [nu,  1, 0       ]
  // [ 0,  0, (1-nu)/2]

  D_voigt.resize(3, vector<double>(3,0.0));
  D_voigt[0][0] = sc * 1;
  D_voigt[1][1] = sc * 1;
  D_voigt[0][1] = sc * nu;
  D_voigt[1][0] = sc * nu;
  D_voigt[2][2] = sc * 0.5 * (1 - nu);
}

FEMSolver::~FEMSolver() {
}

void FEMSolver::compute_elems() {
  int node_i, node_j;

  elems.resize(num_nodes_x - 1, vector<vector<int> >(num_nodes_y - 1, vector<int>(8, 0)));

  for (int i = 0; i < num_nodes_x - 1; i++){
    for (int j = 0; j < num_nodes_y - 1; j++){
      // -1, -1
      elems[i][j][0] = (i + 0) * num_nodes_y * 2 + (j + 0) * 2 + 0;
      elems[i][j][1] = (i + 0) * num_nodes_y * 2 + (j + 0) * 2 + 1;
      //  1, -1
      elems[i][j][2] = (i + 1) * num_nodes_y * 2 + (j + 0) * 2 + 0;
      elems[i][j][3] = (i + 1) * num_nodes_y * 2 + (j + 0) * 2 + 1;
      //  1,  1
      elems[i][j][4] = (i + 1) * num_nodes_y * 2 + (j + 1) * 2 + 0;
      elems[i][j][5] = (i + 1) * num_nodes_y * 2 + (j + 1) * 2 + 1;
      // -1,  1
      elems[i][j][6] = (i + 0) * num_nodes_y * 2 + (j + 1) * 2 + 0;
      elems[i][j][7] = (i + 0) * num_nodes_y * 2 + (j + 1) * 2 + 1;
    }
  }
}

void FEMSolver::compute_nodes() {
  nodes.resize(num_nodes_x, vector<vector<double> >(num_nodes_y, vector<double>(2,0))); // (x,y) for indx. (i,j)
  double sp_x = length_x/(num_nodes_x-1), sp_y = length_y/(num_nodes_y-1); // spacing btw. nodes
  for (int ii = 0; ii < num_nodes_x; ii++){
    for (int jj = 0; jj < num_nodes_y; jj++){
      nodes[ii][jj][0] = ii*sp_x;
      nodes[ii][jj][1] = jj*sp_y;
    }
  }

  // print(nodes);
}

void FEMSolver::compute_Ke(){
  Matrix B;
  double Area = (length_x*length_y)/((num_nodes_x-1)*(num_nodes_y-1));

  Ke.resize(8, vector<double>(8,0.0));
  Ke0.resize(8, vector<double>(8,0.0));
  Ke1.resize(8, vector<double>(8,0.0));
  Ke2.resize(8, vector<double>(8,0.0));
  Ke3.resize(8, vector<double>(8,0.0));

  // vector<double> ri = {-1./sqrt(3), +1./sqrt(3), +1./sqrt(3), -1./sqrt(3)};
  // vector<double> si = {-1./sqrt(3), -1./sqrt(3), +1./sqrt(3), +1./sqrt(3)};
  // vector<double> wi = {1, 1, 1, 1};

  vector<double> ri, si, wi;
  ri.resize(4);
  si.resize(4);
  wi.resize(4);

  ri[0] = -1./sqrt(3);
  ri[1] =  1./sqrt(3);
  ri[2] =  1./sqrt(3);
  ri[3] = -1./sqrt(3);

  si[0] = -1./sqrt(3);
  si[1] = -1./sqrt(3);
  si[2] =  1./sqrt(3);
  si[3] =  1./sqrt(3);

  wi[0] = 1;
  wi[1] = 1;
  wi[2] = 1;
  wi[3] = 1;


  // u = N1 * u1 + N2 * u2 + N3 * u3 + N4 * u4
  // v = N1 * v1 + N2 * v2 + N3 * v3 + N4 * v4

  // N1 = 0.25 * (1 - r) * (1 - s);
  // N2 = 0.25 * (1 + r) * (1 - s);
  // N3 = 0.25 * (1 + r) * (1 + s);
  // N4 = 0.25 * (1 - r) * (1 + s);


  for (int gg = 0; gg < 4; gg++){
    B.resize(3, vector<double>(8,0.0));
    double r = ri[gg], s = si[gg], w = wi[gg];

    // du_dx
    B[0][0] = -0.25 * (1-s) * dr_dx;
    B[0][2] =  0.25 * (1-s) * dr_dx;
    B[0][4] =  0.25 * (1+s) * dr_dx;
    B[0][6] = -0.25 * (1+s) * dr_dx;

    // dv_dy
    B[1][1] = -0.25 * (1-r) * ds_dy;
    B[1][3] = -0.25 * (1+r) * ds_dy;
    B[1][5] =  0.25 * (1+r) * ds_dy;
    B[1][7] =  0.25 * (1-r) * ds_dy;

    // du_dy
    B[2][0] = -0.25 * (1-r) * ds_dy;
    B[2][2] = -0.25 * (1+r) * ds_dy;
    B[2][4] =  0.25 * (1+r) * ds_dy;
    B[2][6] =  0.25 * (1-r) * ds_dy;

    // dv_dx
    B[2][1] = -0.25 * (1-s) * dr_dx;
    B[2][3] =  0.25 * (1-s) * dr_dx;
    B[2][5] =  0.25 * (1+s) * dr_dx;
    B[2][7] = -0.25 * (1+s) * dr_dx;


    Matrix BT = transpose(B);
    Matrix BD_tmp = dot(BT, D_voigt);
    Matrix BDB = dot(BD_tmp, B);
    double detjw = Area/4.0*w;
    BDB = BDB * detjw;

    if (gg == 0) {
      Ke0 = BDB;
    } else if (gg == 1) {
      Ke1 = BDB;
    } else if (gg == 2) {
      Ke2 = BDB;
    } else if (gg == 3) {
      Ke3 = BDB;
    }
    Ke = Ke0 + Ke1;
    Ke = Ke + Ke2;
    Ke = Ke + Ke3;
  }
}

void FEMSolver::get_stiffness_matrix(double* data, int* rows, int* cols) {
  int index = 0;

  for (int ielem_x = 0; ielem_x < num_nodes_x - 1; ielem_x++) {
    for (int ielem_y = 0; ielem_y < num_nodes_y - 1; ielem_y++) {
      for (int imat_x = 0; imat_x < 8; imat_x++) {
        for (int imat_y = 0; imat_y < 8; imat_y++) {
          data[index] = Ke[imat_x][imat_y] * area_fraction[ielem_x][ielem_y];
          rows[index] = elems[ielem_x][ielem_y][imat_x];
          cols[index] = elems[ielem_x][ielem_y][imat_y];
          index += 1;
        }
      }
    }
  }

  // Lagrange multipliers
  int inode_x = 0;
  int idof = 0;
  int num_dofs = num_nodes_x * num_nodes_y * 2;

  for (int k = 0; k < 2; k++) {
    for (int inode_y = 0; inode_y < num_nodes_y; inode_y++) {
      idof = inode_x * num_nodes_y * 2 + inode_y * 2 + k;

      data[index] = 1;
      rows[index] = idof;
      cols[index] = inode_y * 2 + k + num_dofs;
      index += 1;

      data[index] = 1;
      rows[index] = inode_y * 2 + k + num_dofs;
      cols[index] = idof;
      index += 1;
    }
  }
}

void FEMSolver::get_stiffness_matrix(double* multipliers, double* data, int* rows, int* cols) {
  int index = 0;

  for (int ielem_x = 0; ielem_x < num_nodes_x - 1; ielem_x++) {
    for (int ielem_y = 0; ielem_y < num_nodes_y - 1; ielem_y++) {
      for (int imat_x = 0; imat_x < 8; imat_x++) {
        for (int imat_y = 0; imat_y < 8; imat_y++) {
          data[index] = Ke0[imat_x][imat_y] * multipliers[(ielem_x + 0) * num_nodes_y + (ielem_y + 0)];
          rows[index] = elems[ielem_x][ielem_y][imat_x];
          cols[index] = elems[ielem_x][ielem_y][imat_y];
          index += 1;

          data[index] = Ke1[imat_x][imat_y] * multipliers[(ielem_x + 1) * num_nodes_y + (ielem_y + 0)];
          rows[index] = elems[ielem_x][ielem_y][imat_x];
          cols[index] = elems[ielem_x][ielem_y][imat_y];
          index += 1;

          data[index] = Ke2[imat_x][imat_y] * multipliers[(ielem_x + 1) * num_nodes_y + (ielem_y + 1)];
          rows[index] = elems[ielem_x][ielem_y][imat_x];
          cols[index] = elems[ielem_x][ielem_y][imat_y];
          index += 1;

          data[index] = Ke3[imat_x][imat_y] * multipliers[(ielem_x + 0) * num_nodes_y + (ielem_y + 1)];
          rows[index] = elems[ielem_x][ielem_y][imat_x];
          cols[index] = elems[ielem_x][ielem_y][imat_y];
          index += 1;
        }
      }
    }
  }

  // Lagrange multipliers
  int inode_x = 0;
  int idof = 0;
  int num_dofs = num_nodes_x * num_nodes_y * 2;

  for (int k = 0; k < 2; k++) {
    for (int inode_y = 0; inode_y < num_nodes_y; inode_y++) {
      idof = inode_x * num_nodes_y * 2 + inode_y * 2 + k;

      data[index] = 1;
      rows[index] = idof;
      cols[index] = inode_y * 2 + k + num_dofs;
      index += 1;

      data[index] = 1;
      rows[index] = inode_y * 2 + k + num_dofs;
      cols[index] = idof;
      index += 1;
    }
  }
}

void FEMSolver::get_stiffness_matrix_derivs(double* states, double* data, int* rows, int* cols) {
  int index = 0;

  for (int ielem_x = 0; ielem_x < num_nodes_x - 1; ielem_x++) {
    for (int ielem_y = 0; ielem_y < num_nodes_y - 1; ielem_y++) {
      for (int imat_x = 0; imat_x < 8; imat_x++) {
        for (int imat_y = 0; imat_y < 8; imat_y++) {
          data[index] = Ke0[imat_x][imat_y] * states[elems[ielem_x][ielem_y][imat_y]];
          rows[index] = elems[ielem_x][ielem_y][imat_x];
          cols[index] = (ielem_x + 0) * num_nodes_y + (ielem_y + 0);
          index += 1;

          data[index] = Ke1[imat_x][imat_y] * states[elems[ielem_x][ielem_y][imat_y]];
          rows[index] = elems[ielem_x][ielem_y][imat_x];
          cols[index] = (ielem_x + 1) * num_nodes_y + (ielem_y + 0);
          index += 1;

          data[index] = Ke2[imat_x][imat_y] * states[elems[ielem_x][ielem_y][imat_y]];
          rows[index] = elems[ielem_x][ielem_y][imat_x];
          cols[index] = (ielem_x + 1) * num_nodes_y + (ielem_y + 1);
          index += 1;

          data[index] = Ke3[imat_x][imat_y] * states[elems[ielem_x][ielem_y][imat_y]];
          rows[index] = elems[ielem_x][ielem_y][imat_x];
          cols[index] = (ielem_x + 0) * num_nodes_y + (ielem_y + 1);
          index += 1;
        }
      }
    }
  }
}

void FEMSolver::get_sensitivity_LSTO(double* u, double* xpos, double* ypos, double* sens){
  // calculate sensitivities at the gausspoints
  double Area = (length_x*length_y)/((num_nodes_x-1)*(num_nodes_y-1));
  double wj = Area/4.0;

  vector<double> ri = {-1./sqrt(3), +1./sqrt(3), +1./sqrt(3), -1./sqrt(3)};
  vector<double> si = {-1./sqrt(3), -1./sqrt(3), +1./sqrt(3), +1./sqrt(3)};
  Vector x0, y0, u0;
  Vector tmp_;
  double N1, N2, N3, N4;
  // vector<double> wi = {1, 1, 1, 1};
  int index = 0;
  for (int ii = 0; ii < num_nodes_x-1; ii++){
    for (int jj = 0; jj < num_nodes_y-1; jj++){
      for (int gg = 0; gg < 4; gg ++ ){
        double r = ri[gg], s = si[gg];
        x0 = {nodes[ii][jj][0], nodes[ii+1][jj][0], nodes[ii+1][jj+1][0],nodes[ii][jj+1][0]};
        y0 = {nodes[ii][jj][1], nodes[ii+1][jj][1], nodes[ii+1][jj+1][1],nodes[ii][jj+1][1]};
        u0 = {u[elems[ii][jj][0]], u[elems[ii][jj][1]], u[elems[ii][jj][2]], u[elems[ii][jj][3]],
          u[elems[ii][jj][4]], u[elems[ii][jj][5]], u[elems[ii][jj][6]], u[elems[ii][jj][7]]};
        N1 = 0.25 * (1 - r) * (1 - s);
        N2 = 0.25 * (1 + r) * (1 - s);
        N3 = 0.25 * (1 + r) * (1 + s); 
        N4 = 0.25 * (1 - r) * (1 + s);
        double x_ = x0[0]*N1 + x0[1]*N2 + x0[2]*N3 + x0[3]*N4;  
        double y_ = y0[0]*N1 + y0[1]*N2 + y0[2]*N3 + y0[3]*N4;  
        xpos[index] = x_; ypos[index] = y_;
        
        if (gg == 0) {
          tmp_ = dot(Ke0,u0);
          tmp_ = tmp_ * area_fraction[ii][jj];
          sens[index] = dot(u0,tmp_);
          sens[index] /= wj;
        } else if (gg == 1) {
          tmp_ = dot(Ke1,u0);
          tmp_ = tmp_ * area_fraction[ii][jj];
          sens[index] =dot(u0,tmp_);
          sens[index] /= wj;
        } else if (gg == 2) {
          tmp_ = dot(Ke2,u0);
          tmp_ = tmp_ * area_fraction[ii][jj]; 
          sens[index] =dot(u0,tmp_);
          sens[index] /= wj;
        } else if (gg == 3) {
          tmp_ = dot(Ke3,u0);
          tmp_ = tmp_ * area_fraction[ii][jj]; 
          sens[index] =dot(u0,tmp_);
          sens[index] /= wj;
        }
                
        index ++;
      }
    }
  }
  
}

void FEMSolver::set_area_fractions(double* areafraction){
  this->area_fraction.resize(num_nodes_x-1, Vector(num_nodes_y-1, 0.0));
  
  for (int ii = 0; ii < num_nodes_x-1; ii++){
    for (int jj = 0; jj < num_nodes_y-1; jj++){
      this->area_fraction[ii][jj] = areafraction[jj*(num_nodes_x-1)+ii];
    }
  }
}

//int main(){
  //FEMSolver fem2(1,1,1,1,1,0.5);
  // for (int ii = 0; ii < 3; ii++){
  //   cout << "\n";
  //   for (int jj = 0; jj < 3; jj++){
  //     cout << fem2.D_voigt[ii][jj] << "\t";
  //   }
  //   cout << "\n";
  // }

//}
