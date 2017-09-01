// A python-wrapping attribute with Cython

#ifndef M2DO_FEA_MA57_SOLVER_H
#define M2DO_FEA_MA57_SOLVER_H


#include "/home/hac210/Dropbox/packages/lib/eigen3/Eigen/Dense"
#include "/home/hac210/Dropbox/packages/lib/eigen3/Eigen/Sparse"
#include <vector>
#include <iostream>

using namespace std ;
using namespace Eigen ;
/* 
	MA57 Functions
*/

extern "C" void ma57id_(double *cntl, int *icntl) ;

extern "C" void ma57ad_(int *n, int *nz, int *irn, int *jcn, int *lkeep, int *keep, int *iw, int *icntl, int *info, double *rinfo) ;
    
extern "C" void ma57bd_(int *n, int *nz, double *a, double *fact, int *lfact, int *ifact, int *lifact, int *lkeep, int *keep, int *iw, int *icntl, double *cntl, int *info, double *rinfo) ;
    
extern "C" void ma57cd_(int *job, int *n, double *fact, int *lfact, int *ifact, int *lifact, int *nrhs, double *rhs, int *lrhs, double *w, int *lw, int *iw, int *icntl, int *info) ;

class MA57Solver {
	
	private:
		//

	public:
		// Properties:
		VectorXd cntl, fact ;
		VectorXi icntl, ifact, iw, info ;
		int n, lfact, lifact ;

		// Methods:
		MA57Solver (bool use_metis, bool print_output) ;
		void print () ;
		
		// python -> Cython
		void pyCompute (int n, int nz, int* irn, int* jcn, double* a); 
		void pySolve (double* b); 
		
		// check
		void compute (SparseMatrix<double> & A) ;
		VectorXd solve (VectorXd & b) ;
		
		// VectorXd compute_pyCompute(SparseMatrix<double> & A, VectorXd & b);


} ;

// #include "/home/hac210/Dropbox/packages/topOpt_MDO/FEM2D/speed_test_HSL/m2do_ma57Wrapper/ma57_solver.cpp"
#endif
