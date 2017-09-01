#include "LS_toolbox.h"
#include "LinAlg_MDO.cpp"
#include <assert.h>

Levelset2D::Levelset2D(int num_nodes_x, int num_nodes_y, int maxArea)
: mesh(num_nodes_x-1,num_nodes_y-1), levelset(mesh), boundary(levelset)
{
    this->num_nodes_x = num_nodes_x;
    this->num_nodes_y = num_nodes_y;
    this->num_elem_x = num_nodes_x-1;
    this->num_elem_y = num_nodes_y-1;
    this->discretize();
    this->maxArea = maxArea;
    this->meshArea = (double) (num_elem_x*num_elem_y);
    this->Aux.reserve(6);
}

void Levelset2D::get_area_fractions(double* areafraction){
    for (int ii = 0; ii < (num_elem_x * num_elem_y) ; ii++){
        areafraction[ii] = mesh.elements[ii].area;
    }
}

void Levelset2D::set_signedDistance(double* phi){
    for (int ii = 0; ii < (num_nodes_x*num_nodes_y) ; ii++){
        this->levelset.signedDistance[ii] = phi[ii];
    }
    this->levelset.reinitialise();
}

void Levelset2D::discretize(){
    this->boundary.discretise();
    this->boundary.computeAreaFractions();
}


// this computes functions for the sub-optimization 
void Levelset2D::get_suboptim_functions(double* lambdas, double* Aux){
    vector<double> lambdas_internal(2);
    lambdas_internal[0] = lambdas[0];
    lambdas_internal[1] = lambdas[1];

    this->set_suboptim(lambdas_internal);

    Aux = this->Aux.data();    
    // f_obj = Aux[0]; // objective function
    // G_const = Aux[1]; // Constraints function
    // limit0[0] = Aux[2]; // min of lambda[0]
    // limit0[1] = Aux[3]; // max of lambda[0]
    // limit1[0] = Aux[4]; // min of lambda[1]
    // limit1[1] = Aux[5]; // max of lambda[1]  
}

void Levelset2D::set_suboptim(vector<double> lambdas0){
    vector<double> constraintDistances;
    constraintDistances.push_back(meshArea*maxArea - boundary.area);
    double timestep = 0; //-lambdas0[0];
    LSM::Optimise optim(boundarypoints, constraintDistances, lambdas0, timestep, levelset.moveLimit, false);
    optim.computeScaleFactors(); // scaling objective function, constraint (function of sensitivity only)
    for (int i = 0; i<2; i++){
        lambdas0[i] /= optim.scaleFactors[i];
    }
    optim.computeLambdaLimits();
    optim.computeConstraintDistances();
    Vector emptyGrad; 
    Aux[0] = optim.callback(lambdas0,emptyGrad,0); // objective
    Aux[1] = optim.callback(lambdas0,emptyGrad,1); // constraint
    Aux[2] = optim.negativeLambdaLimits[0];
    Aux[3] = optim.positiveLambdaLimits[0];
    Aux[4] = optim.negativeLambdaLimits[1];
    Aux[5] = optim.positiveLambdaLimits[1];
}

void Levelset2D::compute_bndVel(vector<double> lambdas){
    vector<double> constraintDistances;
    constraintDistances.push_back(meshArea*maxArea - boundary.area);
    double timestep = -lambdas[0];
    LSM::Optimise optim(boundarypoints, constraintDistances, lambdas, timestep, levelset.moveLimit,false);
    optim.computeScaleFactors(); // scaling objective function, constraint (function of sensitivity only)
    
    // compute velocities at boundary points
    // Compute the optimum displacement vector.
    optim.computeDisplacements(lambdas);
    
    // Calculate the unscaled lambda values.
    for (unsigned int i=0;i< 2;i++)
        lambdas[i] *= optim.scaleFactors[i];

    // Calculate boundary point velocities.
    for (unsigned int i=0;i<boundarypoints.size();i++)
        boundarypoints[i].velocity = displacements[i] / std::abs(lambdas[0]) / 0.05;
        
}

void get_velocity(double* lambdas, double* velocities){

}
// double* Levelset2D::compute_scaleFactors(){
//     double scaleFactors[2];
//     double maxSens = 0;
//     for (unsigned int i = 0; i < 2; i++){ // objective = 1, constraint = 1;
//         for (unsigned int j = 0; j < nPoints; j ++){
//             if (!boundarypoints[j].isFixed){
//                 double sens = abs(boundarypoints[j].sensitivities[i]);
//                 if (sens > maxSens) maxSens = sens;
//             }
//         }
//     }
//     scaleFactors[i] = (1.0/maxSens);

//     vector<double> lambda(2,0.0);
//     vector<double> gradient(2,0.0);
//     for (unsigned int i = 0; i < 2; i++ ){
//         computeGradients(lambda, gradient, i);
//         scaleFactors[i] *= (1.0 / abs(gradient[i]));
//     }
//     return scaleFactors;
// }




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