// This optim.* refactors optimization.* in slsm namespace in order to be used in openMDAO project
// nlopt callback funciton is used in optimise namespace, but would be deprecated when another solver is used

#ifndef _OPTIM_H
#define _OPTIM_H

#include "Boundary.h"
#include <nlopt.hpp>

namespace Optimise_m
{

struct NLoptWrapper
{
    unsigned int index;
    void* callback;

};

// double callbackWrapper(const std::vector<double>& lambda, std::vector<double>& gradient, void* data);
// NOTE: only one constraint is considered: 

class Optimise
{
    public: 

    Optimise(std::vector<slsm::BoundaryPoint>& boundarypoints_, std::vector<double>& constraintDistance_,
            std::vector<double>& lambdas_, std::vector<double>& Multipliers_, double max_Displacement_ = 0.5,
            bool isMax_ = false, const std::vector<bool>& isEquality_ = {});

    nlopt::algorithm algorithm;

    // double solve();
    // double callback(const std::vector<double>& lambda, std::vector<double>& gradient, unsigned int index);

    // MEMO: reference must be initialized as they are declared
    std::vector<slsm::BoundaryPoint>& boundaryPoints;
    std::vector<double>& constraintDistances;
    std::vector<double> constraintDistancesScaled;
    const std::vector<bool> isEquality;

    unsigned int nConstraints = 1;
    unsigned int nConstraintsInitial;
    std::vector<double>& lambdas;
    std::vector<double>& LagrangeMultipliers;

    std::vector<unsigned int> indexMap; // 1 if obj/constraint is activated
    
    double timeStep;
    double max_Displacement;
    bool isMax;
    
    unsigned int nPoints; // number of boundarypoints

    std::vector<double> negativeLambdaLimits;
    std::vector<double> positiveLambdaLimits;
    std::vector<double> scaleFactors;

   
    // private:

};

class Pre_optimise : public Optimise
{
    public:
    // constraintDistance is a scalar value
    // timeStep = abs(lambdas[0]) hence not used here
    Pre_optimise(std::vector<slsm::BoundaryPoint>& boundarypoints_, std::vector<double>& constraintDistance_,
            std::vector<double>& lambdas_, std::vector<double>& Multipliers_, double max_Displacement_ = 0.5,
            bool isMax_ = false);
        
    void computeScaleFactors();
    // for effective optimization

    void computeConstraintDistances(unsigned int nCurrentConstraints=1);
    /*        if constraints are violated: scale constraint distances
        if constraints are satisfied: remove constraints*/

    void computeLambdaLimits();
    /*computes lamba limits by minimum displacements that violates CFL condition 
        for each function*/


    // double computeFunction(unsigned int index);
    // // compute the change in the objective/constraint function

    void computeGradients(const std::vector<double>& lambda, std::vector<double>& gradients, unsigned int index);
    // comptue the gradient of the objective/constraints function

    // double callback_obj(const std::vector<double>& lambda, std::vector<double>& gradients);
    // double callback_const(const std::vector<double>& lambda, std::vector<double>& gradients);
};

class Post_optimise : public Optimise
{
    public: 
    Post_optimise(std::vector<slsm::BoundaryPoint>& boundarypoints_, std::vector<double>& constraintDistance_,
            std::vector<double>& lambdas_, std::vector<double>& Multipliers_, double max_Displacement_ = 0.5,
            bool isMax_ = false);

    // void comptueDisplacement(std::vector<double> lambda);
    // // displacement = scaleFactor * lambda *boundarySensitivity * length
    
    // void rescaleDisplacements();
    // // check for cfl violation

    // void update_velocity();
    // // boundarypoints.velocity = displacement / timestep
    
};
}

#endif