/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+lsm@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <iostream>

#include "Debug.h"
#include "Optimise_boost.h"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

/*! \file Optimise.cpp
    \brief A class for finding the solution for the optimum velocity vector.
 */

namespace optim
{

double callbackWrapper(const std::vector<double>& lambda, std::vector<double>& gradient, void* data)
{
    NLoptWrapper* wrapperData = reinterpret_cast<NLoptWrapper*>(data);
    return reinterpret_cast<Optimise*>(wrapperData->callback)->callback(lambda, gradient, wrapperData->index);
}

Optimise::Optimise(slsm::LevelSet& levelSet_,
    std::vector<slsm::BoundaryPoint>& boundaryPoints_,
                   const std::vector<double>& constraintDistances_,
                   std::vector<double>& lambdas_,
                   std::vector <double>& Multipliers_,
                   //double& timeStep_,
                   double maxDisplacement_,
                   bool isMax_,
                   const std::vector<bool>& isEquality_,
                   nlopt::algorithm algorithm_) :
                   levelSet(levelSet_),
                   boundaryPoints(boundaryPoints_),
                   constraintDistances(constraintDistances_),
                   lambdas(lambdas_),
                   LagrangeMultipliers (Multipliers_),
                //    timeStep(timeStep_),
                   maxDisplacement(maxDisplacement_),
                   isMax(isMax_),
                   isEquality(isEquality_),
                   algorithm(algorithm_)
{
    errno = EINVAL;
    slsm_check(((maxDisplacement_ > 0) && (maxDisplacement_ < 1)), "Move limit must be between 0 and 1.");
    
    // Store the initial number of constraints.
    nConstraints = lambdas.size() - 1;
    nConstraintsInitial = nConstraints;

    // Check for empty constraint distances vector.
    if (constraintDistances.empty())
    {
        errno = 0;
        slsm_log_warn("Empty constraint distances vector. Assuming unconstrained optimisation.");
        nConstraints = 0;
    }

    // Resize data structures.
    negativeLambdaLimits.resize(nConstraints + 1);
    positiveLambdaLimits.resize(nConstraints + 1);
    scaleFactors.resize(nConstraints + 1);

    // Copy constraint distances.
    constraintDistancesScaled = constraintDistances;

    // Set default as inequality constraints.
    if (isEquality.size() == 0) isEquality.resize(nConstraints, false);

    // WIP
    num_points = levelSet.mesh.nNodes;
    phi.resize(num_points);
    dphi_dt.resize(num_points);

    return;

error:
    exit(EXIT_FAILURE);
}

double Optimise::callback(const std::vector<double>& lambda, std::vector<double>& gradient, unsigned int index)
{
    // Calculate the boundary displacement vector.
    computeDisplacements(lambda);

    // Compute the gradients.
    if (!gradient.empty())
        computeGradients(lambda, gradient, index);

    // Return the function value.
    return computeFunction(index);
}

void Optimise::preprocess(){
    // Store the number of boundary points.
    // This can change between successive optimisation calls.
    nPoints = boundaryPoints.size();

    // Resize boundary point dependent data structures.
    displacements.resize(nPoints);

		// Initializing index map making all constraints to be active.
    indexMap.resize(nConstraints + 1);
    std::iota(indexMap.begin(), indexMap.end(), 0);

		// Compute the scale factors for the objective and constraints.
    computeScaleFactors();

    // Scale inital lambda estimates.
    for (unsigned int i=0;i<nConstraints+1;i++)
    	lambdas[i] /= scaleFactors[i];

    // Compute the lambda limits.
    computeLambdaLimits();

    // Compute scaled constraint change distances and remove inactive
    // inequality constraints.
    if (nConstraints > 0) computeConstraintDistances(nConstraints);
}

void Optimise::compute_dphi_dt(){
    // postprocess(0);

    // narrowbanded properties(other = 0)
    levelSet.computeVelocities(boundaryPoints);
    levelSet.computeGradients();

    std::cout << "lambdas( " << lambdas.size() << ":) " << lambdas[0] << "," << lambdas[1] << std::endl;

    for (unsigned int i = 0; i < levelSet.mesh.nNodes; i++){
        dphi_dt[i] = levelSet.gradient[i]*levelSet.velocity[i];
    }

    //levelSet.reinitialise();    
    phi = levelSet.signedDistance;
}

double Optimise::postprocess(double optObjective_){
    // Compute the optimum displacement vector.
    computeDisplacements(lambdas);

    // Rescale the displacements and lambda values (if necessary).
    optObjective_ *= rescaleDisplacements();

    // Calculate the unscaled lambda values.
    for (unsigned int i=0;i<nConstraints+1;i++)
        lambdas[i] *= scaleFactors[i];

    // Effective time step.
    timeStep = std::abs(lambdas[0]);

    // Calculate boundary point velocities.
    for (unsigned int i=0;i<nPoints;i++)
        boundaryPoints[i].velocity = displacements[i] / timeStep;

    // Remap data if there are inactive constraints.
    if (nConstraints < nConstraintsInitial)
    {
        // Return the lambda vector to its original size.
        lambdas.resize(1 + nConstraintsInitial);

        // Check for constraints.
        if (nConstraints > 0)
        {
            // Copy the lambda values back into position.
            // Loop backwards to avoid overwriting data.
            for (unsigned int i=nConstraints+1;i>0;i--)
                lambdas[indexMap[i]] = lambdas[i];
        }
    }

    // Return the unscaled objective change.
    return (optObjective_ / scaleFactors[0]);
}

double Optimise::nlopt_activate(){
    // Create wrapper for objective.
    NLoptWrapper objectiveWrapper;
    objectiveWrapper.index = 0;
    objectiveWrapper.callback = this;

    // Create wrappers for constraints.
    NLoptWrapper constraintWrappers[nConstraints];
    for (unsigned int i=0;i<nConstraints;i++)
    {
        constraintWrappers[i].index = i + 1;
        constraintWrappers[i].callback = this;
    }

    // Instantiate NLopt optimisation object.
    nlopt::opt opt(algorithm, 1 + nConstraints);

    // Set limits.
    opt.set_lower_bounds(negativeLambdaLimits);
    opt.set_upper_bounds(positiveLambdaLimits);

    // Set a maximum of 500 iterations.
    opt.set_maxeval(500);
    //opt.set_maxeval(1000);

    // Set convergence tolerance.
    opt.set_xtol_rel(1e-8);
    //opt.set_xtol_rel(1e-2);

    // Specify whether we want to minimise or maximise the objective function.
    if (isMax) opt.set_max_objective(callbackWrapper, &objectiveWrapper);
    else       opt.set_min_objective(callbackWrapper, &objectiveWrapper);

    // Add the constraints.
    for (unsigned int i=0;i<nConstraints;i++)
    {
        // Add equality constraint.
        if (isEquality[i])
            opt.add_equality_constraint(callbackWrapper, &constraintWrappers[i], 1e-8);

        // Add inequality constraint.
        else
            opt.add_inequality_constraint(callbackWrapper, &constraintWrappers[i], 1e-8);
    }

    // The optimum value of the objective function.
    double optObjective;

    // Whether a single optimisation has been attempted.
    bool isAttempted = false;

    // Keep trying optmisation until a solution is found.
    while (!isAttempted || (returnCode < 0))
    {
        isAttempted = true;

        // Attempt optimisation.
        try
        {
            returnCode = opt.optimize(lambdas, optObjective);
        }

        // Catch roundoff errors.
        catch (nlopt::roundoff_limited)
        {
            // Reduce the constraint change targets.
            for (unsigned int i=0;i<nConstraints;i++)
                constraintDistancesScaled[i] *= 0.7;
                //constraintDistancesScaled[i] *= 0.2;
        }

        // Catch argument errors.
        catch (std::invalid_argument)
        {
            errno = EINVAL;
            slsm_log_err("Invalid arguments!");
        }
    }
    return optObjective;
}

double Optimise::solve()
{
    // Store the number of boundary points.
    // This can change between successive optimisation calls.
    nPoints = boundaryPoints.size();

    // Resize boundary point dependent data structures.
    displacements.resize(nPoints);

		// Initializing index map making all constraints to be active.
		indexMap.resize(nConstraints + 1);
    std::iota(indexMap.begin(), indexMap.end(), 0);

		// Compute the scale factors for the objective and constraints.
    computeScaleFactors();

    // Scale inital lambda estimates.
    for (unsigned int i=0;i<nConstraints+1;i++)
    	lambdas[i] /= scaleFactors[i];

    // Compute the lambda limits.
    computeLambdaLimits();

    // Compute scaled constraint change distances and remove inactive
    // inequality constraints.
    if (nConstraints > 0) computeConstraintDistances(nConstraints);
    
    // Create wrapper for objective.
    NLoptWrapper objectiveWrapper;
    objectiveWrapper.index = 0;
    objectiveWrapper.callback = this;

    // Create wrappers for constraints.
    NLoptWrapper constraintWrappers[nConstraints];
    for (unsigned int i=0;i<nConstraints;i++)
    {
        constraintWrappers[i].index = i + 1;
        constraintWrappers[i].callback = this;
    }

    // Instantiate NLopt optimisation object.
    nlopt::opt opt(algorithm, 1 + nConstraints);

    // Set limits.
    opt.set_lower_bounds(negativeLambdaLimits);
    opt.set_upper_bounds(positiveLambdaLimits);

    // Set a maximum of 500 iterations.
    opt.set_maxeval(500);
    //opt.set_maxeval(1000);

    // Set convergence tolerance.
    opt.set_xtol_rel(1e-8);
    //opt.set_xtol_rel(1e-2);

    // Specify whether we want to minimise or maximise the objective function.
    if (isMax) opt.set_max_objective(callbackWrapper, &objectiveWrapper);
    else       opt.set_min_objective(callbackWrapper, &objectiveWrapper);

    // Add the constraints.
    for (unsigned int i=0;i<nConstraints;i++)
    {
        // Add equality constraint.
        if (isEquality[i])
            opt.add_equality_constraint(callbackWrapper, &constraintWrappers[i], 1e-8);

        // Add inequality constraint.
        else
            opt.add_inequality_constraint(callbackWrapper, &constraintWrappers[i], 1e-8);
    }

    // The optimum value of the objective function.
    double optObjective;

    // Whether a single optimisation has been attempted.
    bool isAttempted = false;

    // Keep trying optmisation until a solution is found.
    while (!isAttempted || (returnCode < 0))
    {
        isAttempted = true;

        // Attempt optimisation.
        try
        {
            returnCode = opt.optimize(lambdas, optObjective);
        }

        // Catch roundoff errors.
        catch (nlopt::roundoff_limited)
        {
            // Reduce the constraint change targets.
            for (unsigned int i=0;i<nConstraints;i++)
                constraintDistancesScaled[i] *= 0.7;
                //constraintDistancesScaled[i] *= 0.2;
        }

        // Catch argument errors.
        catch (std::invalid_argument)
        {
            errno = EINVAL;
            slsm_log_err("Invalid arguments!");
        }
    }

    // Compute the optimum displacement vector.
    computeDisplacements(lambdas);

    // Rescale the displacements and lambda values (if necessary).
    optObjective *= rescaleDisplacements();

    // Calculate the unscaled lambda values.
    for (unsigned int i=0;i<nConstraints+1;i++)
        lambdas[i] *= scaleFactors[i];

    // Effective time step.
    timeStep = std::abs(lambdas[0]);

    // Calculate boundary point velocities.
    for (unsigned int i=0;i<nPoints;i++)
        boundaryPoints[i].velocity = displacements[i] / timeStep;

    // Remap data if there are inactive constraints.
    if (nConstraints < nConstraintsInitial)
    {
        // Return the lambda vector to its original size.
        lambdas.resize(1 + nConstraintsInitial);

        // Check for constraints.
        if (nConstraints > 0)
        {
            // Copy the lambda values back into position.
            // Loop backwards to avoid overwriting data.
            for (unsigned int i=nConstraints+1;i>0;i--)
                lambdas[indexMap[i]] = lambdas[i];
        }
    }

    // Return the unscaled objective change.
    return (optObjective / scaleFactors[0]);
}

void Optimise::computeScaleFactors()
{
    /* In order for the optimiser to work effectively it is important
       that small changes in the lambdas result in small changes in the
       functions for the objective and constraints and their respective
       gradients. Since we are solving a multi-dimensional optimisation
       problem it is important that all variables are on the same scale.
       This enables us to use a universal convergence tolerance.

       Our scaling protocol is described below. Note that we choose to
       store scale factors for each function, rather than scaling (then
       rescaling) the input data. Sensitivites are scaled down (reduced)
       and the lambda values are scale up (enlarged).

        1) For each function, scale by the largest absolute sensitivity,
           i.e. the maximum magnitude is one.

        2) Scale by the gradient at the origin (all lambdas are zero).
           This ensures that the gradient is close to one.

       There is no need to independently scale the boundary integral
       coefficients since the boundary lengths are independent of the
       function, i.e.

         c^f_i = s^f_i * l_i

       The l_i variables are the same for the objective and constraints
       so, by definition, the c^f's are on the same scale if we simply
       normalise by max(abs(s^f_i)).
     */

    // Loop over all functions: objective first, then constraints.
    for (unsigned int i=0;i<nConstraints+1;i++)
    {
        // Initialise maximum sensitivity.
        double maxSens = 0;

        // Loop over all boundary points.
        for (unsigned int j=0;j<nPoints;j++)
        {
            // Don't consider fixed points.
            if (!boundaryPoints[j].isFixed)
            {
                // Test whether sensitivity magnitude is current maximum.
                double sens = std::abs(boundaryPoints[j].sensitivities[i]); // Later try c instead of s
                if (sens > maxSens) maxSens = sens;
            }
        }

        // Store scale factor.
        scaleFactors[i] = (1.0 / maxSens);
    }

    // Create lambda vector (all zeros, i.e. at the origin).
    std::vector<double> lambda(nConstraints + 1, 0.0);

    // Initialise gradient vector.
    std::vector<double> gradient(nConstraints + 1);

    // Loop over all functions: objective first, then constraints.
    for (unsigned int i=0;i<nConstraints+1;i++)
    {
        // Calculate the gradient.
        computeGradients(lambda, gradient, i);

        // Scale by diagonal gradient entry (absolute value).
        scaleFactors[i] *= (1.0 / std::abs(gradient[i]));
    }
}

void Optimise::computeConstraintDistances(unsigned int nCurrentConstraints)
{
    /* If we are far from satisfying the constraint then we need
       to scale the constraint distance so that it can be "satisfied"
       by simply moving in the correct direction, i.e. moving towards
       satisying the constraint.

       If we are well within the region where an inequality constraint
       is satisfied, then the constraint can be removed from the optimisation
       problem. Here we create a map between indices for the vector of active
       constraints and the original constraints
       vector.
     */


    // Number of active contraints.
    unsigned int nActive = 0;

		// Whether each constraint is active.
    std::vector<bool> isActive(nCurrentConstraints);

    // Loop over all constraints.
    for (unsigned int i=0;i<nCurrentConstraints;i++)
    {
        // Flag constraint as active.
        isActive [i] = true;

        // Min and max constraint changes.
        double min = 0;
        double max = 0;

        // Integrate over boundary points.
        for (unsigned int j=0;j<nPoints;j++)
        {
            // Don't consider fixed points.
            if (!boundaryPoints[j].isFixed)
            {
                if (boundaryPoints[j].sensitivities[indexMap [i+1]] > 0)
                {
                    min += boundaryPoints[j].sensitivities[indexMap [i+1]]
                        * boundaryPoints[j].length
                        * boundaryPoints[j].negativeLimit;

                    max += boundaryPoints[j].sensitivities[indexMap [i+1]]
                        * boundaryPoints[j].length
                        * boundaryPoints[j].positiveLimit;
                }
                else
                {
                    min += boundaryPoints[j].sensitivities[indexMap [i+1]]
                        * boundaryPoints[j].length
                        * boundaryPoints[j].positiveLimit;

                    max += boundaryPoints[j].sensitivities[indexMap [i+1]]
                        * boundaryPoints[j].length
                        * boundaryPoints[j].negativeLimit;
                }
            }
        }

        // Scale (20% is arbitrary, but seems to work well).
        min *= 0.2;
        max *= 0.2;

        // Constraint is violated.
        if (constraintDistances[indexMap [i+1] - 1] < 0)
        {
            if (constraintDistances[indexMap [i+1] - 1] < min)
                constraintDistancesScaled[i] = min;
        }

        // Constraint is satisfied.
        else
        {
            if (constraintDistances[indexMap [i+1] - 1] > max)
            {
                // Flag inequality constraint as inactive.
                if (!isEquality[indexMap [i+1] - 1]) isActive [i] = false;

                else constraintDistancesScaled[i] = max;
            }
        }

    }

		for (unsigned i=0;i<nConstraints;i++)
		{
			// Constraint is active.
        if (isActive [i])
        {
            // Shift initial lambda estimate.
            lambdas[nActive+1] = lambdas[i+1];
            
            // Shift constraint distance.
            constraintDistancesScaled[nActive] = constraintDistancesScaled[i];
            
						// Shift negative lambda limit.
            negativeLambdaLimits[nActive + 1] = negativeLambdaLimits[i + 1];

            // Shift positive lambda limit.
            positiveLambdaLimits[nActive + 1] = positiveLambdaLimits[i + 1];
            
            // Shift equality flag.
            isEquality[nActive] = isEquality[i];
            
            // Map the constraint index: active --> original
            indexMap[nActive+1] = indexMap [i + 1];

            // Incremement the number of active constraints.
            nActive++;
        }
		}

    // Resize vectors if constraints have been removed.
    if (nActive < nConstraints)
    {
        lambdas.resize(nActive + 1);
        negativeLambdaLimits.resize(nActive + 1);
        positiveLambdaLimits.resize(nActive + 1);
        scaleFactors.resize(nActive + 1);
        constraintDistancesScaled.resize(nActive);
        isEquality.resize(nActive);

        // Reduce the number of constraints.
        nConstraints = nActive;
    }
    // If a constraint has been removed, then repeat the process.
        if ((nActive > 0) && (nActive < nCurrentConstraints))
            computeConstraintDistances(nConstraints);
    //  std::cout << isActive [0] << "\t" << isActive [1] << "\n";
}

void Optimise::computeLambdaLimits()
{
    /* The lambda limits are set by computing the minimum displacement
    that violates the CFL condition independently for each function, i.e.
    when setting other lambda values equal to zero.

    In this case the displacement for a given function is simply

        z = lambda x sensitivity x length

    and the largest lambda that doesn't trigger the CFL limit is

        lambda = CFL / (max(abs(sensitivity)) x length)
    */

    // Loop over objective and constraints.
    for (unsigned int i=0;i<nConstraints+1;i++)
    {
        // Remap the sensitivity index: active --> original
        unsigned int k = indexMap[i];

        // Initialise limits.
        negativeLambdaLimits[i] = positiveLambdaLimits[i] = 0;

        // Loop over all boundary points.
        for (unsigned int j=0;j<boundaryPoints.size();j++)
        {
            // Don't consider fixed points.
            if (!boundaryPoints[j].isFixed)
            {
                // Store sensitivity.
                double sens = boundaryPoints[j].sensitivities[k] * boundaryPoints [j].length;

                // Initialise min and max displacements.
                double minDisp, maxDisp;

                // Set limits depending on sign of the sensitivity.
                if (sens > 0)
                {
                    minDisp = boundaryPoints[j].negativeLimit / sens;
                    maxDisp = boundaryPoints[j].positiveLimit / sens;
                }
                else
                {
                    maxDisp = boundaryPoints[j].negativeLimit / sens;
                    minDisp = boundaryPoints[j].positiveLimit / sens;
                }

                // Update limits.

                if (maxDisp > positiveLambdaLimits[i])
                    positiveLambdaLimits[i] = maxDisp;

                if (minDisp < negativeLambdaLimits[i])
                    negativeLambdaLimits[i] = minDisp;
            }
        }

        // Scale limits.
        negativeLambdaLimits[i] /= scaleFactors[i];
        positiveLambdaLimits[i] /= scaleFactors[i];

        /* Rescale the lambda values so that they are in range.

            N.B. Resetting the lambda values to zero can cause spurious
            errors with the optimiser when close to convergence. Solutions
            may be found with lambda = 0 (for the objective) which will
            result in undefined velocities.
        */
        while (lambdas[i] < negativeLambdaLimits[i]) lambdas[i] *= 0.9;
        while (lambdas[i] > positiveLambdaLimits[i]) lambdas[i] *= 0.9;
    }

    /* Check for a zero negative lambda limit for the objective function.
       This can occur for compliance minimisation when the initial domain
       is completely filled.
     */
    if (negativeLambdaLimits[0] == 0) negativeLambdaLimits[0] = -1e-2;
}

void Optimise::computeDisplacements(const std::vector<double>& lambda)
{
    // Loop over all boundary points.
    for (unsigned int i=0;i<nPoints;i++)
    {
        // Don't consider fixed points.
        if (!boundaryPoints[i].isFixed)
        {

            // Initialise component for objective.
            displacements[i] = scaleFactors[0] * lambda[0] * boundaryPoints[i].sensitivities[0] * boundaryPoints [i].length;

            // Add components for active constraints.
            for (unsigned int j=1;j<nConstraints+1;j++)
            {
                // Remap the sensitivity index: active --> original
                unsigned int k = indexMap[j];

                // Update displacement vector.
                displacements[i] += scaleFactors[j] * lambda[j] * boundaryPoints[i].sensitivities[k] * boundaryPoints [i].length;
            }

            // Check side limits if point lies close to domain boundary.
            if (boundaryPoints[i].isDomain)
            {
                // Apply side limit.
                if (displacements[i] < boundaryPoints[i].negativeLimit)
                {
                    displacements[i] = boundaryPoints[i].negativeLimit;
                }
            }
        }
    }
}

double Optimise::computeFunction(unsigned int index)
{
    // This method assumes that displacements have already been calculated.

    // Remap the sensitivity index: active --> original
    unsigned int j = indexMap[index];

    // Initialise function.
    double func = 0;

    // Integrate function over boundary points.
    for (unsigned int i=0;i<nPoints;i++)
    {
        // Don't consider fixed points.
        if (!boundaryPoints[i].isFixed)
            func += (scaleFactors[index] * displacements[i] * boundaryPoints[i].sensitivities[j] * boundaryPoints[i].length);
    }

    if (index == 0) return func;
    else return (func - (scaleFactors[index] * constraintDistancesScaled[index - 1]));
}

void Optimise::computeGradients(const std::vector<double>& lambda, std::vector<double>& gradient, unsigned int index)
{
    std::fill (gradient.begin (), gradient.end (), 0.0);

    // Calculate the derivative with respect to each lambda. ({c1*ci, c2*ci, ..., ci*ci, ..., cn*ci})

    // Loop over all points.
    for (unsigned int i=0;i<nPoints;i++)
    {
        // Don't consider fixed points.
        if (!boundaryPoints[i].isFixed)
        {
            // Loop over all functions (objective, then constraints).
            for (unsigned int j=0;j<nConstraints+1;j++)
            {
                // Remap the sensitivity index: active --> original
                unsigned int k = indexMap[j];

                // Scale factor.
                double scaleFactor = scaleFactors[index] * scaleFactors[k];

                gradient[k] += (boundaryPoints[i].sensitivities[index] * boundaryPoints [i].length
                                * boundaryPoints[i].sensitivities[k] * boundaryPoints [i].length
                                * scaleFactor);
            }
        }
    }
}

double Optimise::rescaleDisplacements()
{
    // Check for CFL violation and rescale the displacments
    // and lambda values if necessary.

    // Maximum displacement magnitude.
    double maxDisp = 0;

    // Displacment scale factor.
    double scale;

    // Loop over all boundary points.
    for (unsigned int i=0;i<nPoints;i++)
    {
        // Don't consider fixed points.
        if (!boundaryPoints[i].isFixed)
        {
            // Absolute displacement.
            double disp = std::abs(displacements[i]);

            // Displacement exceeds the CFL limit.
            if (disp > maxDisplacement)
            {
                // Check if current maximum is exceeded.
                if (disp > maxDisp)
                {
                    // Store maximum displacement and scaling factor.
                    maxDisp = disp;
                    scale = maxDisplacement / maxDisp;
                }
            }
        }
    }

    // CFL condition is violated, rescale displacements.
    if (maxDisp)
    {
        // Scale lambda values.
        for (unsigned int i=0;i<nConstraints+1;i++)
            lambdas[i] *= scale;

        // Recompute the displacement vector.
        computeDisplacements(lambdas);

        // Return the scale factor.
        return scale;
    }

    return 1.0;
}

void Optimise::queryReturnCode()
{
    /* N.B.
       Despite providing an extensive list of return codes NLopt does not
       report when the optmisation exits because a lower or upper lambda
       limit is hit. In this case the return code is 4. However, it's easy
       to test for this situation by comparing the "optimum" lambda values
       to the limits.
     */

    // Success.
    if (returnCode == 1) std::cout << "[INFO] Success: Generic success return value.\n";
    else if (returnCode == 2) std::cout << "[INFO] Success: stopval was reached.\n";
    else if (returnCode == 3) std::cout << "[INFO] Success: ftol_rel or ftol_abs was reached.\n";
    else if (returnCode == 4) std::cout << "[INFO] Success: xtol_rel or xtol_abs was reached.\n";
    else if (returnCode == 5) std::cout << "[INFO] Success: maxeval was reached.\n";
    else if (returnCode == 6) std::cout << "[INFO] Success: maxtime was reached.\n";

    // Error codes.
    else if (returnCode == -1) std::cout << "[INFO] Failed: Generic failure code.\n";
    else if (returnCode == -2) std::cout << "[INFO] Failed: Invalid arguments.\n";
    else if (returnCode == -3) std::cout << "[INFO] Failed: Ran out of memory.\n";
    else if (returnCode == -4) std::cout << "[INFO] Failed: Halted because roundoff errors limited progress.\n";
    else if (returnCode == -5) std::cout << "[INFO] Failed: Halted because of a forced termination.\n";

    // No optimisation.
    else std::cout << "[INFO] No optimisation has been performed.\n";
}

std::vector<double> Optimise::get_lambdas(){
    return lambdas;
}

void Optimise::set_lambdas(const double x, const double y){
    lambdas[0] = x;
    lambdas[1] = y;
}

}

using namespace boost::python;

class dummyC{};


BOOST_PYTHON_MODULE(Opt_Module){


    {scope optim = class_<dummyC>("optim");
    
    class_<optim::Optimise>("Optimise", init<slsm::LevelSet&, std::vector<slsm::BoundaryPoint>&, std::vector<double>&, 
    std::vector<double>&, std::vector<double>&, double, bool>() )
        .def("solve",&optim::Optimise::solve)
        .def("preprocess",&optim::Optimise::preprocess)
        .def("postprocess",&optim::Optimise::postprocess)
        .def("nlopt_activate",&optim::Optimise::nlopt_activate)
        .add_property("negativeLambdaLimits",&optim::Optimise::negativeLambdaLimits)
        .add_property("positiveLambdaLimits",&optim::Optimise::positiveLambdaLimits)
        .add_property("phi",&optim::Optimise::phi)
        .add_property("dphi_dt", &optim::Optimise::dphi_dt)
        .def("compute_dphi_dt", &optim::Optimise::compute_dphi_dt)
        //.add_property('BoundaryPoints',&optim::Optimise::boundaryPoints)
        .def("callback",&optim::Optimise::callback)
        .def("computeDisplacements",&optim::Optimise::computeDisplacements)
        // .def("get_lambdas",&optim::Optimise::get_lambdas)//,return_value_policy<reference_existing_object>())
        // .def("set_lambdas",&optim::Optimise::set_lambdas)//,return_value_policy<reference_existing_object>())
	// .add_property("displacements",&optim::Optimise::displacements)
	// .add_property("dphi_dt",&optim::Optimise::dphi_dt)
        ;
     }
    
}

