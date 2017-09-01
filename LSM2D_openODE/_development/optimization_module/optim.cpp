#include <algorithm>
#include <iostream>

#include "Debug.h"
#include "optim.h"

namespace Optimise_m
{
    
// double callbackWrapper(const std::vector<double>& lambda, std::vector<double>& gradient, void* data)
// {
//     NLoptWrapper* wrapperData = reinterpret_cast<NLoptWrapper*>(data);
//     return reinterpret_cast<Optimise*>(wrapperData->callback)->Optmise::callback(lambda, gradient, wrapperData->index);
// }

Optimise::Optimise(std::vector<slsm::BoundaryPoint>& boundarypoints_, std::vector<double>& constraintDistance_,
            std::vector<double>& lambdas_, std::vector<double>& Multipliers_, 
            double max_Displacement_, bool isMax_, const std::vector<bool>& isEquality_ ): 
            boundaryPoints(boundarypoints_), constraintDistances(constraintDistance_),
            lambdas(lambdas_), LagrangeMultipliers(Multipliers_),
            max_Displacement(max_Displacement_), isMax(isMax_), isEquality(isEquality_)

{
    nlopt::algorithm algorithm = nlopt::LD_SLSQP;
    
    nPoints = boundaryPoints.size();
    nConstraints = lambdas.size() - 1;
    nConstraintsInitial = nConstraints;

    // resize to match lambda size
    positiveLambdaLimits.resize(nConstraints+1);
    negativeLambdaLimits.resize(nConstraints+1);
    scaleFactors.resize(nConstraints + 1);

    if (isEquality.size() == 0) isEquality.resize(nConstraints, false);

}

Pre_optimise::Pre_optimise(std::vector<slsm::BoundaryPoint>& boundarypoints_, std::vector<double>& constraintDistance_,
            std::vector<double>& lambdas_, std::vector<double>& Multipliers_, double max_Displacement_,
            bool isMax_): Optimise(boundarypoints_, constraintDistance_, lambdas, Multipliers_, max_Displacement, isMax_){}
        
void Pre_optimise::computeScaleFactors()
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
    std::vector<double> lambda(nConstraints + 1, 0.0); // temporary vars

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

void Pre_optimise::computeConstraintDistances(unsigned int nCurrentConstraints)
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

void Pre_optimise::computeLambdaLimits()
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

void Pre_optimise::computeGradients(const std::vector<double>& lambda, std::vector<double>& gradient, unsigned int index)
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

}