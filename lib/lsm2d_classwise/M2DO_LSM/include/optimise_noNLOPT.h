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

/* changed a bit to make it exposed to cython: HAC*/

#ifndef _OPTIMISE_H
#define _OPTIMISE_H

#include "./common.h"
#include "./debug.h"

/*! \file Optimise.h
    \brief A class for finding the solution for the optimum velocity vector.
 */

// ASSOCIATED DATA TYPES

class Optimise
{
public:
    Optimise(std::vector<BoundaryPoint>&, const std::vector<double>&,
        std::vector<double>&, double&, double maxDisplacement_ = 0.5, bool isMax_ = false, bool isCapped_ = false,
        const std::vector<bool>& isEquality_ = {});

    /// The number of boundary points.
    unsigned int nPoints;

    /// The number of active constraints.
    unsigned int nConstraints;

    /// The number of initial constraints.
    unsigned int nConstraintsInitial;

    /// A reference to a vector of boundary points.
    std::vector<BoundaryPoint>& boundaryPoints;

    /// A reference to a vector of constraint violation distances.
    const std::vector<double>& constraintDistances;

    /// A reference to a vector of optimum lambda values (to be found by solver).
    std::vector<double>& lambdas;

    /// The effective time step.
    double& timeStep;

    /// The maximum displacement.
    double maxDisplacement;

    /// Whether to maximise the objective function.
    bool isMax;

    /// Whether displacements are capped to the CFL condition.
    bool isCapped;

    /// Whether the constraints are equality constraints.
    std::vector<bool> isEquality;

    /// A map between indices for active constraints.
    std::vector<unsigned int> indexMap;

    /// The boundary point displacement vector.
    std::vector<double> displacements;

    /// Whether the displacement side limit is active for each boundary point.
    std::vector<bool> isSideLimit;

    /// Negative lambda limits.
    std::vector<double> negativeLambdaLimits;

    /// Positive lambda limits.
    std::vector<double> positiveLambdaLimits;

    /// Scale factor for each function.
    std::vector<double> scaleFactors;

    /// Scaled constraint distances.
    std::vector<double> constraintDistancesScaled;

    double callback(const std::vector<double>&, std::vector<double>&, unsigned int);
    
    //! Compute scale factors.
    void computeScaleFactors();

    //! Compute constraint distances.
    void computeConstraintDistances();

    //! Compute lambda limits.
    void computeLambdaLimits();

    //! Compute the boundary movement vector.
    /*! \param lambda
            A vector of lambda values (objective, then constraints).
     */
    void computeDisplacements(const std::vector<double>&);

    //! Compute the change in the objective or constraint functions.
    /*! \param index
            Function index, 0 = objective, 1, 2, 3, ... = constraints.
     */
    double computeFunction(unsigned int);
    double computeFunction(std::vector<double>&, unsigned int);
    

    //! Compute the gradient of the objective or constraint functions.
    /*! \param lambda
            A vector of lambda values (objective, then constraints).

        \param gradient
            The gradient of the change in objective or constraint with respect to each lambda.

        \param index
            Function index, 0 = objective, 1, 2, 3, ... = constraints.
     */
    void computeGradients(const std::vector<double>&, std::vector<double>&, unsigned int index);

    //! Rescale displacements and lambda values if the CFL condition is violated.
    /*! \return
            The scale factor.
     */
    double rescaleDisplacements(std::vector<double>& lambdas);

    //! Apply Laplacian smoothing to computed displacements and rescale them to the CFL condition.
    void LaplacianDisplacementsFilter(int);
};

// #include "../src/optimise.cpp"

#endif  /* _OPTIMISE_H */
