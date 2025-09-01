/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#ifndef nrsolver_h
#define nrsolver_h

#include <set>
#include <vector>

#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "floatarray.h"
#include "linesearch.h"
#include "convergedreason.h"

#include <memory>
#include <map>
#include <cmath>
///@name Input fields for NRSolver
//@{
#define _IFT_NRSolver_Name "nrsolver"
#define _IFT_NRSolver_maxiter "maxiter"
#define _IFT_NRSolver_miniterations "miniter"
#define _IFT_NRSolver_minsteplength "minsteplength"
#define _IFT_NRSolver_manrmsteps "manrmsteps"
#define _IFT_NRSolver_lstype "lstype"
#define _IFT_NRSolver_ddm "ddm"
#define _IFT_NRSolver_ddv "ddv"
#define _IFT_NRSolver_ddfunc "ddltf"
#define _IFT_NRSolver_linesearch "linesearch"
#define _IFT_NRSolver_rtolv "rtolv"
#define _IFT_NRSolver_rtolf "rtolf"
#define _IFT_NRSolver_rtold "rtold"
#define _IFT_NRSolver_calcstiffbeforeres "calcstiffbeforeres"
#define _IFT_NRSolver_constrainedNRalpha "constrainednralpha"
#define _IFT_NRSolver_constrainedNRminiter "constrainednrminiter"
#define _IFT_NRSolver_maxinc "maxinc"
#define _IFT_NRSolver_forceScale "forcescale"
#define _IFT_NRSolver_forceScaleDofs "forcescaledofs"
#define _IFT_NRSolver_solutionDependentExternalForces "soldepextforces"
#define _IFT_NRSolver_broyden "broyden"
//@}

namespace oofem {
class Domain;
class EngngModel;

/**
 * This class implements Newton-Raphson Method, derived from abstract NumericalMethod class
 * for solving non-linear problems.
 * Except the traditional load control, it also provides direct displacement control without
 * requiring BC applied,
 * so that the equation renumbering is not required when combining arc-length and Newton-Raphson
 * solvers within a simulation.
 *
 * The direct displacement control is achieved by adding a large number alpha to the corresponding
 * diagonal member of K and replacing the right-hand side by alpha*prescribed_value.
 * If alpha is very much larger than other stiffness coefficients then this alteration
 * effectively replaces the corresponding equation by
 * alpha*unknown_value = alpha*prescribed_value
 * that is, the required condition, but the whole system remains symmetric and minimal
 * changes are necessary in the computational sequence.
 * The above artifice has been introduced by Payne and Irons.
 */
class OOFEM_EXPORT NRSolver : public SparseNonLinearSystemNM
{
protected:
    enum nrsolver_ModeType { nrsolverModifiedNRM, nrsolverFullNRM, nrsolverAccelNRM };

    int nsmax, minIterations;
    double minStepLength;
    nrsolver_ModeType NR_Mode, NR_OldMode;
    int NR_ModeTick;
    int MANRMSteps;

    /// linear system solver
    std :: unique_ptr< SparseLinearSystemNM >linSolver;
    /// linear system solver ID
    LinSystSolverType solverType;
    /// sparse matrix version, used to control constrains application to stiffness
    SparseMtrx :: SparseMtrxVersionType smConstraintVersion;
    /// number of prescribed displacement
    int numberOfPrescribedDofs;
    /**
     * Flag indicating that some dofs are controlled under displacement control.
     * In parallel mode, numberOfPrescribedDofs is local (related to specific partition)
     * so its nonzero value does not mean that there are no prescribed dofs on
     * other partitions.
     */
    bool prescribedDofsFlag;

    /// Array of pairs identifying prescribed dofs (node, dof)
    IntArray prescribedDofs;
    /// Array of prescribed values
    FloatArray prescribedDofsValues;
    /// Load Time Function of prescribed values
    int prescribedDisplacementTF;
    /// Array of prescribed equations
    IntArray prescribedEqs;
    /// Flag indicating that prescribedEqs were initialized
    bool prescribedEqsInitFlag;
    /// Computed reactions. They are stored in order to print them in printState method.
    FloatArray lastReactions;
    /// Flag indicating whether to use line-search
    bool lsFlag;
    /// Line search solver
    std :: unique_ptr< LineSearchNM >linesearchSolver;
    /// Flag indicating if the stiffness should be evaluated before the residual in the first iteration.
    bool mCalcStiffBeforeRes;
    /// Flag indicating whether to use constrained Newton
    bool constrainedNRFlag;
    /// Scale factor for dX, dX_new = alpha * dX
    double constrainedNRalpha;
    /// Minimum number of iterations before constraint is activated
    int constrainedNRminiter;

    /// Relative unbalanced force tolerance for each group
    FloatArray rtolf;
    /// Relative iterative displacement change tolerance for each group
    FloatArray rtold;

    ///@todo This doesn't check units, it is nonsense and must be corrected / Mikael
    FloatArray forceErrVec;
    FloatArray forceErrVecOld;

    /// Solution dependent external forces - updating then each NR iteration
    bool solutionDependentExternalForcesFlag;



    /// Optional user supplied scale of forces used in convergence check.
    std :: map< int, double >dg_forceScale;

    double maxIncAllowed;


    /// using broyden method
    bool broydenFlag = false;
    
public:
    NRSolver(Domain *d, EngngModel *m);
    virtual ~NRSolver();

    // Overloaded methods:
    ConvergedReason solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                    FloatArray &X, FloatArray &dX, FloatArray &F,
                    const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                    int &nite, TimeStep *) override;
    void printState(FILE *outputStream) override;

    void initializeFrom(InputRecord &ir) override;
    const char *giveClassName() const override { return "NRSolver"; }
    virtual const char *giveInputRecordName() const { return _IFT_NRSolver_Name; }

    void setDomain(Domain *d) override {
        this->domain = d;
        if ( linSolver ) {
            linSolver->setDomain(d);
        }
        if ( linesearchSolver ) {
            linesearchSolver->setDomain(d);
        }
    }
    void reinitialize() override {
        if ( linSolver ) {
            linSolver->reinitialize();
        }
    }

    SparseLinearSystemNM *giveLinearSolver() override;

protected:
    /// Constructs and returns a line search solver.
    LineSearchNM *giveLineSearchSolver();

    /// Initiates prescribed equations
    void initPrescribedEqs();
    void applyConstraintsToStiffness(SparseMtrx &k);
    void applyConstraintsToLoadIncrement(int nite, const SparseMtrx &k, FloatArray &R,
                                         referenceLoadInputModeType rlm, TimeStep *tStep);

    /**
     * Determines whether or not the solution has reached convergence.
     * @return True if solution has converged, otherwise false.
     */
    bool checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs, FloatArray &ddX, FloatArray &X,
                          double RRT, const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange);
};



// Forward declarations are optional if nrsolver.h already includes these.
// class SparseMtrx;
// class SparseLinearSystemNM;
// class FloatArray;

// ---------------------------------------------------------------------
// BroydenUpdater — rank-1 (good) Broyden over a frozen K factorization
// ---------------------------------------------------------------------
class BroydenUpdater {
public:
    BroydenUpdater() = default;

    void reset(SparseMtrx *K_in, SparseLinearSystemNM *solver_in) {
        K = K_in;
        solver = solver_in;
        havePrev = false;
        ddXn.resize(0);
        rn.resize(0);
    }



  
  // dd_out = K^{-1} r_k − K^{-1}Δr + Δd_prev * ( (Δd_prev^T K^{-1} r_k)/(Δd_prev^T K^{-1}Δr) )
  // Uses one new solve: kr = K^{-1} r_k.
  void computeStep(const FloatArray &rhs, FloatArray &dd_out) {
    // kr = K^{-1} rhs
    FloatArray rhs_copy = rhs;   // <- make it non-const for solver API
    FloatArray ikr, ikdr;
    solver->solve(*K, rhs_copy, ikr);
    //  
    if (!havePrev) { dd_out = ikr; return; }

    auto dr = rhs - rn;
    auto dd = ikr - ddXn;
    //
    solver->solve(*K, dr, ikdr);
    ikdr.add(dd);
    //
    const double denom = dd.dotProduct(ikdr);
    const double numer = dd.dotProduct(ikr);
    //
    const double tol  = 1e-14 * sqrt(dd.computeSquaredNorm());

    dd_out = ikr;
    if (std::abs(denom) > tol) {
      dd_out.add(-numer/denom, ikdr);
    }
  }
  
  // Update with actually applied increment and the new residual (after X += ddX)
  void updateHistory(const FloatArray &ddX, const FloatArray &r) {
    ddXn = ddX;
    rn = r;
    havePrev = true;
  }

  void clear() {
    havePrev = false;
    ddXn.resize(0);
    rn.resize(0);
  }

private:
    SparseMtrx            *K       = nullptr; // frozen stiffness (assembled & factorized externally)
    SparseLinearSystemNM  *solver  = nullptr; // linear system interface

    bool        havePrev = false;
    FloatArray  ddXn;  // d_n
    FloatArray  rn;  //   r_n
};

 
} // end namespace oofem

#endif // nrsolver_h

