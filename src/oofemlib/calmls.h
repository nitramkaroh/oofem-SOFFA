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

#ifndef calmls_h
#define calmls_h

#include <set>
#include <vector>
#include <memory>

#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "convergedreason.h"
#include "floatarray.h"
#include "intarray.h"
#include "dofiditem.h"

#include "bifurcationinterface.h"
#include "eigenmtrx.h"
#include "linesearch.h"

///@name Input fields for CylindricalALM
//@{
#define _IFT_CylindricalALM_Name "calm"
#define _IFT_CylindricalALM_psi "psi"
#define _IFT_CylindricalALM_maxiter "maxiter"
#define _IFT_CylindricalALM_maxrestarts "maxrestarts"
#define _IFT_CylindricalALM_minsteplength "minsteplength"
#define _IFT_CylindricalALM_steplength "steplength"
#define _IFT_CylindricalALM_initialsteplength "initialsteplength"
#define _IFT_CylindricalALM_forcedinitialsteplength "forcedinitialsteplength"
#define _IFT_CylindricalALM_reqiterations "reqiterations"
#define _IFT_CylindricalALM_miniterations "miniter"
#define _IFT_CylindricalALM_manrmsteps "manrmsteps"
#define _IFT_CylindricalALM_hpcmode "hpcmode"
#define _IFT_CylindricalALM_hpc "hpc"
#define _IFT_CylindricalALM_hpcw "hpcw"
#define _IFT_CylindricalALM_lstype "lstype"
#define _IFT_CylindricalALM_linesearch "linesearch"
#define _IFT_CylindricalALM_lsearchtol "lsearchtol"
#define _IFT_CylindricalALM_lsearchamp "lsearchamp"
#define _IFT_CylindricalALM_lsearchmaxeta "lsearchmaxeta"
#define _IFT_CylindricalALM_nccdg "nccdg"
#define _IFT_CylindricalALM_ccdg "ccdg"
#define _IFT_CylindricalALM_rtolv "rtolv"
#define _IFT_CylindricalALM_rtolf "rtolf"
#define _IFT_CylindricalALM_rtold "rtold"
#define _IFT_CylindricalALM_rootselectiontype "rootselectiontype"

/// Type of the line search solver object (LineSearchType: 1 LineSearchNM, 2 and 3 ExactLineSearchNM), read with linesearch 1.
#define _IFT_CylindricalALM_linesearchtype "linesearchtype"

//@}



////////////////////////
/// Input record of MaterialCylindricalALM: the load level scales a material parameter.
#define _IFT_MaterialCylindricalALM_Name "matcalm"
////////////////////////


namespace oofem {
class Domain;
class EngngModel;
class ParallelContext;

#define calm_SMALL_NUM 1.e-20
#define calm_SMALL_ERROR_NUM 1.e-6

/**
 * Implementation of sparse nonlinear solver with indirect control.
 * It uses Cylindrical Arc-Length Method algorithm.
 *
 * DESCRIPTION :
 * Perform solution of non-linear problem in the form
 * Kt * deltaR = g
 * where g is defined as g = g(Lambda,X) = Lambda* R - F(X).
 * During iteration process, we want g became zero vector.
 *
 * This method uses Modified Newton Raphson iteration scheme
 *
 * If we solve non-linear static we can interpret symbols as follows:
 *
 * Kt     - tangential stiffness
 * dX     - increment of displacements
 * g      - vector of unbalanced forces (at the end should be zero one)
 * Lambda - Load level
 * Rt     - Quasi Total Load vector Rt = R0 + Lambda*R
 * X      - total displacement vector
 * F(X)   - Nodal representation of (real) internal forces.
 * Psi    - control parameter (=0. means displacement control)
 * Bergan_k0 - variable used to compute bergan's parameters of current stiffness.
 * calm_NR_Mode - variable controlling the mode of NRM (ModifiedNR, Full NRM (stiffness update after each iteration),
 *              Modified Accelerated NRM (we perform iteration with stiffness matrix updated only after calm_MANRMSteps)
 * calm_NR_OldMode - variable containing the old mode of NRM, which will be restored after
 *              calm_NR_ModeTick iterations.
 * calm_NR_ModeTick - see calm_NR_OldMode.
 * calm_MANRMSteps - if calm_NR_Mode == calm_accelNRM, it specifies, that new updated
 *                 stiffness matrix is assembled after calm_MANRMSteps.
 * calm_Control - variable indicating the ALM control.
 * calm_HPCIndirectDofMask - Mask, telling which dofs are used for HPC.
 * calm_HPCWeights - dofs weights in constrain.
 *
 * Tasks:
 * - solving problem: solveYourselfAt.
 * - returning results (increment of displacement, reached level of loading and so on)
 *
 * Variable description  :
 * - K(N,N) = Stiffness matrix (assumed positive definite)
 * - dX(N) = Iterative increment of displacement
 * - Rt = "quasi" total load vector
 * - R0 = Fixed (Initial) load vector
 * - deltaL = Step length
 * - deltaLambda = Increment of load level
 * - Psi = Determines loading or displacement control
 * - dX = Current total increment
 * - F = Nodal representation of (real) internal forces
 * - DeltaLambda = Current total increment of load level
 * - Lambda = Total load level
 * - RTOL= Convergence tolerance
 * - NSMAX = Maximum number of sweeps allowed
 *
 * OUTPUT : (after call solveYourselfAt)
 * - dX = Reached displacement increment
 * - Lambda = Total load level
 * - nite = Number of iterations required to fulfill balance
 * - status = NM_status with flags set to reached state
 */
class OOFEM_EXPORT CylindricalALM : public SparseNonLinearSystemNM
{
protected:
    /// CALM mode type; determines the calm step length control.
    enum calm_ControlType {
        calm_hpc_off = 0, ///< Full ALM with quadratic constrain and all dofs.
        calm_hpc_on, ///< Full ALM with quadratic constrain, taking into account only selected dofs.
        calml_hpc, ///<  Linearized ALM (only displacements), taking into account only selected dofs with given weight.
    };

    /// Controlling mode of Newton-Raphson Method.
    enum calm_NR_ModeType {
        calm_modifiedNRM, ///< Keeps the old tangent.
        calm_fullNRM, ///< Updates the tangent every iteration.
        calm_accelNRM, ///< Updates the tangent after a few steps.
    };

    typedef std :: set< DofIDItem > __DofIDSet;

    int nsmax;
    int maxRestarts;
    double Psi;
    double deltaL, minStepLength, maxStepLength;
    int solved, numberOfRequiredIterations;
    calm_NR_ModeType calm_NR_Mode, calm_NR_OldMode;
    int calm_NR_ModeTick;
    int calm_MANRMSteps;

    /// Minimum hard number of iteration.s
    int minIterations;

    /// Variables for HyperPlaneControl.
    int calm_hpc_init;
    calm_ControlType calm_Control;
    FloatArray calm_HPCWeights;
    /// Array containing equation numbers of dofs under indirect control.
    IntArray calm_HPCIndirectDofMask;
    /// Input array containing dofmanagers and corresponding dof numbers under indirect control.
    IntArray calm_HPCDmanDofSrcArray;
    /// Input array of dofman weights (for hpcmode 2).
    FloatArray calm_HPCDmanWeightSrcArray;

    /// Linear system solver.
    std :: unique_ptr< SparseLinearSystemNM > linSolver;
    /// linear system solver ID.
    LinSystSolverType solverType;

    /// Line search flag.
    int lsFlag;
    /// Line search tolerance.
    double ls_tolerance;
    /// Line search amplification factor.
    double amplifFactor;
    /// Line search parameters (limits).
    double maxEta, minEta;

    // Support for evaluation of error norms for user defined dof-groups.
    /// Number of convergence criteria dof groups.
    int nccdg;
    /// Convergence criteria dof groups.
    std :: vector< __DofIDSet >ccDofGroups;
    /// Relative unbalanced force tolerance for each group.
    FloatArray rtolf;
    /// Relative iterative displacement change tolerance for each group.
    FloatArray rtold;
    /// Parallel context for computing norms, dot products and such.
    ParallelContext *parallel_context;

    /** Type characterizing how the roots of the quadratic equation 
     *  for the arc-length parameter are selected
     *  0 -- cosinus based criterion
     *  1 -- dot product based criterion
     */
    enum RootSelectionType {
      RST_Cos=0,
      RST_Dot=1,
    };

    /// Root selection type
    RootSelectionType rootselectiontype;
    /// previous increment of dX, needed by root selection type 1
    FloatArray old_dX;

    
    //////
    /// Type of linesearchSolver (input linesearchtype).
    LineSearchType LsType = LST_Default;
    /// Line search solver for postbifurcation analysis (MaterialCylindricalALM declares and uses its own).
    std ::unique_ptr<ExactLineSearchNM> linesearchSolverPostBifurcation;
    /**
     * Line search solver of type LsType, created by giveLineSearchSolver. With linesearch 1, initializeFrom
     * passes it the input record; the iterations of this class use do_lineSearch, those of
     * MaterialCylindricalALM solve_linesearch.
     */
    std ::unique_ptr<LineSearchNM> linesearchSolver;
    //////
    //////

public:
    CylindricalALM(Domain * d, EngngModel * m);
    virtual ~CylindricalALM();

    // Overloaded methods:
    ConvergedReason solve(SparseMtrx &K, FloatArray &R, FloatArray *R0,
                    FloatArray &X, FloatArray &dX, FloatArray &F,
                    const FloatArray &internalForcesEBENorm, double &ReachedLambda, referenceLoadInputModeType rlm,
                    int &nite, TimeStep *) override;
    double giveCurrentStepLength() override { return deltaL; }
    void setStepLength(double s) override { deltaL = s; }
    /// Prints the definiteness of the converged tangent stiffness when the linear solver checks it (BifurcationInterface).
    void printState(FILE *outputStream) override;
    void initializeFrom(InputRecord &ir) override;
    bool referenceLoad() const override { return true; }
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    void setDomain(Domain *d) override {
        this->domain = d;
        if ( linSolver ) {
            linSolver->setDomain(d);
        }
    }
    void reinitialize() override {
        calm_hpc_init = 1;
        if ( linSolver ) {
            linSolver->reinitialize();
        }
    }
    const char *giveClassName() const override { return "CylindricalALM"; }
    virtual const char *giveInputRecordName() const { return _IFT_CylindricalALM_Name; }

    SparseLinearSystemNM *giveLinearSolver() override;

protected:
    void convertHPCMap();

    /**
     * @return If 0 then ok, 1 then failure (restart).
     */
    int computeDeltaLambda(double &deltaLambda, const FloatArray &dX, const FloatArray &deltaXt,
                           const FloatArray &deltaX_, const FloatArray &R, double RR, double eta,
                           double deltaL, double DeltaLambda0, int neq);

    void search(int istep, FloatArray &prod, FloatArray &eta, double amp,
                double maxeta, double mineta, int &status);

    /**
     * Evaluates the convergence criteria: the relative norms of the unbalanced forces (tolerance rtolf)
     * and of the iterative displacement increment (rtold), for all DOFs or per DOF group (nccdg).
     * @param errorOutOfRange Set if an error exceeds its tolerance CALM_MAX_REL_ERROR_BOUND times or is NaN.
     * @return True if both errors are within their tolerances.
     */
    virtual bool checkConvergence(const FloatArray &R, const FloatArray *R0, const FloatArray &F,
                          const FloatArray &X, const FloatArray &ddX,
                          double Lambda, double RR0, double RR, double drProduct,
                          const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange);

    /// Perform line search optimization of step length
    void do_lineSearch(FloatArray &X, const FloatArray &XInitial, const FloatArray &deltaX_, const FloatArray &deltaXt,
                       const FloatArray &dXm1, FloatArray &dX, FloatArray &ddX,
                       const FloatArray &R, const FloatArray *R0, const FloatArray &F,
                       double &DeltaLambda, double &DeltaLambdam1, double &deltaLambda,
                       double &Lambda, double &ReachedLambda, double RR, double &drProduct, TimeStep *tStep);

        /// Returns linesearchSolver, created on the first call with the type LsType.
    LineSearchNM *giveLineSearchSolver();
};


////////////////////////
/**
 * Cylindrical arc-length method in which the load level Lambda scales a material parameter instead of
 * the reference load (input record matcalm, otherwise the parameters of calm).
 *
 * Both external load vectors, the reference load R and the fixed load R0, are applied in full in every
 * step (giveReferenceLoadScale returns 1). Every material with an arc-length parameter
 * (ArcLengthMaterialInterface, e.g. arclengthpar of the surface materials) multiplies that parameter by
 * Lambda (setMaterialLambda). The equilibrium is
 *   R + R0 - F(X, Lambda) = 0,
 * and the reference load of CylindricalALM is replaced by -Fint_lam, where Fint_lam are the internal
 * forces of the controlled parameter alone at its input value (EngngModel::updateArcLengthInternalForces),
 * i.e. dF/dLambda when the internal forces are linear in that parameter. The tangent displacement solves
 * K deltaXt = -Fint_lam and the arc-length constraint reads
 *   ||dX||^2 + Psi^2 DeltaLambda^2 ||Fint_lam||^2 = deltaL^2,
 * over the HPC DOFs with hpcmode 1, linearised with hpcmode 2, as in CylindricalALM.
 *
 * With linesearch 1 the corrector uses the exact line search solve_linesearch. When the linear solver
 * implements BifurcationInterface, the definiteness of the converged tangent is checked and a stable
 * branch is sought past a critical point (performBifurcationAnalysis).
 */
class OOFEM_EXPORT MaterialCylindricalALM : public CylindricalALM
{
protected:
    /// Exact line search used after a critical point, see giveSetPostBifurcationLineSearchSolver.
    std ::unique_ptr<ExactLineSearchNM> linesearchSolverPostBifurcation;
    bool dL_update = false; // used only by a commented-out step length control in solve


public:
    MaterialCylindricalALM( Domain *d, EngngModel *m );
    virtual ~MaterialCylindricalALM() {};

    // Overloaded methods:
    ConvergedReason solve( SparseMtrx &K, FloatArray &R, FloatArray *R0,
        FloatArray &X, FloatArray &dX, FloatArray &F,
        const FloatArray &internalForcesEBENorm, double &ReachedLambda, referenceLoadInputModeType rlm,
        int &nite, TimeStep * ) override;
    //void initializeFrom( InputRecord &ir ) override;

    const char *giveClassName() const override { return "MaterialCylindricalALM"; }
    const char *giveInputRecordName() const { return _IFT_MaterialCylindricalALM_Name; }

    /// Lambda controls the material; both external load vectors remain fixed.
    double giveReferenceLoadScale( double ) const override { return 1.; }

    /**
     * Switches every material with an arc-length parameter between all its parameters (returnAll true)
     * and the controlled parameter alone (false). StaticStructural::updateArcLengthInternalForces does
     * the same switching itself.
     */
    void setCALMmaterial( bool returnAll);
    /// Sets the load level Lambda in every material with an arc-length parameter (it multiplies the controlled parameter).
    void setMaterialLambda( double lambdaSet );

    /**
     * Stability treatment after the convergence check of an iteration; called when the linear solver
     * implements BifurcationInterface.
     * - While the deflation is active, convergence is re-evaluated with the residual deflated at the
     *   unstable solution X0, i.e. multiplied by ||X0||^p / ||X - X0||^p + 1, so that the iterations do
     *   not return to X0.
     * - At a converged state the tangent is assembled and its definiteness checked, with the pressure
     *   DOFs condensed out if there are any (findPressureDofIndex). If it is not positive definite and the
     *   deflation or the eigenvector variant is on, converged is reset to false, the current X is stored
     *   as X0 and the solution is moved along the eigenvector of the stability solver: by alphaStability
     *   times it (deflation only), by an exact line search with deflation of X0 (deflation and eigenvector,
     *   the setting of solve), or the stability solver gets a larger eigenvector amplitude (eigenvector
     *   only). The ordinary line search is switched off until a converged tangent is positive definite again.
     * - After a critical point, an iteration that has not converged redoes its last increment with the
     *   exact line search, if postBifurcationLineSearch was set and the deflation is not active.
     * @param rhs Set to the new residual R + R0 - F when the solution is moved here.
     * @param converged Convergence of the iteration; re-evaluated or reset as described above.
     * @param isBifurcationSet Whether the bifurcation treatment is active.
     * @param alphamax Upper bound of the eigenvector amplitude (eigenvector variant).
     * @param alphaStability Perturbation amplitude (deflation only) or amplitude increment (eigenvector
     * variant); halved when the bound would be reached, and when NaN appears.
     * @param bifurcTypes Switches { LDL^T (Cholesky) update, deflation, eigenvector perturbation }.
     * @param postBifurcationLineSearch Whether the exact line search is used after a critical point.
     */
    void performBifurcationAnalysis( SparseMtrx &k, FloatArray &X, FloatArray &dX, FloatArray &F, const FloatArray &internalForcesEBENorm,
        referenceLoadInputModeType rlm, int &nite, TimeStep *tStep, FloatArray &rhs, bool &converged, const FloatArray &R, const FloatArray &R0,
        bool &errorOutOfRangeFlag, BifurcationInterface *stabSolver, bool &isBifurcationSet, double &alphamax, double &alphaStability,
        FloatArray &ddX, std::vector<bool> &bifurcTypes, bool &postBifurcationLineSearch, double Lambda, double RR0, double RR, double drProduct );

    /**
     * Returns the exact line search used after a critical point, created on the first call with 20
     * iterations, tolerance 1e-10, the deflation switch doDeflation and the exponent p_set; later calls
     * only set X0.
     * @param X0 Unstable solution to deflate; an empty array keeps the previous one.
     */
    ExactLineSearchNM *giveSetPostBifurcationLineSearchSolver( FloatArray &X0, bool doDeflation, double p_set = 2. );

    /**
     * Finds the smallest equation number of the pressure DOFs (P_f) of the internal DOF managers of the
     * elements (mixed u-p elements).
     * @return False if there are no such DOFs.
     */
    bool findPressureDofIndex( int &indMin );
    /// Collects the equation numbers of all pressure DOFs (P_f) of the internal DOF managers of the elements.
    bool findPressureDofIndices( IntArray &pinds );
    //std::unique_ptr<SparseMtrx> getCondensedMatrixPressure( SparseMtrx &k );
    //std::unique_ptr<SparseMtrx> getCondensedMatrixPressure( std::unique_ptr<SparseMtrx> k );
    /// Condenses the pressure DOFs out of k in place (EigenMtrx only; they must be numbered last).
    void getCondensedMatrixPressure( SparseMtrx &k );


    /**
     * Exact line search of the corrector: Newton iterations for the step length eta along the fixed
     * direction d = deltaX_ + deltaLambda deltaXt, which make the residual orthogonal to it,
     * (RT - F(Xinitial + dXm1 + eta d)) . d = 0; they stop when the absolute value of this product
     * (deflated if requested) is below ls_tolerance. deltaLambda is not scaled by eta; Lambda follows it.
     * If the search does not converge in 20 iterations or gives NaN, ddX and deltaLambda get back their
     * values on entry and X = Xinitial + dXm1 + ddX (the internal forces are not recomputed then).
     * @param RT Total external load R + R0.
     * @param deflation If true, the residual is deflated at Xdefl (factor 1 / ||X - Xdefl||^2 + 1).
     */
    void solve_linesearch( FloatArray &X, const FloatArray &Xinitial, FloatArray &ddX, FloatArray &F, const FloatArray &RT, TimeStep *tStep, SparseMtrx &k,
        const double &ReachedLambda, double &Lambda, double &DeltaLambda, const double &DeltaLambdam1, double &deltaLambda,
        const FloatArray &deltaX_, const FloatArray &deltaXt, FloatArray &dX, const FloatArray &dXm1, bool deflation = false, const FloatArray &Xdefl = FloatArray( 0 ) );

    /**
     * Convergence criteria of the material arc length, for all DOFs (DOF groups nccdg are not supported):
     * forceErr = sqrt( ||R + R0 - F||^2 / ||R + R0|| ) <= rtolf and dispErr = ||ddX|| / ||X|| <= rtold
     * (absolute measures when the denominators are below calm_SMALL_ERROR_NUM). R is not scaled by Lambda.
     */
    bool checkConvergence( const FloatArray &R, const FloatArray *R0, const FloatArray &F,
        const FloatArray &X, const FloatArray &ddX,
        double Lambda, double RR0, double RR, double drProduct,
        const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange ) override;



};
////////////////////////

////////////////////////
#define _IFT_CylindricalALMStability_Name "calmstability"
/// Flag: the arc-length increment is the time increment of the solution step (see timeAsArcLength).
#define _IFT_CylindricalALMStability_timeAsArcLength "timeasarclength"
/// <node, nodeDof> pairs removed from the tangent before its definiteness check (see pdConstrainedDmanDofSrcArray).
#define _IFT_CylindricalALMStability_pdConstrainedDofs "pdconstraineddofs"
/// Arc length after which the branch is traced in the opposite direction (see reverseAfterTime).
#define _IFT_CylindricalALMStability_reverseAfterTime "reverseaftertime"
/// Obsolete keyword, accepted and ignored: the predictor always uses the increment of R0.
#define _IFT_CylindricalALMStability_r0Increment "r0increment"

/**
 * Cylindrical arc-length method for tracing equilibrium paths through critical points.
 * In addition to CylindricalALM it
 *  - checks the definiteness of the converged tangent stiffness when the linear solver can
 *    (BifurcationInterface), optionally with some DOFs removed (pdconstraineddofs);
 *  - can take the arc-length increment from the solution step (timeasarclength) and reverse
 *    the traversal after a given arc length (reverseaftertime);
 *  - loads the predictor by the increment of the non-reference load since the last converged
 *    step, R0 - R0_old, instead of by the whole R0, see solve.
 */
class OOFEM_EXPORT CylindricalALMStability : public CylindricalALM
{
protected:
    /**
     * If set, the arc-length increment is identified with the solution step time
     * increment: deltaL = tStep->giveTimeIncrement().
     *
     * The step length then stops being solver-owned state and becomes a quantity
     * handed down by the TimeStepController, so that a single lever (the metastep
     * deltaT) controls how far the step advances. Consequences:
     *  - pseudo-time becomes the cumulated arc length by construction, since
     *    targetTime accumulates deltaT;
     *  - @c prescribedtimes becomes an explicit arc-length schedule and
     *    @c requiredtimes pins exact arc-length stations;
     *  - @c finalt becomes the total arc length to be traced;
     *  - the solver must no longer modify deltaL on its own, so the internal
     *    step-length reduction is replaced by a ConvergenceException (letting
     *    MetaStep::reduceTimeStep own the retry) and the post-convergence
     *    adaptation is delegated to MetaStep::adaptTimeStep.
     *
     * The step length keywords @c steplength / @c minsteplength are then unused;
     * @c dtmax / @c dtmin of the metastep reduction strategy take over.
     * Default is false, which keeps the legacy solver-owned behaviour.
     */
    bool timeAsArcLength = false;

    /**
     * Optional set of DOFs to be removed from the tangent stiffness before its
     * definiteness is evaluated, given as <node, nodeDof> pairs.
     *
     * This decouples the stability test from the path-following control. The arc-length
     * constraint (@c hpc) is a numerical device for traversing the branch, whereas removing
     * rows and columns asks a physical question: is the equilibrium stable when these DOFs
     * are held fixed. Keying the latter on the HPC mask forced arc-length control onto the
     * very DOF whose definiteness was of interest, which is degenerate wherever that DOF
     * has a turning point, because there the projection of the tangent onto it vanishes.
     * With this field the constrained check works under any @c hpcmode, in particular
     * under full ALM (@c hpcmode 0).
     *
     * When left unset the HPC mask is used, so existing input behaves as before.
     */
    IntArray pdConstrainedDmanDofSrcArray;
    /// Equation numbers of pdConstrainedDmanDofSrcArray; built on the first solve.
    IntArray pdConstrainedDofMask;
    /// Nonzero while pdConstrainedDofMask still has to be built.
    int pd_init = 0;

    /// Converts the <node, nodeDof> pairs into equation numbers. Mirrors convertHPCMap.
    void convertPDConstrainedMap();

    /**
     * Returns the DOFs to be removed before the definiteness check, empty if the check
     * should not be constrained. Prefers the dedicated mask and falls back to the HPC
     * mask, the latter only in calm_hpc_on mode where it is an actual constraint set.
     */
    const IntArray &givePDConstrainedDofMask() const;

    /**
     * Arc length after which the branch is traced in the opposite direction. Negative
     * (the default) means never, zero reverses from the very first step.
     *
     * Keyed on arc length rather than on a step count on purpose: the step count is an
     * artefact of the stepping, since a time step reduction subdivides a step, whereas
     * the traversed arc length identifies the state on the branch. Note that with
     * @c timeasarclength off the field keys on the load-time-function parameter instead,
     * because that is then the meaning of time.
     *
     * Listing the same value in @c requiredtimes makes a step boundary land exactly on
     * it, so that the turn happens at the intended arc length instead of overshooting.
     */
    double reverseAfterTime = -1.;

    /**
     * Orientation of the traversal, +1 forward and -1 reversed. Multiplies the memoryless
     * reference deltaXt^T R of the predictor, which is what makes one flag serve both root
     * selection types: RST_Cos consults that reference on every step, so the factor acts
     * persistently, whereas RST_Dot consults it only when seeding and then continues from
     * old_dX, so the factor acts once and the continuation carries the new orientation.
     * Applying it persistently to RST_Dot would negate the history every step and oscillate.
     */
    double pathDirection = 1.;

    /// R0 of the last converged step; the predictor is loaded by R0 - R0_old. Empty before the first step.
    FloatArray R0_old;

public:
    CylindricalALMStability( Domain *d, EngngModel *m );

    void initializeFrom( InputRecord &ir ) override;

    ConvergedReason solve(
        SparseMtrx &K,
        FloatArray &R,
        FloatArray *R0,
        FloatArray &X,
        FloatArray &dX,
        FloatArray &F,
        const FloatArray &internalForcesEBENorm,
        double &ReachedLambda,
        referenceLoadInputModeType rlm,
        int &nite,
        TimeStep *tStep ) override;

    /// Returns true if the arc-length increment is taken from the solution step.
    bool giveTimeAsArcLengthFlag() const { return timeAsArcLength; }

    const char *giveClassName() const override
    {
        return "CylindricalALMStability";
    }

    const char *giveInputRecordName() const override
    {
        return _IFT_CylindricalALMStability_Name;
    }
};

} // end namespace oofem
#endif // calmls_h
