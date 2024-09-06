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

#include "nrsolver.h"
#include "verbose.h"
#include "timestep.h"
#include "mathfem.h"
// includes for ddc - not very clean (NumMethod knows what is "node" and "dof")
#include "node.h"
#include "element.h"
#include "generalboundarycondition.h"
#include "dof.h"
#include "function.h"
#include "linesearch.h"
#include "classfactory.h"
#include "exportmodulemanager.h"
#include "engngm.h"
#include "parallelcontext.h"
#include "unknownnumberingscheme.h"


#ifdef __PETSC_MODULE
 #include "petscsolver.h"
 #include "petscsparsemtrx.h"
#endif

#include <cstdio>

namespace oofem {
#define nrsolver_ERROR_NORM_SMALL_NUM 1.e-6
#define NRSOLVER_MAX_REL_ERROR_BOUND 1.e20
#define NRSOLVER_MAX_RESTARTS 4
#define NRSOLVER_RESET_STEP_REDUCE 0.25
#define NRSOLVER_DEFAULT_NRM_TICKS 10

REGISTER_SparseNonLinearSystemNM(NRSolver)

NRSolver :: NRSolver(Domain *d, EngngModel *m) :
    SparseNonLinearSystemNM(d, m), prescribedDofs(), prescribedDofsValues()
{
    //
    // constructor
    //
    nsmax = 60;     // default maximum number of sweeps allowed
    deltaL = 1.0;
    NR_Mode = NR_OldMode = nrsolverModifiedNRM;
    NR_ModeTick = -1; // do not switch to calm_NR_OldMode
    MANRMSteps = 0;
    numberOfPrescribedDofs = 0;
    prescribedDofsFlag = false;
    prescribedEqsInitFlag = false;
    prescribedDisplacementTF = 0;
    lsFlag = 0; // no line-search

    constrainedNRFlag = false;
    this->forceErrVecOld.resize(0);
    this->forceErrVec.resize(0);
    constrainedNRalpha = 0.5; // default

    smConstraintVersion = 0;
    mCalcStiffBeforeRes = true;

    maxIncAllowed = 1.0e20;
}


NRSolver :: ~NRSolver()
{}


void
NRSolver :: initializeFrom(InputRecord &ir)
{
    SparseNonLinearSystemNM :: initializeFrom(ir);

    // Choosing a big "enough" number. (Alternative: Force input of maxinter)
    nsmax = ( int ) 1e8;
    IR_GIVE_OPTIONAL_FIELD(ir, nsmax, _IFT_NRSolver_maxiter);
    if ( nsmax < 0 ) {
        OOFEM_ERROR("nsmax < 0");
    }

    minIterations = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, minIterations, _IFT_NRSolver_miniterations);

    minStepLength = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, _IFT_NRSolver_minsteplength);

    // read if MANRM method is used
    MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, MANRMSteps, _IFT_NRSolver_manrmsteps);
    if ( MANRMSteps > 0 ) {
        NR_Mode = NR_OldMode = nrsolverAccelNRM;
    } else {
        NR_Mode = nrsolverModifiedNRM;
    }

    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_NRSolver_lstype);
    solverType = ( LinSystSolverType ) _val;
    this->giveLinearSolver()->initializeFrom(ir); // make sure that linear solver is initialized from it as well

    // read relative error tolerances of the solver
    // if rtolv provided set to this tolerance both rtolf and rtold
    rtolf.resize(1);
    rtolf.at(1) = 1.e-3; // Default value.
    rtold.resize(1);
    rtold = FloatArray {
        0.0
    };                       // Default off (0.0 or negative values mean that residual is ignored)
    IR_GIVE_OPTIONAL_FIELD(ir, rtolf.at(1), _IFT_NRSolver_rtolv);
    IR_GIVE_OPTIONAL_FIELD(ir, rtold.at(1), _IFT_NRSolver_rtolv);

    // read optional force and displacement tolerances
    IR_GIVE_OPTIONAL_FIELD(ir, rtolf.at(1), _IFT_NRSolver_rtolf);
    IR_GIVE_OPTIONAL_FIELD(ir, rtold.at(1), _IFT_NRSolver_rtold);

    prescribedDofs.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDofs, _IFT_NRSolver_ddm);
    prescribedDofsValues.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDofsValues, _IFT_NRSolver_ddv);
    prescribedDisplacementTF = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDisplacementTF, _IFT_NRSolver_ddfunc);

    numberOfPrescribedDofs = prescribedDofs.giveSize() / 2;
    if ( numberOfPrescribedDofs != prescribedDofsValues.giveSize() ) {
        OOFEM_ERROR("direct displacement mask size mismatch");
    }

    if ( numberOfPrescribedDofs ) {
        prescribedDofsFlag = true;
    } else {
        prescribedDofsFlag = false;
    }

    this->lsFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, lsFlag, _IFT_NRSolver_linesearch);

    ////////
    IR_GIVE_OPTIONAL_FIELD( ir, this->LsType, "linesearchtype" );
    ////////
    

    if ( this->lsFlag ) {
        this->giveLineSearchSolver()->initializeFrom(ir);
    }

    int calcStiffBeforeResFlag = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, calcStiffBeforeResFlag, _IFT_NRSolver_calcstiffbeforeres);
    if ( calcStiffBeforeResFlag == 0 ) {
        mCalcStiffBeforeRes = false;
    }

    solutionDependentExternalForcesFlag = ir.hasField(_IFT_NRSolver_solutionDependentExternalForces);


    this->constrainedNRminiter = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->constrainedNRminiter, _IFT_NRSolver_constrainedNRminiter);
    this->constrainedNRFlag = this->constrainedNRminiter != 0;


    IR_GIVE_OPTIONAL_FIELD(ir, this->maxIncAllowed, _IFT_NRSolver_maxinc);

    dg_forceScale.clear();
    if ( ir.hasField(_IFT_NRSolver_forceScale) ) {
        IntArray dofs;
        FloatArray forces;
        IR_GIVE_FIELD(ir, forces, _IFT_NRSolver_forceScale);
        IR_GIVE_FIELD(ir, dofs, _IFT_NRSolver_forceScaleDofs);
        for ( int i = 0; i < dofs.giveSize(); ++i ) {
            dg_forceScale [ dofs [ i ] ] = forces [ i ];
        }
    }

}


ConvergedReason
NRSolver :: solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                  FloatArray &X, FloatArray &dX, FloatArray &F,
                  const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tStep)
//
// this function solve the problem of the unbalanced equilibrium
// using NR scheme
//
//
{
    // residual, iteration increment of solution, total external force
    FloatArray rhs, ddX, RT;
    double RRT;
    int neq = X.giveSize();
    bool converged, errorOutOfRangeFlag;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );

    if ( engngModel->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("NRSolver: Iteration");
        if ( rtolf.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" ForceError");
        }
        if ( rtold.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" DisplError");
        }
        OOFEM_LOG_INFO("\n----------------------------------------------------------------------------\n");
    }

    l = 1.0;

    ConvergedReason status = CR_UNKNOWN;
    this->giveLinearSolver();

    // compute total load R = R+R0
    RT = R;
    if ( R0 ) {
        RT.add(* R0);
    }

    RRT = parallel_context->localNorm(RT);
    RRT *= RRT;

    ddX.resize(neq);
    ddX.zero();

    // Fetch the matrix before evaluating internal forces.
    // This is intentional, since its a simple way to drastically increase convergence for nonlinear problems.
    // (This old tangent is just used)
    // This improves convergence for many nonlinear problems, but not all. It may actually
    // cause divergence for some nonlinear problems. Therefore a flag is used to determine if
    // the stiffness should be evaluated before the residual (default yes). /ES

    engngModel->updateComponent(tStep, NonLinearLhs, domain);
    if ( this->prescribedDofsFlag ) {
        if ( !prescribedEqsInitFlag ) {
            this->initPrescribedEqs();
        }
        applyConstraintsToStiffness(k);
    }

    nite = 0;

    //////////////////
    // bifurcatiuon modification
    BifurcationInterface *stabSolver = dynamic_cast<BifurcationInterface *>( this->linSolver.get() ); // Check if bifurcation solver is used 
    bool isBifurcationSet = false;
    double alphaStability  = 1e-6; // eigenvector multiplicator
    double alphamax        = 1e8;

    bool LDLTbif                = false; // Only if choelsky bifurcation should be performed, if FALSE, bifurcation using eigenvectors is done
    bool deflationBifurcation   = true; // only if deflation bifurcation if performed
    bool eigenvectorBifurcation = true; // only if eigenvector bifurcation if performed
    bool postBifurcationLineSearch = true;
    std::vector<bool> bifurcTypes( { LDLTbif, deflationBifurcation, eigenvectorBifurcation } );
    bool exactLineSearchBeforeBifurcation = true;
    //this->setExactPreBifurcationLineSearch( exactLineSearchBeforeBifurcation );
    ///// for provisional data extracting 
    //FloatArray detasDil, detasRot;
    //////////////////

    for ( nite = 0; ; ++nite ) {
        // Compute the residual
        engngModel->updateComponent(tStep, InternalRhs, domain);

        rhs.beDifferenceOf(RT, F);

        if ( this->prescribedDofsFlag ) {
            this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tStep);
        }
        
        // convergence check
        converged = this->checkConvergence( RT, F , rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag);
        
        ////////////////////  
        //// linesearch
        //if ( isnan( F.computeNorm() ) ) { // restart step and tunr Linesearch on
        //    X = X - dX;
        //    dX.zero();
        //    this->setExactPreBifurcationLineSearch( true );
        //    engngModel->updateComponent( tStep, InternalRhs, domain );
        //    rhs.beDifferenceOf( RT, F );
        //    if ( this->prescribedDofsFlag )  this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );
        //}


        //if ( nite > 0 && this->doExactPreBifurcationLineSearch) {
        //    if ( !converged ) {
        //        FloatArray direction, X0;
        //        double Eta, alp = 1e-8;
        //        if ( isnan( F.computeNorm() ) ) { // start from skretch
        //            X = X - dX;
        //            dX.zero();
        //        } else {
        //            X  = X - ddX;
        //            dX = dX - ddX;
        //            direction = ddX;
        //        }
        //        this->exactLineSearch( k, X, dX, F, rlm, nite, tStep, rhs, RT, alp, Eta, X0, direction, 0., false );
        //    } else {
        //        this->doExactPreBifurcationLineSearch = false;
        //    }
        //    
        //}

        //////////////////  
        // bifurcatiuon modification
        if ( stabSolver ) this->performBifurcationAnalysis( k, X, dX, F, internalForcesEBENorm, rlm, nite, tStep, rhs, converged, RT, RRT,
            errorOutOfRangeFlag, stabSolver, isBifurcationSet, alphamax, alphaStability, ddX, bifurcTypes, postBifurcationLineSearch );
        
        // for provisional data extracting 
        //if ( false ) {
        //    if ( nite >= 2 && ( stabSolver->getBifurcation() || tStep->giveNumber() >= 2 ) ) {
        //        FloatArray direction;
        //        stabSolver->getEigenVectors().copyColumn( direction, stabSolver->getEigenVectors().giveNumberOfColumns() );

        //        detasRot.push_back( dX.at( 1 ) / direction.at( 1 ) );
        //        detasDil.push_back( -X.at( 2 ) + stabSolver->giveX0Defl().at( 2 ) );
        //        if ( converged || nite >= 20  ) {
        //             detasDil.printYourselfToFile( "detasDil2LS.txt", false );
        //             detasRot.printYourselfToFile( "detasRot2LS.txt", false );
        //        }
        //    }
        //}
        ////////////////// 
        

        if ( errorOutOfRangeFlag ) {
            status = CR_DIVERGED_TOL;
            OOFEM_WARNING("Divergence reached after %d iterations", nite);
            break;
        } else if ( converged && ( nite >= minIterations ) ) {
            status = CR_CONVERGED;
            break;
        } else if ( nite >= nsmax ) {
            OOFEM_LOG_DEBUG("Maximum number of iterations reached\n");
            status = CR_DIVERGED_ITS;
            break;
        }

        if ( nite > 0 || !mCalcStiffBeforeRes ) {
            if ( ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
                engngModel->updateComponent(tStep, NonLinearLhs, domain);
                applyConstraintsToStiffness(k);
            }
        }

        if ( ( nite == 0 ) && ( deltaL < 1.0 ) ) { // deltaL < 1 means no increment applied, only equilibrate current state
            rhs.zero();
            R.zero();
            ddX = rhs;
        } else {
            //            if ( engngModel->giveProblemScale() == macroScale ) {
            //              k.writeToFile("k.txt");
            //            }

            linSolver->solve(k, rhs, ddX);
        }

        //
        // update solution
        //
        //if ( this->lsFlag && ( nite > 0 ) ) { // Why not nite == 0 ?
        if ( this->lsFlag ) { // Why not nite == 0 ?
            // line search
            LineSearchNM :: LS_status LSstatus;
            double eta;
            this->giveLineSearchSolver()->solve(X, ddX, F, R, R0, prescribedEqs, 1.0, eta, LSstatus, tStep, k);
        } else if ( this->constrainedNRFlag && ( nite > this->constrainedNRminiter ) ) {
            ///@todo This doesn't check units, it is nonsense and must be corrected / Mikael
            if ( this->forceErrVec.computeSquaredNorm() > this->forceErrVecOld.computeSquaredNorm() ) {
                OOFEM_LOG_INFO("Constraining increment to be %e times full increment...\n", this->constrainedNRalpha);
                ddX.times(this->constrainedNRalpha);
            }
            //this->giveConstrainedNRSolver()->solve(X, & ddX, this->forceErrVec, this->forceErrVecOld, status, tStep);
        }


        /////////////////////////////////////////

        double maxInc = 0.0;
        for ( double inc : ddX ) {
            if ( fabs(inc) > maxInc ) {
                maxInc = fabs(inc);
            }
        }

        if ( maxInc > maxIncAllowed ) {
            if ( engngModel->giveProblemScale() == macroScale ) {
                printf("Restricting increment. maxInc: %e\n", maxInc);
            }
            ddX.times(maxIncAllowed / maxInc);
        }

        /////////////////////////////////////////


        X.add(ddX);
        dX.add(ddX);

        if ( solutionDependentExternalForcesFlag ) {
            engngModel->updateComponent(tStep, ExternalRhs, domain);
            RT = R;
            if ( R0 ) {
                RT.add(* R0);
            }
        }



        tStep->incrementStateCounter(); // update solution state counter
        tStep->incrementSubStepNumber();

        engngModel->giveExportModuleManager()->doOutput(tStep, true);
    }

    // Modify Load vector to include "quasi reaction"
    if ( R0 ) {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = F.at( prescribedEqs.at(i) ) - R0->at( prescribedEqs.at(i) ) - R.at( prescribedEqs.at(i) );
        }
    } else {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = F.at( prescribedEqs.at(i) ) - R.at( prescribedEqs.at(i) );
        }
    }

    this->lastReactions.resize(numberOfPrescribedDofs);

#ifdef VERBOSE
    if ( numberOfPrescribedDofs ) {
        // print quasi reactions if direct displacement control used
        OOFEM_LOG_INFO("\n");
        OOFEM_LOG_INFO("NRSolver:     Quasi reaction table                                 \n");
        OOFEM_LOG_INFO("NRSolver:     Node            Dof             Displacement    Force\n");
        double reaction;
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            reaction = R.at( prescribedEqs.at(i) );
            if ( R0 ) {
                reaction += R0->at( prescribedEqs.at(i) );
            }
            lastReactions.at(i) = reaction;
            OOFEM_LOG_INFO("NRSolver:     %-15d %-15d %-+15.5e %-+15.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i),
                           X.at( prescribedEqs.at(i) ), reaction);
        }
        OOFEM_LOG_INFO("\n");
    }
#endif

    return status;
}

void NRSolver::performBifurcationAnalysis( SparseMtrx &k, FloatArray &X, FloatArray &dX, FloatArray &F, const FloatArray &internalForcesEBENorm, 
    referenceLoadInputModeType rlm, int &nite, TimeStep *tStep, FloatArray &rhs, bool &converged, FloatArray &RT, double &RRT, 
    bool &errorOutOfRangeFlag, BifurcationInterface *stabSolver, bool &isBifurcationSet, double &alphamax, double &alphaStability, 
    FloatArray &ddX, std::vector<bool> &bifurcTypes, bool &postBifurcationLineSearch )
{
    // Choose bifurcation type
    bool LDLTbif                = bifurcTypes[0]; // Only if choelsky bifurcation should be performed, if FALSE, bifurcation using eigenvectors is done
    bool deflationBifurcation   = bifurcTypes[1]; // only if deflation bifurcation if performed
    bool eigenvectorBifurcation = bifurcTypes[2]; // only if deflation bifurcation if performed
    bool postDeflation          = false;
    bool exactLineSearchBeforeBifurcation = true;
    

    bool performBifurcation = LDLTbif || deflationBifurcation || eigenvectorBifurcation; // Any bifurcation is on, eigenvectors are default

    // For bifurcation using eigenvectors
    double alphaStability0 = alphaStability;
    int maxBifIte = 100, count = 0;

    //if ( exactLineSearchBeforeBifurcation && isnan( F.computeNorm() ) ) { 
    //    X  = X - ddX;
    //    dX = dX - ddX;


    //    FloatArray X0;
    //    double Eta;
    //    this->exactLineSearch( k, X, dX, F, rlm, nite, tStep, rhs, RT, alphaStability, Eta, X0, ddX, 0., false );
    //}

    // bifurcation analysis modifications
    if ( stabSolver && performBifurcation ) {
        if ( deflationBifurcation && isBifurcationSet ) { // If deflation, modify residuum and recheck convergence
            double dx_defl_Norm = ( X - stabSolver->giveX0Defl() ).computeNorm();
            FloatArray rhs_mod  = rhs * ( 1 / ( dx_defl_Norm * dx_defl_Norm ) + 1 );
            converged           = this->checkConvergence( RT, F, rhs_mod, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag ); // Residuum needs to be modified
        } else if ( isnan( F.computeNorm() ) && isBifurcationSet && eigenvectorBifurcation ) { // if eigenvector bif. and not defl. and NaNs appear
            OOFEM_LOG_INFO( "Retarting step due to NaN\n" );
            X -= dX; // go back to start of the step
            dX.zero();
            alphamax = stabSolver->giveAlpha(); // modify amplitude of the perturbation
            stabSolver->setAlpha( alphamax - alphaStability );
            alphaStability = alphaStability / 2;

            engngModel->updateComponent( tStep, InternalRhs, domain );
            rhs.beDifferenceOf( RT, F );
            if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );

            converged = this->checkConvergence( RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag ); // recheck the convergence
        }


        if ( converged ) {
            engngModel->updateComponent( tStep, NonLinearLhs, domain );
            applyConstraintsToStiffness( k );
            bool isPD = stabSolver->checkPD( k );
            if ( !isPD ) { // If matrix is not PD
                OOFEM_LOG_INFO( "Negative eigenvalue\n" );
                stabSolver->setFoundLimitPoint( true );
                stabSolver->setPostBifurcationLineSearchSolver( postBifurcationLineSearch ); // Set the post bifurcation exact linesearch
                this->setExactPreBifurcationLineSearch( false ); // turn off the pre Bifurcational LS

                if ( !isBifurcationSet ) { // If first time not PD
                    stabSolver->setBifurcation( true ); // perform bifurcation anaysis
                    stabSolver->setAlpha( alphaStability0 - alphaStability ); // set eigenvector multiplicator to starting value
                    isBifurcationSet = true; // Not sure if correct
                }

                if ( deflationBifurcation || eigenvectorBifurcation ) {
                    converged = false;

                    FloatArray X0 = X; // Store the instable root
                    stabSolver->setX0Defl( X0 ); // set the instable root to the solver

                    // Find out if eigenvalues should be recomputed
                    if ( stabSolver->getXeigs().giveSize() == 0 || ( X - stabSolver->getXeigs() ).computeNorm() > 1e-8 ) { // if not the same, recompute the eigenvectors
                        engngModel->updateComponent( tStep, NonLinearLhs, domain );
                        applyConstraintsToStiffness( k );
                        stabSolver->setEigenValuesVectors( k, X );
                    }
                    FloatArray direction;// get eigenvectors
                    stabSolver->getEigenVectors().copyColumn( direction, stabSolver->getEigenVectors().giveNumberOfColumns() );

                    if ( !eigenvectorBifurcation ) { // Only deflation (but initial perturbation by eigenvectors)
                        ddX = alphaStability * direction; // perturbation in direction of eigenvector
                        X   = X + ddX; // Perturbe the solution with chosen amplitude
                        dX  = dX + ddX;
                        stabSolver->compute_dx_defl( X ); // compute deflation difference

                        // Update residuum
                        engngModel->updateComponent( tStep, InternalRhs, domain );
                        rhs.beDifferenceOf( RT, F );
                        if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );

                    } else { // Eigenvector bifurcation with/without deflation
                        if ( deflationBifurcation ) { // If combination of deflation and eigenvector bifurcation

                            //if ( false ) this->provisionalOutput( X, F, tStep, RT, direction, X0 ); // PROVISIONAL output residual


                            ddX = direction;
                            this->giveSetPostBifurcationLineSearchSolver( X0, true )->solve( X, ddX, F, RT, prescribedEqs, tStep, k ); // Do linesearch with deflation
                            X += ddX;
                            dX += ddX;
                            engngModel->updateComponent( tStep, InternalRhs, domain );
                            rhs.beDifferenceOf( RT, F );
                            if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );


                            stabSolver->compute_dx_defl( X );
                            if ( !postDeflation ) { // Continue without deflation
                                 isBifurcationSet = false;
                                 stabSolver->setBifurcation( false );
                                 //deflationBifurcation = false;
                                 //eigenvectorBifurcation = false;
                            } 
                        } else if ( count < maxBifIte ) { // if just eigenvectors
                            count++;
                            if ( stabSolver->giveAlpha() + alphaStability >= alphamax ) alphaStability = alphaStability / 2;
                            stabSolver->setAlpha( stabSolver->giveAlpha() + alphaStability ); // increase alpha (eigenvector multiplicator)
                            OOFEM_LOG_INFO( "alpha =  %f\n", stabSolver->giveAlpha() );
                            stabSolver->setBifurcation( true );
                        } else { // if eigenvectors and maximum iterations reached
                            OOFEM_LOG_INFO( "Maximum bifurcation iterations reached\n" );
                        }
                    }
                }
            } else { // is PD,
                // Compute distance from the instable solution
                if ( stabSolver->getBifurcation() ) {
                    FloatArray solDif;
                    solDif.beDifferenceOf( X, stabSolver->giveX0Defl() );
                    double distX0 = solDif.computeNorm();
                }
                stabSolver->setBifurcation( false ); // only needed for cholesky bifurcation
                isBifurcationSet = false; // Not sure if correct

                // Turn the linesearch off after the limit point is found (turn it on only if NAN appears)
                if ( stabSolver->giveFoundLimitPoint() && stabSolver->givePostBifurcationLineSearchSolver() ) stabSolver->setPostBifurcationLineSearchSolver( false );
            }
        } else if ( nite == 0 ) { // not converged, check if the bifurcation analysis should be performed
            stabSolver->setCholesky( LDLTbif ); // This should be done only once
            stabSolver->setDeflation( deflationBifurcation ); // This should be done only once
            isBifurcationSet = stabSolver->getBifurcation();
        } else if ( stabSolver->getBifurcation() && deflationBifurcation ) {
            stabSolver->compute_dx_defl( X );
        //} else if ( stabSolver->giveFoundLimitPoint() && (stabSolver->givePostBifurcationLineSearchSolver( ) ||  isnan( F.computeNorm()))) { // if limit point found and exact LS should be done
        } else if ( stabSolver->giveFoundLimitPoint() && ( stabSolver->givePostBifurcationLineSearchSolver() ) ) { // if limit point found and exact LS should be done
            FloatArray directionLS;
            if ( !stabSolver->givePostBifurcationLineSearchSolver() && isnan( F.computeNorm() ) ) { // If nan, perform linesearch (if enabled)
                stabSolver->setPostBifurcationLineSearchSolver( postBifurcationLineSearch ); // turn ON if nan
                X  = X - dX;
                dX.zero();

                // update components
                engngModel->updateComponent( tStep, InternalRhs, domain );
                rhs.beDifferenceOf( RT, F );
                if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );
                engngModel->updateComponent( tStep, NonLinearLhs, domain );
                applyConstraintsToStiffness( k );

                this->linSolver->solve( k, rhs, ddX ); // Direction must be recomputed
            } else {
                X  = X - ddX;
                dX = dX - ddX;
            }


            this->giveSetPostBifurcationLineSearchSolver( stabSolver->giveX0Defl(), true )->solve( X, ddX, F, RT, prescribedEqs, tStep, k );
            X += ddX;
            dX += ddX;

            engngModel->updateComponent( tStep, InternalRhs, domain );
            rhs.beDifferenceOf( RT, F );
            if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );

            double Eta = directionLS.computeNorm();
            //this->exactLineSearch( k, X, dX, F, rlm, nite, tStep, rhs, RT, alphaStability, Eta, X0, directionLS, 0., true );

            OOFEM_LOG_INFO( "Eta = %.7e\n", Eta );
            
            //if ( Eta == alphaStability ) {
            //    stabSolver->setPostBifurcationLineSearchSolver( false ); // dont line search anymore
            //}            
        }
    }
    bifurcTypes[0] = LDLTbif;
    bifurcTypes[1] = deflationBifurcation; 
    bifurcTypes[2] = eigenvectorBifurcation;
}

void NRSolver::exactLineSearch(SparseMtrx& k, FloatArray& X, FloatArray& dX, FloatArray& F, referenceLoadInputModeType rlm, 
                                int &nite, TimeStep *tStep, FloatArray &rhs, FloatArray &RT, double &alphaStability, double &Eta,
                                FloatArray &X0, FloatArray &directionOrig, double mult, bool deflation )
{
    // Normalize the direction
    FloatArray direction = directionOrig;
    direction.normalize();
    // do exact linesearch including delfation
    int maxIteLS = 100, p = 2;
    double RHS, tolLS = 1e-12, ResNormLS, RHSdefl, deta = 0., alph_defl = 1., dx_norm, gamma, dx_normSq, kdefl;
    FloatArray Ku( X.giveSize() ), dx_Defl, rhsDefl;

    /*deta = alphaStability;*/
    deta           = 1e-6;
    Eta = deta; // start guess
    FloatArray ddX = deta * direction;
    dX  = dX + ddX;
    X   = X + ddX; // Initial moification of the solution
    for ( size_t niteLS = 0; niteLS < maxIteLS; niteLS++ ) { // inner Newton loop to find local minimum in the direction
        engngModel->updateComponent( tStep, InternalRhs, domain );
        rhs.beDifferenceOf( RT, F );
        //if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );

        RHS = rhs.dotProduct( direction );
        RHSdefl = RHS;
        if ( deflation ) { // Compute deflated residuum
            dx_Defl          = X - X0;
            dx_normSq = dx_Defl.dotProduct( dx_Defl );
            RHSdefl *= 1. / dx_normSq + 1.; 
        }  

        // Check convergence
        OOFEM_LOG_INFO( "|resLS| = %.3e\n", abs( RHSdefl ) );
        if ( abs( RHSdefl ) < tolLS ) break;

        engngModel->updateComponent( tStep, NonLinearLhs, domain );
        //applyConstraintsToStiffness( k );
        k.times( direction, Ku );
        kdefl = Ku.dotProduct( direction );

        if ( deflation ) {
            gamma = p / ( dx_normSq + alph_defl * pow( dx_normSq, p ) );
            kdefl +=  gamma * RHS * dx_Defl.dotProduct( direction );
        }
        
        deta = RHS / kdefl;
        Eta += deta;
        ddX = deta * direction;
        X   = X + ddX;
        dX  = dX + ddX;
    } 

    //// In case eta should be modified
    //if ( abs(mult) > 1e-10 ) {
    //    ddX = mult* Eta * direction;
    //    dX  = dX + ddX;
    //    X   = X + ddX; // Increase the step to get out of the valley
    //    // Update residuum
    //    engngModel->updateComponent( tStep, InternalRhs, domain );
    //    rhs.beDifferenceOf( RT, F );
    //    if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );
    //}
 }



void NRSolver::provisionalOutput( FloatArray &X, FloatArray &F, TimeStep *tStep, FloatArray &RT, FloatArray &direction,FloatArray &X0 )
 {
    int Np = 301;
    // double etaRotMin = 0., etaRotMax = 0.0008, etaRotj; // for exact limit point
    // double etaDilMin = 0., etaDilMax = 1e-7, etaDilj;
    // double etaRotMin = 0., etaRotMax = 0.19, etaRotj;
    // double etaDilMin = 0., etaDilMax = 0.015, etaDilj;
    //  double etaRotMin = 0., etaRotMax = 0.5, etaRotj; // for load 10
    // double etaDilMin = -0.01, etaDilMax = 0.03, etaDilj;
    // double etaRotMin = 0.02, etaRotMax = 0.04, etaRotj; // for load 10 detail
    // double etaDilMin = -0.001, etaDilMax = 0.001, etaDilj;
    double etaRotMin = 0., etaRotMax = 4., etaRotj; // for load 10 large
    double etaDilMin = 0., etaDilMax = 0.2, etaDilj;
    // double etaRotMin = 0., etaRotMax = 0.25, etaRotj;
    // double etaDilMin = 0., etaDilMax = 0.07, etaDilj;
    FloatArray dilDir( { 0., -1., 0. } ), rhs;
    FloatMatrix dEn( Np, Np ), etasRot( Np, Np ), etasDil( Np, Np ), rhs1( Np, Np ), rhs2( Np, Np ), RHSdefl;
    for ( size_t jj = 0; jj < Np; jj++ ) {
        etaRotj = etaRotMin + jj * ( etaRotMax - etaRotMin ) / ( Np - 1 );
        for ( size_t kk = 0; kk < Np; kk++ ) {
            etaDilj = etaDilMin + kk * ( etaDilMax - etaDilMin ) / ( Np - 1 );

            etasDil.at( jj + 1, kk + 1 ) = etaDilj;
            etasRot.at( jj + 1, kk + 1 ) = etaRotj;

            X = X0 + direction * etaRotj + dilDir * etaDilj;
            engngModel->updateComponent( tStep, InternalRhs, domain );
            rhs.beDifferenceOf( RT, F );
            // rhs = rhs * ( 1. / ( X - X0 ).dotProduct( X - X0 ) + 1. );
            rhs1.at( jj + 1, kk + 1 ) = rhs.at( 1 );
            rhs2.at( jj + 1, kk + 1 ) = rhs.at( 2 );
        }
    }
    // etasRot.printYourselfToFile( "etasRot2.txt", false );
    // etasDil.printYourselfToFile( "etasDil2.txt", false );
    // rhs1.printYourselfToFile( "rhs1m2.txt", false );
    // rhs2.printYourselfToFile( "rhs2m2.txt", false );
    etasRot.printYourselfToFile( "etasRot2lar.txt", false );
    etasDil.printYourselfToFile( "etasDil2lar.txt", false );
    rhs1.printYourselfToFile( "rhs1m2lar.txt", false );
    rhs2.printYourselfToFile( "rhs2m2lar.txt", false ); 

 }




SparseLinearSystemNM *
NRSolver :: giveLinearSolver()
{
    if ( linSolver ) {
        if ( linSolver->giveLinSystSolverType() == solverType ) {
            return linSolver.get();
        }
    }

    linSolver = classFactory.createSparseLinSolver(solverType, domain, engngModel);
    if ( !linSolver ) {
        OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
    }

    return linSolver.get();
}


LineSearchNM *
NRSolver :: giveLineSearchSolver()
{
    if ( !linesearchSolver ) {
        switch ( this->LsType ) {
        case LST_Default:
            linesearchSolver = std ::make_unique<LineSearchNM>( domain, engngModel );
            break;
        case LST_Exact:
            linesearchSolver = std ::make_unique<ExactLineSearchNM>( domain, engngModel );
            break;
        }
    }

    return linesearchSolver.get();
}

ExactLineSearchNM *NRSolver ::giveSetPostBifurcationLineSearchSolver(  FloatArray & X0, bool doDeflation )
{
    if ( !linesearchSolverPostBifurcation ) {
        linesearchSolverPostBifurcation = std ::make_unique<ExactLineSearchNM>( domain, engngModel );
    }
    this->linesearchSolverPostBifurcation->setDeflation( doDeflation );
    this->linesearchSolverPostBifurcation->setX0defl( X0 );
    return linesearchSolverPostBifurcation.get();
}

void
NRSolver :: initPrescribedEqs()
{
    EModelDefaultEquationNumbering dn;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );
    int count = 0, ndofman = domain->giveNumberOfDofManagers();
    IntArray localPrescribedEqs(numberOfPrescribedDofs);

    for ( int j = 1; j <= ndofman; j++ ) {
        if ( !parallel_context->isLocal( domain->giveNode(j) ) ) {
            continue;
        }
        int jglobnum = domain->giveNode(j)->giveGlobalNumber();
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            int inode = prescribedDofs.at(2 * i - 1);
            int idofid = prescribedDofs.at(2 * i);
            if ( inode == jglobnum ) {
                localPrescribedEqs.at(++count) = domain->giveNode(j)->giveDofWithID(idofid)->giveEquationNumber(dn);
                continue;
            }
        }
    }

    prescribedEqs.resize(count);
    for ( int i = 1; i <= count; i++ ) {
        prescribedEqs.at(i) = localPrescribedEqs.at(i);
    }

    numberOfPrescribedDofs = count;

    this->prescribedEqsInitFlag = true;
}


void
NRSolver :: applyConstraintsToStiffness(SparseMtrx &k)
{
    if ( this->smConstraintVersion == k.giveVersion() ) {
        return;
    }

#ifdef __PETSC_MODULE
    PetscSparseMtrx *lhs = dynamic_cast< PetscSparseMtrx * >(& k);
    if ( lhs ) {
        Vec diag;
        PetscScalar *ptr;
        int eq;

        lhs->createVecGlobal(& diag);
        MatGetDiagonal(* lhs->giveMtrx(), diag);
        VecGetArray(diag, & ptr);
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            eq = prescribedEqs.at(i) - 1;
            MatSetValue(* ( lhs->giveMtrx() ), eq, eq, ptr [ eq ] * 1.e6, INSERT_VALUES);
        }

        MatAssemblyBegin(* lhs->giveMtrx(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(* lhs->giveMtrx(), MAT_FINAL_ASSEMBLY);
        VecRestoreArray(diag, & ptr);
        VecDestroy(& diag);
        if ( numberOfPrescribedDofs ) {
            this->smConstraintVersion = k.giveVersion();
        }

        return;
    }

#endif // __PETSC_MODULE
    for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
        k.at( prescribedEqs.at(i), prescribedEqs.at(i) ) *= 1.e6;
    }

    if ( numberOfPrescribedDofs ) {
        this->smConstraintVersion = k.giveVersion();
    }
}


void
NRSolver :: applyConstraintsToLoadIncrement(int nite, const SparseMtrx &k, FloatArray &R,
                                            referenceLoadInputModeType rlm, TimeStep *tStep)
{
    double factor = engngModel->giveDomain(1)->giveFunction(prescribedDisplacementTF)->evaluateAtTime( tStep->giveTargetTime() );
    if ( ( rlm == rlm_total ) && ( !tStep->isTheFirstStep() ) ) {
        //factor -= engngModel->giveDomain(1)->giveFunction(prescribedDisplacementTF)->
        // at(tStep->givePreviousStep()->giveTime()) ;
        factor -= engngModel->giveDomain(1)->giveFunction(prescribedDisplacementTF)->
                  evaluateAtTime( tStep->giveTargetTime() - tStep->giveTimeIncrement() );
    }

    if ( nite == 0 ) {
#if 0
 #ifdef __PETSC_MODULE
        if ( solverType == ST_Petsc ) {
            //Natural2LocalOrdering* n2lpm = engngModel->giveParallelContext(1)->giveN2Lmap();
            //IntArray* map = n2lpm->giveN2Lmap();
            for ( i = 1; i <= prescribedEqs.giveSize(); i++ ) {
                eq = prescribedEqs.at(i);
                R.at(eq) = prescribedDofsValues.at(i) * factor; // local eq
            }

            return;
        }

 #endif
#else
 #ifdef __PETSC_MODULE
        const PetscSparseMtrx *lhs = dynamic_cast< const PetscSparseMtrx * >(& k);
        if ( lhs ) {
            Vec diag;
            PetscScalar *ptr;
            lhs->createVecGlobal(& diag);
            MatGetDiagonal(* ( const_cast< PetscSparseMtrx * >(lhs)->giveMtrx() ), diag);
            VecGetArray(diag, & ptr);

            for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
                int eq = prescribedEqs.at(i) - 1;
                R.at(eq + 1) = ptr [ eq ] * prescribedDofsValues.at(i) * factor;
            }

            VecRestoreArray(diag, & ptr);
            VecDestroy(& diag);
            return;
        }
 #endif
#endif
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            int eq = prescribedEqs.at(i);
            R.at(eq) = k.at(eq, eq) * prescribedDofsValues.at(i) * factor;
        }
    } else {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = 0.0;
        }
    }
}


void
NRSolver :: printState(FILE *outputStream)
{
#ifdef VERBOSE
    // print quasi reactions if direct displacement control used
    fprintf(outputStream, "\nQuasi reaction table:\n\n");
    fprintf(outputStream, "  node  dof            force\n");
    fprintf(outputStream, "============================\n");
    if ( lastReactions.giveSize() == 0 ) {
        return;
    }

    double reaction;
    for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
        reaction = lastReactions.at(i);
        fprintf(outputStream, "%6d  %3d   %+11.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i), reaction);
    }
    fprintf(outputStream, "============================\n\n");
#endif
}


bool
NRSolver :: checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs,  FloatArray &ddX, FloatArray &X,
                             double RRT, const FloatArray &internalForcesEBENorm,
                             int nite, bool &errorOutOfRange)
{
    double forceErr, dispErr;
    FloatArray dg_forceErr, dg_dispErr, dg_totalLoadLevel, dg_totalDisp;
    bool answer;
    EModelDefaultEquationNumbering dn;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );

    /*
     * The force errors are (if possible) evaluated as relative errors.
     * If the norm of applied load vector is zero (one may load by temperature, etc)
     * then the norm of reaction forces is used in relative norm evaluation.
     *
     * Note: This is done only when all dofs are included (nccdg = 0). Not implemented if
     * multiple convergence criteria are used.
     *
     */

    answer = true;
    errorOutOfRange = false;

    // Store the errors associated with the dof groups
    if ( this->constrainedNRFlag ) {
        this->forceErrVecOld = this->forceErrVec; // copy the old values
        this->forceErrVec.resize( internalForcesEBENorm.giveSize() );
        forceErrVec.zero();
    }

    if ( internalForcesEBENorm.giveSize() > 1 ) { // Special treatment when just one norm is given; No grouping
        int nccdg = this->domain->giveMaxDofID();
        // Keeps tracks of which dof IDs are actually in use;
        IntArray idsInUse(nccdg);
        idsInUse.zero();
        // zero error norms per group
        dg_forceErr.resize(nccdg);
        dg_forceErr.zero();
        dg_dispErr.resize(nccdg);
        dg_dispErr.zero();
        dg_totalLoadLevel.resize(nccdg);
        dg_totalLoadLevel.zero();
        dg_totalDisp.resize(nccdg);
        dg_totalDisp.zero();
        // loop over dof managers
        for ( auto &dofman : domain->giveDofManagers() ) {
            if ( !parallel_context->isLocal( dofman.get() ) ) {
                continue;
            }

            // loop over individual dofs
            for ( Dof *dof : *dofman ) {
                if ( !dof->isPrimaryDof() ) {
                    continue;
                }
                int eq = dof->giveEquationNumber(dn);
                int dofid = dof->giveDofID();
                if ( !eq ) {
                    continue;
                }

                dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                idsInUse.at(dofid)++;
            } // end loop over DOFs
        } // end loop over dof managers

        // loop over elements and their DOFs
        for ( auto &elem : domain->giveElements() ) {
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

            // loop over element internal Dofs
            for ( int idofman = 1; idofman <= elem->giveNumberOfInternalDofManagers(); idofman++ ) {
                DofManager *dofman = elem->giveInternalDofManager(idofman);
                // loop over individual dofs
                for ( Dof *dof : *dofman ) {
                    if ( !dof->isPrimaryDof() ) {
                        continue;
                    }
                    int eq = dof->giveEquationNumber(dn);
                    int dofid = dof->giveDofID();

                    if ( !eq ) {
                        continue;
                    }

                    dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                    dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                    dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                    dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                    idsInUse.at(dofid)++;
                } // end loop over DOFs
            } // end loop over element internal dofmans
        } // end loop over elements

        // loop over boundary conditions and their internal DOFs
        for ( auto &bc : domain->giveBcs() ) {
            // loop over element internal Dofs
            for ( int idofman = 1; idofman <= bc->giveNumberOfInternalDofManagers(); idofman++ ) {
                DofManager *dofman = bc->giveInternalDofManager(idofman);
                // loop over individual dofs
                for ( Dof *dof : *dofman ) {
                    if ( !dof->isPrimaryDof() ) {
                        continue;
                    }
                    int eq = dof->giveEquationNumber(dn);
                    int dofid = dof->giveDofID();

                    if ( !eq ) {
                        continue;
                    }

                    dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                    dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                    dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                    dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                    idsInUse.at(dofid)++;
                } // end loop over DOFs
            } // end loop over element internal dofmans
        } // end loop over elements

        // exchange individual partition contributions (simultaneously for all groups)
        FloatArray collectiveErr(nccdg);
        parallel_context->accumulate(dg_forceErr,       collectiveErr);
        dg_forceErr       = collectiveErr;
        parallel_context->accumulate(dg_dispErr,        collectiveErr);
        dg_dispErr        = collectiveErr;
        parallel_context->accumulate(dg_totalLoadLevel, collectiveErr);
        dg_totalLoadLevel = collectiveErr;
        parallel_context->accumulate(dg_totalDisp,      collectiveErr);
        dg_totalDisp      = collectiveErr;

        if ( engngModel->giveProblemScale() == macroScale ) {
            OOFEM_LOG_INFO("NRSolver: %-5d", nite);
        }

        int maxNumPrintouts = 6;
        int numPrintouts = 0;

        //bool zeroNorm = false;
        // loop over dof groups and check convergence individually
        for ( int dg = 1; dg <= nccdg; dg++ ) {
            bool zeroFNorm = false, zeroDNorm = false;
            // Skips the ones which aren't used in this problem (the residual will be zero for these anyway, but it is annoying to print them all)
            if ( !idsInUse.at(dg) ) {
                continue;
            }

            numPrintouts++;

            if ( engngModel->giveProblemScale() == macroScale && numPrintouts <= maxNumPrintouts ) {
                OOFEM_LOG_INFO( "  %s:", __DofIDItemToString( ( DofIDItem ) dg ).c_str() );
            }

            if ( rtolf.at(1) > 0.0 ) {
                //  compute a relative error norm
                if ( dg_forceScale.find(dg) != dg_forceScale.end() ) {
                    forceErr = sqrt( dg_forceErr.at(dg) / ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) +
                                                            idsInUse.at(dg) * dg_forceScale [ dg ] * dg_forceScale [ dg ] ) );
                } else if ( ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) ) >= nrsolver_ERROR_NORM_SMALL_NUM ) {
                    forceErr = sqrt( dg_forceErr.at(dg) / ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) ) );
                } else {
                    // If both external forces and internal ebe norms are zero, then the residual must be zero.
                    //zeroNorm = true; // Warning about this afterwards.
                    zeroFNorm = true;
                    forceErr = sqrt( dg_forceErr.at(dg) );
                }

                if ( forceErr > rtolf.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                    errorOutOfRange = true;
                }
                if ( forceErr > rtolf.at(1) ) {
                    answer = false;
                }
		if ( (forceErr > 0.0) && (nite == 0)) {
		  // if forceError > 0 and first iteration we report no convergence to actually
		  // apply the loading
		  answer = false;
		}

                if ( engngModel->giveProblemScale() == macroScale  && numPrintouts <= maxNumPrintouts ) {
                    OOFEM_LOG_INFO(zeroFNorm ? " *%.3e" : "  %.3e", forceErr);
                }

                // Store the errors from the current iteration
                if ( this->constrainedNRFlag ) {
                    forceErrVec.at(dg) = forceErr;
                }
            }

            if ( rtold.at(1) > 0.0 ) {
                // compute displacement error
                if ( dg_totalDisp.at(dg) >  nrsolver_ERROR_NORM_SMALL_NUM ) {
                    dispErr = sqrt( dg_dispErr.at(dg) / dg_totalDisp.at(dg) );
                } else {
                    ///@todo This is almost always the case for displacement error. nrsolveR_ERROR_NORM_SMALL_NUM is no good.
                    //zeroNorm = true; // Warning about this afterwards.
                    //zeroDNorm = true;
                    dispErr = sqrt( dg_dispErr.at(dg) );
                }
                if ( dispErr  > rtold.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                    errorOutOfRange = true;
                }
                if ( dispErr > rtold.at(1) ) {
                    answer = false;
                }

                if ( engngModel->giveProblemScale() == macroScale  && numPrintouts <= maxNumPrintouts ) {
                    OOFEM_LOG_INFO(zeroDNorm ? " *%.3e" : "  %.3e", dispErr);
                }
            }
        }


        if ( engngModel->giveProblemScale() == macroScale ) {
            OOFEM_LOG_INFO("\n");
        }

        //if ( zeroNorm ) OOFEM_WARNING("Had to resort to absolute error measure (marked by *)");
    } else { // No dof grouping
        double dXX, dXdX;

        if ( engngModel->giveProblemScale() == macroScale ) {
            OOFEM_LOG_INFO("NRSolver:     %-15d", nite);
        } else {
            //            OOFEM_LOG_INFO("  NRSolver:     %-15d", nite);
        }


        forceErr = parallel_context->localNorm(rhs);
        forceErr *= forceErr;
        dXX = parallel_context->localNorm(X);
        dXX *= dXX;                                       // Note: Solutions are always total global values (natural distribution makes little sense for the solution)
        dXdX = parallel_context->localNorm(ddX);
        dXdX *= dXdX;

        if ( rtolf.at(1) > 0.0 ) {
            // we compute a relative error norm
            if ( ( RRT + internalForcesEBENorm.at(1) ) > nrsolver_ERROR_NORM_SMALL_NUM ) {
                forceErr = sqrt( forceErr / ( RRT + internalForcesEBENorm.at(1) ) );
            } else {
                forceErr = sqrt(forceErr);   // absolute norm as last resort
            }
            if ( fabs(forceErr) > rtolf.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                errorOutOfRange = true;
            }
            if ( fabs(forceErr) > rtolf.at(1) ) {
                answer = false;
            }

            if ( engngModel->giveProblemScale() == macroScale ) {
                OOFEM_LOG_INFO(" %-15e", forceErr);
            }

            if ( this->constrainedNRFlag ) {
                // store the errors from the current iteration for use in the next
                forceErrVec.at(1) = forceErr;
            }
        }

        if ( rtold.at(1) > 0.0 ) {
            // compute displacement error
            // err is relative displacement change
            if ( dXX > nrsolver_ERROR_NORM_SMALL_NUM ) {
                dispErr = sqrt(dXdX / dXX);
            } else {
                dispErr = sqrt(dXdX);
            }
            if ( fabs(dispErr)  > rtold.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                errorOutOfRange = true;
            }
            if ( fabs(dispErr)  > rtold.at(1) ) {
                answer = false;
            }

            if ( engngModel->giveProblemScale() == macroScale ) {
                OOFEM_LOG_INFO(" %-15e", dispErr);
            }
        }

        if ( engngModel->giveProblemScale() == macroScale ) {
            OOFEM_LOG_INFO("\n");
        }
    } // end default case (all dofs contributing)

    return answer;
}
} // end namespace oofem
