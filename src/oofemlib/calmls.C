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

#include "calmls.h"
#include "verbose.h"
#include "timestep.h"
#include "floatmatrix.h"
#include "datastream.h"
#include "mathfem.h"
#include "element.h"
#include "classfactory.h"
#include "engngm.h"
// includes for HPC - not very clean (NumMethod knows what is "node" and "dof")
#include "node.h"
#include "dof.h"
#include "contextioerr.h"
#include "exportmodulemanager.h"
#include "parallelcontext.h"
#include "unknownnumberingscheme.h"
#include "convergenceexception.h"


#include "sm/Materials/arclengthmaterialinterface.h"
#include "material.h"

#include "sm/Elements/nlstructuralelement.h"
#include "sm/Elements/structural3delement.h"
#include "initialcondition.h"
#include "igaclampboundaryconditionlagrange.h"
#include "linesearch.h"
//#include "linesearch.h"

namespace oofem {
#define CALM_RESET_STEP_REDUCE 0.25
//#define CALM_RESET_STEP_REDUCE 0.1
#define CALM_TANGENT_STIFF_TRESHOLD 0.1
#define CALM_DEFAULT_NRM_TICKS 2
#define CALM_MAX_REL_ERROR_BOUND 1.e10

REGISTER_SparseNonLinearSystemNM(CylindricalALM)

CylindricalALM :: CylindricalALM(Domain *d, EngngModel *m) :
    SparseNonLinearSystemNM(d, m), calm_HPCWeights(), calm_HPCIndirectDofMask(), calm_HPCDmanDofSrcArray(), ccDofGroups(), old_dX()
{
    nsmax  = 60;       // default maximum number of sweeps allowed
    numberOfRequiredIterations = 3;
    //rtol   = 10.E-3  ;   // convergence tolerance
    //Psi    = 0.1;       // displacement control on
    Psi    = 1.0;     // load control on
    solved = 0;
    calm_NR_Mode = calm_NR_OldMode = calm_modifiedNRM;
    calm_NR_ModeTick = -1; // do not swith to calm_NR_OldMode
    calm_MANRMSteps = 0;

    //Bergan_k0 = 0.;    // value used for computing Bergan's parameter
    // of current stiffness.

    deltaL    = -1.0;
    // TangenStiffnessTreshold = 0.10;

    // Variables for Hyper Plane Control
    calm_Control = calm_hpc_off; // HPControl is not default
    // linesearch default off
    lsFlag = 0;

    ls_tolerance = 1e-6;
    amplifFactor = 2.5;
    maxEta = 8.5;
    minEta = 0.2;

    // number of convergence_criteria dof groups, set to 0 (default behavior)
    nccdg = 0;

    minIterations = 0;

    calm_hpc_init = 0;

    // Maximum number of restarts when convergence not reached during maxiter
    maxRestarts = 4;

    parallel_context = engngModel->giveParallelContext( d->giveNumber() );
}


CylindricalALM :: ~CylindricalALM()
{}


ConvergedReason
CylindricalALM :: solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                        FloatArray &X, FloatArray &dX, FloatArray &F,
                        const FloatArray &internalForcesEBENorm, double &ReachedLambda, referenceLoadInputModeType rlm,
                        int &nite, TimeStep *tStep)
{
    FloatArray rhs, deltaXt, deltaX_, dXm1, XInitial;
    FloatArray ddX; // total increment of displacements in iteration
    //double Bergan_k0 = 1.0, bk;
    double XX, RR, RR0, XR, p = 0.0;
    double deltaLambda, Lambda, eta, DeltaLambdam1, DeltaLambda = 0.0;
    double drProduct = 0.0;
    int neq = R.giveSize();
    int irest = 0;
    int HPsize, i, ind;
    double _RR, _XX;
    ConvergedReason status = CR_UNKNOWN;
    bool converged, errorOutOfRangeFlag;
    // print iteration header

    OOFEM_LOG_INFO("CALMLS:       Initial step length: %-15e\n", deltaL);
    if ( nccdg == 0 ) {
        OOFEM_LOG_INFO("CALMLS:       Iteration       LoadLevel       ForceError      DisplError    \n");
        OOFEM_LOG_INFO("----------------------------------------------------------------------------\n");
    } else {
        OOFEM_LOG_INFO("Iter  LoadLevel       ");
        for ( i = 1; i <= nccdg; i++ ) {
            OOFEM_LOG_INFO("ForceError(%02d)  DisplError(%02d)  ", i, i);
        }

        OOFEM_LOG_INFO("\n__________________________________________________________\n");
    }

    //
    // Now smarter method is default:
    // after convergence troubles for (calm_NR_ModeTick subsequent
    // steps newly set calm_NR_Mode will be used. After the calm_NR_OldMode will be restored
    //
    if ( calm_NR_ModeTick == 0 ) {
        calm_NR_Mode = calm_NR_OldMode;
    }

    if ( calm_NR_ModeTick > 0 ) {
        calm_NR_ModeTick--;
    }

    XInitial = X;
    ddX.resize(neq);
    ddX.zero();
    deltaXt.resize(neq);
    deltaXt.zero();

    status = CR_UNKNOWN;
    this->giveLinearSolver();

    // create HPC Map if needed
    if ( calm_hpc_init ) {
        this->convertHPCMap();
        calm_hpc_init = 0;
    }

    if ( R0 ) {
        RR0 = parallel_context->localNorm(* R0);
        RR0 *= RR0;
    } else {
        RR0 = 0.0;
    }

    //
    // A  initial step (predictor)
    //
    //

    //
    // A.2.   We assume positive-definite (0)Kt (tangent stiffness mtrx).
    //
    RR = parallel_context->localNorm(R);
    RR *= RR;

restart:
    //
    // A.1. calculation of (0)VarRt
    //
    dX.zero();
    //engngModel->updateComponent(tStep, InternalRhs, domain); // By not updating this, one obtains the old equilibrated tangent.
    engngModel->updateComponent(tStep, NonLinearLhs, domain);
    linSolver->solve(k, R, deltaXt);

    // If desired by the user, the solution is (slightly) perturbed, so that various symmetries can be broken.
    // This is useful e.g. to trigger localization in a homogeneous material under uniform stress without
    // the need to introduce material imperfections. The problem itself remains symmetric but the iterative
    // solution is brought to a nonsymmetric state and it gets a chance to converge to a nonsymmetric solution.
    // Parameters of the perturbation technique are specified by the user and by default no perturbation is done.
    // Milan Jirasek
    SparseNonLinearSystemNM :: applyPerturbation(& deltaXt);

    if ( calm_Control == calm_hpc_off ) {
        XX = parallel_context->localNorm(deltaXt);
        XX *= XX;
        p = sqrt(XX + Psi * Psi * RR);
    } else if ( calm_Control == calm_hpc_on ) {
        HPsize = calm_HPCIndirectDofMask.giveSize();
        _XX = 0;
        _RR = 0;
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                _XX += deltaXt.at(ind) * deltaXt.at(ind);
                _RR += R.at(ind) * R.at(ind);
            }
        }

        // In case of paralllel analysis:
        FloatArray collected_XXRR;
        parallel_context->accumulate({ _XX, _RR }, collected_XXRR);
        _XX = collected_XXRR(0);
        _RR = collected_XXRR(1);

        p = sqrt(_XX + Psi * Psi * _RR);
    } else if ( calm_Control == calml_hpc ) {
        HPsize = calm_HPCIndirectDofMask.giveSize();
        p = 0.;
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                p += deltaXt.at(ind) * calm_HPCWeights.at(i);
            }
        }

        // In case of paralllel analysis:
        p = parallel_context->accumulate(p);
    }

    Lambda = ReachedLambda;
    if(rootselectiontype == RST_Cos) {
      /* XR is unscaled Bergan's param of current stiffness XR = deltaXt^T k deltaXt
       * this is used to test whether k has negative or positive slope */
      XR = parallel_context->localDotProduct(deltaXt, R);
    } else {
      if(old_dX.giveSize()) {
        XR = parallel_context->localDotProduct(deltaXt, old_dX);
      } else {
	XR = parallel_context->localDotProduct(deltaXt, R);
      }
    }
    DeltaLambda = deltaLambda = sgn(XR) * deltaL / p;
    Lambda += DeltaLambda;
    //
    // A.3.
    //

    rhs = R;
    rhs.times(DeltaLambda);
    if ( R0 ) {
        rhs.add(* R0);
    }
    linSolver->solve(k, rhs, dX);
    X.add(dX);

    nite = 0;

    // update solution state counter
    tStep->incrementStateCounter();
    engngModel->updateComponent(tStep, InternalRhs, domain);

    do {
        nite++;
        tStep->incrementSubStepNumber();

        dXm1 = dX;
        DeltaLambdam1 = DeltaLambda;

        //
        // B  - iteration MNRM is used
        //
        // B.1. is ommited because MNRM is used instead of NRM.
        //
        if ( ( calm_NR_Mode == calm_fullNRM ) || ( ( calm_NR_Mode == calm_accelNRM ) && ( nite % calm_MANRMSteps == 0 ) ) ) {
            //
            // ALM with full NRM
            //
            // we assemble new tangent stiffness and compute new deltaXt
            // internal state of elements is updated by previous calling
            // InternalRhs
            //
            engngModel->updateComponent(tStep, NonLinearLhs, domain);
            //
            // compute deltaXt for i-th iteration
            //
            linSolver->solve(k, R, deltaXt);
        }

        // B.2.
        //

        rhs =  R;
        rhs.times(Lambda);
        if ( R0 ) {
            rhs.add(* R0);
        }

        rhs.subtract(F);
        deltaX_.resize(neq);
        linSolver->solve(k, rhs, deltaX_);
        eta = 1.0;
        //
        // B.3.
        //
        if ( this->computeDeltaLambda(deltaLambda, dX, deltaXt, deltaX_, R, RR, eta, deltaL, DeltaLambda, neq) ) {
            irest++;
            if ( irest <= maxRestarts ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length
                deltaL =  deltaL * CALM_RESET_STEP_REDUCE;
                if ( deltaL < minStepLength ) {
                    deltaL = minStepLength;
                }

                // restore previous total displacement vector
                X = XInitial;
                // reset all changes fro previous equilibrium state
                dX.zero();
                // restore initial stiffness
                engngModel->updateComponent(tStep, NonLinearLhs, domain);

                OOFEM_LOG_INFO("CALMLS:       Iteration Reset ...\n");

                calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = CR_DIVERGED_ITS;
                OOFEM_ERROR("can't continue further");
            }
        }

        //
        // B.4.+ B.5.
        //


        if ( this->lsFlag && ( nite != 1 ) ) {
            //
            //  LINE SEARCH
            //
            this->do_lineSearch(X, XInitial, deltaX_, deltaXt,
                                dXm1, dX, ddX,
                                R, R0, F,
                                DeltaLambda, DeltaLambdam1, deltaLambda, Lambda, ReachedLambda,
                                RR, drProduct, tStep);
        } else { // no line search
            //
            // update solution vectors
            //
            ddX.clear();
            ddX.add(eta * deltaLambda, deltaXt);
            ddX.add(eta, deltaX_);
            dX = dXm1;
            dX.add(ddX);
            X = XInitial;
            X.add(dX);

            drProduct = parallel_context->localNorm(ddX);
            drProduct *= drProduct;

            tStep->incrementStateCounter();     // update solution state counter
            //
            // B.6.
            //
            DeltaLambda = DeltaLambdam1 + deltaLambda;
            Lambda = ReachedLambda + DeltaLambda;

            tStep->incrementStateCounter();     // update solution state counter
            engngModel->updateComponent(tStep, InternalRhs, domain);
        }

        //
        // B.7.
        //
        // convergence check
        //

        converged = this->checkConvergence(R, R0, F, X, ddX, Lambda, RR0, RR, drProduct,
                                           internalForcesEBENorm, nite, errorOutOfRangeFlag);
        if ( ( nite >= nsmax ) || errorOutOfRangeFlag ) {
            irest++;
            if ( irest <= maxRestarts ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length
                deltaL =  deltaL * CALM_RESET_STEP_REDUCE;
                if ( deltaL < minStepLength ) {
                    deltaL = minStepLength;
                }

                // restore previous total displacement vector
                X = XInitial;
                // reset all changes from previous equilibrium state
                engngModel->initStepIncrements();
                dX.zero();
                // restore initial stiffness
                engngModel->updateComponent(tStep, NonLinearLhs, domain);

                OOFEM_LOG_INFO("CALMLS:       Iteration Reset ...\n");

                calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = CR_DIVERGED_ITS;
                OOFEM_WARNING("Convergence not reached after %d iterations", nsmax);
                // exit(1);
                break;
            }
        }

        // output of per iteration data
        engngModel->giveExportModuleManager()->doOutput(tStep, true);
    } while ( !converged || ( nite < minIterations ) );

    //
    // update dofs,nodes,Elemms and print result
    //
#ifdef VERBOSE
    // printf ("\nCALM - step iteration finished") ;
#endif

    // Seems to be completely unused / Mikael
#if 0
    /* compute Bergan's parameter of current stiffness */
    double RDR = parallel_context->localDotProduct(* R, * dX);
    double DR = parallel_context->localNorm(* dX);
    bk = DeltaLambda * RDR / ( DR * DR );

    if ( tStep->giveNumber() == 1 ) {
        // compute Bergan_k0
        Bergan_k0 = bk;
    } else {
        bk = bk / Bergan_k0;
        //if (fabs(bk) < CALM_TANGENT_STIFF_TRESHOLD) status |= NM_KeepTangent;
    }
#endif

    /* we have computed Bergan's  parameter -> you can adjust deltaL according
     * to this value */
    //
    // there has been restart already - set nite to maxiter
    //
    // if (irest > 0) nite = nsmax;

    if ( nite > numberOfRequiredIterations ) {
        deltaL =  deltaL * numberOfRequiredIterations / nite;
    } else {
        deltaL =  deltaL * sqrt( sqrt( ( double ) numberOfRequiredIterations / ( double ) nite ) );
    }

    if ( deltaL > maxStepLength ) {
        deltaL = maxStepLength;
    }

    if ( deltaL < minStepLength ) {
        deltaL = minStepLength;
    }

    OOFEM_LOG_INFO("CALMLS:       Adjusted step length: %-15e\n", deltaL);

    status = CR_CONVERGED;
    solved = 1;
    ReachedLambda = Lambda;
    old_dX = dX;
    return status;
}

bool
CylindricalALM :: checkConvergence(const FloatArray &R, const FloatArray *R0, const FloatArray &F,
                                   const FloatArray &X, const FloatArray &ddX,
                                   double Lambda, double RR0, double RR, double drProduct,
                                   const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange)
{
    /*
     * typedef std::set<DofID> __DofIDSet;
     * std::list<__DofIDSet> __ccDofGroups;
     * int nccdg; // number of Convergence Criteria Dof Groups
     */
    int _ng = nccdg;
    double forceErr, dispErr;
    FloatArray rhs; // residual of momentum balance eq (unbalanced nodal forces)
    FloatArray dg_forceErr(nccdg), dg_dispErr(nccdg), dg_totalLoadLevel(nccdg), dg_totalDisp(nccdg);
    bool answer;
    EModelDefaultEquationNumbering dn;

    answer = true;
    errorOutOfRange = false;

    // compute residual vector
    rhs =  R;
    rhs.times(Lambda);
    if ( R0 ) {
        rhs.add(* R0);
    }

    rhs.subtract(F);

    if ( _ng > 0 ) {
        forceErr = dispErr = 0.0;
        // zero error norms per group
        dg_forceErr.zero();
        dg_dispErr.zero();
        dg_totalLoadLevel.zero();
        dg_totalDisp.zero();
        // loop over dof managers
        for ( auto &dman : domain->giveDofManagers() ) {
            if ( !dman->isLocal() ) {
                continue;
            }

            // loop over individual dofs
            for ( Dof *_idofptr : *dman ) {
                // loop over dof groups
                for ( int _dg = 1; _dg <= _ng; _dg++ ) {
                    // test if dof ID is in active set
                    if ( ccDofGroups.at(_dg - 1).find( _idofptr->giveDofID() ) != ccDofGroups.at(_dg - 1).end() ) {
                        int _eq = _idofptr->giveEquationNumber(dn);

                        if ( _eq ) {
                            continue;
                        }
                        dg_forceErr.at(_dg) += rhs.at(_eq) * rhs.at(_eq);
                        dg_dispErr.at(_dg)  += ddX.at(_eq) * ddX.at(_eq);
                        // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
                        if ( R0 ) {
                            dg_totalLoadLevel.at(_dg) += R0->at(_eq) * R0->at(_eq);
                        }

                        dg_totalLoadLevel.at(_dg) += R.at(_eq) * R.at(_eq) * Lambda * Lambda;
                        dg_totalDisp.at(_dg) += X.at(_eq) * X.at(_eq);
                    }
                } // end loop over dof groups
            } // end loop over DOFs
        } // end loop over dof managers

        // loop over elements and their DOFs
        for ( auto &elem : domain->giveElements() ) {
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

            // loop over element internal Dofs
            for ( int _idofman = 1; _idofman <= elem->giveNumberOfInternalDofManagers(); _idofman++ ) {
                // loop over individual dofs
                for ( Dof *_idofptr : *elem->giveInternalDofManager(_idofman) ) {
                    // loop over dof groups
                    for ( int _dg = 1; _dg <= _ng; _dg++ ) {
                        // test if dof ID is in active set
                        if ( ccDofGroups.at(_dg - 1).find( _idofptr->giveDofID() ) != ccDofGroups.at(_dg - 1).end() ) {
                            int _eq = _idofptr->giveEquationNumber(dn);

                            if ( _eq ) {
                                continue;
                            }
                            dg_forceErr.at(_dg) += rhs.at(_eq) * rhs.at(_eq);
                            dg_dispErr.at(_dg)  += ddX.at(_eq) * ddX.at(_eq);
                            // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
                            if ( R0 ) {
                                dg_totalLoadLevel.at(_dg) += R0->at(_eq) * R0->at(_eq);
                            }

                            dg_totalLoadLevel.at(_dg) += R.at(_eq) * R.at(_eq) * Lambda * Lambda;
                            dg_totalDisp.at(_dg) += X.at(_eq) * X.at(_eq);
                        }
                    } // end loop over dof groups
                } // end loop over DOFs
            } // end loop over internal element dofmans
        } // end loop over elements

        // exchange individual partition contributions (simultaneously for all groups)
        FloatArray collectiveErr;
        parallel_context->accumulate(dg_forceErr, collectiveErr);
        dg_forceErr = collectiveErr;
        parallel_context->accumulate(dg_dispErr, collectiveErr);
        dg_dispErr = collectiveErr;
        parallel_context->accumulate(dg_totalLoadLevel, collectiveErr);
        dg_totalLoadLevel = collectiveErr;
        parallel_context->accumulate(dg_totalDisp, collectiveErr);
        dg_totalDisp = collectiveErr;

        OOFEM_LOG_INFO("CALMLS:       %-15d %-15e ", nite, Lambda);
        // loop over dof groups
        for ( int _dg = 1; _dg <= _ng; _dg++ ) {
            //  compute a relative error norm
            if ( ( dg_totalLoadLevel.at(_dg) ) < calm_SMALL_ERROR_NUM ) {
                dg_forceErr.at(_dg) = sqrt( dg_forceErr.at(_dg) );
            } else {
                dg_forceErr.at(_dg) = sqrt( dg_forceErr.at(_dg) / dg_totalLoadLevel.at(_dg) );
            }

            //
            // compute displacement error
            //
            if ( dg_totalDisp.at(_dg) < calm_SMALL_ERROR_NUM ) {
                dg_dispErr.at(_dg) = sqrt( dg_dispErr.at(_dg) );
            } else {
                dg_dispErr.at(_dg) = sqrt( dg_dispErr.at(_dg) / dg_totalDisp.at(_dg) );
            }

            if ( ( fabs( dg_forceErr.at(_dg) ) > rtolf.at(_dg) * CALM_MAX_REL_ERROR_BOUND ) ||
                 ( fabs( dg_dispErr.at(_dg) )  > rtold.at(_dg) * CALM_MAX_REL_ERROR_BOUND ) || std :: isnan( dg_forceErr.at(_dg) ) || std :: isnan( dg_dispErr.at(_dg) ) ) {
                errorOutOfRange = true;
            }

            if ( ( fabs( dg_forceErr.at(_dg) ) > rtolf.at(_dg) ) || ( fabs( dg_dispErr.at(_dg) ) > rtold.at(_dg) ) ) {
                answer = false;
            }


            OOFEM_LOG_INFO( "%-15e %-15e ", dg_forceErr.at(_dg), dg_dispErr.at(_dg) );
        }

        OOFEM_LOG_INFO("\n");
    } else {
        //
        // _ng==0 (errors computed for all dofs - this is the default)
        //

        //
        // compute force error(s)
        //
        double dXX;
        forceErr = parallel_context->localNorm(rhs);
        forceErr *= forceErr;
        dXX = parallel_context->localNorm(X);
        dXX *= dXX;

        double eNorm = internalForcesEBENorm.sum();
        // we compute a relative error norm
        if ( ( RR0 + RR * Lambda * Lambda ) > calm_SMALL_ERROR_NUM ) {
            forceErr = sqrt( forceErr / ( RR0 + RR * Lambda * Lambda ) );
        } else if ( eNorm > calm_SMALL_ERROR_NUM ) {
            forceErr = sqrt(forceErr / eNorm);
        } else {
            forceErr = sqrt(forceErr);
        }

        //
        // compute displacement error
        //
        if ( dXX < calm_SMALL_ERROR_NUM ) {
            dispErr = drProduct;
        } else {
            dispErr = drProduct / dXX;
            dispErr = sqrt(dispErr);
        }

        if ( ( fabs(forceErr) > rtolf.at(1) * CALM_MAX_REL_ERROR_BOUND ) ||
             ( fabs(dispErr)  > rtold.at(1) * CALM_MAX_REL_ERROR_BOUND ) || ( std :: isnan(forceErr) ) || ( std :: isnan(dispErr) ) ) {
            errorOutOfRange = true;
        }

        if ( ( fabs(forceErr) > rtolf.at(1) ) || ( fabs(dispErr) > rtold.at(1) ) ) {
            answer = false;
        }

        OOFEM_LOG_INFO("CALMLS:       %-15d %-15e %-15e %-15e\n", nite, Lambda, forceErr, dispErr);
    } // end default case (all dofs conributing)

    return answer;
}



void
CylindricalALM :: initializeFrom(InputRecord &ir)
{
    SparseNonLinearSystemNM :: initializeFrom(ir);

    double oldPsi =  Psi; // default from constructor
    double initialStepLength, forcedInitialStepLength;
    int hpcMode;

    IR_GIVE_OPTIONAL_FIELD(ir, Psi, _IFT_CylindricalALM_psi);
    if ( Psi < 0. ) {
        Psi = oldPsi;
    }

    nsmax = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, nsmax, _IFT_CylindricalALM_maxiter);
    //if ( nsmax < 30 ) {
    //    nsmax = 30;
    //}

    IR_GIVE_OPTIONAL_FIELD(ir, maxRestarts, _IFT_CylindricalALM_maxrestarts);

    minStepLength = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, _IFT_CylindricalALM_minsteplength);

    IR_GIVE_FIELD(ir, maxStepLength, _IFT_CylindricalALM_steplength);
    initialStepLength = maxStepLength;
    IR_GIVE_OPTIONAL_FIELD(ir, initialStepLength, _IFT_CylindricalALM_initialsteplength);

    //if (deltaL <= 0.0)  deltaL=maxStepLength;
    // This method (instanciate) is called not only at the beginning but also
    // after restart from engngModel updateAttributes -> and in this case
    // we want to keep restored deltaL)
    if ( ( deltaL <= 0.0 ) || ( deltaL > maxStepLength ) ) {
        deltaL = initialStepLength;
    }

    // using this, one can enforce deltaL from the input file after restart
    forcedInitialStepLength = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, forcedInitialStepLength, _IFT_CylindricalALM_forcedinitialsteplength);
    if ( forcedInitialStepLength > 0. ) {
        deltaL = forcedInitialStepLength;
    }

    numberOfRequiredIterations = 3;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfRequiredIterations, _IFT_CylindricalALM_reqiterations);
    if ( numberOfRequiredIterations < 3 ) {
        numberOfRequiredIterations = 3;
    }

    if ( numberOfRequiredIterations > 1000 ) {
        numberOfRequiredIterations = 1000;
    }

    try {
        IR_GIVE_OPTIONAL_FIELD(ir, minIterations, _IFT_CylindricalALM_miniterations);
    } catch ( InputException & ) {
        if ( minIterations > 3 && minIterations < 1000 ) {
            numberOfRequiredIterations = minIterations;
        }

        if ( nsmax <= minIterations ) {
            nsmax = minIterations + 1;
        }
    }

    // read if MANRM method is used
    calm_MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, calm_MANRMSteps, _IFT_CylindricalALM_manrmsteps);
    if ( calm_MANRMSteps > 0 ) {
        calm_NR_Mode = calm_NR_OldMode = calm_accelNRM;
    } else {
        calm_NR_Mode = calm_modifiedNRM;
    }

    // read if HPC is requested
    hpcMode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, hpcMode, _IFT_CylindricalALM_hpcmode);

    calm_HPCDmanDofSrcArray.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, calm_HPCDmanDofSrcArray, _IFT_CylindricalALM_hpc);

    calm_HPCDmanWeightSrcArray.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, calm_HPCDmanWeightSrcArray, _IFT_CylindricalALM_hpcw);
    // in calm_HPCIndirectDofMask are stored pairs with following meaning:
    // inode idof
    // example HPC 4 1 2 6 1
    // will yield to HPC for node 4 dof 1 and node 6 dof 1
    // calm_HPCIndirectDofMask must be converted to indirect map
    // -> because dof eqs. are not known now, we derefer this to
    // solveYourselfAt() subroutine. The need for converting is indicated by
    // calm_HPControl = hpc_init
    if ( calm_HPCDmanDofSrcArray.giveSize() != 0 ) {
        int nsize;
        if ( hpcMode == 1 ) {
            calm_Control = calm_hpc_on;
        } else if ( hpcMode == 2 ) {
            calm_Control = calml_hpc;
        } else {
            calm_Control = calm_hpc_on; // default is to use hpc_on
        }

        if ( ( calm_HPCDmanDofSrcArray.giveSize() % 2 ) != 0 ) {
            OOFEM_ERROR("HPC Map size must be even number, it contains pairs <node, nodeDof>");
        }

        nsize = calm_HPCDmanDofSrcArray.giveSize() / 2;
        if ( calm_HPCWeights.giveSize() == 0 ) {
            // no weights -> set to 1.0 by default
            calm_HPCWeights.resize(nsize);
            for ( int i = 1; i <= nsize; i++ ) {
                calm_HPCWeights.at(i) = 1.0;
            }
        } else if ( nsize != calm_HPCWeights.giveSize() ) {
            OOFEM_ERROR("HPC map size and weight array size mismatch");
        }

        calm_hpc_init = 1;
    } else {
        if ( hpcMode ) {
            OOFEM_ERROR("HPC Map must be specified");
        }
    }

    int _value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _value, _IFT_CylindricalALM_lstype);
    solverType = ( LinSystSolverType ) _value;

    this->lsFlag = 0; // linesearch default is off
    IR_GIVE_OPTIONAL_FIELD(ir, lsFlag, _IFT_CylindricalALM_linesearch);

    ////////
    int valLS = 1;
    if ( this->lsFlag ) IR_GIVE_OPTIONAL_FIELD( ir, valLS, _IFT_CylindricalALM_linesearchtype );
    this->LsType = (LineSearchType)valLS;
    ////////


    if ( this->lsFlag ) {
        this->giveLineSearchSolver()->initializeFrom( ir );
    }

    IR_GIVE_OPTIONAL_FIELD(ir, ls_tolerance, _IFT_CylindricalALM_lsearchtol);
    //if ( ls_tolerance < 0.6 ) {
    //    ls_tolerance = 0.6;
    //}

    //if ( ls_tolerance > 0.95 ) {
    //    ls_tolerance = 0.95;
    //}

    IR_GIVE_OPTIONAL_FIELD(ir, amplifFactor, _IFT_CylindricalALM_lsearchamp);
    if ( amplifFactor < 1.0 ) {
        amplifFactor = 1.0;
    }

    if ( amplifFactor > 10.0 ) {
        amplifFactor = 10.0;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, maxEta, _IFT_CylindricalALM_lsearchmaxeta);
    if ( maxEta < 1.5 ) {
        maxEta = 1.5;
    }

    if ( maxEta > 15.0 ) {
        maxEta = 15.0;
    }

    /* initialize optional dof groups for convergence criteria evaluation */
    this->nccdg = 0; // default, no dof cc group, all norms evaluated for all dofs
    IR_GIVE_OPTIONAL_FIELD(ir, nccdg, _IFT_CylindricalALM_nccdg);

    if ( nccdg >= 1 ) {
        IntArray _val;
        char name [ 16 ];
        // create an empty set
        __DofIDSet _set;
        // resize dof group vector
        this->ccDofGroups.resize(nccdg, _set);
        for ( int _i = 0; _i < nccdg; _i++ ) {
            sprintf(name, "%s%d", _IFT_CylindricalALM_ccdg, _i + 1);
            // read dof group as int array under ccdg# keyword
            IR_GIVE_FIELD(ir, _val, name);
            // convert aray into set
            for ( int _j = 1; _j <= _val.giveSize(); _j++ ) {
                ccDofGroups.at(_i).insert( ( DofIDItem ) _val.at(_j) );
            }
        }

        // read relative error tolerances of the solver fo each cc
        // if common rtolv provided, set to this tolerace both rtolf and rtold
        IR_GIVE_OPTIONAL_FIELD(ir, rtolf, _IFT_CylindricalALM_rtolv);
        rtold = rtolf;
        // read optional force and displacement tolerances
        IR_GIVE_OPTIONAL_FIELD(ir, rtolf, _IFT_CylindricalALM_rtolf);
        IR_GIVE_OPTIONAL_FIELD(ir, rtold, _IFT_CylindricalALM_rtold);

        if ( ( rtolf.giveSize() != nccdg ) || ( rtold.giveSize() != nccdg ) ) {
            OOFEM_ERROR("Incompatible size of rtolf or rtold params, expected size %d (nccdg)", nccdg);
        }
    } else {
        nccdg = 0;
        double _rtol = 1.e-3; // default tolerance
        rtolf.resize(1);
        rtold.resize(1);
        // read relative error tolerances of the solver
        // if common rtolv provided, set to this tolerace both rtolf and rtold
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, _IFT_CylindricalALM_rtolv);
        rtolf.at(1) = rtold.at(1) = _rtol;
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, _IFT_CylindricalALM_rtolf);
        rtolf.at(1) = _rtol;
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, _IFT_CylindricalALM_rtold);
        rtold.at(1) = _rtol;
    }

    this->giveLinearSolver()->initializeFrom(ir);
    
    rootselectiontype = RST_Cos;
    if ( ( calm_Control == calm_hpc_off ) || ( calm_Control == calm_hpc_on ) ) {
      int rst = 0;
      IR_GIVE_OPTIONAL_FIELD(ir, rst, _IFT_CylindricalALM_rootselectiontype);
      if(rst == 0) {
	rootselectiontype = RST_Cos;
      } else if(rst == 1) {
	rootselectiontype = RST_Dot;
      }	  
    }
}


void CylindricalALM :: convertHPCMap()
//
// converts HPC map from user input map to HPC indirect Map
//
// indirect map:
// size of indirect map is number of controlled DOFs
// on i-th position this map contain equation number of i-th controlled dof
//
// user input map;
// user input map has size 2*controlled DOFs
// and contain pairs (node, nodeDof);
//
// This is used in order to hide equation numbering from user
//
{
    IntArray indirectMap;
    FloatArray weights;
    int size;
    EModelDefaultEquationNumbering dn;

    int count = 0;
    size = calm_HPCDmanDofSrcArray.giveSize() / 2;
    indirectMap.resize(size);
    weights.resize(size);
    for ( auto &dman : domain->giveDofManagers() ) {
        int jglobnum = dman->giveLabel();
        for ( int i = 1; i <= size; i++ ) {
            int inode = calm_HPCDmanDofSrcArray.at(2 * i - 1);
            int idofid = calm_HPCDmanDofSrcArray.at(2 * i);
            if ( inode == jglobnum ) {
                if ( parallel_context->isLocal( dman.get() ) ) {
                    indirectMap.at(++count) = dman->giveDofWithID(idofid)->giveEquationNumber(dn);
                    if ( calm_Control == calml_hpc ) {
                        weights.at(count) = calm_HPCDmanWeightSrcArray.at(i);
                    }
                }

                continue;
            }
        }
    }

    if ( count != size ) {
        OOFEM_WARNING("some dofmans/Dofs in HPCarray not recognized");
    }

    calm_HPCIndirectDofMask.resize(count);
    if ( calm_Control == calml_hpc ) {
        calm_HPCWeights.resize(count);
    }

    for ( int i = 1; i <= count; i++ ) {
        calm_HPCIndirectDofMask.at(i) = indirectMap.at(i);
        calm_HPCWeights.at(i) = weights.at(i);
    }
}


void
CylindricalALM :: saveContext(DataStream &stream, ContextMode mode)
{
    if ( !stream.write(deltaL) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
CylindricalALM :: restoreContext(DataStream &stream, ContextMode mode)
{
    if ( !stream.read(deltaL) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


SparseLinearSystemNM *
CylindricalALM :: giveLinearSolver()
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


void
CylindricalALM :: printState(FILE *outputStream)
{
    // A stability-capable linear solver has recorded the definiteness of the converged
    // tangent stiffness. Nothing is printed when no such solver is used, or when no
    // converged solution has been checked in this step.
    if ( auto *stabSolver = dynamic_cast< BifurcationInterface * >( this->linSolver.get() ) ) {
        stabSolver->printPDStatus(outputStream);
    }
}


int
CylindricalALM :: computeDeltaLambda(double &deltaLambda, const FloatArray &dX, const FloatArray &deltaXt,
                                     const FloatArray &deltaX_, const FloatArray &R, double RR, double eta,
                                     double deltaL, double DeltaLambda0, int neq)
{
    double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4, a5;
    //
    // B.3.
    //
    if ( ( calm_Control == calm_hpc_off ) || ( calm_Control == calm_hpc_on ) ) {
        if ( calm_Control == calm_hpc_off ) {
            // this two lines are necessary if NRM is used
            // (for MNRM they can be computed at startup A1).
            double XX = parallel_context->localNorm(deltaXt);
            XX *= XX;
            double XXt = parallel_context->localDotProduct(dX, deltaXt);
            double X_Xt = parallel_context->localDotProduct(deltaX_, deltaXt);
            double dXdX = parallel_context->localNorm(dX);
            dXdX *= dXdX;
            double dXX_ = parallel_context->localDotProduct(dX, deltaX_);
            double X_X_ = parallel_context->localNorm(deltaX_);
            X_X_ *= X_X_;

            a1 = eta * eta * XX + Psi * Psi * RR;
            a2 = RR * Psi * Psi * DeltaLambda0 * 2.0
                 + 2.0 * eta * XXt
                 + 2.0 * eta * eta * X_Xt;
            a3 = dXdX
                 + 2.0 * eta * dXX_
                 + eta * eta * X_X_
                 + DeltaLambda0 * DeltaLambda0 * RR * Psi * Psi - deltaL * deltaL;
        } else if ( calm_Control == calm_hpc_on ) {
            double _rr = 0.;
            double _RR = 0.;
            double _a2 = 0.;
            double _a3 = 0.;
            for ( int i = 1; i <= calm_HPCIndirectDofMask.giveSize(); i++ ) {
                int ind;
                if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                    double _pr   = ( dX.at(ind) + eta * deltaX_.at(ind) );
                    _rr += deltaXt.at(ind) * deltaXt.at(ind);
                    _RR += R.at(ind) * R.at(ind);
                    _a2 += eta * deltaXt.at(ind)  * _pr;
                    _a3 += _pr * _pr;
                }
            }

            FloatArray col_;
            parallel_context->accumulate({ _rr, _RR, _a2, _a3 }, col_);
            a1 = eta * eta * col_(0) + Psi *Psi *col_(1);
            a2 = col_(1) * Psi * Psi * DeltaLambda0 * 2.0;
            a2 += 2.0 * col_(2);
            a3 = col_(3) - deltaL * deltaL + DeltaLambda0 *DeltaLambda0 *col_(1) * Psi * Psi;
            //printf ("\na1=%e, a2=%e, a3=%e", a1,a2,a3);
        }

        // solution of quadratic eqn.
        double discr = a2 * a2 - 4.0 * a1 * a3;
        if ( discr < 0.0 ) {
            OOFEM_WARNING("discriminant is negative, restarting the step");
            deltaLambda = 0;
            return 0;
            //OOFEM_ERROR("discriminant is negative, solution failed");
        }

        discr = sqrt(discr);
        double lam1 = ( -a2 + discr ) / 2. / a1;
        double lam2 = ( -a2 - discr ) / 2. / a1;

        // select better lam (according to angle between deltar0 and deltar1(2).
        //
        // up to there rewritten
        //
        if ( calm_Control == calm_hpc_off ) {
            double tmp = parallel_context->localNorm(dX);
            tmp *= tmp;
            a4 = eta * parallel_context->localDotProduct(dX, deltaX_) + tmp;
            a5 = eta * parallel_context->localDotProduct(dX, deltaXt);
        } else {
            a4 = 0.;
            a5 = 0.;
            for ( int i = 1; i <= calm_HPCIndirectDofMask.giveSize(); i++ ) {
                int ind;
                if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                    a4 += eta * dX.at(ind) * deltaX_.at(ind) + dX.at(ind) * dX.at(ind);
                    a5 += eta * dX.at(ind) * deltaXt.at(ind);
                }
            }

            // In case of parallel simulations (equiv to no-op on seq sim):
            FloatArray cola;
            parallel_context->accumulate({ a4, a5 }, cola);
            a4 = cola(0);
            a5 = cola(1);
        }
	if(rootselectiontype == RST_Cos) {
	  double cos1 = ( a4 + a5 * lam1 ) / deltaL / deltaL;
	  double cos2 = ( a4 + a5 * lam2 ) / deltaL / deltaL;
	  if ( cos1 > cos2 ) {
            deltaLambda = lam1;
	  } else {
            deltaLambda = lam2;
	  }
	} else {
	    double XX = parallel_context->localNorm(deltaXt);
	    XX *= XX;
	    double XXt = parallel_context->localDotProduct(dX, deltaXt);
            double dXdX = parallel_context->localNorm(dX);
            dXdX *= dXdX;
            double dXX_ = parallel_context->localDotProduct(dX, deltaX_);
            double X_X_ = parallel_context->localNorm(deltaX_);
            X_X_ *= X_X_;
	  
	    double c1 = (dXdX + dXX_ + lam1 * XXt) + (Psi * Psi  * ((DeltaLambda0 + deltaLambda) * (DeltaLambda0 + deltaLambda) + (DeltaLambda0 + deltaLambda) * lam1)) ;
	    double c2 = (dXdX + dXX_ + lam2 * XXt) + (Psi * Psi  * ((DeltaLambda0 + deltaLambda) * (DeltaLambda0 + deltaLambda) + (DeltaLambda0 + deltaLambda) * lam2));
	    
	    if ( c1 > c2 ) {
	      deltaLambda = lam1;
	    } else {
	      deltaLambda = lam2;
	    }
	    
	}

        //printf ("eta=%e, lam1=%e, lam2=%e", eta, lam1, lam2);
    } else if ( calm_Control == calml_hpc ) {
        // linearized control
        double nom = 0., denom = 0.;

        for ( int i = 1; i <= calm_HPCIndirectDofMask.giveSize(); i++ ) {
            int ind;
            if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                nom += ( dX.at(ind) + eta * deltaX_.at(ind) ) * calm_HPCWeights.at(i);
                denom += eta * deltaXt.at(ind) * calm_HPCWeights.at(i);
            }
        }

        // In case of parallel simulations (equiv to no-op on seq sim):
        FloatArray colv;
        parallel_context->accumulate({ nom, denom }, colv);
        nom = colv(0);
        denom = colv(1);

        if ( fabs(denom) < calm_SMALL_NUM ) {
            OOFEM_ERROR("zero denominator in linearized control");
        }

        deltaLambda = ( deltaL - nom ) / denom;
    }

    return 0;
}


void
CylindricalALM :: search(int istep, FloatArray &prod, FloatArray &eta, double amp,
                         double maxetalim, double minetalim, int &ico)
{
    int ineg = 0;
    double etaneg = 1.0;
    double etamax = 0.0;


    // obtain ineg (number of previous line search iteration with negative ratio nearest to origin)
    // as well as max previous step length, etamax

    for ( int i = 1; i <= istep; i++ ) {
        etamax = max( etamax, eta.at(i) );
        if ( prod.at(i) >= 0.0 ) {
            continue;
        }

        if ( eta.at(i) >= etaneg ) {
            continue;
        }

        etaneg = eta.at(i);
        ineg = i;
    }

    if ( ineg ) {
        // allow interpolation
        // first find ipos (position of previous s-l with positive ratio that is
        // closest to ineg (but with smaller s-l)
        int ipos = 1;
        for ( int i = 1; i <= istep; i++ ) {
            if ( prod.at(i) <= 0.0 ) {
                continue;
            }

            if ( eta.at(i) > eta.at(ineg) ) {
                continue;
            }

            if ( eta.at(i) < eta.at(ipos) ) {
                continue;
            }

            ipos = i;
        }

        // interpolate to get step-length
        double etaint = ( prod.at(ineg) * eta.at(ipos) - prod.at(ipos) * eta.at(ineg) ) / ( prod.at(ineg) - prod.at(ipos) );
        // alternatively get eta ensuring reasonable change
        double etaalt = eta.at(ipos) + 0.2 * ( eta.at(ineg) - eta.at(ipos) );
        etaint = max(etaint, etaalt);
        if ( etaint < minetalim ) {
            etaint = minetalim;
            if ( ico == 1 ) {
                ico = 2;
            } else {
                ico = 1;
            }
        }

        eta.at(istep + 1) = etaint;
    } else { // ineq == 0
        // allow extrapolation
        double etamaxstep = amp * etamax;
        // extrapolate between current and previous
        double etaextrap = ( prod.at(istep) * eta.at(istep - 1) - prod.at(istep - 1) * eta.at(istep) ) /
                           ( prod.at(istep) - prod.at(istep - 1) );
        eta.at(istep + 1) = etaextrap;
        // check if in limits
        if ( ( etaextrap <= 0.0 ) || ( etaextrap > etamaxstep ) ) {
            eta.at(istep + 1) = etamaxstep;
        }

        if ( ( eta.at(istep + 1) > maxetalim ) && ( ico == 1 ) ) {
            ico = 2;
            return;
        }

        if ( ( eta.at(istep + 1) > maxetalim ) ) {
            ico = 1;
            eta.at(istep + 1) = maxetalim;
        }
    }
}

void
CylindricalALM :: do_lineSearch(FloatArray &r, const FloatArray &rInitial, const FloatArray &deltaX_, const FloatArray &deltaXt,
                                const FloatArray &dXm1, FloatArray &dX, FloatArray &ddX,
                                const FloatArray &R, const FloatArray *R0, const FloatArray &F,
                                double &DeltaLambda, double &DeltaLambdam1, double &deltaLambda,
                                double &Lambda, double &ReachedLambda, double RR, double &drProduct, TimeStep *tStep)
{
    //
    //  LINE SEARCH
    //

    int neq = r.giveSize();
    int ls_failed, dl_failed = 0;
    int _iter = 0;
    int ico;
    int ls_maxiter = 10;

    DeltaLambda = DeltaLambdam1 + deltaLambda;
    Lambda = ReachedLambda + DeltaLambda;
    double deltaLambdaForEta1 = deltaLambda;

    double d6, d7, d8, d9;

    d6 = parallel_context->localDotProduct(deltaX_, F);
    d7 = parallel_context->localDotProduct(deltaXt, F);
    d8 = -1.0 * parallel_context->localDotProduct(deltaX_, R);
    d9 = -1.0 * parallel_context->localDotProduct(deltaXt, R);

    double e1, e2, d10 = 0.0, d11 = 0.0;
    double s0, si;
    double prevEta, currEta;

    FloatArray eta(ls_maxiter + 1), prod(ls_maxiter + 1);


    if ( R0 ) {
        d10 = -1.0 * parallel_context->localDotProduct(deltaX_, * R0);
        d11 = -1.0 * parallel_context->localDotProduct(deltaXt, * R0);
    }

    // prepare starting product ratios and step lengths
    prod.at(1) = 1.0;
    eta.at(1) = 0.0;
    currEta = eta.at(2) = 1.0;
    // following counter shows how many times the max or min step length has been reached
    ico = 0;

    //
    // begin line search loop
    //
    ls_failed = 1;
    for ( int ils = 2; ils <= ls_maxiter; ils++ ) {
        // update displacements
        drProduct = 0.0; // dotproduct of iterative displacement increment vector

        ddX.clear();
        ddX.add(eta.at(ils) * deltaLambda, deltaXt);
        ddX.add(eta.at(ils), deltaX_);
        dX = dXm1;
        dX.add(ddX);

        drProduct = parallel_context->localNorm(ddX);
        drProduct *= drProduct;

        tStep->incrementStateCounter(); // update solution state counter
        // update internal forces according to new state
        engngModel->updateComponent(tStep, InternalRhs, domain);

        e1 = parallel_context->localDotProduct(deltaX_, F);
        e2 = parallel_context->localDotProduct(deltaXt, F);

        s0 = d6 + deltaLambda * d7 + Lambda * d8 + deltaLambda * Lambda * d9 + d10 + deltaLambda * d11;
        si = e1 + deltaLambda * e2 + Lambda * d8 + deltaLambda * Lambda * d9 + d10 + deltaLambda * d11;
        prod.at(ils) = si / s0;

        //printf ("\ns0=%e, si=%e, prod=%e", s0, si, prod.at(ils));
        if ( s0 >= 0.0 ) {
            //printf ("solve starting inner product uphill, val=%e",s0);
            ls_failed = 3;
            currEta = 1.0;
            break;
        }

        if ( fabs(si / s0) < ls_tolerance ) {
            ls_failed = 0;
            currEta = eta.at(ils);
            break;
        }

        _iter = 0;

        currEta = eta.at(ils);
        //printf ("\n_ite=%d, deltaLambda=%e, eta=%e", _iter, deltaLambda, currEta);
        do { // solve simultaneously the equations for eta and lambda
            _iter++;
            prevEta = currEta;
            s0 = d6 + deltaLambda * d7 + Lambda * d8 + deltaLambda * Lambda * d9 + d10 + deltaLambda * d11;
            si = e1 + deltaLambda * e2 + Lambda * d8 + deltaLambda * Lambda * d9 + d10 + deltaLambda * d11;
            prod.at(ils) = si / s0;

            // call line-search routine to get new estimate of eta.at(ils+1)
            this->search(ils, prod, eta, amplifFactor, maxEta, minEta, ico);
            if ( ico == 2 ) {
                ls_failed = 2;
                break; // exit the loop
            }

            currEta = eta.at(ils + 1);
            // solve for deltaLambda
            dl_failed = this->computeDeltaLambda(deltaLambda, dXm1, deltaXt, deltaX_, R, RR, currEta, deltaL, DeltaLambdam1, neq);
            if ( dl_failed ) {
                eta.at(ils + 1) = 1.0;
                deltaLambda = deltaLambdaForEta1;
                break;
            }

            DeltaLambda = DeltaLambdam1 + deltaLambda;
            Lambda = ReachedLambda + DeltaLambda;
            //printf ("\n_ite=%d, deltaLambda=%e, eta=%e", _iter, deltaLambda, currEta);
        } while ( ( _iter < 10 ) && ( fabs( ( currEta - prevEta ) / prevEta ) > 0.01 ) );

        if ( ( ls_failed > 1 ) || dl_failed ) {
            break;
        }

        //printf ("\ncalm ls...");
        //printf ("eta = %e, err=%d, ", currEta,ls_failed);
        //printf ("dLambda=%e, lerr=%d ", deltaLambda, dl_failed);
    } // end of line search loop

    if ( ls_failed || dl_failed ) {
        // last resort
        deltaLambda = deltaLambdaForEta1;
        drProduct = 0.0; // dotproduct of iterative displacement increment vector

        ddX.clear();
        ddX.add(deltaLambda, deltaXt);
        ddX.add(deltaX_);
        dX = dXm1;
        dX.add(ddX);

        drProduct = parallel_context->localNorm(ddX);
        drProduct *= drProduct;

        tStep->incrementStateCounter(); // update solution state counter
        engngModel->updateComponent(tStep, InternalRhs, domain);
        DeltaLambda = DeltaLambdam1 + deltaLambda;
        Lambda = ReachedLambda + DeltaLambda;

        OOFEM_LOG_INFO("CALMLS:       Line search - err_id=%d, eta=%e, dlambda=%e\n", ls_failed, 1.0, deltaLambda);
    } else {
        OOFEM_LOG_INFO("CALMLS:       Line search - err_id=%d, eta=%e, dlambda=%e\n", ls_failed, currEta, deltaLambda);
    }
}


LineSearchNM *
CylindricalALM ::giveLineSearchSolver()
{
    if ( !linesearchSolver ) {
        switch ( this->LsType ) {
        case LST_Default:
            linesearchSolver = std ::make_unique<LineSearchNM>( domain, engngModel );
            break;
        case LST_Exact:
        case LST_Exact_Adaptive:
            linesearchSolver = std ::make_unique<ExactLineSearchNM>( domain, engngModel );
            break;
        }
    }

    // if ( this->LsType == LST_Exact_Adaptive ) this->lsFlag = 0;
    return linesearchSolver.get();
}



///////////////////////////////////
// material cylindrical arc-length
REGISTER_SparseNonLinearSystemNM( MaterialCylindricalALM )

MaterialCylindricalALM ::MaterialCylindricalALM( Domain *d, EngngModel *m ) :
    CylindricalALM( d, m )
{
}

void MaterialCylindricalALM::setCALMmaterial( bool returnAll ){ // maybe not used
    auto& matList = this->domain->giveMaterials();
    for ( auto &imat : matList ) {
        auto iface = dynamic_cast<ArcLengthMaterialInterface*>( imat->giveInterface( ArcLengthMaterialInterfaceType ));
        if ( iface && iface->hasArcLengthParameter() ) {
            iface->setReturnAll( returnAll );
        }
    }

}

void MaterialCylindricalALM::setMaterialLambda( double lambdaSet )
{
    auto &matList = this->domain->giveMaterials();
    for ( auto &imat : matList ) {
        auto iface = dynamic_cast<ArcLengthMaterialInterface *>( imat->giveInterface( ArcLengthMaterialInterfaceType ) );
        if ( iface && iface->hasArcLengthParameter() ) {
            iface->setLambda( lambdaSet );
        }
    }
}


ConvergedReason
MaterialCylindricalALM ::solve( SparseMtrx &k, FloatArray &R, FloatArray *R0,
    FloatArray &X, FloatArray &dX, FloatArray &F,
    const FloatArray &internalForcesEBENorm, double &ReachedLambda, referenceLoadInputModeType rlm,
    int &nite, TimeStep *tStep )
{
    FloatArray rhs, deltaXt, deltaX_, dXm1, XInitial;
    FloatArray ddX; // total increment of displacements in iteration
    // double Bergan_k0 = 1.0, bk;
    double XX, RR, RR0, XR, p                                   = 0.0;
    double deltaLambda, Lambda, eta, DeltaLambdam1, DeltaLambda = 0.0;
    double drProduct = 0.0;
    int neq          = R.giveSize();
    int irest        = 0;
    int HPsize, i, ind;
    double _RR, _XX;
    ConvergedReason status = CR_UNKNOWN;
    bool converged, errorOutOfRangeFlag;
    // print iteration header

    ///////////
    FloatArray RT = R + *R0;
    ///////////

    OOFEM_LOG_INFO( "CALMLS:       Initial step length: %-15e\n", deltaL );
    if ( nccdg == 0 ) {
        OOFEM_LOG_INFO( "CALMLS:       Iteration       LoadLevel       ForceError      DisplError    \n" );
        OOFEM_LOG_INFO( "----------------------------------------------------------------------------\n" );
    } else {
        OOFEM_LOG_INFO( "Iter  LoadLevel       " );
        for ( i = 1; i <= nccdg; i++ ) {
            OOFEM_LOG_INFO( "ForceError(%02d)  DisplError(%02d)  ", i, i );
        }

        OOFEM_LOG_INFO( "\n__________________________________________________________\n" );
    }

    //
    // Now smarter method is default:
    // after convergence troubles for (calm_NR_ModeTick subsequent
    // steps newly set calm_NR_Mode will be used. After the calm_NR_OldMode will be restored
    //
    if ( calm_NR_ModeTick == 0 ) {
        calm_NR_Mode = calm_NR_OldMode;
    }

    if ( calm_NR_ModeTick > 0 ) {
        calm_NR_ModeTick--;
    }

    XInitial = X;
    ddX.resize( neq );
    ddX.zero();
    deltaXt.resize( neq );
    deltaXt.zero();

    status = CR_UNKNOWN;
    this->giveLinearSolver();

    // create HPC Map if needed
    if ( calm_hpc_init ) {
        this->convertHPCMap();
        calm_hpc_init = 0;
    }

    if ( R0 ) {
        RR0 = parallel_context->localNorm( *R0 );
        RR0 *= RR0;
    } else {
        RR0 = 0.0;
    }

    //
    // A  initial step (predictor)
    //
    //

    //
    // A.2.   We assume positive-definite (0)Kt (tangent stiffness mtrx).
    //
    
    /////////////
    this->setMaterialLambda( ReachedLambda ); // set material to current lambda
    FloatArray Fint_lam;
    Fint_lam.resize( neq );
    // array of vector of lambda dependent internal forces, role of -Fext in classical arc-length
    engngModel->updateArcLengthInternalForces( Fint_lam, tStep, InternalRhs, domain );

    RR = parallel_context->localNorm( Fint_lam );
    /////////////
    // RR = parallel_context->localNorm( R );
    RR *= RR;

restart:
    //
    // A.1. calculation of (0)VarRt
    //
    dX.zero();
    // engngModel->updateComponent(tStep, InternalRhs, domain); // By not updating this, one obtains the old equilibrated tangent.

    engngModel->updateComponent( tStep, NonLinearLhs, domain );
    /*linSolver->solve( k, R, deltaXt );*/
    linSolver->solve( k, -1.*Fint_lam, deltaXt );

    // If desired by the user, the solution is (slightly) perturbed, so that various symmetries can be broken.
    // This is useful e.g. to trigger localization in a homogeneous material under uniform stress without
    // the need to introduce material imperfections. The problem itself remains symmetric but the iterative
    // solution is brought to a nonsymmetric state and it gets a chance to converge to a nonsymmetric solution.
    // Parameters of the perturbation technique are specified by the user and by default no perturbation is done.
    // Milan Jirasek
    SparseNonLinearSystemNM ::applyPerturbation( &deltaXt );

    if ( calm_Control == calm_hpc_off ) {
        XX = parallel_context->localNorm( deltaXt );
        XX *= XX;
        p = sqrt( XX + Psi * Psi * RR );
    } 
    else if ( calm_Control == calm_hpc_on ) {
        HPsize = calm_HPCIndirectDofMask.giveSize();
        _XX    = 0;
        _RR    = 0;
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at( i ) ) != 0 ) {
                _XX += deltaXt.at( ind ) * deltaXt.at( ind );
                _RR += R.at( ind ) * R.at( ind );
            }
        }

        // In case of paralllel analysis:
        FloatArray collected_XXRR;
        parallel_context->accumulate( { _XX, _RR }, collected_XXRR );
        _XX = collected_XXRR( 0 );
        _RR = collected_XXRR( 1 );

        p = sqrt( _XX + Psi * Psi * _RR );
    } else if ( calm_Control == calml_hpc ) {
        HPsize = calm_HPCIndirectDofMask.giveSize();
        p      = 0.;
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at( i ) ) != 0 ) {
                p += deltaXt.at( ind ) * calm_HPCWeights.at( i );
            }
        }

        // In case of paralllel analysis:
        p = parallel_context->accumulate( p );
    }

    Lambda = ReachedLambda;
    if ( rootselectiontype == RST_Cos ) {
        /* XR is unscaled Bergan's param of current stiffness XR = deltaXt^T k deltaXt
         * this is used to test whether k has negative or positive slope */

        //XR = parallel_context->localDotProduct( deltaXt, R );
        XR = parallel_context->localDotProduct( deltaXt, -1. * Fint_lam );
    } else {
        if ( old_dX.giveSize() ) {
            XR = parallel_context->localDotProduct( deltaXt, old_dX );
        } else {
            //XR = parallel_context->localDotProduct( deltaXt, R );
            XR = parallel_context->localDotProduct( deltaXt, -1. * Fint_lam );
        }
    }
    ///////////////
    //if ( abs( p ) < 1e-8) {
    //    DeltaLambda = 0.;
    //} else {
        DeltaLambda = deltaLambda = sgn( XR ) * deltaL / p;
    //}
    ////////////////////

    Lambda += DeltaLambda;
    //
    // A.3.
    //

    //////////
    this->setMaterialLambda( Lambda ); // set material to current lambda
    //////////
    engngModel->updateComponent( tStep, InternalRhs, domain ); // since lambda has changed

    rhs = R;
    /*rhs.times( DeltaLambda );*/
    if ( R0 ) {
        rhs.add( *R0 );
    }
    /////////////
    rhs.subtract( F );
    /////////////

    linSolver->solve( k, rhs, dX );
    X.add( dX );

    nite = 0;

    // update solution state counter
    tStep->incrementStateCounter();
    engngModel->updateComponent( tStep, InternalRhs, domain );

        //////////////////
    // bifurcatiuon modification
    BifurcationInterface *stabSolver = dynamic_cast<BifurcationInterface *>( this->linSolver.get() ); // Check if bifurcation solver is used
    if ( stabSolver ) {
        // Discard the definiteness recorded for the previous step, or for a previous
        // attempt of this step, so that it is never reported for the current one.
        stabSolver->resetPDStatus();
    }
    bool isBifurcationSet            = false;
    double alphaStability = 1e-6, alphamax = 1e8; // eigenvector multiplicator

    bool LDLTbif = false; // Only if choelsky bifurcation should be performed, if FALSE, bifurcation using eigenvectors is done
     bool deflationBifurcation = true, eigenvectorBifurcation = true, postBifurcationLineSearch = true; // deflation/eigenvectorpostbif. LS options
    //bool deflationBifurcation = false, eigenvectorBifurcation = false, postBifurcationLineSearch = false; // deflation/eigenvector/postbif. LS options
    std::vector<bool> bifurcTypes( { LDLTbif, deflationBifurcation, eigenvectorBifurcation } );

    //if ( this->LsType == LST_Exact_Adaptive ) this->lsFlag = false; // For adaptive LS
    bool afterBeforeBifurcation = true;

    //double tfun = this->engngModel->giveDomain( 1 )->giveFunction( 4 )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //OOFEM_LOG_INFO( "t=%f\n", tfun );
    //////////////////
    double deltaLold = deltaL;
    do {
        nite++;
        tStep->incrementSubStepNumber();

        dXm1          = dX;
        DeltaLambdam1 = DeltaLambda;

        //
        // B  - iteration MNRM is used
        //
        // B.1. is ommited because MNRM is used instead of NRM.
        //
        

        if ( ( calm_NR_Mode == calm_fullNRM ) || ( ( calm_NR_Mode == calm_accelNRM ) && ( nite % calm_MANRMSteps == 0 ) ) ) {
            //
            // ALM with full NRM
            //
            // we assemble new tangent stiffness and compute new deltaXt
            // internal state of elements is updated by previous calling
            // InternalRhs
            //
            ///////////
            engngModel->updateArcLengthInternalForces( Fint_lam, tStep, InternalRhs, domain );
            RR = parallel_context->localNorm( Fint_lam );
            RR *= RR; // RR now changes
            ///////////

            engngModel->updateComponent( tStep, NonLinearLhs, domain );

            ////////////////
            //EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &k );
            //auto mat      = Ae->giveMatrix();
            //k.printYourself();
            ////////////////
            //
            // compute deltaXt for i-th iteration
            //
            //linSolver->solve( k, R, deltaXt );
            linSolver->solve( k, -1. * Fint_lam, deltaXt );
        }

        // B.2.
        //

        rhs = R;
        //rhs.times( Lambda );
        if ( R0 ) {
            rhs.add( *R0 );
        }

        rhs.subtract( F );
        deltaX_.resize( neq );
        linSolver->solve( k, rhs, deltaX_ );
        eta = 1.0;
        //
        // B.3.
        //
        if ( this->computeDeltaLambda( deltaLambda, dX, deltaXt, deltaX_, R, RR, eta, deltaL, DeltaLambda, neq ) ) {
            irest++;
            if ( irest <= maxRestarts ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length
                deltaL = deltaL * CALM_RESET_STEP_REDUCE;
                if ( deltaL < minStepLength ) {
                    deltaL = minStepLength;
                }

                // restore previous total displacement vector
                X = XInitial;
                // reset all changes fro previous equilibrium state
                dX.zero();
                // restore initial stiffness
                engngModel->updateComponent( tStep, NonLinearLhs, domain );

                OOFEM_LOG_INFO( "CALMLS:       Iteration Reset ...\n" );

                calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = CR_DIVERGED_ITS;
                OOFEM_ERROR( "can't continue further" );
            }
        }

        //
        // B.4.+ B.5.
        //


        if ( this->lsFlag && ( nite != 1 ) ) {
            //
            //  LINE SEARCH
            //
            // this->do_lineSearch( X, XInitial, deltaX_, deltaXt,
            //    dXm1, dX, ddX,
            //    R, R0, F,
            //    DeltaLambda, DeltaLambdam1, deltaLambda, Lambda, ReachedLambda,
            //    RR, drProduct, tStep );

            ////////////////////////
            //// line search
            // this->solve_linesearch( X, XInitial, ddX, F, RT, tStep, k, ReachedLambda, Lambda, DeltaLambda, DeltaLambdam1, deltaLambda, deltaX_, deltaXt, dX, dXm1 );
            //if ( stabSolver && stabSolver->giveX0Defl().giveSize() > 0 ) {
            //    this->solve_linesearch( X, XInitial, ddX, F, RT, tStep, k, ReachedLambda, Lambda, DeltaLambda, DeltaLambdam1, deltaLambda,
            //        deltaX_, deltaXt, dX, dXm1, true, stabSolver->giveX0Defl() );
            //} else {
                this->solve_linesearch( X, XInitial, ddX, F, RT, tStep, k, ReachedLambda, Lambda, DeltaLambda, DeltaLambdam1, deltaLambda, deltaX_, deltaXt, dX, dXm1 );
            //}


            drProduct = parallel_context->localNorm( ddX );
            drProduct *= drProduct;
            ////////////////////////
        //} else if ( this->lsFlag && nite == 1 && stabSolver && stabSolver->giveX0Defl().giveSize() > 0 ) {
        //    this->solve_linesearch( X, XInitial, ddX, F, RT, tStep, k, ReachedLambda, Lambda, DeltaLambda, DeltaLambdam1, deltaLambda,
        //        deltaX_, deltaXt, dX, dXm1, true, stabSolver->giveX0Defl() );
        } else { // no line search
            //
            // update solution vectors
            //
            ddX.clear();
            ddX.add( eta * deltaLambda, deltaXt );
            ddX.add( eta, deltaX_ );
            dX = dXm1;
            dX.add( ddX );
            /////////////
            double temp = dX.computeNorm();
            OOFEM_LOG_INFO( "dX_norm = %f \n", temp );
            OOFEM_LOG_INFO( "dX(10) = %.3e \n", dX.at( calm_HPCIndirectDofMask.at(1) ) );
            ////////////
            X = XInitial;
            X.add( dX );

            drProduct = parallel_context->localNorm( ddX );
            drProduct *= drProduct;

            tStep->incrementStateCounter(); // update solution state counter
            //
            // B.6.
            //
            DeltaLambda = DeltaLambdam1 + deltaLambda;
            Lambda      = ReachedLambda + DeltaLambda;

            //////////
            this->setMaterialLambda( Lambda ); // set material to current lambda
            //OOFEM_LOG_INFO( "Lambda = %-15e\n", Lambda );
            //////////

            tStep->incrementStateCounter(); // update solution state counter
            engngModel->updateComponent( tStep, InternalRhs, domain );
        }

        //
        // B.7.
        //
        // convergence check
        //

        converged = this->checkConvergence( R, R0, F, X, ddX, Lambda, RR0, RR, drProduct,
            internalForcesEBENorm, nite, errorOutOfRangeFlag );
        //converged = CylindricalALM::checkConvergence( R, R0, F, X, ddX, Lambda, RR0, RR, drProduct,
        //    internalForcesEBENorm, nite, errorOutOfRangeFlag );

        // bifurcatiuon modification
        if ( stabSolver ) this->performBifurcationAnalysis( k, X, dX, F, internalForcesEBENorm, rlm, nite, tStep, rhs, converged, R, *R0,
            errorOutOfRangeFlag, stabSolver, isBifurcationSet, alphamax, alphaStability, ddX, bifurcTypes, postBifurcationLineSearch, Lambda, RR0, RR, drProduct );
        ////////////////// 

        if ( ( nite >= nsmax ) || errorOutOfRangeFlag ) {
            irest++;
            if ( irest <= maxRestarts ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length
                deltaL = deltaL * CALM_RESET_STEP_REDUCE;
                if ( deltaL < minStepLength ) {
                    deltaL = minStepLength;
                }

                // restore previous total displacement vector
                X = XInitial;
                // reset all changes from previous equilibrium state
                engngModel->initStepIncrements();
                dX.zero();
                // restore initial stiffness
                engngModel->updateComponent( tStep, NonLinearLhs, domain );

                OOFEM_LOG_INFO( "CALMLS:       Iteration Reset ...\n" );

                calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = CR_DIVERGED_ITS;
                OOFEM_WARNING( "Convergence not reached after %d iterations", nsmax );
                // exit(1);
                break;
            }
        }

        // output of per iteration data
        engngModel->giveExportModuleManager()->doOutput( tStep, true );
    } while ( !converged || ( nite < minIterations ) );

    //
    // update dofs,nodes,Elemms and print result
    //
#ifdef VERBOSE
    // printf ("\nCALM - step iteration finished") ;
#endif

    // Seems to be completely unused / Mikael
#if 0
    /* compute Bergan's parameter of current stiffness */
    double RDR = parallel_context->localDotProduct(* R, * dX);
    double DR = parallel_context->localNorm(* dX);
    bk = DeltaLambda * RDR / ( DR * DR );

    if ( tStep->giveNumber() == 1 ) {
        // compute Bergan_k0
        Bergan_k0 = bk;
    } else {
        bk = bk / Bergan_k0;
        //if (fabs(bk) < CALM_TANGENT_STIFF_TRESHOLD) status |= NM_KeepTangent;
    }
#endif

    /* we have computed Bergan's  parameter -> you can adjust deltaL according
     * to this value */
    //
    // there has been restart already - set nite to maxiter
    //
    // if (irest > 0) nite = nsmax;

    if ( nite > numberOfRequiredIterations ) {
        deltaL = deltaL * numberOfRequiredIterations / nite;
    } else {
        deltaL = deltaL * sqrt( sqrt( (double)numberOfRequiredIterations / (double)nite ) );
    }


    ///////////////
    //if ( stabSolver ) {
    //    double val = 6.6;
    //    if ( Lambda > val && !this->dL_update ) {
    //        deltaLold       = deltaL;
    //        deltaL          = minStepLength;
    //        //this->dL_update = true;
    //    }
    //    // else if ( stabSolver ->giveFoundLimitPoint()) {
    //    //     deltaL = deltaLold;
    //    // }
    //}
    /////////////////

    if ( deltaL > maxStepLength ) {
        deltaL = maxStepLength;
    }

    if ( deltaL < minStepLength ) {
        deltaL = minStepLength;
    }



    OOFEM_LOG_INFO( "CALMLS:       Adjusted step length: %-15e\n", deltaL );

    status        = CR_CONVERGED;
    solved        = 1;
    ReachedLambda = Lambda;
    old_dX        = dX;
    return status;
}


void MaterialCylindricalALM::performBifurcationAnalysis( SparseMtrx &k, FloatArray &X, FloatArray &dX, FloatArray &F, const FloatArray &internalForcesEBENorm,
    referenceLoadInputModeType rlm, int &nite, TimeStep *tStep, FloatArray &rhs, bool &converged, const FloatArray &R, const FloatArray &R0,
    bool &errorOutOfRangeFlag, BifurcationInterface *stabSolver, bool &isBifurcationSet, double &alphamax, double &alphaStability,
    FloatArray &ddX, std::vector<bool> &bifurcTypes, bool &postBifurcationLineSearch, double Lambda, double RR0, double RR, double drProduct )
{
    // Choose bifurcation type
    bool LDLTbif                          = bifurcTypes[0]; // Only if choelsky bifurcation should be performed, if FALSE, bifurcation using eigenvectors is done
    bool deflationBifurcation             = bifurcTypes[1]; // only if deflation bifurcation if performed
    bool eigenvectorBifurcation           = bifurcTypes[2]; // only if deflation bifurcation if performed
    //bool postDeflation                    = false;
    bool postDeflation                    = false;
    bool exactLineSearchBeforeBifurcation = true;


    bool performBifurcation = LDLTbif || deflationBifurcation || eigenvectorBifurcation; // Any bifurcation is on, eigenvectors are default
    // For bifurcation using eigenvectors
    double alphaStability0 = alphaStability;
    int maxBifIte = 100, count = 0;


    FloatArray RT = R + R0;

    // bifurcation analysis modifications
    if ( stabSolver && performBifurcation ) {
        if ( deflationBifurcation && isBifurcationSet ) { // If deflation, modify residuum and recheck convergence
            double dx_defl_Norm = ( X - stabSolver->giveX0Defl() ).computeNorm();
            double p            = stabSolver->give_p();
            //double factor       = ( 1. / pow( dx_defl_Norm, p ) + 1 );
            double xd_norm = stabSolver->giveX0Defl().computeNorm();
            double factor  =  pow( xd_norm, p ) / pow( dx_defl_Norm, p ) + 1. ;
            //FloatArray rhs_mod  = rhs * ( 1. / ( dx_defl_Norm * dx_defl_Norm ) + 1 );
            //FloatArray rhs_mod = rhs * ( 1. / pow( dx_defl_Norm, p ) + 1 );
            //converged           = this->checkConvergence( RT, F, rhs_mod, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag ); // Residuum needs to be modified
            FloatArray R0mod = R0 * factor;
            FloatArray Rmod = R * factor;
            FloatArray Fmod = F * factor;
            //FloatArray R0mod = R0 ;
            //FloatArray Rmod  = R ;
            //FloatArray Fmod  = F;
            converged        = this->checkConvergence( Rmod, &R0mod, Fmod, X, ddX, Lambda, RR0, RR, drProduct, internalForcesEBENorm, nite, errorOutOfRangeFlag );
            //converged = CylindricalALM::checkConvergence( R, &R0, F, X, ddX, Lambda, RR0, RR, drProduct, internalForcesEBENorm, nite, errorOutOfRangeFlag );

        } else if ( isnan( F.computeNorm() ) && isBifurcationSet && eigenvectorBifurcation ) { // if eigenvector bif. and not defl. and NaNs appear
            OOFEM_LOG_INFO( "Retarting step due to NaN\n" );
            X -= dX; // go back to start of the step
            dX.zero();
            alphamax = stabSolver->giveAlpha(); // modify amplitude of the perturbation
            stabSolver->setAlpha( alphamax - alphaStability );
            alphaStability = alphaStability / 2;

            engngModel->updateComponent( tStep, InternalRhs, domain );
            //rhs.beDifferenceOf( RT, F );
            //if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );

            //converged = this->checkConvergence( RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag ); // recheck the convergence
            converged = this->checkConvergence( R, &R0, F, X, ddX, Lambda, RR0, RR, drProduct, internalForcesEBENorm, nite, errorOutOfRangeFlag );
        }


        if ( converged ) {
            engngModel->updateComponent( tStep, NonLinearLhs, domain );
            //applyConstraintsToStiffness( k );
            bool isPD = stabSolver->checkPD( k );


            ///////////////////////
            //// // condense lagrange multipliers from bc
            // EigenMtrx *k_eigmat_lam = dynamic_cast<EigenMtrx *>( &k );
            // for ( int i = 0; i < domain->giveBcs().size(); i++ ) {
            //     auto &bc = domain->giveBcs()[i];
            //     auto abc = dynamic_cast<IGAClampBoundaryConditionLagrange *>( bc.get() );
            //     if ( abc ) {
            //         IntArray lambdaLoc( 1 );
            //         EModelDefaultEquationNumbering dn;
            //         abc->giveLagrangianMultiplierLocationArray( dn, lambdaLoc );

            //        // now this just takes the K11 block
            //        std::unique_ptr<EigenMtrx> k_condensed = k_eigmat_lam->doStaticCondensationLagrange( lambdaLoc.at(1) );
            //        // k_condensed->printYourself();
            //        isPD = stabSolver->checkPD( *k_condensed );

            //        break; // we just find the first lambda and assume the multipliers are ordered, this shouold be made more general
            //    }
            //}
            ///////////////////////

            ///////////////
            // condense the pressure DOFs
            int indPressure = 0;
            bool exVal      = this->findPressureDofIndex( indPressure ); // thios goves starting index of pressure dofs

            IntArray pinds;
            bool exVal2 = findPressureDofIndices( pinds );

            EigenMtrx *k_eigmat = dynamic_cast<EigenMtrx *>( &k );
            if ( exVal && k_eigmat ) {
                std::unique_ptr<EigenMtrx> k_condensed = k_eigmat->doStaticCondensation( indPressure ); // Doesn work if pressure dofs are not at end!!
                // k_condensed->printYourself();
                isPD = stabSolver->checkPD( *k_condensed );
                //////////////////
                
                /////////////////
            }
            //if ( !isPD ) {
            //    std::unique_ptr<EigenMtrx> k_condensed2 = k_eigmat->doStaticCondensation( indPressure );
            //}
            //////////////////////
            ///////////////////
            //if ( !isPD && stabSolver->giveNumFoundSols() == 1 ) {
            //    double norm = ( X - stabSolver->giveX0Defl() ).computeNorm();
            //    OOFEM_LOG_INFO( "dXdefl =  %f\n", norm );
            //}
            ///////////////////

             //if ( !isPD && stabSolver->giveNumFoundSols() < 1 ) { // If matrix is not PD
            if ( !isPD ) { // If matrix is not PD

                OOFEM_LOG_INFO( "Negative eigenvalue\n" );
                /////////

                /////////////////////
                // just a check, should be commented
                //if ( this->calm_HPCWeights.giveSize() == 2 ) {
                //    calm_HPCWeights.at( 2 ) = -calm_HPCWeights.at( 1 );
                //    converged               = false;
                //    return;
                //}
                // 
                //this->deltaL = this->deltaL * 0.1;
                //converged               = false;
                //return;

                if (false && calm_Control == calm_hpc_on && calm_HPCIndirectDofMask.giveSize() == 1 ) { // just testing
                    converged = false;
                    ///
                    // do deflation test
                    FloatArray X0 = X; // Store the instable root
                    stabSolver->setX0Defl( X0 ); // set the instable root to the solver
                    stabSolver->setBifurcation( true ); // perform bifurcation anaysis
                    stabSolver->setDeflation( deflationBifurcation );
                    ///

                    if ( stabSolver->getXeigs().giveSize() == 0 || ( X - stabSolver->getXeigs() ).computeNorm() > 1e-8 ) { // if not the same, recompute the eigenvectors
                        engngModel->updateComponent( tStep, NonLinearLhs, domain );
                        // applyConstraintsToStiffness( k );

                        if ( exVal && k_eigmat ) {
                            // EigenMtrx *k_eigmat2 = dynamic_cast<EigenMtrx *>( &k );
                            // k.printYourself();

                            std::unique_ptr<EigenMtrx> k_condensed = k_eigmat->doStaticCondensation( indPressure );
                            stabSolver->setEigenValuesVectors( *k_condensed, X );
                        } else {
                            stabSolver->setEigenValuesVectors( k, X );
                        }
                    }
                    FloatArray direction; // get eigenvectors
                    stabSolver->getEigenVectors().copyColumn( direction, stabSolver->getEigenVectors().giveNumberOfColumns() );
                    // scale direction
                    int HPsize = calm_HPCIndirectDofMask.giveSize();
                    int ind;
                    double dDir = 0., dDX = 0.;
                    for ( int i = 1; i <= HPsize; i++ ) {
                        if ( ( ind = calm_HPCIndirectDofMask.at( i ) ) != 0 ) {
                            dDir += direction.at( ind );
                            // dDX += dX.at( ind ) * calm_HPCWeights.at( i );
                        }
                    }
                    direction = -1. * direction * deltaL / dDir;

                    X -= dX;
                    if ( exVal ) {
                        for ( int i = 1; i < indPressure; i++ ) {
                            dX.at( i ) = direction.at( i );
                        }
                    } else {
                        dX = direction;
                    }
                    ddX = dX;
                    X += dX;

                    ////// for deflation
                    stabSolver->compute_dx_defl( X );
                    isBifurcationSet = true; // Not sure if correct
                    ///////

                    OOFEM_LOG_INFO( "dX(10) = %.3e \n", dX.at( calm_HPCIndirectDofMask.at( 1 ) ) );

                    engngModel->updateComponent( tStep, InternalRhs, domain );

                    ////////////////
                    // engngModel->updateComponent( tStep, NonLinearLhs, domain );
                    // EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &k );
                    // auto mat      = Ae->giveMatrix();
                    // k.printYourself();
                    // EigenMtrx *k_eigmat = dynamic_cast<EigenMtrx *>( &k );
                    //
                    // std::unique_ptr<EigenMtrx> k_condensed = k_eigmat->doStaticCondensation( indPressure ); // Doesn work if pressure dofs are not at end!!
                    //// k_condensed->printYourself();
                    // isPD = stabSolver->checkPD( *k_condensed );
                    ////////////////
                    // rhs.beDifferenceOf( RT, F );
                    // isBifurcationSet = false;
                    // stabSolver->setBifurcation( false );

                    
                    return;
                } ///////////////////////
               
                /////////////////////

                if ( stabSolver->giveX0Defl().giveSize() > 0 ) {
                    double norm = ( X - stabSolver->giveX0Defl() ).computeNorm();
                    OOFEM_LOG_INFO( "dXdefl =  %f\n", norm );
                }
                
                /////////
                //double tfun = this->engngModel->giveDomain( 1 )->giveFunction( 3 )->evaluateAtTime( tStep->giveIntrinsicTime() );
                //OOFEM_LOG_INFO( "gamma=%f\n", tfun );
                stabSolver->setFoundLimitPoint( true );
                stabSolver->setPostBifurcationLineSearchSolver( postBifurcationLineSearch ); // Set the post bifurcation exact linesearch
                // this->setExactPreBifurcationLineSearch( false ); // turn off the pre Bifurcational LS
                stabSolver->setLineSearchState( this->lsFlag ); // remember linesearch status
                /////////////////////////
                this->lsFlag = false; // turn off the pre Bifurcational LS /// THIS IS CHANGED!
                //////////////////////////

                if ( !isBifurcationSet ) { // If first time not PD
                    stabSolver->setBifurcation( true ); // perform bifurcation anaysis
                    stabSolver->setAlpha( alphaStability0 - alphaStability ); // set eigenvector multiplicator to starting value
                    isBifurcationSet = true; // Not sure if correct
                }

                if ( deflationBifurcation || eigenvectorBifurcation ) {
                    converged = false;
                    ////
                    stabSolver->incrementNumFoundSols();
                    ////
                    FloatArray X0 = X; // Store the instable root
                    stabSolver->setX0Defl( X0 ); // set the instable root to the solver

                    // Find out if eigenvalues should be recomputed
                    if ( stabSolver->getXeigs().giveSize() == 0 || ( X - stabSolver->getXeigs() ).computeNorm() > 1e-8 ) { // if not the same, recompute the eigenvectors
                        engngModel->updateComponent( tStep, NonLinearLhs, domain );
                        //applyConstraintsToStiffness( k );

                        if ( exVal && k_eigmat ) {
                            //EigenMtrx *k_eigmat2 = dynamic_cast<EigenMtrx *>( &k );
                            //k.printYourself();

                            std::unique_ptr<EigenMtrx> k_condensed = k_eigmat->doStaticCondensation( indPressure );
                            stabSolver->setEigenValuesVectors( *k_condensed, X );
                        } else {
                            stabSolver->setEigenValuesVectors( k, X );
                        }
                        
                    }
                    FloatArray direction; // get eigenvectors
                    stabSolver->getEigenVectors().copyColumn( direction, stabSolver->getEigenVectors().giveNumberOfColumns() );

                    /////////////////////
                    //// just a check, should be commented
                    //if ( this->calm_HPCWeights.giveSize() == 2 ) {
                    //    calm_HPCWeights.at( 2 ) = -calm_HPCWeights.at( 1 );
                    //    converged               = false;
                    //    // scale direction
                    //    int HPsize = calm_HPCIndirectDofMask.giveSize();
                    //    int ind;
                    //    double dDir = 0., dDX = 0.;
                    //    for ( int i = 1; i <= HPsize; i++ ) {
                    //        if ( ( ind = calm_HPCIndirectDofMask.at( i ) ) != 0 ) {
                    //            dDir += direction.at( ind ) * calm_HPCWeights.at( i );
                    //            //dDX += dX.at( ind ) * calm_HPCWeights.at( i );
                    //        }
                    //    }
                    //    direction = direction * deltaL / dDir;
                    //}
                    /////////////////////
                    
                    
                    //if ( calm_Control == calm_hpc_on && calm_HPCIndirectDofMask.giveSize() == 1 ) { // just testing
                    //    converged = false;
                    //    // scale direction
                    //    int HPsize = calm_HPCIndirectDofMask.giveSize();
                    //    int ind;
                    //    double dDir = 0., dDX = 0.;
                    //    for ( int i = 1; i <= HPsize; i++ ) {
                    //        if ( ( ind = calm_HPCIndirectDofMask.at( i ) ) != 0 ) {
                    //            dDir += direction.at( ind );
                    //            // dDX += dX.at( ind ) * calm_HPCWeights.at( i );
                    //        }
                    //    }
                    //    direction = -1.* direction * deltaL / dDir;

                    //    X -= dX;
                    //    if ( exVal ) {
                    //        for ( int i = 1; i < indPressure; i++ ) {
                    //            dX.at( i ) = direction.at( i );
                    //        }
                    //    }
                    //    ddX  = dX;
                    //    X += dX;

                    //    OOFEM_LOG_INFO( "dX(10) = %.3e \n", dX.at( calm_HPCIndirectDofMask.at(1) ) );

                    //    engngModel->updateComponent( tStep, InternalRhs, domain );

                    //    ////////////////
                    //    //engngModel->updateComponent( tStep, NonLinearLhs, domain );
                    //    //EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &k );
                    //    //auto mat      = Ae->giveMatrix();
                    //    //k.printYourself();
                    //    //EigenMtrx *k_eigmat = dynamic_cast<EigenMtrx *>( &k );
                    //    //
                    //    //std::unique_ptr<EigenMtrx> k_condensed = k_eigmat->doStaticCondensation( indPressure ); // Doesn work if pressure dofs are not at end!!
                    //    //// k_condensed->printYourself();
                    //    //isPD = stabSolver->checkPD( *k_condensed );
                    //    ////////////////
                    //    //rhs.beDifferenceOf( RT, F );
                    //    //isBifurcationSet = false;
                    //    //stabSolver->setBifurcation( false );
                    //    //return;
                    //}///////////////////////

                    if ( !eigenvectorBifurcation ) { // Only deflation (but initial perturbation by eigenvectors)
                        ddX = alphaStability * direction; // perturbation in direction of eigenvector
                        X   = X + ddX; // Perturbe the solution with chosen amplitude
                        dX  = dX + ddX;
                        stabSolver->compute_dx_defl( X ); // compute deflation difference

                        // Update residuum
                        engngModel->updateComponent( tStep, InternalRhs, domain );
                        rhs.beDifferenceOf( RT, F );
                        //if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );

                    } else { // Eigenvector bifurcation with/without deflation
                        if ( deflationBifurcation ) { // If combination of deflation and eigenvector bifurcation
                            // if ( false ) this->provisionalOutput( X, F, tStep, RT, direction, X0 ); // PROVISIONAL output residual
                            ddX = direction;
                            IntArray dummyInt( 0 );
                            this->giveSetPostBifurcationLineSearchSolver( X0, true, stabSolver->give_p() )->solve( X, ddX, F, RT, dummyInt, tStep, k ); // Do linesearch with deflation
                            X += ddX;
                            dX += ddX;
                            ///
                            //this->solve_linesearch( X, XInitial, ddX, F, RT, tStep, k, ReachedLambda, Lambda, DeltaLambda, DeltaLambdam1, deltaLambda, deltaX_, deltaXt, dX, dXm1 );
                            ///
                            engngModel->updateComponent( tStep, InternalRhs, domain );
                            rhs.beDifferenceOf( RT, F );
                            //if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );

                            stabSolver->compute_dx_defl( X );
                            if ( !postDeflation ) { // Continue without deflation
                                isBifurcationSet = false;
                                stabSolver->setBifurcation( false );
                            }
                            /////////////////////// just a test
                            else {
                                stabSolver->setDeflation( deflationBifurcation ); // This should be done only once
                            }
                            /////////////////////
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
                ////////////////
                stabSolver->nullNumFoundSols();
                ////////////////
                // Compute distance from the instable solution
                if ( stabSolver->getBifurcation() ) {
                    FloatArray solDif;
                    solDif.beDifferenceOf( X, stabSolver->giveX0Defl() );
                    double distX0 = solDif.computeNorm();
                }
                stabSolver->setBifurcation( false ); // only needed for cholesky bifurcation
                isBifurcationSet = false; // Not sure if correct

                // Turn the post bifurcation linesearch off after the limit point is found (turn it on only if NAN appears)
                if ( stabSolver->giveFoundLimitPoint() ) {
                    stabSolver->setPostBifurcationLineSearchSolver( false );
                    this->lsFlag = stabSolver->giveLineSearchState(); // Set the lsflag to the value it was before bifurcation
                }
            }
        } else if ( nite == 0 ) { // not converged, check if the bifurcation analysis should be performed
            stabSolver->setCholesky( LDLTbif ); // This should be done only once
            stabSolver->setDeflation( deflationBifurcation ); // This should be done only once
            isBifurcationSet = stabSolver->getBifurcation();
        } else if ( stabSolver->getBifurcation() && deflationBifurcation ) {
            stabSolver->compute_dx_defl( X );

        } else if ( stabSolver->giveFoundLimitPoint() && ( stabSolver->givePostBifurcationLineSearchSolver() ) ) { // if limit point found and exact post bif LS should be done
            X  = X - ddX;
            dX = dX - ddX;
            FloatArray dummy( 0 );
            IntArray dummyInt( 0 );
            this->giveSetPostBifurcationLineSearchSolver( dummy, true, stabSolver->give_p() )->solve( X, ddX, F, RT, dummyInt, tStep, k );
            X += ddX;
            dX += ddX;

            engngModel->updateComponent( tStep, InternalRhs, domain );
            rhs.beDifferenceOf( RT, F );
            //if ( this->prescribedDofsFlag ) this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );
        }
    }
    bifurcTypes[0] = LDLTbif;
    bifurcTypes[1] = deflationBifurcation;
    bifurcTypes[2] = eigenvectorBifurcation;
}

ExactLineSearchNM *MaterialCylindricalALM ::giveSetPostBifurcationLineSearchSolver( FloatArray &X0, bool doDeflation, double p_set )
{
    if ( !this->linesearchSolverPostBifurcation ) {
        this->linesearchSolverPostBifurcation = std ::make_unique<ExactLineSearchNM>( domain, engngModel );
        this->linesearchSolverPostBifurcation->setMaxIter( 20 );
        this->linesearchSolverPostBifurcation->setTolerance( 1e-10 );
        this->linesearchSolverPostBifurcation->setDeflation( doDeflation );
        this->linesearchSolverPostBifurcation->set_p( p_set );
    }

    if ( X0.giveSize() > 0 ) this->linesearchSolverPostBifurcation->setX0defl( X0 );
    return this->linesearchSolverPostBifurcation.get();
}

bool MaterialCylindricalALM::findPressureDofIndex( int &indMin )
{
    EModelDefaultEquationNumbering dn;
    indMin         = INT_MAX;
    bool ReturnVal = false;

    for ( auto &elem : domain->giveElements() ) {
        for ( int idofman = 1; idofman <= elem->giveNumberOfInternalDofManagers(); idofman++ ) {
            DofManager *dofman = elem->giveInternalDofManager( idofman );
            for ( Dof *dof : *dofman ) {
                auto dofIdType = dof->giveDofID();
                if ( dofIdType == P_f ) {
                    int eq = dof->giveEquationNumber( dn );
                    if ( eq < indMin ) {
                        indMin    = eq;
                        ReturnVal = true;
                    }
                }
            }
        }
    }

    /////////////////////

    // EModelDefaultEquationNumbering dn;
    // indMin = 1e32;
    // bool ReturnVal = false;
    // for ( auto &dman : this->domain->giveDofManagers() ) {
    //     //if ( dman->isNull() ) continue;
    //     for ( auto &dof : *dman ) {
    //         auto dofIdType = dof->giveDofID();
    //         if ( dofIdType == P_f ) {
    //             int eq = dof->giveEquationNumber( dn );
    //             if ( eq < indMin ) {
    //                 indMin = eq;
    //                 ReturnVal = true;
    //             }
    //         }
    //     }
    // }
    return ReturnVal;
}


bool MaterialCylindricalALM::findPressureDofIndices( IntArray &pinds )
{
    pinds.resize(0);
    EModelDefaultEquationNumbering dn;
    bool ReturnVal = false;

    for ( auto &elem : domain->giveElements() ) {
        for ( int idofman = 1; idofman <= elem->giveNumberOfInternalDofManagers(); idofman++ ) {
            DofManager *dofman = elem->giveInternalDofManager( idofman );
            for ( Dof *dof : *dofman ) {
                auto dofIdType = dof->giveDofID();
                if ( dofIdType == P_f ) {
                    int eq = dof->giveEquationNumber( dn );
                    pinds.followedBy( eq );
                    ReturnVal = true;
                }
            }
        }
    }

    return ReturnVal;
}

void MaterialCylindricalALM::getCondensedMatrixPressure( SparseMtrx &k )
{
    int indPressure = 0;
    bool exVal      = this->findPressureDofIndex( indPressure );

    EigenMtrx *k_eigmat = dynamic_cast<EigenMtrx *>( &k );
    if ( exVal && k_eigmat ) {
        k_eigmat->doStaticCondensationOwn( indPressure );

    } 
}

//std::unique_ptr<SparseMtrx> MaterialCylindricalALM::getCondensedMatrixPressure(std::unique_ptr<SparseMtrx> k )
//{
//    int indPressure = 0;
//    bool exVal      = this->findPressureDofIndex( indPressure );
//
//    // cast the raw pointer inside the unique_ptr
//    if ( exVal ) {
//        if ( auto *k_eigmat = dynamic_cast<EigenMtrx *>( k.get() ) ) {
//            return k_eigmat->doStaticCondensation( indPressure );
//        }
//    }
//
//    // no condensation; return the original matrix
//    return k; // moved implicitly
//}


void MaterialCylindricalALM ::solve_linesearch( FloatArray &X, const FloatArray &Xinitial, FloatArray &ddX, FloatArray &F, const FloatArray &RT, TimeStep *tStep, SparseMtrx &k, 
    const double &ReachedLambda, double &Lambda, double &DeltaLambda, const double &DeltaLambdam1, double &deltaLambda, 
    const FloatArray &deltaX_, const FloatArray &deltaXt, FloatArray &dX, const FloatArray &dXm1, bool deflation , const FloatArray &Xdefl)
{
    int max_iter = 20;
    // Normalize the direction
    int neq                  = X.giveSize();
    FloatArray direction, rhs, Ku( neq ), dx_Defl, rhsDefl, Fint_lam( neq );

    // save old arrays
    FloatArray ddXold( ddX );
    double  deltaLambda0( deltaLambda );

    //// if residuum is small enough already
    //if ( fabs( ( RT - F ).dotProduct( direction ) ) < this->ls_tolerance ) {
    //    return;
    //}

    //double deltaX__d, deltaXt_d, RR;
    double RR;

    // do exact linesearch including delfation
    int p = 2;
    double RHS, ResNormLS, RHSdefl, deta, alph_defl = 1., dx_norm, gamma, dx_normSq, kdefl;

    double Fint_lam_d;
    double Eta = 1.;

    double rhs_0;
    double resMult = 0.25;
    for ( size_t niteLS = 0; niteLS < max_iter; niteLS++ ) { // inner Newton loop to find local minimum in the direction
        //OOFEM_LOG_INFO( "deltaLambda = %.7e\n", deltaLambda );
        direction = deltaX_ + deltaXt*deltaLambda;
        ddX       = Eta * direction;
        dX = dXm1 + ddX;
        X  = Xinitial + dX;

        // Update lambda
        DeltaLambda = DeltaLambdam1 + deltaLambda;
        Lambda      = ReachedLambda + DeltaLambda;
        this->setMaterialLambda( Lambda ); // set material to current lambda

        // Update internal RHS
        engngModel->updateComponent( tStep, InternalRhs, domain );
        rhs     = RT - F;
        RHS     = rhs.dotProduct( direction );
        RHSdefl = RHS;

        if ( deflation ) { // Compute deflated residuum
        dx_Defl   = X - Xdefl;
        dx_normSq = dx_Defl.dotProduct( dx_Defl );
        RHSdefl *= 1. / dx_normSq + 1.;
        }

        if ( niteLS  == 0) {
            rhs_0 = RHS;
        }


        // Update F_lambda
        engngModel->updateArcLengthInternalForces( Fint_lam, tStep, InternalRhs, domain );
        RR = parallel_context->localNorm( Fint_lam );
        RR *= RR; // RR now changes
        //Fint_lam_d = Fint_lam.dotProduct( direction );


        // Check convergence
        OOFEM_LOG_INFO( "|resLS| = %.3e\n", fabs( RHSdefl ) );
        /*if ( fabs( RHSdefl ) < this->ls_tolerance || fabs( RHSdefl ) < resMult*fabs( rhs_0 ) ) {*/
        if ( fabs( RHSdefl ) < this->ls_tolerance ) {
            // X, ddX already updated
            OOFEM_LOG_INFO( "Eta = %.7e\n", Eta );
            return;
        } else if ( isnan( RHSdefl ) ) {
            break;
        }

        engngModel->updateComponent( tStep, NonLinearLhs, domain );
        k.times( direction, Ku );
        kdefl = Ku.dotProduct( direction );

         if ( deflation ) {
             gamma = p / ( dx_normSq + alph_defl * pow( dx_normSq, p ) );
             kdefl += gamma * RHS * dx_Defl.dotProduct( direction );
         }

        deta = RHS / kdefl;
        Eta += deta;
        //deltaXt_d = -Fint_lam_d / kdefl;

        //Eta = deltaLambda*deltaXt_d + deltaX__d;

        //this->computeDeltaLambda( deltaLambda, dX, deltaXt, deltaX_, RT, RR, Eta, deltaL, DeltaLambda, neq );

        

        //deta = deltaXt_d + deltaX__d * deltaLambda;
        //if ( isnan( deta ) ) {
        //    break;
        //}

        /*Eta += deta;*/

        //ddX = deta * direction;
        //dX += ddX;
        //X  += ddX;
    }
    // No update of X, dX, ddX, and lambdas
    ddX = ddXold;
    dX  = dXm1 + ddX;
    X   = Xinitial + dX;
    
    deltaLambda = deltaLambda0;
    DeltaLambda = DeltaLambdam1 + deltaLambda;
    Lambda      = ReachedLambda + DeltaLambda;
    

    OOFEM_LOG_INFO( "Eta = %.7e\n", 1. );

}

bool MaterialCylindricalALM ::checkConvergence( const FloatArray &R, const FloatArray *R0, const FloatArray &F,
    const FloatArray &X, const FloatArray &ddX,
    double Lambda, double RR0, double RR, double drProduct,
    const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange )
{
    /*
     * typedef std::set<DofID> __DofIDSet;
     * std::list<__DofIDSet> __ccDofGroups;
     * int nccdg; // number of Convergence Criteria Dof Groups
     */
    int _ng = nccdg;
    double forceErr, dispErr;
    FloatArray rhs; // residual of momentum balance eq (unbalanced nodal forces)
    FloatArray dg_forceErr( nccdg ), dg_dispErr( nccdg ), dg_totalLoadLevel( nccdg ), dg_totalDisp( nccdg );
    bool answer;
    EModelDefaultEquationNumbering dn;

    answer          = true;
    errorOutOfRange = false;

    // compute residual vector
    rhs = R;
    //rhs.times( Lambda );
    if ( R0 ) {
        rhs.add( *R0 );
    }
    //
    double normExternal = rhs.computeNorm();
    //
    rhs.subtract( F );

    //if ( _ng > 0 ) { // This is not treated in the modified version
    //    forceErr = dispErr = 0.0;
    //    // zero error norms per group
    //    dg_forceErr.zero();
    //    dg_dispErr.zero();
    //    dg_totalLoadLevel.zero();
    //    dg_totalDisp.zero();
    //    // loop over dof managers
    //    for ( auto &dman : domain->giveDofManagers() ) {
    //        if ( !dman->isLocal() ) {
    //            continue;
    //        }

    //        // loop over individual dofs
    //        for ( Dof *_idofptr : *dman ) {
    //            // loop over dof groups
    //            for ( int _dg = 1; _dg <= _ng; _dg++ ) {
    //                // test if dof ID is in active set
    //                if ( ccDofGroups.at( _dg - 1 ).find( _idofptr->giveDofID() ) != ccDofGroups.at( _dg - 1 ).end() ) {
    //                    int _eq = _idofptr->giveEquationNumber( dn );

    //                    if ( _eq ) {
    //                        continue;
    //                    }
    //                    dg_forceErr.at( _dg ) += rhs.at( _eq ) * rhs.at( _eq );
    //                    dg_dispErr.at( _dg ) += ddX.at( _eq ) * ddX.at( _eq );
    //                    // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
    //                    if ( R0 ) {
    //                        dg_totalLoadLevel.at( _dg ) += R0->at( _eq ) * R0->at( _eq );
    //                    }

    //                    dg_totalLoadLevel.at( _dg ) += R.at( _eq ) * R.at( _eq ) * Lambda * Lambda;
    //                    dg_totalDisp.at( _dg ) += X.at( _eq ) * X.at( _eq );
    //                }
    //            } // end loop over dof groups
    //        } // end loop over DOFs
    //    } // end loop over dof managers

    //    // loop over elements and their DOFs
    //    for ( auto &elem : domain->giveElements() ) {
    //        if ( elem->giveParallelMode() != Element_local ) {
    //            continue;
    //        }

    //        // loop over element internal Dofs
    //        for ( int _idofman = 1; _idofman <= elem->giveNumberOfInternalDofManagers(); _idofman++ ) {
    //            // loop over individual dofs
    //            for ( Dof *_idofptr : *elem->giveInternalDofManager( _idofman ) ) {
    //                // loop over dof groups
    //                for ( int _dg = 1; _dg <= _ng; _dg++ ) {
    //                    // test if dof ID is in active set
    //                    if ( ccDofGroups.at( _dg - 1 ).find( _idofptr->giveDofID() ) != ccDofGroups.at( _dg - 1 ).end() ) {
    //                        int _eq = _idofptr->giveEquationNumber( dn );

    //                        if ( _eq ) {
    //                            continue;
    //                        }
    //                        dg_forceErr.at( _dg ) += rhs.at( _eq ) * rhs.at( _eq );
    //                        dg_dispErr.at( _dg ) += ddX.at( _eq ) * ddX.at( _eq );
    //                        // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
    //                        if ( R0 ) {
    //                            dg_totalLoadLevel.at( _dg ) += R0->at( _eq ) * R0->at( _eq );
    //                        }

    //                        dg_totalLoadLevel.at( _dg ) += R.at( _eq ) * R.at( _eq ) * Lambda * Lambda;
    //                        dg_totalDisp.at( _dg ) += X.at( _eq ) * X.at( _eq );
    //                    }
    //                } // end loop over dof groups
    //            } // end loop over DOFs
    //        } // end loop over internal element dofmans
    //    } // end loop over elements

    //    // exchange individual partition contributions (simultaneously for all groups)
    //    FloatArray collectiveErr;
    //    parallel_context->accumulate( dg_forceErr, collectiveErr );
    //    dg_forceErr = collectiveErr;
    //    parallel_context->accumulate( dg_dispErr, collectiveErr );
    //    dg_dispErr = collectiveErr;
    //    parallel_context->accumulate( dg_totalLoadLevel, collectiveErr );
    //    dg_totalLoadLevel = collectiveErr;
    //    parallel_context->accumulate( dg_totalDisp, collectiveErr );
    //    dg_totalDisp = collectiveErr;

    //    OOFEM_LOG_INFO( "CALMLS:       %-15d %-15e ", nite, Lambda );
    //    // loop over dof groups
    //    for ( int _dg = 1; _dg <= _ng; _dg++ ) {
    //        //  compute a relative error norm
    //        if ( ( dg_totalLoadLevel.at( _dg ) ) < calm_SMALL_ERROR_NUM ) {
    //            dg_forceErr.at( _dg ) = sqrt( dg_forceErr.at( _dg ) );
    //        } else {
    //            dg_forceErr.at( _dg ) = sqrt( dg_forceErr.at( _dg ) / dg_totalLoadLevel.at( _dg ) );
    //        }

    //        //
    //        // compute displacement error
    //        //
    //        if ( dg_totalDisp.at( _dg ) < calm_SMALL_ERROR_NUM ) {
    //            dg_dispErr.at( _dg ) = sqrt( dg_dispErr.at( _dg ) );
    //        } else {
    //            dg_dispErr.at( _dg ) = sqrt( dg_dispErr.at( _dg ) / dg_totalDisp.at( _dg ) );
    //        }

    //        if ( ( fabs( dg_forceErr.at( _dg ) ) > rtolf.at( _dg ) * CALM_MAX_REL_ERROR_BOUND ) || ( fabs( dg_dispErr.at( _dg ) ) > rtold.at( _dg ) * CALM_MAX_REL_ERROR_BOUND ) || std ::isnan( dg_forceErr.at( _dg ) ) || std ::isnan( dg_dispErr.at( _dg ) ) ) {
    //            errorOutOfRange = true;
    //        }

    //        if ( ( fabs( dg_forceErr.at( _dg ) ) > rtolf.at( _dg ) ) || ( fabs( dg_dispErr.at( _dg ) ) > rtold.at( _dg ) ) ) {
    //            answer = false;
    //        }


    //        OOFEM_LOG_INFO( "%-15e %-15e ", dg_forceErr.at( _dg ), dg_dispErr.at( _dg ) );
    //    }

    //    OOFEM_LOG_INFO( "\n" );
    //} else {
        //
        // _ng==0 (errors computed for all dofs - this is the default)
        //

        //
        // compute force error(s)
        //
        double dXX;
        forceErr = parallel_context->localNorm( rhs );
        forceErr *= forceErr;
        dXX = parallel_context->localNorm( X );
        dXX *= dXX;

        double eNorm = internalForcesEBENorm.sum();
        // we compute a relative error norm
        //if ( ( RR0 + RR * Lambda * Lambda ) > calm_SMALL_ERROR_NUM ) {
        //    forceErr = sqrt( forceErr / ( RR0 + RR * Lambda * Lambda ) );
        //} else if ( eNorm > calm_SMALL_ERROR_NUM ) {
        //    forceErr = sqrt( forceErr / eNorm );
        //} else {
        //    forceErr = sqrt( forceErr );
        //}
        //////////////////////
        if ( ( normExternal ) > calm_SMALL_ERROR_NUM ) {
            forceErr = sqrt( forceErr / normExternal );
        //} else if ( eNorm > calm_SMALL_ERROR_NUM ) { // dont know this, so not considering
        //    forceErr = sqrt( forceErr / eNorm );
        } else {
            forceErr = sqrt( forceErr );
        }
        //////////////////////

        //
        // compute displacement error
        //
        if ( dXX < calm_SMALL_ERROR_NUM ) {
            dispErr = drProduct;
        } else {
            dispErr = drProduct / dXX;
            dispErr = sqrt( dispErr );
        }

        if ( ( fabs( forceErr ) > rtolf.at( 1 ) * CALM_MAX_REL_ERROR_BOUND ) || ( fabs( dispErr ) > rtold.at( 1 ) * CALM_MAX_REL_ERROR_BOUND ) || ( std ::isnan( forceErr ) ) || ( std ::isnan( dispErr ) ) ) {
            errorOutOfRange = true;
        }

        if ( ( fabs( forceErr ) > rtolf.at( 1 ) ) || ( fabs( dispErr ) > rtold.at( 1 ) ) ) {
            answer = false;
        }

        OOFEM_LOG_INFO( "CALMLS:       %-15d %-15e %-15e %-15e\n", nite, Lambda, forceErr, dispErr );
    //} // end default case (all dofs conributing)

    return answer;
}

////////////////////////////////////
// CylindricalALMStability
REGISTER_SparseNonLinearSystemNM( CylindricalALMStability )

CylindricalALMStability::CylindricalALMStability( Domain *d, EngngModel *m ) :
    CylindricalALM( d, m )
{
}

void
CylindricalALMStability::initializeFrom( InputRecord &ir )
{
    CylindricalALM::initializeFrom( ir );

    this->timeAsArcLength = ir.hasField( _IFT_CylindricalALMStability_timeAsArcLength );

    if ( this->timeAsArcLength ) {
        // deltaL is handed down by the TimeStepController from now on, so the
        // solver-owned step length read by the base class is not used. Say so,
        // because the base class still requires "steplength" on input.
        OOFEM_LOG_INFO( "CALMLS:       Arc-length increment taken from the solution step (deltaL = deltaT);\n"
                        "CALMLS:       'steplength'/'minsteplength' are ignored, use metastep 'deltat' with 'dtmax'/'dtmin'\n" );
    }

    // Optional DOF set for the constrained definiteness check. Deliberately independent
    // of hpcmode, so the branch can be traced with full ALM while the definiteness is
    // still evaluated on the constrained tangent.
    this->pdConstrainedDmanDofSrcArray.clear();
    IR_GIVE_OPTIONAL_FIELD( ir, this->pdConstrainedDmanDofSrcArray, _IFT_CylindricalALMStability_pdConstrainedDofs );
    if ( this->pdConstrainedDmanDofSrcArray.giveSize() ) {
        if ( ( this->pdConstrainedDmanDofSrcArray.giveSize() % 2 ) != 0 ) {
            OOFEM_ERROR( "'%s' must contain <node, nodeDof> pairs", _IFT_CylindricalALMStability_pdConstrainedDofs );
        }
        // Equation numbers are unknown at this point, defer the conversion to solve.
        this->pd_init = 1;
    }

    // Arc length at which the traversal direction is reversed; negative disables it.
    this->reverseAfterTime = -1.;
    IR_GIVE_OPTIONAL_FIELD( ir, this->reverseAfterTime, _IFT_CylindricalALMStability_reverseAfterTime );
    if ( this->reverseAfterTime >= 0. ) {
        OOFEM_LOG_INFO( "CALMLS:       Traversal direction will be reversed after arc length %e\n",
            this->reverseAfterTime );
    }
}


void CylindricalALMStability::convertPDConstrainedMap()
//
// Converts the user <node, nodeDof> pairs into equation numbers, in the same spirit as
// CylindricalALM::convertHPCMap, so that equation numbering stays hidden from the user.
//
{
    EModelDefaultEquationNumbering dn;
    int size = this->pdConstrainedDmanDofSrcArray.giveSize() / 2;
    IntArray found( size );
    int count = 0;

    for ( auto &dman : domain->giveDofManagers() ) {
        int jglobnum = dman->giveLabel();
        for ( int i = 1; i <= size; i++ ) {
            if ( this->pdConstrainedDmanDofSrcArray.at( 2 * i - 1 ) != jglobnum ) {
                continue;
            }
            if ( !parallel_context->isLocal( dman.get() ) ) {
                continue;
            }

            int idofid = this->pdConstrainedDmanDofSrcArray.at( 2 * i );
            int eq     = dman->giveDofWithID( idofid )->giveEquationNumber( dn );
            // A prescribed DOF carries no equation. Row/column 0 does not exist and would
            // be rejected by EigenMtrx::removeRowsAndColumns, so drop it with a warning
            // instead of aborting the analysis.
            if ( eq < 1 ) {
                OOFEM_WARNING( "node %d dof %d given in '%s' has no equation (prescribed?), ignored",
                    jglobnum, idofid, _IFT_CylindricalALMStability_pdConstrainedDofs );
                continue;
            }

            found.at( ++count ) = eq;
        }
    }

    this->pdConstrainedDofMask.resize( count );
    for ( int i = 1; i <= count; i++ ) {
        this->pdConstrainedDofMask.at( i ) = found.at( i );
    }

    if ( count != size ) {
        OOFEM_WARNING( "'%s': only %d of %d requested DOFs resolved to equations",
            _IFT_CylindricalALMStability_pdConstrainedDofs, count, size );
    }

    OOFEM_LOG_INFO( "CALMLS:       Definiteness check constrained on %d DOF(s) from '%s'\n",
        count, _IFT_CylindricalALMStability_pdConstrainedDofs );
}


const IntArray &
CylindricalALMStability::givePDConstrainedDofMask() const
{
    if ( this->pdConstrainedDofMask.giveSize() ) {
        return this->pdConstrainedDofMask;
    }

    // Fall back to the HPC mask, but only when the arc-length control actually uses it as a
    // constraint set. Reproduces the behaviour of input predating pdconstraineddofs.
    if ( calm_Control == calm_hpc_on ) {
        return calm_HPCIndirectDofMask;
    }

    static const IntArray emptyMask;
    return emptyMask;
}

ConvergedReason
CylindricalALMStability ::solve( SparseMtrx &k, FloatArray &R, FloatArray *R0,
    FloatArray &X, FloatArray &dX, FloatArray &F,
    const FloatArray &internalForcesEBENorm, double &ReachedLambda, referenceLoadInputModeType rlm,
    int &nite, TimeStep *tStep )
{
    const double basePerturbationAmplitude = 1.0e-4;
    double perturbationAmplitude           = basePerturbationAmplitude;
    int numtrials                          = 0;

    FloatArray rhs, deltaXt, deltaX_, dXm1, XInitial;
    FloatArray ddX; // total increment of displacements in iteration
    // double Bergan_k0 = 1.0, bk;
    double XX, RR, RR0, XR, p                                   = 0.0;
    double deltaLambda, Lambda, eta, DeltaLambdam1, DeltaLambda = 0.0;
    double drProduct = 0.0;
    int neq          = R.giveSize();
    int irest        = 0;
    int HPsize, i, ind;
    double _RR, _XX;
    ConvergedReason status = CR_UNKNOWN;
    bool converged, errorOutOfRangeFlag;

    // The arc-length increment is the step time increment, so that the metastep
    // deltaT is the single lever controlling how far this step advances and
    // pseudo-time is the cumulated arc length. See timeAsArcLength.
    if ( this->timeAsArcLength ) {
        double dT = tStep->giveTimeIncrement();
        if ( dT <= 0. ) {
            OOFEM_ERROR( "non-positive time increment (%e) cannot be used as arc-length increment", dT );
        }

        deltaL = dT;
    }

    // print iteration header

    OOFEM_LOG_INFO( "CALMLS:       Initial step length: %-15e\n", deltaL );
    if ( nccdg == 0 ) {
        OOFEM_LOG_INFO( "CALMLS:       Iteration       LoadLevel       ForceError      DisplError    \n" );
        OOFEM_LOG_INFO( "----------------------------------------------------------------------------\n" );
    } else {
        OOFEM_LOG_INFO( "Iter  LoadLevel       " );
        for ( i = 1; i <= nccdg; i++ ) {
            OOFEM_LOG_INFO( "ForceError(%02d)  DisplError(%02d)  ", i, i );
        }

        OOFEM_LOG_INFO( "\n__________________________________________________________\n" );
    }

    //
    // Now smarter method is default:
    // after convergence troubles for (calm_NR_ModeTick subsequent
    // steps newly set calm_NR_Mode will be used. After the calm_NR_OldMode will be restored
    //
    if ( calm_NR_ModeTick == 0 ) {
        calm_NR_Mode = calm_NR_OldMode;
    }

    if ( calm_NR_ModeTick > 0 ) {
        calm_NR_ModeTick--;
    }

    XInitial = X;
    ddX.resize( neq );
    ddX.zero();
    deltaXt.resize( neq );
    deltaXt.zero();

    status = CR_UNKNOWN;
    this->giveLinearSolver();

    // A stability-capable linear solver can check the definiteness of the
    // converged tangent stiffness matrix.
    BifurcationInterface *stabSolver = dynamic_cast<BifurcationInterface *>( this->linSolver.get() );
    if ( stabSolver ) {
        // Discard the definiteness recorded for the previous step, or for a previous
        // attempt of this step, so that it is never reported for the current one.
        stabSolver->resetPDStatus();
    }

    // create HPC Map if needed
    if ( calm_hpc_init ) {
        this->convertHPCMap();
        calm_hpc_init = 0;
    }

    // create the constrained definiteness-check map if needed
    if ( this->pd_init ) {
        this->convertPDConstrainedMap();
        this->pd_init = 0;
    }

    if ( R0 ) {
        RR0 = parallel_context->localNorm( *R0 );
        RR0 *= RR0;
    } else {
        RR0 = 0.0;
    }

    // Decide the orientation of this step. The state we start from is
    // targetTime - timeIncrement, which TimeStepController::reduceTimeStep keeps
    // invariant, so the trigger is stable when a step is retried with a smaller
    // increment. Deliberately placed before the restart label: the direction is a
    // per-step decision and must not be re-evaluated by a mid-step goto restart.
    {
        double sNow = tStep->giveTargetTime() - tStep->giveTimeIncrement();
        double aNow = fabs( sNow );
        // sNow is an accumulated sum, so compare with a relative tolerance, otherwise the
        // threshold can be missed by an ULP.
        double tol    = 1.e-10 * ( aNow > 1. ? aNow : 1. );
        double newDir = ( this->reverseAfterTime >= 0. && sNow >= this->reverseAfterTime - tol ) ? -1. : 1.;

        if ( newDir != this->pathDirection ) {
            this->pathDirection = newDir;
            // Drop the accumulated direction so that the predictor is re-seeded and the new
            // orientation takes effect on this step rather than one step later.
            old_dX.clear();
            OOFEM_LOG_INFO( "CALMLS:       Reversing traversal direction at arc length %e\n", sNow );
        }
    }

    //
    // A  initial step (predictor)
    //
    //

    //
    // A.2.   We assume positive-definite (0)Kt (tangent stiffness mtrx).
    //
    RR = parallel_context->localNorm( R );
    RR *= RR;

restart:
    //
    // A.1. calculation of (0)VarRt
    //
    dX.zero();
    // engngModel->updateComponent(tStep, InternalRhs, domain); // By not updating this, one obtains the old equilibrated tangent.
    engngModel->updateComponent( tStep, NonLinearLhs, domain );
    linSolver->solve( k, R, deltaXt );

    // If desired by the user, the solution is (slightly) perturbed, so that various symmetries can be broken.
    // This is useful e.g. to trigger localization in a homogeneous material under uniform stress without
    // the need to introduce material imperfections. The problem itself remains symmetric but the iterative
    // solution is brought to a nonsymmetric state and it gets a chance to converge to a nonsymmetric solution.
    // Parameters of the perturbation technique are specified by the user and by default no perturbation is done.
    // Milan Jirasek
    SparseNonLinearSystemNM ::applyPerturbation( &deltaXt );

    if ( calm_Control == calm_hpc_off ) {
        XX = parallel_context->localNorm( deltaXt );
        XX *= XX;
        p = sqrt( XX + Psi * Psi * RR );
    } else if ( calm_Control == calm_hpc_on ) {
        HPsize = calm_HPCIndirectDofMask.giveSize();
        _XX    = 0;
        _RR    = 0;
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at( i ) ) != 0 ) {
                _XX += deltaXt.at( ind ) * deltaXt.at( ind );
                _RR += R.at( ind ) * R.at( ind );
            }
        }

        // In case of paralllel analysis:
        FloatArray collected_XXRR;
        parallel_context->accumulate( { _XX, _RR }, collected_XXRR );
        _XX = collected_XXRR( 0 );
        _RR = collected_XXRR( 1 );

        p = sqrt( _XX + Psi * Psi * _RR );
    } else if ( calm_Control == calml_hpc ) {
        HPsize = calm_HPCIndirectDofMask.giveSize();
        p      = 0.;
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at( i ) ) != 0 ) {
                p += deltaXt.at( ind ) * calm_HPCWeights.at( i );
            }
        }

        // In case of paralllel analysis:
        p = parallel_context->accumulate( p );
    }

    Lambda = ReachedLambda;
    if ( rootselectiontype == RST_Cos ) {
        /* XR is unscaled Bergan's param of current stiffness XR = deltaXt^T k deltaXt
         * this is used to test whether k has negative or positive slope */
        XR = this->pathDirection * parallel_context->localDotProduct( deltaXt, R );
    } else {
        if ( old_dX.giveSize() ) {
            XR = parallel_context->localDotProduct( deltaXt, old_dX );
        } else {
            XR = this->pathDirection * parallel_context->localDotProduct( deltaXt, R );
        }
    }
    DeltaLambda = deltaLambda = sgn( XR ) * deltaL / p;
    Lambda += DeltaLambda;
    //
    // A.3.
    //

    rhs = R;
    rhs.times( DeltaLambda );
    if ( R0 ) {
        rhs.add( *R0 );
    }
    linSolver->solve( k, rhs, dX );
    X.add( dX );

    nite = 0;

    // update solution state counter
    tStep->incrementStateCounter();
    engngModel->updateComponent( tStep, InternalRhs, domain );

    do {
        nite++;
        tStep->incrementSubStepNumber();

        dXm1          = dX;
        DeltaLambdam1 = DeltaLambda;

        //
        // B  - iteration MNRM is used
        //
        // B.1. is ommited because MNRM is used instead of NRM.
        //
        if ( ( calm_NR_Mode == calm_fullNRM ) || ( ( calm_NR_Mode == calm_accelNRM ) && ( nite % calm_MANRMSteps == 0 ) ) ) {
            //
            // ALM with full NRM
            //
            // we assemble new tangent stiffness and compute new deltaXt
            // internal state of elements is updated by previous calling
            // InternalRhs
            //
            engngModel->updateComponent( tStep, NonLinearLhs, domain );
            //
            // compute deltaXt for i-th iteration
            //
            linSolver->solve( k, R, deltaXt );
        }

        // B.2.
        //

        rhs = R;
        rhs.times( Lambda );
        if ( R0 ) {
            rhs.add( *R0 );
        }

        rhs.subtract( F );
        deltaX_.resize( neq );
        linSolver->solve( k, rhs, deltaX_ );
        eta = 1.0;
        //
        // B.3.
        //
        if ( this->computeDeltaLambda( deltaLambda, dX, deltaXt, deltaX_, R, RR, eta, deltaL, DeltaLambda, neq ) ) {
            irest++;
            if ( irest <= maxRestarts ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length
                deltaL = deltaL * CALM_RESET_STEP_REDUCE;
                if ( deltaL < minStepLength ) {
                    deltaL = minStepLength;
                }

                // restore previous total displacement vector
                X = XInitial;
                // reset all changes fro previous equilibrium state
                dX.zero();
                // restore initial stiffness
                engngModel->updateComponent( tStep, NonLinearLhs, domain );

                // the perturbation trial count/amplitude is tied to the (now reset) displacement state
                numtrials             = 0;
                perturbationAmplitude = basePerturbationAmplitude;

                OOFEM_LOG_INFO( "CALMLS:       Iteration Reset ...\n" );

                calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = CR_DIVERGED_ITS;
                OOFEM_ERROR( "can't continue further" );
            }
        }

        //
        // B.4.+ B.5.
        //


        if ( this->lsFlag && ( nite != 1 ) ) {
            //
            //  LINE SEARCH
            //
            this->do_lineSearch( X, XInitial, deltaX_, deltaXt,
                dXm1, dX, ddX,
                R, R0, F,
                DeltaLambda, DeltaLambdam1, deltaLambda, Lambda, ReachedLambda,
                RR, drProduct, tStep );
        } else { // no line search
            //
            // update solution vectors
            //
            ddX.clear();
            ddX.add( eta * deltaLambda, deltaXt );
            ddX.add( eta, deltaX_ );
            dX = dXm1;
            dX.add( ddX );
            X = XInitial;
            X.add( dX );

            drProduct = parallel_context->localNorm( ddX );
            drProduct *= drProduct;

            tStep->incrementStateCounter(); // update solution state counter
            //
            // B.6.
            //
            DeltaLambda = DeltaLambdam1 + deltaLambda;
            Lambda      = ReachedLambda + DeltaLambda;

            tStep->incrementStateCounter(); // update solution state counter
            engngModel->updateComponent( tStep, InternalRhs, domain );
        }

        //
        // B.7.
        //
        // convergence check
        //

        converged = this->checkConvergence( R, R0, F, X, ddX, Lambda, RR0, RR, drProduct,
            internalForcesEBENorm, nite, errorOutOfRangeFlag );

        // Check the tangent stiffness at the converged equilibrium state. Do
        // not perturb the solution here; that will be handled separately.
        bool doPerturbation = false;
        if ( converged && stabSolver && stabSolver->canComputeEigenVectors() ) {
            engngModel->updateComponent( tStep, NonLinearLhs, domain );
            bool isPD = stabSolver->checkPD( k );

            ////////////////////////
            //// Independently verify the negative result with the smallest
            //// eigenpair computed by the stability solver.
            //FloatArray eigenvalues;
            //int nconv = stabSolver->computeEigenValues( k, eigenvalues );
            //if ( nconv > 0 ) {
            //    bool isPD_check = eigenvalues.at( 1 ) > 0.0;
            //    if ( isPD_check != isPD ) {
            //        OOFEM_WARNING( "Positive-definiteness checks disagree" );
            //    }
            //} else {
            //    OOFEM_WARNING( "Could not compute eigenpairs for positive-definiteness check" );
            //}
            ////////////////////////

            if ( !isPD ) {
                //// Independently verify the negative result with the smallest
                //// eigenvalue computed by the stability solver.
                //FloatArray eigenvalues;
                //int nconv = stabSolver->computeEigenValues( k, eigenvalues );
                //if ( nconv > 0 ) {
                //    if ( eigenvalues.at( 1 ) > 0.0 ) {
                //        OOFEM_WARNING( "Positive-definiteness checks disagree" );
                //    }
                //} else {
                //    OOFEM_WARNING( "Could not compute eigenpairs for positive-definiteness check" );
                //}

                OOFEM_LOG_INFO( "Negative eigenvalue\n" );
            }
        }
        auto *eigenK = dynamic_cast<EigenMtrx *>( &k );
        if ( converged && stabSolver && eigenK ) {
            engngModel->updateComponent( tStep, NonLinearLhs, domain );
            // Make a copy of the stiffness matrix
            EigenMtrx Kcopy( *eigenK );
            // Remove the prescribed DOFs corresponding to the HPC mask
            // Constrain the check whenever a DOF set is available, see givePDConstrainedDofMask.
            const IntArray &pdMask = this->givePDConstrainedDofMask();
            if ( pdMask.giveSize() ) {
                std::unique_ptr<EigenMtrx> Kdel = Kcopy.removeRowsAndColumns( pdMask );

                ////////////////////
                // std::string filename = "C:\\Users\\msmej\\OneDrive - �esk� vysok� u�en� technick� v Praze\\cvut\\Doktor�t\\ERC project\\OOFEM\\vypocty\\K.txt";
                // std::string filenameDel = "C:\\Users\\msmej\\OneDrive - �esk� vysok� u�en� technick� v Praze\\cvut\\Doktor�t\\ERC project\\OOFEM\\vypocty\\Kdel.txt";
                //eigenK->printToFile( filename );
                //Kdel->printToFile( filenameDel );
                ////////////////////

                bool isPD2 = stabSolver->checkPD( *Kdel );
                if ( !isPD2 ) {

                    OOFEM_LOG_INFO( "Negative eigenvalue of the constrained matrix\n" );

                    if (doPerturbation) {
                        FloatArray eigenvalues;
                        FloatMatrix eigenvectors;
                        int nconv = stabSolver->computeEigenValuesVectors( *Kdel, eigenvalues, eigenvectors );
                        FloatArray direction; // get eigenvectors
                        

                        if ( nconv > 0 ) {
                            eigenvectors.copyColumn( direction, 1 );

                            if ( eigenvalues.at( 1 ) > 0.0 ) {
                                OOFEM_WARNING( "Positive-definiteness checks disagree" );
                            } else {
                                // Reconstruct the full eigenvector by inserting zeros
                                FloatArray fullEigenvector;
                                int totalNrOfEquations = X.giveSize();
                                fullEigenvector.resize( totalNrOfEquations );
                                fullEigenvector.zero();

                                std::vector<bool> removed( totalNrOfEquations + 1, false );
                                for ( int i = 1; i <= pdMask.giveSize(); ++i ) {
                                    removed[pdMask.at( i )] = true;
                                }

                                int redEq = 1;
                                for ( int eq = 1; eq <= totalNrOfEquations; ++eq ) {
                                    if ( !removed[eq] ) {
                                        fullEigenvector.at( eq ) = direction.at( redEq++ );
                                    }
                                }

                                // Normalize
                                fullEigenvector.normalize();

                                // Perturb current displacement vector
                                perturbationAmplitude = basePerturbationAmplitude * ( 1.5 * numtrials + 1 );
                                numtrials++;

                                OOFEM_LOG_INFO(
                                    "Perturbing along eigenmode: "
                                    "lambda_min = %.8e, "
                                    "amplitude = %.8e\n",
                                    eigenvalues.at( 1 ),
                                    perturbationAmplitude );


                                X.add( perturbationAmplitude, fullEigenvector );

                                OOFEM_LOG_INFO( "Perturbed solution along the unstable eigenmode; "
                                                 "restarting as a new arc-length step from the perturbed state.\n" );

                                // Correcting at fixed Lambda would solve against the (near-)singular
                                // tangent right at the bifurcation point, and folding the (large) kick
                                // into the current step's dX would blow past the arc-length radius deltaL
                                // and trip the negative-discriminant restart below, undoing the perturbation.
                                // Instead, treat the perturbed configuration as the baseline for a fresh
                                // arc-length predictor step, so the arc-length constraint stays in control.
                                XInitial = X;
                                dX.zero();
                                ddX.zero();
                                engngModel->updateComponent( tStep, InternalRhs, domain );

                                calm_NR_OldMode  = calm_NR_Mode;
                                calm_NR_Mode     = calm_fullNRM;
                                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                                goto restart;
                            }
                        } else {
                            OOFEM_WARNING( "Could not compute eigenpairs for positive-definiteness check" );
                        }

                    }
                    
                }
            }
        }

        if ( ( nite >= nsmax ) || errorOutOfRangeFlag ) {
            if ( this->timeAsArcLength ) {
                // deltaL is owned by the TimeStepController, so the solver must not
                // shrink it behind the controller's back - that would desynchronise
                // deltaL from tStep->giveTimeIncrement() and hence corrupt the
                // cumulated arc length carried by targetTime. Hand control back
                // instead: MetaStep::reduceTimeStep will shorten deltaT and
                // EngngModel::solveYourself will re-run the step.
                X = XInitial;
                dX.zero();

                throw ConvergenceException( errorOutOfRangeFlag
                        ? "Error out of range in arc-length step"
                        : "Maximum number of iterations reached without convergence in arc-length step" );
            }

            irest++;
            if ( irest <= maxRestarts ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length
                deltaL = deltaL * CALM_RESET_STEP_REDUCE;
                if ( deltaL < minStepLength ) {
                    deltaL = minStepLength;
                }

                // restore previous total displacement vector
                X = XInitial;
                // reset all changes from previous equilibrium state
                engngModel->initStepIncrements();
                dX.zero();
                // restore initial stiffness
                engngModel->updateComponent( tStep, NonLinearLhs, domain );

                // the perturbation trial count/amplitude is tied to the (now reset) displacement state
                numtrials             = 0;
                perturbationAmplitude = basePerturbationAmplitude;

                OOFEM_LOG_INFO( "CALMLS:       Iteration Reset ...\n" );

                calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = CR_DIVERGED_ITS;
                OOFEM_WARNING( "Convergence not reached after %d iterations", nsmax );
                // exit(1);
                break;
            }
        }

        // output of per iteration data
        engngModel->giveExportModuleManager()->doOutput( tStep, true );
    } while ( !converged || ( nite < minIterations ) );

    //
    // update dofs,nodes,Elemms and print result
    //
#ifdef VERBOSE
    // printf ("\nCALM - step iteration finished") ;
#endif

    // Seems to be completely unused / Mikael
#if 0
    /* compute Bergan's parameter of current stiffness */
    double RDR = parallel_context->localDotProduct(* R, * dX);
    double DR = parallel_context->localNorm(* dX);
    bk = DeltaLambda * RDR / ( DR * DR );

    if ( tStep->giveNumber() == 1 ) {
        // compute Bergan_k0
        Bergan_k0 = bk;
    } else {
        bk = bk / Bergan_k0;
        //if (fabs(bk) < CALM_TANGENT_STIFF_TRESHOLD) status |= NM_KeepTangent;
    }
#endif

    /* we have computed Bergan's  parameter -> you can adjust deltaL according
     * to this value */
    //
    // there has been restart already - set nite to maxiter
    //
    // if (irest > 0) nite = nsmax;

    // Post-convergence adaptation of the step length. When deltaL is the step time
    // increment this is MetaStep::adaptTimeStep's job (driven by minrequirediter /
    // maxrequirediter and clamped by dtmin / dtmax), so the solver keeps its hands off.
    if ( !this->timeAsArcLength ) {
        if ( nite > numberOfRequiredIterations ) {
            deltaL = deltaL * numberOfRequiredIterations / nite;
        } else {
            deltaL = deltaL * sqrt( sqrt( (double)numberOfRequiredIterations / (double)nite ) );
        }

        if ( deltaL > maxStepLength ) {
            deltaL = maxStepLength;
        }

        if ( deltaL < minStepLength ) {
            deltaL = minStepLength;
        }

        OOFEM_LOG_INFO( "CALMLS:       Adjusted step length: %-15e\n", deltaL );
    }

    // Do not report a step that gave up as converged: the caller stores this in
    // tStep->convergedReason and would otherwise accept a non-equilibrated state.
    if ( status != CR_DIVERGED_ITS ) {
        status = CR_CONVERGED;
        solved = 1;
    }

    ReachedLambda = Lambda;
    old_dX        = dX;
    return status;
}





} // end namespace oofem
