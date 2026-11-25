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

#include "tm/EngineeringModels/nltransienttransportproblem.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "verbose.h"
#include "tm/Elements/transportelement.h"
#include "classfactory.h"
#include "mathfem.h"
#include "assemblercallback.h"
#include "unknownnumberingscheme.h"
#include <algorithm>
#include <cmath>

namespace oofem {
REGISTER_EngngModel(NLTransientTransportProblem);

namespace {
class NewmarkWaveLhsAssembler : public MatrixAssembler
{
    double beta;
    double gamma;
    double damping;
public:
    NewmarkWaveLhsAssembler(double beta, double gamma, double damping) : beta(beta), gamma(gamma), damping(damping) { }

    void matrixFromElement(FloatMatrix &answer, Element &element, TimeStep *tStep) const override
    {
        FloatMatrix capacity;
        element.giveCharacteristicMatrix(answer, TangentStiffnessMatrix, tStep);
        element.giveCharacteristicMatrix(capacity, MassMatrix, tStep);
        const double dt = tStep->giveTimeIncrement();
        const bool dynamicBoundary = element.giveGeometryType() == EGT_line_1;
        const double factor = dynamicBoundary
            ? gamma / ( beta * dt )
            : 1.0 / ( beta * dt * dt ) + damping * gamma / ( beta * dt );
        answer.add(factor, capacity);
    }
};

class BulkWaveMassAssembler : public MatrixAssembler
{
public:
    void matrixFromElement(FloatMatrix &answer, Element &element, TimeStep *tStep) const override
    {
        if ( element.giveGeometryType() == EGT_line_1 ) {
            answer.clear(); // boundary equation is first order, not inertial
        } else {
            element.giveCharacteristicMatrix(answer, MassMatrix, tStep);
        }
    }
};
}

NLTransientTransportProblem :: NLTransientTransportProblem(int i, EngngModel *_master = nullptr) : NonStationaryTransportProblem(i, _master)
{
}

void
NLTransientTransportProblem :: initializeFrom(InputRecord &ir)
{
     NonStationaryTransportProblem :: initializeFrom(ir);

    int val = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_NLTransientTransportProblem_nsmax);
    nsmax = val;

    IR_GIVE_FIELD(ir, rtol, _IFT_NLTransientTransportProblem_rtol);

    NR_Mode = nrsolverModifiedNRM;
    MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, MANRMSteps, _IFT_NLTransientTransportProblem_manrmsteps);
    IR_GIVE_OPTIONAL_FIELD(ir, bdfOrder, _IFT_NLTransientTransportProblem_bdforder);
    IR_GIVE_OPTIONAL_FIELD(ir, waveMode, _IFT_NLTransientTransportProblem_wave);
    IR_GIVE_OPTIONAL_FIELD(ir, newmarkBeta, _IFT_NLTransientTransportProblem_newmarkbeta);
    IR_GIVE_OPTIONAL_FIELD(ir, newmarkGamma, _IFT_NLTransientTransportProblem_newmarkgamma);
    IR_GIVE_OPTIONAL_FIELD(ir, initialVelocity, _IFT_NLTransientTransportProblem_initialvelocity);
    IR_GIVE_OPTIONAL_FIELD(ir, initialVelocityValues, _IFT_NLTransientTransportProblem_initialvelocityvalues);
    IR_GIVE_OPTIONAL_FIELD(ir, waveDamping, _IFT_NLTransientTransportProblem_wavedamping);
    if ( bdfOrder < 1 || bdfOrder > 5 ) {
        throw ValueInputException(ir, _IFT_NLTransientTransportProblem_bdforder,
                                  "BDF order must be between 1 and 5");
    }
    if ( bdfOrder > 1 && alpha != 1.0 ) {
        throw ValueInputException(ir, _IFT_NonStationaryTransportProblem_alpha,
                                  "BDF integration requires alpha 1.0");
    }
    if ( waveMode ) {
        if ( alpha != 1.0 ) {
            throw ValueInputException(ir, _IFT_NonStationaryTransportProblem_alpha,
                                      "Newmark wave integration requires alpha 1.0");
        }
        if ( bdfOrder > 1 ) {
            throw ValueInputException(ir, _IFT_NLTransientTransportProblem_bdforder,
                                      "BDF integration cannot be combined with Newmark wave integration");
        }
        if ( newmarkBeta <= 0.0 || newmarkGamma <= 0.0 ) {
            throw ValueInputException(ir, _IFT_NLTransientTransportProblem_newmarkbeta,
                                      "Newmark beta and gamma must be positive");
        }
    }
    if ( bdfOrder > 1 && ( dtFunction != 0 || discreteTimes.giveSize() > 0 ) ) {
        throw ValueInputException(ir, _IFT_NonStationaryTransportProblem_deltat,
                                  "BDF integration requires a fixed time increment");
    }
    if ( MANRMSteps > 0 ) {
        NR_Mode = nrsolverAccelNRM;
    } else {
        NR_Mode = nrsolverModifiedNRM;
    }
    if ( waveMode ) {
        // The Newmark effective tangent depends on the current nonlinear
        // state through the Allen-Cahn reaction term, so refresh it each iteration.
        NR_Mode = nrsolverFullNRM;
    }
}

void
NLTransientTransportProblem :: giveBdfCoefficients(int order, std::vector<double> &coefficients)
{
    static const std::vector< std::vector<double> > table = {
        {},
        { 1.0, -1.0 },
        { 1.5, -2.0, 0.5 },
        { 11.0 / 6.0, -3.0, 1.5, -1.0 / 3.0 },
        { 25.0 / 12.0, -4.0, 3.0, -4.0 / 3.0, 0.25 },
        { 137.0 / 60.0, -5.0, 5.0, -10.0 / 3.0, 1.25, -0.2 }
    };
    if ( order < 1 || order > 5 ) {
        OOFEM_ERROR("invalid BDF order %d", order);
    }
    coefficients = table[order];
}

void
NLTransientTransportProblem :: captureBdfElementState(TimeStep *tStep)
{
    std::vector<FloatArray> state;
    Domain *domain = this->giveDomain(1);
    state.reserve(domain->giveNumberOfElements());
    for ( auto &element : domain->giveElements() ) {
        FloatArray local;
        element->computeVectorOf(VM_Total, tStep, local);
        state.emplace_back(std::move(local));
    }
    bdfElementHistory.push_front(std::move(state));
    while ( bdfElementHistory.size() > 5 ) {
        bdfElementHistory.pop_back();
    }
}

TimeStep *
NLTransientTransportProblem :: giveNextStep()
{
    double intrinsicTime;
    NonStationaryTransportProblem :: giveNextStep();
    if ( bdfOrder > 1 && bdfPreviousIncrement > 0.0 ) {
        const double currentIncrement = currentStep->giveTimeIncrement();
        if ( std::abs(currentIncrement - bdfPreviousIncrement) >
             1.e-12 * std::max(currentIncrement, bdfPreviousIncrement) ) {
            OOFEM_ERROR("BDF integration requires a uniform time step; previous dt=%g, current dt=%g",
                        bdfPreviousIncrement, currentIncrement);
        }
    }
    intrinsicTime = previousStep->giveTargetTime() + this->alpha*(currentStep->giveTargetTime()-previousStep->giveTargetTime());
    currentStep->setIntrinsicTime(intrinsicTime);
    return currentStep.get();
}


void NLTransientTransportProblem :: solveYourselfAt(TimeStep *tStep)
{
    // creates system of governing eq's and solves them at given time step
    // first assemble problem at current time step

    // Right hand side
    FloatArray rhs;
    double solutionErr, incrementErr;
    int neq =  this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif
    //Delete lhs matrix and create a new one. This is necessary due to growing/decreasing number of equations.
    if ( tStep->isTheFirstStep() || this->changingProblemSize ) {

        conductivityMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        if ( !conductivityMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        conductivityMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling conductivity and capacity matrices\n");
#endif
    }

    //create previous solution from IC or from previous tStep
    if ( tStep->isTheFirstStep() ) {
        if ( !stepWhenIcApply ) {
            stepWhenIcApply = std::make_unique<TimeStep>( *tStep->givePreviousStep() );
        }
        this->applyIC(stepWhenIcApply.get()); //insert solution to hash=1(previous), if changes in equation numbering
    }

    double dTTau = tStep->giveTimeIncrement();
    double Tau = tStep->giveTargetTime() - ( 1. - alpha ) * tStep->giveTimeIncrement();
    //Time step in which material laws are taken into account
    TimeStep TauStep(tStep->giveNumber(), this, tStep->giveMetaStepNumber(), Tau, dTTau, tStep->giveSolutionStateCounter() + 1);

    //Predictor
    FloatArray *solutionVector;
    UnknownsField->advanceSolution(tStep);
    solutionVector = UnknownsField->giveSolutionVector(tStep);

    //Initialize and give solutionVector from previous solution
    if ( changingProblemSize ) {
        if ( !tStep->isTheFirstStep() ) {
            //copy recent solution to previous position, copy from hash=0 to hash=1(previous)
            copyUnknownsInDictionary( VM_Total, tStep, tStep->givePreviousStep() );
        }

        UnknownsField->initialize( VM_Total, tStep->givePreviousStep(), *solutionVector, EModelDefaultEquationNumbering() );
    } else {
        //copy previous solution vector to actual
        *solutionVector = *UnknownsField->giveSolutionVector( tStep->givePreviousStep() );
    }

    if ( waveMode ) {
        const FloatArray *previousSolution = UnknownsField->giveSolutionVector(tStep->givePreviousStep());
        const double dt = tStep->giveTimeIncrement();
        const double predictorAccelerationFactor = dt * dt * ( 0.5 - newmarkBeta );
        for ( int eq = 1; eq <= neq; ++eq ) {
            solutionVector->at(eq) = previousSolution->at(eq) + dt * waveVelocity.at(eq)
                + predictorAccelerationFactor * waveAcceleration.at(eq);
        }
    }

    this->updateInternalState(& TauStep); //insert to hash=0(current), if changes in equation numbering

    FloatArray solutionVectorIncrement(neq);
    int nite = 0;

    OOFEM_LOG_INFO("Time            Iter       ResidNorm       IncrNorm\n__________________________________________________________\n");


    do {
        nite++;

        // Corrector
#ifdef VERBOSE
        // printf("\nAssembling conductivity and capacity matrices");
#endif

        if ( ( nite == 1 ) || ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
            conductivityMatrix->zero();
            if ( waveMode ) {
                this->assemble(*conductivityMatrix, &TauStep,
                               NewmarkWaveLhsAssembler(newmarkBeta, newmarkGamma, waveDamping),
                               EModelDefaultEquationNumbering(), this->giveDomain(1));
                this->assemble(*conductivityMatrix, &TauStep, IntSourceLHSAssembler(),
                               EModelDefaultEquationNumbering(), this->giveDomain(1));
            } else {
                //Assembling left hand side - start with conductivity matrix
                this->assemble( *conductivityMatrix, & TauStep, IntSourceLHSAssembler(),
                               EModelDefaultEquationNumbering(), this->giveDomain(1) );
                conductivityMatrix->times(alpha);
                //Add capacity matrix
                std::vector<double> coefficients;
                const int activeOrder = bdfOrder == 1 ? 1 :
                    std::min(bdfOrder, static_cast<int>(bdfElementHistory.size()));
                giveBdfCoefficients(activeOrder, coefficients);
                this->assemble( *conductivityMatrix, & TauStep,
                               MidpointLhsAssembler(lumpedCapacityStab, alpha, coefficients[0]),
                               EModelDefaultEquationNumbering(), this->giveDomain(1) );
            }
        }

        rhs.resize(neq);
        rhs.zero();
        //edge or surface load on element
        //add internal source vector on elements
        this->assembleVectorFromElements( rhs, tStep, TransportExternalForceAssembler(), VM_Total,
                                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //add nodal load
        this->assembleVectorFromDofManagers( rhs, tStep, ExternalForceAssembler(), VM_Total,
                                            EModelDefaultEquationNumbering(), this->giveDomain(1) );

        // subtract the rhs part depending on previous solution
        assembleAlgorithmicPartOfRhs(rhs, EModelDefaultEquationNumbering(), tStep);
        // set-up numerical model
        this->giveNumericalMethod( this->giveCurrentMetaStep() );

        // call numerical model to solve arised problem
#ifdef VERBOSE
        //OOFEM_LOG_INFO("Solving ...\n");
#endif

        // compute norm of residuals from balance equations
        solutionErr = rhs.computeNorm();

        if ( !std::isfinite(solutionErr) ) {
            OOFEM_ERROR("non-finite nonlinear transport residual at t=%g, iteration %d", tStep->giveTargetTime(), nite);
        }

        linSolver->solve(*conductivityMatrix, rhs, solutionVectorIncrement);
        for ( double value : solutionVectorIncrement ) {
            if ( !std::isfinite(value) ) {
                OOFEM_ERROR("non-finite nonlinear transport increment at t=%g, iteration %d", tStep->giveTargetTime(), nite);
            }
        }
        solutionVector->add(solutionVectorIncrement);
        this->updateInternalState(tStep); //insert to hash=0(current), if changes in equation numbering
        // compute error in the solutionvector increment
        incrementErr = solutionVectorIncrement.computeNorm();
        if ( !std::isfinite(incrementErr) ) {
            OOFEM_ERROR("non-finite nonlinear transport increment norm at t=%g, iteration %d", tStep->giveTargetTime(), nite);
        }

        // update solution state counter
        TauStep.incrementStateCounter();
        tStep->incrementStateCounter();

        OOFEM_LOG_INFO("%-15e %-10d %-15e %-15e\n", tStep->giveTargetTime(), nite, solutionErr, incrementErr);

        currentIterations = nite;

        if ( nite >= nsmax ) {
            OOFEM_ERROR("convergence not reached after %d iterations", nsmax);
        }
    } while ( ( fabs(solutionErr) > rtol ) || ( fabs(incrementErr) > rtol ) );

    // The time-step controller adapts the next increment using this count.
    tStep->numberOfIterations = currentIterations;
}


double NLTransientTransportProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation
// This function translates this request to numerical method language
{
    if ( waveMode && ( mode == VM_Velocity || mode == VM_Acceleration ) ) {
        const int eq = dof->__giveEquationNumber();
        if ( eq > 0 && eq <= waveVelocity.giveSize() ) {
            return mode == VM_Velocity ? waveVelocity.at(eq) : waveAcceleration.at(eq);
        }
        if ( eq > 0 ) {
            return mode == VM_Velocity ? initialVelocity : 0.0;
        }
        return dof->hasBc(tStep) ? dof->giveBcValue(mode, tStep) : 0.0;
    }

    if ( this->requiresUnknownsDictionaryUpdate() ) {
        if ( mode == VM_Incremental ) { //get difference between current and previous time variable
            return dof->giveUnknowns()->at(0) - dof->giveUnknowns()->at(1);
        } else if ( mode == VM_TotalIntrinsic ) { // intrinsic value only for current step
            return this->alpha * dof->giveUnknowns()->at(0) + (1.-this->alpha) * dof->giveUnknowns()->at(1);
        }
        int hash = this->giveUnknownDictHashIndx(mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else {
            OOFEM_ERROR("Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    }

    double t = tStep->giveTargetTime();
    TimeStep *previousStep = this->givePreviousStep(), *currentStep = this->giveCurrentStep();

    if ( dof->__giveEquationNumber() == 0 ) {
        OOFEM_ERROR("invalid equation number on DoF %d", dof->giveDofID() );
    }

    if ( ( t >= previousStep->giveTargetTime() ) && ( t <= currentStep->giveTargetTime() ) ) {
        ///@todo Shouldn't it be enough to just run this?
        //UnknownsField->giveUnknownValue(dof, mode, currentStep);
        double rtdt = UnknownsField->giveUnknownValue(dof, VM_Total, currentStep);
        double rt   = UnknownsField->giveUnknownValue(dof, VM_Total, previousStep);
        double psi = ( t - previousStep->giveTargetTime() ) / currentStep->giveTimeIncrement();
        if ( mode == VM_Velocity ) {
            return ( rtdt - rt ) / currentStep->giveTimeIncrement();
        } else if ( mode == VM_TotalIntrinsic ) {
            // only supported for current step
            return this->alpha * rtdt + ( 1. - this->alpha ) * rt; 
        } else if ( mode == VM_Total ) {
            return psi * rtdt + ( 1. - psi ) * rt;
        } else if ( mode == VM_Incremental ) {
            if ( previousStep->isIcApply() ) {
                return 0;
            } else {
                return ( rtdt - rt );
            }
        } else {
            OOFEM_ERROR("Unknown mode %s is undefined for this problem", __ValueModeTypeToString(mode) );
        }
    } else {
        OOFEM_ERROR("time value %f not within bounds %f and %f", t, previousStep->giveTargetTime(), currentStep->giveTargetTime() );
    }

    return 0.; // to make compiler happy;
}


void
NLTransientTransportProblem :: applyIC(TimeStep *stepWhenIcApply)
{
    Domain *domain = this->giveDomain(1);

    NonStationaryTransportProblem :: applyIC(stepWhenIcApply);

    if ( waveMode && waveVelocity.giveSize() == 0 ) {
        const int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());
        if ( initialVelocityValues.giveSize() != 0 && initialVelocityValues.giveSize() != neq ) {
            OOFEM_ERROR("Initial velocity array has %d values, but the model has %d equations",
                        initialVelocityValues.giveSize(), neq);
        }
        waveVelocity.resize(neq);
        if ( initialVelocityValues.giveSize() == neq ) {
            waveVelocity = initialVelocityValues;
        } else {
            waveVelocity.zero();
        }
        waveAcceleration.resize(neq);
        waveAcceleration.zero();
        if ( initialVelocityValues.giveSize() == 0 ) {
            for ( int eq = 1; eq <= neq; ++eq ) {
                waveVelocity.at(eq) = initialVelocity;
            }
        }
    }

    // update element state according to given ic
    for ( auto &elem : domain->giveElements() ) {
        TransportElement *element = static_cast< TransportElement * >( elem.get() );
        element->updateInternalState(stepWhenIcApply);
        element->updateYourself(stepWhenIcApply);
    }

    if ( waveMode && waveAcceleration.giveSize() > 0 ) {
        // Newmark requires acceleration consistent with the initial equilibrium.
        // Starting with a=0 introduces an O(dt) startup error even for a smooth
        // solution. Solve M_bulk*a0 = source - K*u0 - M_surface*v0; the line
        // capacity is the first-order dynamic-boundary term and contributes via
        // the prescribed initial velocity, not to the inertia matrix.
        const auto &numbering = EModelDefaultEquationNumbering();
        FloatArray initialRhs(waveAcceleration.giveSize());
        initialRhs.zero();
        this->assembleVectorFromElements(initialRhs, stepWhenIcApply,
                                         TransportExternalForceAssembler(), VM_Total,
                                         numbering, domain);
        this->assembleVectorFromDofManagers(initialRhs, stepWhenIcApply,
                                           ExternalForceAssembler(), VM_Total,
                                           numbering, domain);

        IntArray loc;
        FloatMatrix stiffness, capacity;
        FloatArray localU, localVelocity, contribution;
        for ( auto &element : domain->giveElements() ) {
            if ( element->giveParallelMode() == Element_remote || !element->isActivated(stepWhenIcApply) ) {
                continue;
            }
            element->giveLocationArray(loc, numbering);
            element->giveCharacteristicMatrix(stiffness, TangentStiffnessMatrix, stepWhenIcApply);
            element->computeVectorOf(VM_Total, stepWhenIcApply, localU);
            contribution.beProductOf(stiffness, localU);
            contribution.negated();
            initialRhs.assemble(contribution, loc);

            if ( element->giveGeometryType() == EGT_line_1 ) {
                element->giveCharacteristicMatrix(capacity, MassMatrix, stepWhenIcApply);
                localVelocity.resize(loc.giveSize());
                for ( int j = 1; j <= loc.giveSize(); ++j ) {
                    const int eq = loc.at(j);
                    localVelocity.at(j) = eq > 0 ? waveVelocity.at(eq) : initialVelocity;
                }
                contribution.beProductOf(capacity, localVelocity);
                contribution.negated();
                initialRhs.assemble(contribution, loc);
            } else if ( waveDamping != 0.0 ) {
                element->giveCharacteristicMatrix(capacity, MassMatrix, stepWhenIcApply);
                localVelocity.resize(loc.giveSize());
                for ( int j = 1; j <= loc.giveSize(); ++j ) {
                    const int eq = loc.at(j);
                    localVelocity.at(j) = eq > 0 ? waveVelocity.at(eq) : initialVelocity;
                }
                contribution.beProductOf(capacity, localVelocity);
                contribution.times(-waveDamping);
                initialRhs.assemble(contribution, loc);
            }
        }

        conductivityMatrix->zero();
        this->assemble(*conductivityMatrix, stepWhenIcApply, BulkWaveMassAssembler(),
                       numbering, domain);
        this->giveNumericalMethod(this->giveCurrentMetaStep());
        linSolver->solve(*conductivityMatrix, initialRhs, waveAcceleration);
    }
    if ( bdfOrder > 1 && bdfElementHistory.empty() ) {
        captureBdfElementState(stepWhenIcApply);
    }
}

void
NLTransientTransportProblem :: createPreviousSolutionInDofUnknownsDictionary(TimeStep *tStep)
{
    //Copy the last known temperature to be a previous solution
    for ( auto &domain: domainList ) {
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( auto &node : domain->giveDofManagers() ) {
                for ( Dof *dof: *node ) {
                    double val = dof->giveUnknown(VM_Total, tStep); //get number on hash=0(current)
                    dof->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, val);
                }
            }
        }
    }
}


void
NLTransientTransportProblem :: updateYourself(TimeStep *tStep)
{
    if ( waveMode ) {
        this->updateWaveKinematics(tStep);
    }
    if ( bdfOrder > 1 ) {
        captureBdfElementState(tStep);
        bdfPreviousIncrement = tStep->giveTimeIncrement();
    }
    //this->updateInternalState(tStep);
    //Set intrinsic time for a staggered problem here. This is important for materials such as hydratingconcretemat, who keep history of intrinsic times.
    double intrinsicTime = tStep->giveIntrinsicTime();
    double Tau = tStep->giveTargetTime() - ( 1. - alpha ) * tStep->giveTimeIncrement();
    tStep->setIntrinsicTime(Tau);
    NonStationaryTransportProblem :: updateYourself(tStep);
    //return back to original intrinsic time
    tStep->setIntrinsicTime(intrinsicTime);
}

void
NLTransientTransportProblem :: updateWaveKinematics(TimeStep *tStep)
{
    const FloatArray *current = UnknownsField->giveSolutionVector(tStep);
    const FloatArray *previous = UnknownsField->giveSolutionVector(tStep->givePreviousStep());
    const double dt = tStep->giveTimeIncrement();
    const double accelerationFactor = 1.0 / ( newmarkBeta * dt * dt );
    const double velocityPredictorFactor = dt * ( 1.0 - newmarkGamma );
    const double displacementPredictorAccelerationFactor = dt * dt * ( 0.5 - newmarkBeta );
    const int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());
    for ( int eq = 1; eq <= neq; ++eq ) {
        const double displacementPredictor = previous->at(eq) + dt * waveVelocity.at(eq)
            + displacementPredictorAccelerationFactor * waveAcceleration.at(eq);
        const double acceleration = accelerationFactor * ( current->at(eq) - displacementPredictor );
        waveVelocity.at(eq) += velocityPredictorFactor * waveAcceleration.at(eq)
            + newmarkGamma * dt * acceleration;
        waveAcceleration.at(eq) = acceleration;
    }
}

int
NLTransientTransportProblem :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    if ( mode == VM_Total ) { //Nodal temperature
        if ( tStep->giveNumber() == this->giveCurrentStep()->giveNumber() ) { //current time
            return 0;
        } else if ( tStep->giveNumber() == this->giveCurrentStep()->giveNumber() - 1 ) { //previous time
            return 1;
        } else {
            OOFEM_ERROR("No history available at TimeStep %d = %f, called from TimeStep %d = %f", tStep->giveNumber(), tStep->giveTargetTime(), this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveTargetTime() );
        }
    } else {
        OOFEM_ERROR("ValueModeType %s undefined", __ValueModeTypeToString(mode));
    }

    return 0;
}

void
NLTransientTransportProblem :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    // update DoF unknowns dictionary. Store the last and previous temperature only, see giveUnknownDictHashIndx
    for ( Dof *dof: *inode ) {
        int eqNum = dof->__giveEquationNumber();
        double val;
        if ( dof->hasBc(tStep) ) { // boundary condition
            val = dof->giveBcValue(VM_Total, tStep);
        } else {
            FloatArray *vect = this->UnknownsField->giveSolutionVector(tStep);
            val = vect->at(eqNum);
        }

        //update temperature, which is present in every node
        dof->updateUnknownsDictionary(tStep, VM_Total, val);
    }
}


void
NLTransientTransportProblem :: updateInternalState(TimeStep *tStep)
{
    for ( auto &domain: domainList ) {
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( auto &dman : domain->giveDofManagers() ) {
                //update dictionary entry or add a new pair if the position is missing
                this->updateDofUnknownsDictionary(dman.get(), tStep);
            }
        }

        for ( auto &elem : domain->giveElements() ) {
            elem->updateInternalState(tStep);
        }
    }
}

void
NLTransientTransportProblem :: assembleAlgorithmicPartOfRhs(FloatArray &answer,
                                                            const UnknownNumberingScheme &ns, TimeStep *tStep)
{
    //
    // Computes right hand side on all nodes
    //
    double t = tStep->giveTargetTime();
    IntArray loc;
    FloatMatrix charMtrxCond, charMtrxCap;
    FloatArray r, drdt, contrib, help;
    Element *element;
    TimeStep *previousStep = this->givePreviousStep(); //r_t
    TimeStep *currentStep = this->giveCurrentStep(); //r_{t+\Delta t}. Note that *tStep is a Tau step between r_t and r_{t+\Delta t}

    Domain *domain = this->giveDomain(1);
    int nelem = domain->giveNumberOfElements();

    if ( waveMode ) {
        const double dt = tStep->giveTimeIncrement();
        const double accelerationFactor = 1.0 / ( newmarkBeta * dt * dt );
        const double velocityAccelerationFactor = newmarkGamma / ( newmarkBeta * dt );
        const double predictorVelocityFactor = dt * ( 1.0 - newmarkGamma );
        const double predictorDisplacementAccelerationFactor = dt * dt * ( 0.5 - newmarkBeta );
        for ( int i = 1; i <= nelem; i++ ) {
            element = domain->giveElement(i);
            if ( element->giveParallelMode() == Element_remote || !element->isActivated(tStep) ) {
                continue;
            }
            element->giveLocationArray(loc, ns);
            element->giveCharacteristicMatrix(charMtrxCond, TangentStiffnessMatrix, tStep);
            element->giveCharacteristicMatrix(charMtrxCap, MassMatrix, tStep);
            element->computeVectorOf(VM_Total, currentStep, r);
            element->computeVectorOf(VM_Total, previousStep, drdt);

            const bool dynamicBoundary = element->giveGeometryType() == EGT_line_1;
            FloatArray dynamicTerm(r.giveSize());
            for ( int j = 1; j <= loc.giveSize(); ++j ) {
                const int eq = loc.at(j);
                double velocity = initialVelocity;
                double acceleration = 0.0;
                if ( eq > 0 ) {
                    velocity = waveVelocity.at(eq);
                    acceleration = waveAcceleration.at(eq);
                }
                const double displacementPredictor = drdt.at(j) + dt * velocity
                    + predictorDisplacementAccelerationFactor * acceleration;
                const double currentAcceleration = accelerationFactor * ( r.at(j) - displacementPredictor );
                const double currentVelocity = velocity + predictorVelocityFactor * acceleration
                    + velocityAccelerationFactor * ( r.at(j) - displacementPredictor );
                dynamicTerm.at(j) = dynamicBoundary
                    ? currentVelocity
                    : currentAcceleration + waveDamping * currentVelocity;
            }

            help.beProductOf(charMtrxCap, dynamicTerm);
            contrib.beProductOf(charMtrxCond, r);
            contrib.add(help);
            contrib.negated();
            answer.assemble(contrib, loc);
        }
        return;
    }

    if ( bdfOrder > 1 ) {
        const int activeOrder = std::min(bdfOrder, static_cast<int>(bdfElementHistory.size()));
        std::vector<double> coefficients;
        giveBdfCoefficients(activeOrder, coefficients);
        const double inverseDt = 1.0 / tStep->giveTimeIncrement();
        for ( int i = 1; i <= nelem; i++ ) {
            element = domain->giveElement(i);
            if ( element->giveParallelMode() == Element_remote || !element->isActivated(tStep) ) {
                continue;
            }
            element->giveLocationArray(loc, ns);
            element->giveCharacteristicMatrix(charMtrxCond, TangentStiffnessMatrix, tStep);
            element->giveCharacteristicMatrix(charMtrxCap,
                lumpedCapacityStab ? LumpedMassMatrix : MassMatrix, tStep);

            if ( lumpedCapacityStab ) {
                const int size = charMtrxCap.giveNumberOfRows();
                for ( int j = 1; j <= size; j++ ) {
                    double rowSum = 0.0;
                    for ( int k = 1; k <= size; k++ ) {
                        rowSum += charMtrxCap.at(j, k);
                        charMtrxCap.at(j, k) = 0.0;
                    }
                    charMtrxCap.at(j, j) = rowSum;
                }
            }

            element->computeVectorOf(VM_Total, currentStep, r);
            drdt = r;
            drdt.times(coefficients[0]);
            for ( int j = 1; j <= activeOrder; j++ ) {
                if ( bdfElementHistory[j - 1].size() != static_cast<size_t>(nelem) ) {
                    OOFEM_ERROR("BDF solution history is inconsistent with the element numbering");
                }
                drdt.add(coefficients[j], bdfElementHistory[j - 1][i - 1]);
            }
            drdt.times(inverseDt);
            help.beProductOf(charMtrxCap, drdt);
            contrib.beProductOf(charMtrxCond, r);
            contrib.add(help);
            contrib.negated();
            answer.assemble(contrib, loc);
        }
        return;
    }

    for ( int i = 1; i <= nelem; i++ ) {
        element = domain->giveElement(i);
        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

        if ( !element->isActivated(tStep) ) {
            continue;
        }

        element->giveLocationArray(loc, ns);

        element->giveCharacteristicMatrix(charMtrxCond, TangentStiffnessMatrix, tStep);
        element->giveCharacteristicMatrix(charMtrxCap, CapacityMatrix, tStep);


        /*
         *  element -> computeVectorOf (VM_Total, tStep, r);
         *  element -> computeVectorOf (VM_Velocity, tStep, drdt);
         */

        if ( ( t >= previousStep->giveTargetTime() ) && ( t <= currentStep->giveTargetTime() ) ) {
            FloatArray rp, rc;
            element->computeVectorOf(VM_Total, currentStep, rc);
            element->computeVectorOf(VM_Total, previousStep, rp);

            //approximate derivative with a difference
            drdt.beDifferenceOf(rc, rp);
            drdt.times( 1. / currentStep->giveTimeIncrement() );
            //approximate current solution from linear interpolation
            rp.times(1 - alpha);
            rc.times(alpha);
            r = rc;
            r.add(rp);
        } else {
            OOFEM_ERROR("unsupported time value");
        }
	// bp: r can be computed simply as element->computeVectorOf(VM_TotalIntrinsic, currentStep, r);

        if ( lumpedCapacityStab ) {
            int size = charMtrxCap.giveNumberOfRows();
            double s;
            for ( int j = 1; j <= size; j++ ) {
                s = 0.0;
                for ( int k = 1; k <= size; k++ ) {
                    s += charMtrxCap.at(j, k);
                    charMtrxCap.at(j, k) = 0.0;
                }

                charMtrxCap.at(j, j) = s;
            }
        }

        help.beProductOf(charMtrxCap, drdt);

        contrib.beProductOf(charMtrxCond, r);
        contrib.add(help);
        contrib.negated();

        answer.assemble(contrib, loc);
    }
}
} // end namespace oofem
