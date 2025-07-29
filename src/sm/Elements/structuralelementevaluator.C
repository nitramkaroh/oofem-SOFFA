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

#include "sm/Elements/structuralelementevaluator.h"
#include "sm/Elements/structuralelement.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "domain.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "mathfem.h"
#include "iga/iga.h"

namespace oofem {
StructuralElementEvaluator :: StructuralElementEvaluator()
{
    this->rotationMatrix.clear();
}

/*
 * int StructuralElementEvaluator :: giveIntegrationElementCodeNumbers(IntArray &answer, Element *elem,
 *                                                                  IntegrationRule *ie) {
 *  int i;
 *  IntArray mask, nodeDofIDMask, nodalArray;
 *
 *  // first evaluate nonzero basis function mask
 *  if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
 *      IGAIntegrationElement *ee = static_cast< IGAIntegrationElement * >( ie );
 *      elem->giveInterpolation()->giveKnotSpanBasisFuncMask(* ee->giveKnotSpan(), mask);
 *      // loop over nonzero shape functions and assemble localization array
 *      answer.clear();
 *      for ( i = 1; i <= mask.giveSize(); i++ ) {
 *          elem->giveDofManDofIDMask(mask.at(i), nodeDofIDMask);
 *          elem->giveDofManager( mask.at(i) )->giveLocationArray(nodeDofIDMask, nodalArray);
 *          answer.followedBy(nodalArray);
 *      }
 *
 *      return 1;
 *  } else {
 *      return 0;
 *  }
 * }
 */

int StructuralElementEvaluator :: giveIntegrationElementLocalCodeNumbers(IntArray &answer, Element *elem,
                                                                         IntegrationRule *ie)
{
    int nsd;
    IntArray mask, nodeDofIDMask, nodalArray;
    int dofmandof;

    // get number of dofs in node
    elem->giveDofManDofIDMask(1, nodeDofIDMask);
    dofmandof = nodeDofIDMask.giveSize();

    nsd = elem->giveInterpolation()->giveNsd();

    // first evaluate nonzero basis function mask
    if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
        IGAIntegrationElement *ee = static_cast< IGAIntegrationElement * >(ie);
        elem->giveInterpolation()->giveKnotSpanBasisFuncMask(* ee->giveKnotSpan(), mask);
        // loop over nonzero shape functions and assemble localization array
        answer.clear();
        for ( int i = 1; i <= mask.giveSize(); i++ ) {
            nodalArray.resize( nodeDofIDMask.giveSize() );
            for ( int j = 1; j <= nsd; j++ ) {
                nodalArray.at(j) = dofmandof * ( mask.at(i) - 1 ) + j;
            }

            answer.followedBy(nodalArray);
        }

        return 1;
    } else {
        return 0;
    }
}


void StructuralElementEvaluator :: giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type == InternalForcesVector ) {
        this->giveInternalForcesVector(answer, tStep, false); /// @todo Only for total value mode type (?)
    } else if ( type == LastEquilibratedInternalForcesVector ) {
        this->giveInternalForcesVector(answer, tStep, true); /// @todo Only for total value mode type (?)
    } else {
        answer.clear();
    }
}


void StructuralElementEvaluator :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                            CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    if ( mtrx == TangentStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else if ( mtrx == MassMatrix ) {
        double mass;
        this->computeConsistentMassMatrix(answer, tStep, mass);
    } else if ( mtrx == LumpedMassMatrix ) {
        this->computeLumpedMassMatrix(answer, tStep);
    } else {
        OOFEM_ERROR( "Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
    }
}


void StructuralElementEvaluator :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    double mass;
    Element *elem = this->giveElement();
    int numberOfDofMans = elem->giveNumberOfDofManagers();
    IntArray nodeDofIDMask, dimFlag(3);
    int indx = 0, ldofs, dim;
    double summ;

    if ( !this->isActivated(tStep) ) {
        int ndofs = elem->computeNumberOfDofs();
        answer.resize(ndofs, ndofs);
        answer.zero();
        return;
    }

    this->computeConsistentMassMatrix(answer, tStep, mass);
    ldofs = answer.giveNumberOfRows();

    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        elem->giveDofManDofIDMask(i, nodeDofIDMask);
        for ( int j = 1; j <= nodeDofIDMask.giveSize(); j++ ) {
            indx++;
            // zero all off-diagonal terms
            for ( int k = 1; k <= ldofs; k++ ) {
                if ( k != indx ) {
                    answer.at(indx, k) = 0.;
                    answer.at(k, indx) = 0.;
                }
            }

            if ( ( nodeDofIDMask.at(j) != D_u ) && ( nodeDofIDMask.at(j) != D_v ) && ( nodeDofIDMask.at(j) != D_w ) ) {
                // zero corresponding diagonal member too <= no displacement dof
                answer.at(indx, indx) = 0.;
            } else if ( nodeDofIDMask.at(j) == D_u ) {
                dimFlag.at(1) = 1;
            } else if ( nodeDofIDMask.at(j) == D_v ) {
                dimFlag.at(2) = 1;
            } else if ( nodeDofIDMask.at(j) == D_w ) {
                dimFlag.at(3) = 1;
            }
        }
    }

    if ( indx != ldofs ) {
        OOFEM_ERROR("internal consistency check failed");
    }

    dim = dimFlag.at(1) + dimFlag.at(2) + dimFlag.at(3);
    summ = 0.;
    for ( int k = 1; k <= ldofs; k++ ) {
        summ += answer.at(k, k);
    }

    answer.times(dim * mass / summ);
}


void StructuralElementEvaluator :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass)
// Computes numerically the consistent (full) mass matrix of the receiver.
{
    //Element *elem = this->giveElement();
    StructuralElement *elem = static_cast< StructuralElement * >( this->giveElement() );
    int ndofs = elem->computeNumberOfDofs();
    FloatMatrix n;
    IntegrationRule *iRule;
    IntArray mask;

    answer.resize(ndofs, ndofs);
    answer.zero();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    if ( ( iRule = this->giveMassMtrxIntegrationRule() ) ) {
        OOFEM_ERROR("no integration rule available");
    }

    this->giveMassMtrxIntegrationMask(mask);

    mass = 0.;

    for ( GaussPoint *gp: *iRule ) {
        double density = elem->giveStructuralCrossSection()->give('d', gp);
        double dV      = this->computeVolumeAround(gp);
        mass   += density * dV;
        this->computeNMatrixAt(n, gp);

        if ( mask.isEmpty() ) {
            answer.plusProductSymmUpper(n, n, density * dV);
        } else {
            for ( int i = 1; i <= ndofs; i++ ) {
                for ( int j = i; j <= ndofs; j++ ) {
                    double summ = 0.;
                    for ( int k = 1; k <= n.giveNumberOfRows(); k++ ) {
                        if ( mask.at(k) == 0 ) {
                            continue;
                        }

                        summ += n.at(k, i) * n.at(k, j);
                    }

                    answer.at(i, j) += summ * density * dV;
                }
            }
        }
    }

    answer.symmetrized();
}


void StructuralElementEvaluator :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, bool useUpdatedGpRecord)
{
    Element *elem = this->giveElement();
    int ndofs = elem->computeNumberOfDofs();
    FloatMatrix b, Gmat;
    FloatArray strain, stress, u, temp, secondOrderStress;
    IntArray irlocnum;

    elem->computeVectorOf(VM_Total, tStep, u);

    answer.resize(ndofs);
    answer.zero();
    FloatArray *m = & answer;
    if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    int numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
    // loop over individual integration rules
    for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        m->clear();
        IntegrationRule *iRule = elem->giveIntegrationRule(ir);
        for ( GaussPoint *gp: *iRule ) {
            StructuralMaterialStatus *matStat = static_cast<StructuralMaterialStatus *>( gp->giveMaterialStatus() );

            if ( nlGeometry == 0 ) {
                this->computeBMatrixAt(b, gp);
                if ( useUpdatedGpRecord ) {
                    stress = matStat->giveStressVector();
                } else {
                    this->computeStrainVector(strain, gp, tStep, u); ///@todo This part computes the B matrix again; Inefficient.
                    this->computeStressVector(stress, strain, gp, tStep);
                }
            } else if ( nlGeometry == 1 ) { // First Piola-Kirchhoff stress
                if ( this->hasSecondGradient ) {
                    this->computeGMatrixAt( gp, Gmat );
                    this->computeSecondOrderStressVector( secondOrderStress, gp, tStep );
                }

                if ( useUpdatedGpRecord == 1 ) {
                    stress = matStat->givePVector();
                } else {
                    this->computeFirstPKStressVector( stress, gp, tStep );

                    ///@todo This is actaully inefficient since it constructs B and twice and collects the nodal unknowns over and over.
                }
                this->computeBHmatrixAt( gp, b );

            }

            if ( stress.giveSize() == 0 ) {
                break;
            }

            // compute nodal representation of internal forces using f = B^T*Sigma dV
            double dV = this->computeVolumeAround(gp);

            if ( nlGeometry == 1 ) { // First Piola-Kirchhoff stress
                if ( stress.giveSize() == 9 ) {
                    FloatArray stressTemp;
                    StructuralMaterial::giveReducedVectorForm( stressTemp, stress, gp->giveMaterialMode() );
                    m->plusProduct( b, stressTemp, dV );
                } else {
                    m->plusProduct( b, stress, dV );
                }

                if ( this->hasSecondGradient ) {
                    m->plusProduct( Gmat, secondOrderStress, dV );
                }
            } else {
                m->plusProduct( b, stress, dV );
            }
        }
        // localize irule contribution into element matrix
        if ( this->giveIntegrationElementLocalCodeNumbers(irlocnum, elem, iRule) ) {
            answer.assemble(* m, irlocnum);
        }
    } // end loop over irules

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
    }
}


void StructuralElementEvaluator :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, FloatArray &u)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray ur;
    Element *elem = this->giveElement();

    if ( !this->isActivated(tStep) ) {
        answer.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        answer.zero();
        return;
    }

    this->computeBMatrixAt(b, gp);

    // get local code numbers corresponding to ir
    IntArray lc;
    this->giveIntegrationElementLocalCodeNumbers(lc, elem, gp->giveIntegrationRule());
    ur.resize( b.giveNumberOfColumns() );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        ur.at(i) = u.at( lc.at(i) );
    }

    answer.beProductOf(b, ur);
}

void StructuralElementEvaluator :: updateInternalState(TimeStep *tStep)
// Updates the receiver at end of step.
{
    FloatArray u;
    Element *elem = this->giveElement();

    elem->computeVectorOf(VM_Total, tStep, u);

#if 0
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(initialDisplacements);
    }
#endif

    FloatArray strain, stress;

    // force updating strains & stresses
    for ( int i = 0; i < elem->giveNumberOfIntegrationRules(); i++ ) {
#ifdef __PARALLEL_MODE
        if ( this->giveElement()->giveKnotSpanParallelMode(i) == Element_remote ) {
            continue;
        }
#endif
        for ( GaussPoint *gp: *elem->giveIntegrationRule(i) ) {
            this->computeStrainVector(strain, gp, tStep, u);
            this->computeStressVector(stress, strain, gp, tStep);
        }
    }
}


void StructuralElementEvaluator :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    int numberOfIntegrationRules;
    FloatMatrix temp, bj, d, dbj, dAddF, dAdF, Gmat, dAG, dAB;
    Element *elem = this->giveElement();
    StructuralCrossSection *cs = static_cast< StructuralCrossSection * >( elem->giveCrossSection() );
    int ndofs = elem->computeNumberOfDofs();
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);
    IntArray irlocnum;

    answer.resize(ndofs, ndofs);
    answer.zero();

    FloatMatrix *m = & answer;
    if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
    // loop over individual integration rules
    for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
#ifdef __PARALLEL_MODE
        if ( this->giveElement()->giveKnotSpanParallelMode(ir) == Element_remote ) {
            continue;
        }
        //fprintf (stderr, "[%d] Computing element.knotspan %d.%d\n", elem->giveDomain()->giveEngngModel()->giveRank(), elem->giveNumber(), ir);
#endif
        m->clear();
        IntegrationRule *iRule = elem->giveIntegrationRule(ir);
        // loop over individual integration points
        for ( GaussPoint *gp: *iRule ) {
            double dV = this->computeVolumeAround(gp);
            if ( nlGeometry == 0 ) {
                this->computeBMatrixAt(bj, gp);
                this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
            } else if ( nlGeometry == 1 ) {
                if ( this->hasSecondGradient ) { // needs to be first
                    this->computeGMatrixAt( gp, Gmat );
                    this->computeConstitutiveMatrix_dAddF_At( dAddF, rMode, gp, tStep ); 
                    this->computeConstitutiveMatrix_dAdF_At( dAdF, rMode, gp, tStep ); 
                }
                this->computeBHmatrixAt( gp, bj );
                this->computeConstitutiveMatrix_dPdF_At( d, rMode, gp, tStep );
            }

            dbj.beProductOf(d, bj);
            if ( matStiffSymmFlag ) {
                m->plusProductSymmUpper(bj, dbj, dV);
            } else {
                m->plusProductUnsym(bj, dbj, dV);
            }


            if ( this->hasSecondGradient ) {

                /////////////
                //FloatMatrix temp2;
                //temp2.beTProductOf( bj, dbj );
                //temp2.printYourself();
                /////////////

                dAG.beProductOf( dAddF, Gmat ); // term with derivative wrt dF
                m->plusProductSymmUpper( Gmat, dAG, dV );
                /////////////
                //FloatMatrix temp;
                //temp.beTProductOf( Gmat, dAG );
                //temp.printYourself();
                /////////////

                dAB.beProductOf( dAdF, bj ); // term with derivative wrt F
                m->plusProductUnsym( Gmat, dAB, dV );

                ///////////////
                // FloatMatrix temp;
                // temp.beTProductOf( Gmat, dAB );
                // temp.printYourself();
                ///////////////

                FloatMatrix dAG2, dAdFt;
                dAdFt.beTranspositionOf( dAdF );
                dAG2.beProductOf( dAdFt, Gmat ); // term with derivative wrt F
                m->plusProductUnsym( bj, dAG2, dV );
            }
        }

        if ( matStiffSymmFlag ) {
            m->symmetrized();
        }

        // localize irule contribution into element matrix
        if ( this->giveIntegrationElementLocalCodeNumbers(irlocnum, elem, iRule) ) {
            answer.assemble(* m, irlocnum);
        }
    } // end loop over irules
}
void StructuralElementEvaluator::initializeFrom( InputRecord &ir )
{
    nlGeometry = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, nlGeometry, _IFT_StructuralElementEvaluator_nlgeoflag );
}

void StructuralElementEvaluator::computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    Element *elem = this->giveElement();
    StructuralCrossSection *cs = static_cast<StructuralCrossSection *>( elem->giveCrossSection() );

    FloatArray vF;
    this->computeDeformationGradientVector( vF, gp, tStep );
    answer = cs->giveFirstPKStresses( vF, gp, tStep );
}

void StructuralElementEvaluator::computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    Element *elem = this->giveElement();
    elem->computeVectorOf( VM_Total, tStep, u );

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt( gp, B );
    //answer.beProductOf( B, u );


    // get local code numbers corresponding to ir
    IntArray lc;
    FloatArray ur;
    this->giveIntegrationElementLocalCodeNumbers( lc, elem, gp->giveIntegrationRule() );
    ur.resize( B.giveNumberOfColumns() );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        ur.at( i ) = u.at( lc.at( i ) );
    }

    answer.beProductOf( B, ur );


    // Deformation gradient F = H + I
    MaterialMode matMode = gp->giveMaterialMode();
    if ( matMode == _3dMat || matMode == _PlaneStrain ) {
        answer.at( 1 ) += 1.0;
        answer.at( 2 ) += 1.0;
        answer.at( 3 ) += 1.0;
    } else if ( matMode == _PlaneStress ) {
        answer.at( 1 ) += 1.0;
        answer.at( 2 ) += 1.0;
    } else if ( matMode == _1dMat ) {
        answer.at( 1 ) += 1.0;
    } else {
        OOFEM_ERROR( "MaterialMode is not supported yet (%s)", __MaterialModeToString( matMode ) );
    }
}

} // end namespace oofem
