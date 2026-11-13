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

#include "sm/Elements/structuralsurfaceelementevaluator.h"
#include "sm/Elements/structuralelement.h"
#include "sm/CrossSections/simplesurfacecrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralsurfacems.h"

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
StructuralSurfaceElementEvaluator ::StructuralSurfaceElementEvaluator():
    StructuralElementEvaluator()
{}


void StructuralSurfaceElementEvaluator::initializeFrom( InputRecord &ir )
{
    StructuralElementEvaluator::initializeFrom( ir );
    if ( this->nlGeometry != 1 ) {
        throw ValueInputException( ir, _IFT_StructuralElementEvaluator_nlgeoflag, "only nlgeo 1 is implemented for surface elements" );
    }
}


void StructuralSurfaceElementEvaluator::computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    IGAElement *elem              = static_cast<IGAElement *>( this->giveElement() );
    SimpleSurfaceCrossSection *cs = static_cast<SimpleSurfaceCrossSection *>( elem->giveCrossSection() );
    FloatArray normal = elem->giveNormal( gp ); 
    FloatArray vF;
    this->computeDeformationGradientVector( vF, gp, tStep );
    answer = cs->giveFirstPKSurfaceStresses( vF, normal, gp, tStep );
}

void StructuralSurfaceElementEvaluator::computeSecondOrderStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    IGAElement *elem              = static_cast<IGAElement *>( this->giveElement() );
    SimpleSurfaceCrossSection *cs = static_cast<SimpleSurfaceCrossSection *>( elem->giveCrossSection() );
    FloatArray normal             = elem->giveNormal( gp );
    FloatArray vF, vG, vDI;
    this->computeDeformationGradientVector( vF, gp, tStep );
    this->computeGradientOfDeformationGradientVector( vG, gp, tStep );
    this->computeSurfaceIdentityGradient( vDI, gp, tStep );

    answer = cs->giveSecondOrderSurfaceStresses( vF, vG, vDI, normal, gp, tStep );
}


void StructuralSurfaceElementEvaluator ::computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    this->initializeReferenceGradients( gp, tStep );
    IGAElement *elem              = static_cast<IGAElement *>( this->giveElement() );
    SimpleSurfaceCrossSection *cs = static_cast<SimpleSurfaceCrossSection *>( elem->giveCrossSection() );
    answer = cs->giveSurfaceStiffnessMatrix_dPdF_3d( rMode, elem->giveNormal( gp ), gp, tStep );
}

void StructuralSurfaceElementEvaluator ::computeConstitutiveMatrix_dAddF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    this->initializeReferenceGradients( gp, tStep );
    IGAElement *elem              = static_cast<IGAElement *>( this->giveElement() );
    SimpleSurfaceCrossSection *cs = static_cast<SimpleSurfaceCrossSection *>( elem->giveCrossSection() );
    answer = cs->giveSurfaceStiffnessMatrix_dAddF_3d( rMode, elem->giveNormal( gp ), gp, tStep );
}

void StructuralSurfaceElementEvaluator ::computeConstitutiveMatrix_dAdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    this->initializeReferenceGradients( gp, tStep );
    IGAElement *elem              = static_cast<IGAElement *>( this->giveElement() );
    SimpleSurfaceCrossSection *cs = static_cast<SimpleSurfaceCrossSection *>( elem->giveCrossSection() );
    answer = cs->giveSurfaceStiffnessMatrix_dAdF_3d( rMode, elem->giveNormal( gp ), gp, tStep );
}

void StructuralSurfaceElementEvaluator::initializeReferenceGradients( GaussPoint *gp, TimeStep *tStep )
{
    if ( !this->hasSecondGradient ) {
        return;
    }
    auto cs = static_cast<SimpleSurfaceCrossSection *>( this->giveElement()->giveCrossSection() );
    auto status = static_cast<StructuralSurfaceMaterialStatus *>( cs->giveMaterial( gp )->giveStatus( gp ) );
    if ( status->giveTempdIVector().giveSize() == 0 ) {
        // A fresh status represents the reference configuration F = I_s. Keep the usual
        // status-based tangent policy after this first initialization.
        FloatArray gradI;
        this->computeSurfaceIdentityGradient( gradI, gp, tStep );
        status->letTempGVectorBe( gradI );
        status->letTempdIVectorBe( gradI );
    }
}

void StructuralSurfaceElementEvaluator::computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the surface deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Displacement gradient H = du/dX from the control points of the knot span
    FloatMatrix B;
    FloatArray ur;
    this->computeBHmatrixAt( gp, B );
    this->giveSpanDisplacements( ur, gp, tStep );
    answer.beProductOf( B, ur );

    // F = H + I_s with the reference surface identity I_s = I - N x N
    IGAElement *elem  = static_cast<IGAElement *>( this->giveElement() );
    FloatArray normal = elem->giveNormal( gp );
    FloatArrayF<9> du( answer );
    Tensor2_3d du_tens( du ), Ftens;

    FloatArrayF<3> UDn( normal );
    Tensor1_3d UDnormalT( UDn );
    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Identity2d, Identity3d( indent3d );
    Identity2d( i_3, j_3 ) = Identity3d( i_3, j_3 ) - UDnormalT( i_3 ) * UDnormalT( j_3 );

    Ftens( i_3, j_3 ) = du_tens( i_3, j_3 ) + Identity2d( i_3, j_3 );
    
    // reduce back to voigt
    answer = Ftens.to_voigt_form();
}

void StructuralSurfaceElementEvaluator::computeSurfaceIdentityGradient( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // grad I = G X (27 components) at the Gauss point gp, as grad F with the reference positions
    // of the control points of the knot span
    FloatMatrix G;
    FloatArray X;
    this->computeGMatrixAt( gp, G );
    this->giveSpanCoordinates( X, gp );
    answer.beProductOf( G, X );
}


} // end namespace oofem
