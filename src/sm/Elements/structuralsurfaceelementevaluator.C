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
//#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/CrossSections/simplesurfacecrosssection.h"
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
StructuralSurfaceElementEvaluator ::StructuralSurfaceElementEvaluator():
    StructuralElementEvaluator()
{}


void StructuralSurfaceElementEvaluator::initializeFrom( InputRecord &ir )
{
    StructuralElementEvaluator::initializeFrom( ir );
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
    IGAElement *elem              = static_cast<IGAElement *>( this->giveElement() );
    SimpleSurfaceCrossSection *cs = static_cast<SimpleSurfaceCrossSection *>( elem->giveCrossSection() );
    FloatArray normal = elem->giveNormal( gp );
    answer = static_cast<SimpleSurfaceCrossSection *>( this->giveElement()->giveCrossSection() )->giveSurfaceStiffnessMatrix_dPdF_3d( rMode, normal, gp, tStep );
}

void StructuralSurfaceElementEvaluator ::computeConstitutiveMatrix_dAddF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    IGAElement *elem              = static_cast<IGAElement *>( this->giveElement() );
    SimpleSurfaceCrossSection *cs = static_cast<SimpleSurfaceCrossSection *>( elem->giveCrossSection() );
    FloatArray normal             = elem->giveNormal( gp );
    answer                        = static_cast<SimpleSurfaceCrossSection *>( this->giveElement()->giveCrossSection() )->giveSurfaceStiffnessMatrix_dAddF_3d( rMode, normal, gp, tStep );
}

void StructuralSurfaceElementEvaluator ::computeConstitutiveMatrix_dAdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    IGAElement *elem              = static_cast<IGAElement *>( this->giveElement() );
    SimpleSurfaceCrossSection *cs = static_cast<SimpleSurfaceCrossSection *>( elem->giveCrossSection() );
    FloatArray normal             = elem->giveNormal( gp );
    answer                        = static_cast<SimpleSurfaceCrossSection *>( this->giveElement()->giveCrossSection() )->giveSurfaceStiffnessMatrix_dAdF_3d( rMode, normal, gp, tStep );
}

void StructuralSurfaceElementEvaluator::computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
    elem->computeVectorOf( VM_Total, tStep, u );

    ///////////////////////////
    // // Alternative approach
    //FloatArray nodesCoords;
    //elem->giveNodalCoordinates( nodesCoords );
    //FloatArray xDef = nodesCoords + u;
    ///////////////////////////

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt( gp, B );

    // get local code numbers corresponding to ir
    IntArray lc;
    FloatArray ur,xr;
    this->giveIntegrationElementLocalCodeNumbers( lc, elem, gp->giveIntegrationRule() );
    ur.resize( B.giveNumberOfColumns() );

    //xr.resize( B.giveNumberOfColumns() );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        ur.at( i ) = u.at( lc.at( i ) );
        //xr.at( i ) = xDef.at( lc.at( i ) );
    }

    answer.beProductOf( B, ur );

    // Construct undeformed surface identity
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

void StructuralSurfaceElementEvaluator::computeGradientOfDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
    elem->computeVectorOf( VM_Total, tStep, u );

    //// Get nodal coordinates
    FloatArray nodesCoords;
    elem->giveNodalCoordinates( nodesCoords );

    // deformed coordinates
    FloatArray xDef = nodesCoords + u;

    // Displacement gradient H = du/dX
    FloatMatrix G;
    this->computeGMatrixAt( gp, G );

    // get local code numbers corresponding to ir
    IntArray lc;
    FloatArray xr;
    this->giveIntegrationElementLocalCodeNumbers( lc, elem, gp->giveIntegrationRule() );
    xr.resize( G.giveNumberOfColumns() );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        xr.at( i ) = xDef.at( lc.at( i ) );
    }

    answer.beProductOf( G, xr );
}

void StructuralSurfaceElementEvaluator::computeSurfaceIdentityGradient( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );

    //// Get nodal coordinates
    FloatArray nodesCoords;
    elem->giveNodalCoordinates( nodesCoords );

    // Displacement gradient H = du/dX
    FloatMatrix G;
    this->computeGMatrixAt( gp, G );

    // get local code numbers corresponding to ir
    IntArray lc;
    FloatArray xr;
    this->giveIntegrationElementLocalCodeNumbers( lc, elem, gp->giveIntegrationRule() );
    xr.resize( G.giveNumberOfColumns() );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        xr.at( i ) = nodesCoords.at( lc.at( i ) );
    }

    answer.beProductOf( G, xr );
}


} // end namespace oofem
