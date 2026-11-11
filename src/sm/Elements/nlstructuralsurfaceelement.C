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

#include "sm/Elements/nlstructuralsurfaceelement.h"
#include "sm/Elements/Structural3DSurfaceElement.h"

#include "sm/Materials/structuralms.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "sm/CrossSections/simplesurfacecrosssection.h"
#include "integrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"
#include "engngm.h"
#include "mathfem.h"

namespace oofem {
NLStructuralSurfaceElement::NLStructuralSurfaceElement( int n, Domain *aDomain ) :
    NLStructuralElement( n, aDomain )
// Constructor. Creates an element with number n, belonging to aDomain.
{}

void NLStructuralSurfaceElement::computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray normal = this->giveNormal( gp );
    FloatArray vF;
    this->computeDeformationGradientVector( vF, gp, tStep );

    // Include prestrain Prestrain
    Tensor2_3d F0( 1., 0., 0., 0., 1., 0., 0., 0., 0. );
    FloatArray vFmod;
    this->PrestrainDeformationGradient( vF, vFmod, F0, gp, tStep );

    FloatArray answerTmp = cs->giveFirstPKSurfaceStresses( vFmod, normal, gp, tStep ); // Compute stress with updated F
    

    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( cs->giveMaterial( gp ) ); // Get material
    double J0 = mat->compute_surface_determinant( F0 ); // Compute surface jacobian

    FloatArrayF<9> vP( answerTmp );
    Tensor2_3d P( vP ), Pmod;
    Pmod( i_3, j_3 ) = 1. / J0 * P( i_3, k_3 ) * F0( j_3, k_3 ); // Modify the answert 
    answer = Pmod.to_voigt_form();

    // update status 
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( mat->giveStatus( gp ) );
    status->letTempFVectorBe( vF);
    status->letTempPVectorBe( answer );
    
}


void NLStructuralSurfaceElement::computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    auto cs  = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( cs->giveMaterial( gp ) );
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( mat->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );

    Tensor2_3d F0( 1., 0., 0., 0., 1., 0., 0., 0., 0. );
    FloatArray vFmod;
    this->PrestrainDeformationGradient( vF, vFmod, F0, gp, tStep ); // Get the prestrain

    status->letTempFVectorBe( vFmod ); // Update F in status

    FloatArray normal = this->giveNormal( gp );
    FloatMatrix answerTmp = cs->giveSurfaceStiffnessMatrix_dPdF_3d( rMode, normal, gp, tStep );// Compute stiffness with updated F

    status->letTempFVectorBe( vF ); // get status to the previous state

    // Modify the stiffness
    FloatMatrixF<9, 9> mA( answerTmp );
    Tensor4_3d A( mA ), Amod;
    double J0 = mat->compute_surface_determinant( F0 );
    Amod( i_3, j_3, k_3, l_3 ) = 1. / J0 * F0( j_3, m_3 ) * A( i_3, m_3, k_3, n_3 ) * F0( l_3, n_3 );
    answer = Amod.to_voigt_form();
}


void NLStructuralSurfaceElement ::computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep )
{
    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray normal = this->giveNormal( gp );   
    answer = cs->giveRealSurfaceStresses( strain, normal, gp, tStep );
}

void NLStructuralSurfaceElement::computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present).
    // The element's own dof mask is used, so planar (u, v) and spatial (u, v, w) surfaces share this code.
    FloatArray u;
    this->computeVectorOf( VM_Total, tStep, u ); // solution vector, in the element local frame
    if ( initialDisplacements ) {
        u.subtract( *initialDisplacements );
    }

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt( gp, B );
 
    answer.beProductOf( B, u );

    // Deformation gradient F = H + I
    answer.at( 1 ) += 1.0;
    answer.at( 2 ) += 1.0;

}



} // end namespace oofem
