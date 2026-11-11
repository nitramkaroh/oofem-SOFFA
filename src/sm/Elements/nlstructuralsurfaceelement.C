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

#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralmaterial.h"
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

#include <cmath>

namespace oofem {
NLStructuralSurfaceElement::NLStructuralSurfaceElement( int n, Domain *aDomain ) :
    NLStructuralElement( n, aDomain )
// Constructor. Creates an element with number n, belonging to aDomain.
{}

void NLStructuralSurfaceElement::giveMaterialKinematics( const FloatArray &vF, GaussPoint *gp, TimeStep *tStep, FloatArray &vFmat, FloatArray &normal, Tensor2_3d &F0s, double &J0 )
{
    Tensor2_3d F0( 1., 0., 0., 0., 1., 0., 0., 0., 0. ); // no prestrain
    this->PrestrainDeformationGradient( vF, vFmat, F0, gp, tStep ); // vFmat = F * F0

    // F maps e3 to zero, so only the part of F0 in the element plane counts; this part maps the plane
    // of the original configuration onto the element plane and its normal N to zero: cof(F0s) = J0 e3 x N.
    F0s = Tensor2_3d( F0( 0, 0 ), F0( 0, 1 ), F0( 0, 2 ), F0( 1, 0 ), F0( 1, 1 ), F0( 1, 2 ), 0., 0., 0. );
    Tensor2_3d cofF0 = F0s.compute_cofactor();
    normal           = { cofF0( 2, 0 ), cofF0( 2, 1 ), cofF0( 2, 2 ) };
    J0               = normal.computeNorm();
    if ( !std::isfinite( J0 ) || !( J0 > 0. ) ) {
        OOFEM_ERROR( "element %d: the prestrain deformation gradient is singular in the element plane or has a non-finite area ratio", this->giveNumber() );
    }
    normal.times( 1. / J0 );
}


void NLStructuralSurfaceElement::computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray vF, vFmod, normal;
    this->computeDeformationGradientVector( vF, gp, tStep );

    // The material sees Fmod = F * F0 with the normal of its reference configuration; the returned
    // stress is pulled back to the kinematic F.
    Tensor2_3d F0;
    double J0;
    this->giveMaterialKinematics( vF, gp, tStep, vFmod, normal, F0, J0 );

    FloatArray answerTmp = cs->giveFirstPKSurfaceStresses( vFmod, normal, gp, tStep );

    FloatArrayF<9> vP( answerTmp );
    Tensor2_3d P( vP ), Pmod;
    Pmod( i_3, j_3 ) = 1. / J0 * P( i_3, k_3 ) * F0( j_3, k_3 ); // P_element = P_material * F0^T / J0.
    answer = Pmod.to_voigt_form();

    // Store work-conjugate F and P in the same frame as BH, including for cached reactions.
    auto mat = cs->giveMaterial( gp );
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( mat->giveStatus( gp ) );
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( answer );
}


void NLStructuralSurfaceElement::computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    auto cs  = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    auto mat = cs->giveMaterial( gp );
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( mat->giveStatus( gp ) );
    FloatArray vF( status->giveTempFVector() ), vFmod, normal;

    Tensor2_3d F0;
    double J0;
    this->giveMaterialKinematics( vF, gp, tStep, vFmod, normal, F0, J0 );

    status->letTempFVectorBe( vFmod ); // the material takes F from the status
    FloatMatrix answerTmp = cs->giveSurfaceStiffnessMatrix_dPdF_3d( rMode, normal, gp, tStep );
    status->letTempFVectorBe( vF ); // Restore kinematic F after the constitutive evaluation.

    // Chain rule for W(F * F0) / J0: apply F0 to both reference indices.
    FloatMatrixF<9, 9> mA( answerTmp );
    Tensor4_3d A( mA ), Amod;
    Amod( i_3, j_3, k_3, l_3 ) = 1. / J0 * F0( j_3, m_3 ) * A( i_3, m_3, k_3, n_3 ) * F0( l_3, n_3 );
    answer = Amod.to_voigt_form();
}


void NLStructuralSurfaceElement ::computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep )
{
    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray normal = { 0., 0., 1. }; // the reference normal in the frame of the element
    answer = cs->giveRealSurfaceStresses( strain, normal, gp, tStep );
}


int NLStructuralSurfaceElement::checkConsistency()
{
    int result = NLStructuralElement::checkConsistency();
    // C0 interpolation gives no second gradient, which a material with bending stiffness needs
    for ( auto &iRule : this->integrationRulesArray ) {
        for ( GaussPoint *gp : *iRule ) {
            auto mat = dynamic_cast<StructuralMaterial *>( this->giveCrossSection()->giveMaterial( gp ) );
            if ( mat && mat->supportsSecondGradient() ) {
                OOFEM_ERROR( "element %d (%s): material %d (%s) needs the second gradient (bending stiffness), which only IGA surface elements of degree >= 2 provide",
                    this->giveNumber(), this->giveClassName(), mat->giveNumber(), mat->giveClassName() );
            }
        }
    }
    return result;
}

void NLStructuralSurfaceElement::computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present).
    // The element's own dof mask is used, so planar (u, v) and spatial (u, v, w) surfaces share this code.
    FloatArray u;
    this->computeVectorOf( VM_Total, tStep, u ); // Element unknowns, rotated to the local frame when provided.
    if ( initialDisplacements ) {
        u.subtract( *initialDisplacements );
    }

    // Surface displacement gradient in the element's reference surface frame.
    FloatMatrix B;
    this->computeBHmatrixAt( gp, B );
 
    answer.beProductOf( B, u );

    // Surface deformation gradient F = H + I_s, where I_s = diag(1, 1, 0).
    answer.at( 1 ) += 1.0;
    answer.at( 2 ) += 1.0;

}



} // end namespace oofem
