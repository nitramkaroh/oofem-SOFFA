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
 *               Copyright (C) 1993 - 2020   Borek Patzak
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

#include "gradientmooneyrivlin.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "domain.h"
#include "structuralelement.h"
#include "oofemlib/tensor/tensor3.h"
#include "oofemlib/tensor/tensor5.h"
#include "oofemlib/tensor/tensor6.h"

namespace oofem {
REGISTER_Material( GradientMooneyRivlinCompressibleMaterial );

GradientMooneyRivlinCompressibleMaterial::GradientMooneyRivlinCompressibleMaterial( int n, Domain *d ) :
    MooneyRivlinCompressibleMaterial( n, d ),
    omega( 0. )
{
}

void GradientMooneyRivlinCompressibleMaterial::initializeFrom( InputRecord &ir )
{
    // Initialize standard Mooney-Rivlin parameters
    MooneyRivlinCompressibleMaterial::initializeFrom( ir );

    // Read the gradient penalty/stiffness parameter
    IR_GIVE_OPTIONAL_FIELD( ir, omega, _IFT_GradientMooneyRivlinCompressibleMaterial_omega );
}

FloatArrayF<9>
GradientMooneyRivlinCompressibleMaterial::giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
{
    // The standard bulk PK1 stress is independent of the gradient part in this simple formulation.
    // We delegate to the base Mooney-Rivlin class.
    FloatArrayF<9> vP = MooneyRivlinCompressibleMaterial::giveFirstPKStressVector_3d( vF, gp, tStep );

    //// Update status variables for consistent output/state tracking
    //StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    //if ( status ) {
    //    status->letTempFVectorBe( vF );
    //    status->letTempPVectorBe( vP );
    //}

    return vP;
}

FloatMatrixF<9, 9>
GradientMooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dPdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const
{
    // Standard hyperelastic tangent (dP/dF) from the base class
    return MooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dPdF( mode, gp, tStep );
}

FloatArrayF<27>
GradientMooneyRivlinCompressibleMaterial::giveSecondOrderStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vG, GaussPoint *gp, TimeStep *tStep ) const
{
    // vG is the gradient of F (nabla F).
    // For W_grad = 0.5 * omega * (nabla F : nabla F), the conjugate stress A is simply omega * nabla F.

    Tensor3_3d G( vG );
    Tensor3_3d A;

    A( i_3, j_3, k_3 ) = this->omega * G( i_3, j_3, k_3 );

    auto vA = A.to_voigt_form_27();

    // Store in status if elements rely on it (e.g., for numerical tangent checks)
    GradientStructuralMaterialStatus *status = static_cast<GradientStructuralMaterialStatus *>( this->giveStatus( gp ) );
    if ( status ) {
        status->letTempGVectorBe( vG );
        // If you have a custom status that holds higher order stress, it can be saved here.
    }

    return vA;
}

FloatMatrixF<27, 27>
GradientMooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dAddF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const
{
    // The second-order stiffness dA/d(nabla F)
    // Since A = omega * G, the derivative is an isotropic 6th-order tensor scaled by omega.

    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Identity3d( indent3d );

    Tensor6_3d dA;

    // dA_{ijkpqr} = omega * delta_{ip} * delta_{jq} * delta_{kr}
    dA( i_3, j_3, k_3, p_3, q_3, r_3 ) = this->omega * Identity3d( i_3, p_3 ) * Identity3d( j_3, q_3 ) * Identity3d( k_3, r_3 );

    return dA.to_voigt_form();
}

FloatMatrixF<27, 9>
GradientMooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dAdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const
{
    // Mixed tangent stiffness dA/dF.
    // In this purely additive volumetric energy split, A does not depend on F.
    // Return a strictly zero 27x9 matrix.

    Tensor5_3d dA;
    return dA.to_voigt_form_27x9();

}

MaterialStatus *
GradientMooneyRivlinCompressibleMaterial::CreateStatus( GaussPoint *gp ) const
{
    return new GradientStructuralMaterialStatus( gp );
}

} // end namespace oofem