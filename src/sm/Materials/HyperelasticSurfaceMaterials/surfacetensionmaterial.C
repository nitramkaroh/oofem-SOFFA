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

#include "surfacetensionmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "domain.h"
#include "function.h"
namespace oofem {
REGISTER_Material( SurfaceTensionMaterial );

SurfaceTensionMaterial::SurfaceTensionMaterial( int n, Domain *d ) :
    HyperElasticSurfaceMaterial(n,d),
    gamma(0.)
{
}

FloatArrayF<9>
SurfaceTensionMaterial::giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    Tensor2_3d F( vF ), P, S;
    F( 2, 2 ) = 0.00;
    double gamma_t;
    if(this->gamma_ltf == 0) {
      gamma_t = this->gamma;
    } else {
      gamma_t = this->gamma * domain->giveFunction(gamma_ltf)->evaluateAtTime(tStep->giveIntrinsicTime());
    }
    // compute the first Piola-Kirchhoff
    //P( i_3, j_3 ) = gamma_t * F.compute_cofactor()( i_3, j_3 );
    P( i_3, j_3 ) = gamma_t * this->compute_surface_cofactor( F )( i_3, j_3 );
    // compute Cauchy stress vector
    //S( i_3, j_3 ) = ( 1 / F.compute_determinant() ) * P( i_3, k_3 ) * F( j_3, k_3 );
    //
    auto vP = P.to_voigt_form();
    //auto vS = S.to_voigt_form();
    // update gp
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );
    //status->letTempCVectorBe( vS );
    //
    return vP;
}


FloatMatrixF<9, 9>
SurfaceTensionMaterial::give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );

    Tensor2_3d F( vF );
    F( 2, 2 ) = 0.00;
    Tensor4_3d A;
    double gamma_t;
    //
    if(this->gamma_ltf == 0) {
      gamma_t = this->gamma;
    } else {
      gamma_t = this->gamma * domain->giveFunction(gamma_ltf)->evaluateAtTime(tStep->giveIntrinsicTime());
    }
    //
    //A( i_3, j_3, k_3, l_3 ) =  gamma_t *  F.compute_tensor_cross_product()( i_3, j_3, k_3, l_3 ) ;
    A( i_3, j_3, k_3, l_3 ) = gamma_t * this->compute_surface_dCof_dF( F )( i_3, j_3, k_3, l_3 );
    //
    return A.to_voigt_form();
}




MaterialStatus *
SurfaceTensionMaterial::CreateStatus( GaussPoint *gp ) const
{
    return new StructuralMaterialStatus( gp );
}


void SurfaceTensionMaterial::initializeFrom( InputRecord &ir )
{
    HyperElasticSurfaceMaterial::initializeFrom( ir );
    IR_GIVE_FIELD( ir, gamma, _IFT_SurfaceTensionMaterial_gamma );
    IR_GIVE_OPTIONAL_FIELD( ir, gamma_ltf, _IFT_SurfaceTensionMaterial_gammaLTF );
}
} // end namespace oofem
