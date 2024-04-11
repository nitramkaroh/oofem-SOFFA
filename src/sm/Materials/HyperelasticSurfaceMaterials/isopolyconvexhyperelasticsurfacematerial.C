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

#include "isopolyconvexhyperelasticsurfacematerial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "domain.h"
#include "function.h"
namespace oofem {
REGISTER_Material( IsotropicPolyconvexHyperelasticSurfaceMaterial );

IsotropicPolyconvexHyperelasticSurfaceMaterial::IsotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d ) :
    HyperElasticSurfaceMaterial( n, d ),
    gamma1( 0. ),
    gamma2( 0. ),
    alpha1( 0. ),
    alpha2( 0. ),
    delta( 0. )
{
}

FloatArrayF<9>
IsotropicPolyconvexHyperelasticSurfaceMaterial::giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> indent = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d F( vF ), P, S, Identity( indent ), Pgamma1, Pgamma2, Palpha1, Palpha2, Pdelta;

    double gamma_t1, gamma_t2, alpha_t1, alpha_t2, delta_t;
    if(this->gamma_ltf == 0) {
        gamma_t1 = this->gamma1;
        gamma_t2  = this->gamma2;
        alpha_t1 = this->alpha1;
        alpha_t2 = this->alpha2;
        delta_t   = this->delta;
    } else {
        gamma_t1 = this->gamma1 * domain->giveFunction(gamma_ltf)->evaluateAtTime(tStep->giveIntrinsicTime());
        gamma_t2  = this->gamma2 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        alpha_t1 = this->alpha1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        alpha_t2 = this->alpha2 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        delta_t  = this->delta * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    // compute the first Piola-Kirchhoff
    Pgamma1( i_3, j_3 ) = gamma_t1 * this->compute_surface_cofactor( F )( i_3, j_3 ); 
    Pgamma2( i_3, j_3 ) = -gamma_t2 * Identity( i_3, j_3 ); 

    Palpha1( i_3, j_3 ) =  alpha_t1 * this->compute_surface_d_normF_dF( F )( i_3, j_3 ) ;
    Palpha2( i_3, j_3 ) = - alpha_t2 / sqrt( 2 ) * Identity( i_3, j_3 );

    Pdelta( i_3, j_3 ) = delta_t * this->compute_surface_d_traceU_dF( F )( i_3, j_3 );
    //if (1 ) {
    //    Tensor2_3d Pdelta2;
    //    Pdelta2( i_3, j_3 ) = Pdelta( i_3, j_3 ) / delta_t;
    //    Pdelta2.printYourself();
    //    F.printYourself();
    //}
    
    P( i_3, j_3 ) = Pgamma1( i_3, j_3 ) + Pgamma2( i_3, j_3 ) + Palpha1( i_3, j_3 ) + Palpha2( i_3, j_3 ) + Pdelta( i_3, j_3 );
    
    // compute Cauchy stress vector
    S( i_3, j_3 ) = ( 1 / this->compute_surface_determinant(F) ) * P( i_3, k_3 ) * F( j_3, k_3 );

    auto vP = P.to_voigt_form();
    auto vS = S.to_voigt_form();
    // update gp
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );
    status->letTempCVectorBe( vS );

    return vP;
}


FloatMatrixF<9, 9>
IsotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );

    Tensor2_3d F( vF );
    Tensor4_3d A;
    double gamma_t1, alpha_t1, delta_t;
    //
    if ( this->gamma_ltf == 0 ) {
        gamma_t1 = this->gamma1;
        alpha_t1 = this->alpha1;
        delta_t  = this->delta;
    } else {
        gamma_t1 = this->gamma1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        alpha_t1 = this->alpha1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        delta_t  = this->delta * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }

    A( i_3, j_3, k_3, l_3 ) = gamma_t1 * this->compute_surface_dCof_dF( F )( i_3, j_3, k_3, l_3 ) + 
                              alpha_t1 * this->compute_surface_d2_normF_dF2( F )( i_3, j_3, k_3, l_3 ) + 
                              delta_t * this->compute_surface_d2_traceU_dF2( F )( i_3, j_3, k_3, l_3 );


    return A.to_voigt_form();
}




MaterialStatus *
IsotropicPolyconvexHyperelasticSurfaceMaterial::CreateStatus( GaussPoint *gp ) const
{
    return new StructuralMaterialStatus( gp );
}


void IsotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( InputRecord &ir )
{
    HyperElasticSurfaceMaterial::initializeFrom( ir );
    IR_GIVE_OPTIONAL_FIELD( ir, gamma1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gamma1 );
    IR_GIVE_OPTIONAL_FIELD( ir, gamma2, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gamma2 );
    IR_GIVE_OPTIONAL_FIELD( ir, alpha1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha1 );
    IR_GIVE_OPTIONAL_FIELD( ir, alpha2, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha2 );
    IR_GIVE_OPTIONAL_FIELD( ir, delta, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_delta );
    IR_GIVE_OPTIONAL_FIELD( ir, gamma_ltf, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gammaLTF );

    // In case of old input fiule with only one alpha parameter
    double alphaTemp, gammaTemp;
    IR_GIVE_OPTIONAL_FIELD( ir, alphaTemp, "alpha" );
    IR_GIVE_OPTIONAL_FIELD( ir, gammaTemp, "gamma" );
    if ( alphaTemp != 0 && alpha1 == 0 && alpha2 == 0 ) {
        alpha1 = alphaTemp;
    };
    if ( gammaTemp != 0 && gamma1 == 0 && gamma2 == 0 ) {
        gamma1 = gammaTemp;
    };


}
} // end namespace oofem
