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
#include "structuralelement.h"
#include "sm/Elements/SurfaceElements/triangularsurface.h"
#include "sm/Elements/igasurfaceelements.h"

namespace oofem {
REGISTER_Material( IsotropicPolyconvexHyperelasticSurfaceMaterial );

IsotropicPolyconvexHyperelasticSurfaceMaterial::IsotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d ) :
    HyperElasticSurfaceMaterial( n, d ),
    gamma1( 0. ),
    gamma2( 0. ),
    alpha1( 0. ),
    alpha2( 0. ),
    delta( 0. ),
    gaminit( 0. ),
    omega1(0.)
{
}

FloatArrayF<9>
IsotropicPolyconvexHyperelasticSurfaceMaterial::giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> indent            = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    FloatArrayF<9> indent3d          = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d F, Fcurrent( vF ), P, Pfull, S, Identity3d( indent3d ), Identity( indent ), Pgamma1, Pgamma2, Palpha1, Palpha2, Pdelta, Pnh, Pnh2, dJdF;
    Tensor2_3d IdentityUD( indent );

    FloatArrayF<3> UDn( normal );
    Tensor1_3d Nt( UDn );
    // if it is IGA element
    if ( dynamic_cast<IGAElement *>(gp->giveElement()  )){
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
    }

    // get element normal in the original undeformed configuration to contruct surface identity
    FloatArray UDnormal;
    if ( gp->giveElement()->giveIPValue( UDnormal, gp, IST_ElementReferenceNormal, tStep ) ) {
        FloatArrayF<3> UDn( UDnormal );
        Tensor1_3d UDnormalT( UDn );
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - UDnormalT( i_3 ) * UDnormalT( j_3 );
    }

    //// prestrain
    //FloatArray vF0_temp;
    //Tensor2_3d F0( indent );
    //if ( gp->giveElement()->giveIPValue( vF0_temp, gp, IST_PrestrainDeformationGradient, tStep ) ) {
    //    FloatArrayF<9> vF0( vF0_temp );
    //    F0  = Tensor2_3d::Tensor2_3d( vF0 ); 
    //}
    //F( i_3, j_3 ) = Fcurrent( i_3, k_3 ) * F0( k_3, j_3 );
    F( i_3, j_3 ) = Fcurrent( i_3, j_3 );

    double gamma_t1, alpha_t1, omega_t1, lambda_t, mu_t, alpha_s_t1;
    if(this->gamma_ltf == 0) {
        gamma_t1 = this->gamma1;
        alpha_t1 = this->alpha1;
        alpha_s_t1 = this->alpha_s1;
        omega_t1 = this->omega1;
        lambda_t = this->lambda;
        mu_t     = this->mu;
    } else {
        gamma_t1 = this->gaminit + this->gamma1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        alpha_t1 = this->alpha1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        alpha_s_t1 = this->alpha_s1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        lambda_t = this->lambda * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        mu_t = this->mu * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    //alpha_t1 = this->alpha1; // Alpha not increased gradually
    //gamma_t1 = this->gamma1; // Gamma not increased gradually

    if ( fixAlpha_s ) { // this should be done for other parameters as well, but for now we only want to fix alpha_s for testing purpose
        alpha_s_t1 = this->alpha_s1;
    }
        

    // compute the first Piola-Kirchhoff
    dJdF( i_3, j_3 )    = this->compute_surface_cofactor( F, IdentityUD )( i_3, j_3 );
    Pgamma1( i_3, j_3 ) = gamma_t1 * dJdF( i_3, j_3 ); 
    Palpha1( i_3, j_3 ) =  alpha_t1 * this->compute_surface_d_normF_dF( F )( i_3, j_3 ) ;

    double J        = this->compute_surface_determinant( F );
    double logJ      = log( J );
    Pnh( i_3, j_3 )  = dJdF( i_3, j_3 ) / J * ( lambda_t * logJ - mu_t ) + mu_t * F( i_3, j_3 ); // Neo-Hookean part
    //Pnh2( i_3, j_3 )  = ( dJdF( i_3, j_3 ) * lambda_t * logJ ) / J + mu_t * F( i_3, j_3 ) - mu_t*dJdF( i_3, j_3 ) / J; // Neo-Hookean part
    //Pnh2( i_3, j_3 )  = lambda_t * logJ / J * dJdF( i_3, j_3 ) + mu_t * F( i_3, j_3 ) - mu_t * dJdF( i_3, j_3 ) / J; // Neo-Hookean part
    //Pnh2( i_3, j_3 ) = lambda_t * logJ / J * dJdF( i_3, j_3 ) + mu_t * ( F( i_3, j_3 ) - dJdF( i_3, j_3 ) / J );

    //Tensor2_3d dP;
    //dP( i_3, j_3 ) = Pnh( i_3, j_3 ) - Pnh2( i_3, j_3 );
    //double ndp     = sqrt( dP( i_3, j_3 ) * dP( i_3, j_3 ) );

    //if ( ndp > 1e-12 ) {

    //    OOFEM_WARNING( "not same! ndp = ", ndp );
    //}

        

    /////////////////////////////
    //// stress free in ref configuration
    //Pgamma1( i_3, j_3 ) = gamma_t1 * ( this->compute_surface_cofactor( F, IdentityUD )( i_3, j_3 ) - this->compute_surface_cofactor( IdentityUD, IdentityUD )( i_3, j_3 ) ); 
    //Palpha1( i_3, j_3 ) = alpha_t1 * ( this->compute_surface_d_normF_dF( F )( i_3, j_3 ) - this->compute_surface_d_normF_dF( IdentityUD )( i_3, j_3 ) );

    //Tensor2_3d FminI; // second option for stress free ref. config. 
    //FminI( i_3, j_3 )   = F( i_3, j_3 ) - IdentityUD( i_3, j_3 );
    //Palpha1( i_3, j_3 ) = alpha_t1 * this->compute_surface_d_normF_dF( FminI )( i_3, j_3 );
    ///////////////////////////////

    // part (J-1)^2
    Tensor2_3d Palpha_s;
    Palpha_s( i_3, j_3 ) = alpha_s_t1 * ( J - 1. ) * dJdF( i_3, j_3 );

    P( i_3, j_3 ) = Pgamma1( i_3, j_3 ) + Palpha1( i_3, j_3 ) + Pnh( i_3, j_3 ) + Palpha_s( i_3, j_3 );
    
    //double J0     = this->compute_surface_determinant( F0 );
    //Pfull( i_3, j_3 ) = 1. / J0 * P( i_3, k_3 ) * F0( j_3, k_3 ); 

    /////////////////////////////////////////
    //// Get the second order part
    //if ( omega_t1 > 1e-6 ) {
    //    // load objhect from status
    //    FloatArrayF<27> vG( status->giveTempGVector() );
    //    FloatArrayF<27> dI( status->giveTempdIVector() );
    //    FloatArrayF<9> vM( status->giveTempMVector() );
    //    Tensor3_3d G( vG ), gradIt( dI ); // grad. F and grad. of surface identity

    //    Tensor2_3d MK( vM ); // second order stress conjugate to kappa
    //    Tensor2_3d P2 = this->giveFirstPKSurfaceStressBending( F, IdentityUD, Nt, G, MK );
    //    ////////////////////////////////////////

    //    Pfull( i_3, j_3 ) = P( i_3, j_3 ) + P2( i_3, j_3 );
    //} else {
        Pfull( i_3, j_3 ) = P( i_3, j_3 );
    //}
    
    // compute Cauchy stress vector
    S( i_3, j_3 ) = ( 1 / this->compute_surface_determinant( F ) ) * Pfull( i_3, k_3 ) * F( j_3, k_3 );

    auto vP = Pfull.to_voigt_form();
    auto vS = S.to_voigt_form();
    // update gp
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );
    status->letTempCVectorBe( vS );

    return vP;
}

Tensor2_3d
IsotropicPolyconvexHyperelasticSurfaceMaterial::giveFirstPKSurfaceStressBending( const Tensor2_3d &F, Tensor2_3d &IdentityUD, const Tensor1_3d &normal, const Tensor3_3d &G, const Tensor2_3d &MK ) const
{
    Tensor2_3d Finv = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
    Tensor1_3d n_def = compute_surface_deformed_normal( F, normal ); // compute deformed normal
    Tensor2_3d P2; 
    P2( p_3, q_3 ) = n_def( p_3 ) * Finv( q_3, k_3 ) * G( k_3, i_3, j_3 ) * MK( i_3, j_3 );
    return P2;
}



FloatMatrixF<9, 9>
IsotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    FloatArrayF<9> indent   = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d IdentityUD( indent ), Identity3d( indent3d );

    FloatArrayF<3> UDn( normal );
    Tensor1_3d Nt( UDn );
    // if it is IGA element
    if ( dynamic_cast<IGAElement *>( gp->giveElement() ) ) {
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
    }

    // get element normal in the original undeformed configuration
    FloatArray UDnormal;
    if ( gp->giveElement()->giveIPValue( UDnormal, gp, IST_ElementReferenceNormal, tStep ) ) {
        FloatArrayF<3> UDn( UDnormal );
        Tensor1_3d UDnormalT( UDn );
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - UDnormalT( i_3 ) * UDnormalT( j_3 );
    }

    //// prestrain
    // FloatArray vF0_temp;
    Tensor2_3d Fcurrent( vF ), F, dJdF;
    // Tensor2_3d F0( indent );
    // if ( gp->giveElement()->giveIPValue( vF0_temp, gp, IST_PrestrainDeformationGradient, tStep ) ) {
    //     FloatArrayF<9> vF0( vF0_temp );
    //     F0  = Tensor2_3d::Tensor2_3d( vF0 );
    // }
    // F( i_3, j_3 ) = Fcurrent( i_3, k_3 ) * F0( k_3, j_3 );
    F( i_3, j_3 ) = Fcurrent( i_3, j_3 );

    // Tensor2_3d F( vF );
    Tensor4_3d A, Afull, Anh;
    double gamma_t1, alpha_t1, omega_t1, lambda_t, mu_t, alpha_s_t1;
    //
    if ( this->gamma_ltf == 0 ) {
        gamma_t1 = this->gamma1;
        alpha_t1 = this->alpha1;
        alpha_s_t1 = this->alpha_s1;
        omega_t1 = this->omega1;
        lambda_t = this->lambda;
        mu_t     = this->mu;
    } else {
        gamma_t1 = this->gaminit + this->gamma1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        alpha_t1 = this->alpha1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        alpha_s_t1 = this->alpha_s1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        lambda_t = this->lambda * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        mu_t     = this->mu * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    // alpha_t1 = this->alpha1; // Alpha not increased gradually
    // gamma_t1 = this->gamma1; // Gamma not increased gradually
    if ( fixAlpha_s ) { // this should be done for other parameters as well, but for now we only want to fix alpha_s for testing purpose
        alpha_s_t1 = this->alpha_s1;
    }

    // // stressed reference configuration
    Tensor4_3d d2JdFdF;
    d2JdFdF( i_3, j_3, k_3, l_3 ) = this->compute_surface_dCof_dF( F, IdentityUD )( i_3, j_3, k_3, l_3 );
   
    // neo-Hookean part
    double J        = this->compute_surface_determinant( F );
    double logJ               = log( J );
    dJdF( i_3, j_3 )          = this->compute_surface_cofactor( F, IdentityUD )( i_3, j_3 );
    
    Anh( i_3, j_3, k_3, l_3 ) = lambda_t*( logJ / J * d2JdFdF( i_3, j_3, k_3, l_3 ) + (1-logJ) / J / J * dJdF( i_3, j_3 ) * dJdF( k_3, l_3 ) ) + 
        mu_t * ( Identity3d( i_3, k_3 ) * IdentityUD( j_3, l_3 ) - 1 / J * d2JdFdF( i_3, j_3, k_3, l_3 ) + 1/ J / J * dJdF( i_3, j_3 ) * dJdF( k_3, l_3 ) ); // Neo-Hookean part

    // (J-1)^2 part
    Tensor4_3d A_alpha_s;
    A_alpha_s( i_3, j_3, k_3, l_3 ) = alpha_s_t1 * ( J - 1. ) * d2JdFdF( i_3, j_3, k_3, l_3 ) + alpha_s_t1 * dJdF( i_3, j_3 ) * dJdF( k_3, l_3 );

    // all together
    A( i_3, j_3, k_3, l_3 ) = gamma_t1 * d2JdFdF( i_3, j_3, k_3, l_3 ) + alpha_t1 * this->compute_surface_d2_normF_dF2( F, IdentityUD )( i_3, j_3, k_3, l_3 )+ 
                                Anh( i_3, j_3, k_3, l_3 ) + A_alpha_s( i_3, j_3, k_3, l_3 );


    /////////////////// stress-free ref. configuration
    //Tensor2_3d FminI; // second option for stress free ref. config.
    //FminI( i_3, j_3 )   = F( i_3, j_3 ) - IdentityUD( i_3, j_3 );
    //A( i_3, j_3, k_3, l_3 ) = gamma_t1 * this->compute_surface_dCof_dF( F, IdentityUD )( i_3, j_3, k_3, l_3 ) + alpha_t1 * this->compute_surface_d2_normF_dF2( FminI, IdentityUD )( i_3, j_3, k_3, l_3 );
    ////////////////////

    /////////////////////////////////////////
    //// Get the second order part 
    //if ( omega_t1 > 1e-6 ) {
    //    // construct some useful objects
    //    Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
    //    Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
    //    FloatArrayF<27> vG( status->giveTempGVector() );
    //    Tensor3_3d G( vG ); // Gradient of F
    //    FloatArrayF<9> vM( status->giveTempMVector() );
    //    Tensor2_3d Mk( vM ); // Second order stress conjugate to bending strain

    //    // Tensor2_3d idef;
    //    // idef( j_3, k_3 ) = F( j_3, m_3 ) * Finv( m_3, k_3 );


    //    // Compute stiffnesses
    //    FloatMatrixF<9, 9> vdAK( status->giveTempD2eDK2_Matrix() ); // second gradient
    //    Tensor4_3d dAK( vdAK ), A2, A3, A4, temp, dFinvdF;

    //    temp( p_3, q_3, i_3, j_3 ) = ( n_def( p_3 ) * Finv( q_3, k_3 ) ) * G( k_3, i_3, j_3 );
    //    A2( p_3, q_3, a_3, b_3 )   = temp( p_3, q_3, i_3, j_3 ) * dAK( i_3, j_3, r_3, s_3 ) * temp( a_3, b_3, r_3, s_3 ); // classical part

    //    Tensor3_3d dndF;
    //    dndF( p_3, a_3, b_3 ) = -Finv( b_3, p_3 ) * n_def( a_3 ); // derivative of deformed normal

    //    double Js = this->compute_surface_determinant( F );
    //    // derivative of pseudoinverse
    //    dFinvdF( i_3, j_3, k_3, l_3 ) = 1. / Js * this->compute_surface_dCof_dF( F, IdentityUD )( j_3, i_3, k_3, l_3 ) - 1. / Js / Js * this->compute_surface_cofactor( F, IdentityUD )( j_3, i_3 ) * this->compute_surface_cofactor( F, IdentityUD )( k_3, l_3 );

    //    A3( p_3, q_3, a_3, b_3 ) = dFinvdF( q_3, k_3, a_3, b_3 ) * ( G( k_3, i_3, j_3 ) * Mk( i_3, j_3 ) ) * n_def( p_3 ); // part from second variation of kappa
    //    A4( p_3, q_3, a_3, b_3 ) = dndF( p_3, a_3, b_3 ) * ( Finv( q_3, k_3 ) * ( G( k_3, i_3, j_3 ) * Mk( i_3, j_3 ) ) ); // part from second variation of kappa

    //    // Tensor4_3d A34;
    //    // A34( i_3, j_3, k_3, l_3 ) = A3( i_3, j_3, k_3, l_3 ) + A4( i_3, j_3, k_3, l_3 );
    //    // auto A34m = FloatMatrix( A34.to_voigt_form() );
    //    ////  A34m.printYourself();

    //    ////////////////////////////////////////

    //    ////// Transform to to original reference frame, just for testing
    //    //// Tensor4_3d A_rotated;
    //    //// Tensor2_3d F_rotated;
    //    //// this->TransformTensor2( F, F_rotated, gp );
    //    //// this->TransformTensor4( A, A_rotated, gp );

    //    //// double J0 = this->compute_surface_determinant( F0 );
    //    //// Afull( i_3, j_3, k_3, l_3 ) = 1./J0* F0( j_3, m_3 ) * A( i_3, m_3, k_3, n_3 ) * F0( l_3, n_3 );

    //    Afull( i_3, j_3, k_3, l_3 ) = A( i_3, j_3, k_3, l_3 ) + A2( i_3, j_3, k_3, l_3 ) + A3( i_3, j_3, k_3, l_3 ) + A4( i_3, j_3, k_3, l_3 );
    //} else {
        Afull( i_3, j_3, k_3, l_3 ) = A( i_3, j_3, k_3, l_3 );
    //}


    return Afull.to_voigt_form();
}

/////////
// Second order stress and stiffnesses

//FloatArrayF<27>
//IsotropicPolyconvexHyperelasticSurfaceMaterial::giveSecondOrderSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//// returns 27 components of the second order stress conjugate to gradient of F
//{    
//    // create identities
//    FloatArrayF<9> indent            = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
//    FloatArrayF<9> indent3d          = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
//    Tensor2_3d Identity3d( indent3d );
//    Tensor2_3d IdentityUD( indent );
//
//    FloatArrayF<3> UDn( normal );
//    Tensor1_3d Nt( UDn ); // normal tensor
//    // if it is IGA element
//    auto iga_el = dynamic_cast<IGAElement *>( gp->giveElement() );
//    if ( iga_el ) {
//        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
//    }
//
//    // Get material parameter
//    double omega_t1;
//    if ( this->gamma_ltf == 0 ) {
//        omega_t1 = this->omega1;
//    } else {
//        omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
//    }
//    if ( abs( omega_t1 ) > 1e-8 ) { // check if the interpolation is spline and at least quadratic
//        int intOrd = iga_el->giveInterpolation()->giveInterpolationOrder();
//        if ( !iga_el || iga_el->giveInterpolation()->giveInterpolationOrder() <= 1) {
//            OOFEM_ERROR( "For second order surface elements the interpolation must be C1" );
//        }
//    } 
//
//    // construct some useful objects
//    Tensor2_3d F( vF ); // deformation gradient
//    Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
//    Tensor1_3d n_def   = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
//    Tensor3_3d A, G( vG ), gradIt( gradI );
//
//    // Compute curvature measures
//    // formulation based on Dort or Tomassetti using bending strain
//    Tensor2_3d K, K0; // gradient of deformed normal nad bending strain
//    K( j_3, k_3 ) =  - n_def( i_3 ) * G( i_3, j_3, k_3 ) ; // Deformed curvature in referential frame
//    K0( j_3, k_3 ) = - Nt( i_3 ) * gradIt( i_3, j_3, k_3 ); // gradient of undeformed normal
//
//
//    ///////////////////
//    //double trK = K( l_3, l_3 );
//    //Tensor2_3d C;
//    //C( k_3, l_3 ) = F( j_3, k_3 ) * F( j_3, l_3 );
//    //////////////////
//
//    // compute the second order tensor conjugate to K
//    Tensor2_3d MK;
//    MK( i_3, j_3 ) = omega_t1 * ( K( l_3, l_3 ) - K0( k_3, k_3 ) ) * IdentityUD( i_3, j_3 );
//    //MK( i_3, j_3 ) = omega_t1 * ( K( l_3, l_3 ) ) * IdentityUD( i_3, j_3 );
//
//    // obtain the third order tensor from Mk conjugarte to grad. F
//    A( i_3, j_3, k_3 ) = -n_def( i_3 ) * MK( j_3, k_3 ); 
//    auto vA = A.to_voigt_form_27();
//
//    // Save object to status, note that this function needs to be called before the standar surface PK stress is evaluated
//    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
//    status->letTempFVectorBe( vF );
//    status->letTempGVectorBe( vG );
//    status->letTempdIVectorBe( gradI );
//    status->letTempMVectorBe( MK.to_voigt_form() );
//    status->letTempKVectorBe( K.to_voigt_form() );
//
//    return vA;
//}


//FloatMatrixF<27, 27>
//IsotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dAddF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//// returns the 27x27 second order tangent stiffness matrix - dA/ddF (derivative of second order stress wrt to gradient of F)
//{
//    // load objects from status
//    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
//    FloatArrayF<9> vF( status->giveTempFVector() );
//    FloatArrayF<27> vG( status->giveTempGVector() ); // second gradient
//
//    // create identities
//    FloatArrayF<9> indent   = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
//    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
//    Tensor2_3d IdentityUD( indent ), Identity3d( indent3d );
//
//    // if it is IGA element
//    FloatArrayF<3> UDn( normal );
//    Tensor1_3d Nt( UDn );
//    if ( dynamic_cast<IGAElement *>( gp->giveElement() ) ) {
//        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
//    }
//
//    // Get material parameter
//    double omega_t1;
//    if ( this->gamma_ltf == 0 ) {
//        omega_t1 = this->omega1;
//    } else {
//        omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
//    }
//
//    // stiffnes conjugate to MK
//    Tensor4_3d dAK;
//    dAK( i_3, j_3, k_3, l_3 ) = omega_t1 * IdentityUD( i_3, j_3 ) * IdentityUD( k_3, l_3 );
//    //dAK( i_3, j_3, k_3, l_3 ) = omega_t1 * Identity3d( i_3, j_3 ) * Identity3d( k_3, l_3 );
//
//    // compute deformed normal
//    Tensor2_3d F( vF );
//    Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
//    Tensor3_3d G( vG ); // Get gradient of F
//
//    // Compute second order stiffness
//    Tensor6_3d dA;
//    // created in a stupid way by contraction
//    dA( i_3, j_3, k_3, l_3, m_3, n_3 ) = ( (n_def( i_3 ) * n_def( l_3 )) * Identity3d( j_3, s_3 ) ) * dAK( s_3, k_3, m_3, n_3 ); 
//    /////////////////////
//
//    status->lettempD2eDK2_MatrixBe( dAK.to_voigt_form() );
//    return dA.to_voigt_form();
//}



//FloatMatrixF<27, 9>
//IsotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dAdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//// returns the 27x9 mixed tangent stiffness matrix - dA/dF (derivative of second order stress wrt to F)
//{
//    // load objects from status
//    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
//    FloatArrayF<9> vF( status->giveTempFVector() );
//    FloatArrayF<9> vM( status->giveTempMVector() ); // second order stress conjugate to curvature tensor from Silhous
//    FloatArrayF<27> vG( status->giveTempGVector() ); // second gradient
//    FloatMatrixF<9,9> vdAK( status->giveTempD2eDK2_Matrix() ); // second gradient
//
//    // create identities
//    FloatArrayF<9> indent   = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
//    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
//    Tensor2_3d IdentityUD( indent ), Identity3d( indent3d );
//
//    // if it is IGA element
//    FloatArrayF<3> UDn( normal );
//    Tensor1_3d Nt( UDn );
//    if ( dynamic_cast<IGAElement *>( gp->giveElement() ) ) {
//        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
//    }
//
//    // Get material parameter
//    double omega_t1; 
//    if ( this->gamma_ltf == 0 ) {
//        omega_t1 = this->omega1;
//    } else {
//        omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
//    }
//
//    // stiffnes conjugate to MK
//    Tensor4_3d dAK( vdAK );
//    //dAK( i_3, j_3, k_3, l_3 ) = omega_t1 * IdentityUD( i_3, j_3 ) * IdentityUD( k_3, l_3 );
//
//    // construct some useful objects
//    Tensor2_3d F( vF ), Mk( vM ); // get second order stress tensor and F
//    Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
//    Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
//    Tensor3_3d G( vG ); // Get gradient of F
//
//    // Compute the mixed stiffness
//    // creted in a stupid way with the extra identity!
//    Tensor5_3d dA, dA1, dA2, dA3, dA4;
//    dA1( k_3, i_3, j_3, p_3, q_3 ) = 1. * ( Mk( i_3, j_3 ) * n_def( r_3 ) ) * ( Identity3d( r_3, p_3 ) * Finv( q_3, k_3 ) ); 
//    dA2( k_3, i_3, j_3, a_3, b_3 ) = -1. * ( dAK( i_3, j_3, r_3, s_3 ) * G( c_3, r_3, s_3 ) ) * ( n_def( k_3 ) * n_def( a_3 ) * Finv( b_3, c_3 ) );
//
//    dA( p_3, q_3, l_3, i_3, j_3 ) = dA1( p_3, q_3, l_3, i_3, j_3 ) + dA2( p_3, q_3, l_3, i_3, j_3 );
//    return dA.to_voigt_form_27x9();
//}



MaterialStatus *
IsotropicPolyconvexHyperelasticSurfaceMaterial::CreateStatus( GaussPoint *gp ) const
{
    //return new StructuralMaterialStatus( gp );
    return new StructuralSurfaceMaterialStatus( gp );
}


void IsotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( InputRecord &ir )
{
    HyperElasticSurfaceMaterial::initializeFrom( ir );
    IR_GIVE_OPTIONAL_FIELD( ir, gamma1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gamma1 );
    IR_GIVE_OPTIONAL_FIELD( ir, gamma2, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gamma2 );
    IR_GIVE_OPTIONAL_FIELD( ir, alpha1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha1 );
    IR_GIVE_OPTIONAL_FIELD( ir, alpha2, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha2 );
    IR_GIVE_OPTIONAL_FIELD( ir, alpha_s1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha_s1 );
    IR_GIVE_OPTIONAL_FIELD( ir, delta, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_delta );
    IR_GIVE_OPTIONAL_FIELD( ir, gamma_ltf, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gammaLTF );
    IR_GIVE_OPTIONAL_FIELD( ir, gaminit, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gaminit );
    IR_GIVE_OPTIONAL_FIELD( ir, omega1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_omega1 );
    IR_GIVE_OPTIONAL_FIELD( ir, lambda, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_lambda );
    IR_GIVE_OPTIONAL_FIELD( ir, mu, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_mu );

    IR_GIVE_OPTIONAL_FIELD( ir, fixGamma, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_FixGamma );
    IR_GIVE_OPTIONAL_FIELD( ir, fixOmega, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_FixOmega );
    IR_GIVE_OPTIONAL_FIELD( ir, fixAlpha_s, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_FixAlpha_s );

    // In case of old input fiule with only one alpha parameter
    double alphaTemp, gammaTemp;
    alphaTemp = 0.;
    gammaTemp = 0.;
    IR_GIVE_OPTIONAL_FIELD( ir, alphaTemp, "alpha" );
    IR_GIVE_OPTIONAL_FIELD( ir, gammaTemp, "gamma" );
    if ( alphaTemp != 0 && alpha1 == 0 && alpha2 == 0 ) {
        alpha1 = alphaTemp;
    };
    if ( gammaTemp != 0 && gamma1 == 0 && gamma2 == 0 ) {
        gamma1 = gammaTemp;
    };
}

double IsotropicPolyconvexHyperelasticSurfaceMaterial::computeEnergy( const Tensor2_3d &F, GaussPoint *gp, TimeStep *tStep ) const
{
    double gamma_t1, gamma_t2, alpha_t1, alpha_t2, delta_t, beta_t;
    if ( this->gamma_ltf == 0 ) {
        gamma_t1 = this->gamma1;
        alpha_t1 = this->alpha1;
        delta_t  = this->delta;
    } else {
        gamma_t1 = this->gamma1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        alpha_t1 = this->alpha1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        delta_t  = this->delta * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }

    Tensor2_3d C;
    C( j_3, k_3 ) = F( i_3, j_3 ) * F( i_3, k_3 );

    FloatMatrix Cmat, Uc, Sc, Vc;
    C.toFloatMatrix( Cmat );
    Cmat.resizeWithData( 2, 2 );
    Cmat.computeSVD2x2( Uc, Sc, Vc );
    double trU = sqrt( Sc.at( 1, 1 ) ) + sqrt( Sc.at( 2, 2 ) );

    double energyIso = gamma_t1 * this->compute_surface_determinant( F ) + alpha_t1 * this->compute_surface_NormF( F ) + delta_t * trU;
    return energyIso;
}


Tensor2_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::DoRTRtTransform( const Tensor2_3d &T, const FloatMatrix &R ) const
{
    FloatMatrix Tm, temp1, Tmr;
    T.toFloatMatrix( Tm );
    temp1.beProductTOf( Tm, R );
    Tmr.beProductOf( R, temp1 );

    return Tensor2_3d::Tensor2_3d( Tmr );
}

void IsotropicPolyconvexHyperelasticSurfaceMaterial::TransformTensor1( const Tensor1_3d &T, Tensor1_3d &Ttrans, GaussPoint *gp ) const
{
    TriangularSurface *SurfEl = dynamic_cast<TriangularSurface *>( gp->giveElement() );
    FloatMatrix RefRotation, Fmr;
    Tensor2_3d F_rot;
    if ( SurfEl && SurfEl->givePrescribedRotation( RefRotation ) ) {
        FloatMatrix R;
        SurfEl->giveRotationMatrix( R );
        IntArray pos( { 1, 2, 3 } );
        FloatMatrix R3x3, Fmat, FR, RtFR;
        R3x3.beSubMatrixOf( R, pos, pos );

        // Get the rotation from current ref. configuration to old ref. configuration
        FloatMatrix Rotation_Ref_2_Ref_Old;
        Rotation_Ref_2_Ref_Old.beProductTOf( RefRotation, R3x3 );

        Tensor2_3d Rot( Rotation_Ref_2_Ref_Old );
        Ttrans( i_3) = T( m_3 ) * Rot( i_3, m_3 ) ;
    } else {
        Ttrans = T;
    }
}


void IsotropicPolyconvexHyperelasticSurfaceMaterial::TransformTensor2( const Tensor2_3d &T, Tensor2_3d &Ttrans, GaussPoint *gp ) const
{
    TriangularSurface *SurfEl = dynamic_cast<TriangularSurface *>( gp->giveElement() );
    FloatMatrix RefRotation, Fmr;
    Tensor2_3d F_rot;
    if ( SurfEl && SurfEl->givePrescribedRotation( RefRotation ) ) {
        FloatMatrix R;
        SurfEl->giveRotationMatrix( R );
        IntArray pos( { 1, 2, 3 } );
        FloatMatrix R3x3, Fmat, FR, RtFR;
        R3x3.beSubMatrixOf( R, pos, pos );

        // Get the rotation from current ref. configuration to old ref. configuration
        FloatMatrix Rotation_Ref_2_Ref_Old;
        Rotation_Ref_2_Ref_Old.beProductTOf( RefRotation, R3x3 );
        Ttrans    = this->DoRTRtTransform( T, Rotation_Ref_2_Ref_Old );
    } else {
        Ttrans = T;
    }

}

void IsotropicPolyconvexHyperelasticSurfaceMaterial::TransformTensor4( const Tensor4_3d &T, Tensor4_3d &Ttrans, GaussPoint *gp ) const
{
    TriangularSurface *SurfEl = dynamic_cast<TriangularSurface *>( gp->giveElement() );
    FloatMatrix RefRotation, Fmr;
    Tensor2_3d F_rot;
    if ( SurfEl && SurfEl->givePrescribedRotation( RefRotation ) ) {
        FloatMatrix R;
        SurfEl->giveRotationMatrix( R );
        IntArray pos( { 1, 2, 3 } );
        FloatMatrix R3x3, Fmat, FR, RtFR;
        R3x3.beSubMatrixOf( R, pos, pos );

        // Get the rotation from current ref. configuration to old ref. configuration
        FloatMatrix Rotation_Ref_2_Ref_Old;
        Rotation_Ref_2_Ref_Old.beProductTOf( RefRotation, R3x3 );

        Tensor2_3d Rot( Rotation_Ref_2_Ref_Old );
        Ttrans( i_3, j_3, k_3, l_3 ) = T( m_3, n_3, o_3, p_3 ) * Rot( i_3, m_3 ) * Rot( j_3, n_3 ) * Rot( k_3, o_3 ) * Rot( l_3, p_3 );

    } else {
        Ttrans = T;
    }
}

Tensor4_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::computeTangentPK1Num( const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, double EPSILON, const Tensor2_3d &P, const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    Tensor4_3d tangent;
    Tensor2_3d Fplus, F( vF );
    FloatArrayF<9> vPplus;

    int dim = 2;
    for ( int k = 0; k <= dim; ++k ) {
        for ( int l = 0; l <= dim; ++l ) {
            Tensor2_3d dF;
            dF( k, l )        = EPSILON;
            Fplus( i_3, j_3 ) = F( i_3, j_3 ) + dF( i_3, j_3 );

            // compute second order stress from modified F so that it gets updated in the status
            auto M = giveSecondOrderSurfaceStressVector_3d( Fplus.to_voigt_form(), vG, gradI, normal, gp, tStep );
            vPplus = giveFirstPKSurfaceStressVector_3d( Fplus.to_voigt_form(), normal, gp, tStep );
            Tensor2_3d Pplus( vPplus );

            for ( int i = 0; i <= 2; ++i ) {
                for ( int j = 0; j <= 2; ++j ) {
                    double value          = ( Pplus( i, j ) - P( i, j ) ) / EPSILON;
                    tangent( i, j, k, l ) = value;
                }
            }
        }
    }
    return tangent;
}

Tensor3_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::dndFnum( double EPSILON, const Tensor2_3d &F, const Tensor1_3d &normal ) const
{
    Tensor3_3d tangent;
    Tensor2_3d Fplus;

    int dim = 2;
    for ( int k = 0; k <= dim; ++k ) {
        for ( int l = 0; l <= dim; ++l ) {
            Tensor2_3d dF;
            dF( k, l )        = EPSILON;
            Fplus( i_3, j_3 ) = F( i_3, j_3 ) + dF( i_3, j_3 );

            Tensor1_3d nplus = compute_surface_deformed_normal( Fplus, normal ); // compute deformed normal;
            Tensor1_3d n     = compute_surface_deformed_normal( F, normal );

            for ( int i = 0; i <= 2; ++i ) {
                double value       = ( nplus( i ) - n( i ) ) / EPSILON;
                tangent( i, k, l ) = value;
            }
        }
    }

    return tangent;
}

Tensor4_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::dFinvdFnum( double EPSILON, const Tensor2_3d &F, const Tensor2_3d &IdentityUD ) const
{
    Tensor4_3d tangent;
    Tensor2_3d Fplus;

    int dim = 2;
    for ( int k = 0; k <= dim; ++k ) {
        for ( int l = 0; l <= dim; ++l ) {
            Tensor2_3d dF;
            dF( k, l )        = EPSILON;
            Fplus( i_3, j_3 ) = F( i_3, j_3 ) + dF( i_3, j_3 );

            Tensor2_3d Finvplus = this->compute_surface_pseudoinverse( Fplus, IdentityUD );
            Tensor2_3d Finv     = this->compute_surface_pseudoinverse( F, IdentityUD );

            for ( int i = 0; i <= 2; ++i ) {
                for ( int j = 0; j <= 2; ++j ) {
                    double value          = ( Finvplus( i, j ) - Finv( i, j ) ) / EPSILON;
                    tangent( i, j, k, l ) = value;
                }
            }
        }
    }

    return tangent;
}

Tensor4_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::computeTangentPK1NumB( double EPSILON, const Tensor2_3d &F, Tensor2_3d &IdentityUD, const Tensor1_3d &normal, const Tensor3_3d &G, const Tensor2_3d &MK ) const
{
    Tensor4_3d tangent;
    Tensor2_3d Fplus;

    int dim = 2;
    for ( int k = 0; k <= dim; ++k ) {
        for ( int l = 0; l <= dim; ++l ) {
            Tensor2_3d dF;
            dF( k, l )        = EPSILON;
            Fplus( i_3, j_3 ) = F( i_3, j_3 ) + dF( i_3, j_3 );

            // compute second order stress from modified F so that it gets updated in the status
            Tensor2_3d Pplus = this->giveFirstPKSurfaceStressBending( Fplus, IdentityUD, normal, G, MK );
            Tensor2_3d P     = this->giveFirstPKSurfaceStressBending( F, IdentityUD, normal, G, MK );

            for ( int i = 0; i <= 2; ++i ) {
                for ( int j = 0; j <= 2; ++j ) {
                    double value          = ( Pplus( i, j ) - P( i, j ) ) / EPSILON;
                    tangent( i, j, k, l ) = value;
                }
            }
        }
    }

    return tangent;
}

Tensor6_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::computeTangentSecondOrderNum( const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, double EPSILON, const Tensor2_3d &P, const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    Tensor6_3d tangent;
    Tensor3_3d Gplus, G( vG );

    int dim = 2;
    for ( int k = 0; k <= dim; ++k ) {
        for ( int l = 0; l <= dim; ++l ) {
            for ( int m = 0; m <= dim; ++m ) {
                Tensor3_3d dG;
                dG( k, l, m )          = EPSILON;
                Gplus( i_3, j_3, k_3 ) = G( i_3, j_3, k_3 ) + dG( i_3, j_3, k_3 );

                // compute second order stress from modified F so that it gets updated in the status
                auto vAplus = giveSecondOrderSurfaceStressVector_3d( vF, Gplus.to_voigt_form_27(), gradI, normal, gp, tStep );

                auto vA = giveSecondOrderSurfaceStressVector_3d( vF, vG, gradI, normal, gp, tStep ); // compute the original stress
                Tensor3_3d Aplus( vAplus ), A( vA );

                for ( int i = 0; i <= 2; ++i ) {
                    for ( int j = 0; j <= 2; ++j ) {
                        for ( int a = 0; a <= 2; ++a ) {
                            double value                = ( Aplus( i, j, a ) - A( i, j, a ) ) / EPSILON;
                            tangent( i, j, a, k, l, m ) = value;
                        }
                    }
                }
            }
        }
    }

    return tangent;
}

Tensor5_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::computeTangentMixedNum( const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, double EPSILON, const Tensor2_3d &P, const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    Tensor5_3d tangent;
    Tensor2_3d Fplus, F( vF );

    int dim = 2;
    for ( int k = 0; k <= dim; ++k ) {
        for ( int l = 0; l <= dim; ++l ) {
            Tensor2_3d dF;
            dF( k, l )        = EPSILON;
            Fplus( i_3, j_3 ) = F( i_3, j_3 ) + dF( i_3, j_3 );

            // compute second order stress from modified F so that it gets updated in the status
            auto vAplus = giveSecondOrderSurfaceStressVector_3d( Fplus.to_voigt_form(), vG, gradI, normal, gp, tStep );

            auto vA = giveSecondOrderSurfaceStressVector_3d( vF, vG, gradI, normal, gp, tStep ); // compute the original stress
            Tensor3_3d Aplus( vAplus ), A( vA );

            for ( int i = 0; i <= 2; ++i ) {
                for ( int j = 0; j <= 2; ++j ) {
                    for ( int a = 0; a <= 2; ++a ) {
                        double value             = ( Aplus( i, j, a ) - A( i, j, a ) ) / EPSILON;
                        tangent( i, j, a, k, l ) = value;
                    }
                }
            }
        }
    }

    return tangent;
}

Tensor5_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::computeTangentMixedNum2( const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, double EPSILON, const Tensor2_3d &P, const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    Tensor5_3d tangent;
    Tensor3_3d Gplus, G( vG );

    int dim = 2;
    for ( int k = 0; k <= dim; ++k ) {
        for ( int l = 0; l <= dim; ++l ) {
            for ( int m = 0; m <= dim; ++m ) {
                Tensor3_3d dG;
                dG( k, l, m )          = EPSILON;
                Gplus( i_3, j_3, m_3 ) = G( i_3, j_3, m_3 ) + dG( i_3, j_3, m_3 );

                // compute second order stress from modified F so that it gets updated in the status
                auto M      = giveSecondOrderSurfaceStressVector_3d( vF, Gplus.to_voigt_form_27(), gradI, normal, gp, tStep );
                auto vPplus = giveFirstPKSurfaceStressVector_3d( vF, normal, gp, tStep );

                // old versions
                auto Mold  = giveSecondOrderSurfaceStressVector_3d( vF, vG, gradI, normal, gp, tStep );
                auto vPold = giveFirstPKSurfaceStressVector_3d( vF, normal, gp, tStep );

                Tensor2_3d Pplus( vPplus );

                for ( int i = 0; i <= 2; ++i ) {
                    for ( int j = 0; j <= 2; ++j ) {
                        double value             = ( Pplus( i, j ) - P( i, j ) ) / EPSILON;
                        tangent( k, l, m, i, j ) = value;
                    }
                }
            }
        }
    }

    return tangent;
}


/////////////////////////////////////
// Actual curvature model


REGISTER_Material( IsotropicPolyconvexHyperelasticSurfaceMaterial2 );

IsotropicPolyconvexHyperelasticSurfaceMaterial2::IsotropicPolyconvexHyperelasticSurfaceMaterial2( int n, Domain *d ) :
    IsotropicPolyconvexHyperelasticSurfaceMaterial( n, d ),
    ArcLengthMaterialInterface()
{
}

FloatArrayF<9>
IsotropicPolyconvexHyperelasticSurfaceMaterial2::giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> indent                   = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    FloatArrayF<9> indent3d                 = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d F, Fcurrent( vF ), P, Pfull, S, Identity3d( indent3d ), Identity( indent ), Pgamma1, Pgamma2, Palpha1, Palpha2, Pdelta, Pnh, Pnh2, dJdF;
    Tensor2_3d IdentityUD( indent );

    FloatArrayF<3> UDn( normal );
    Tensor1_3d Nt( UDn );
    // if it is IGA element
    if ( dynamic_cast<IGAElement *>( gp->giveElement() ) ) {
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
    }

    // get element normal in the original undeformed configuration to contruct surface identity
    FloatArray UDnormal;
    if ( gp->giveElement()->giveIPValue( UDnormal, gp, IST_ElementReferenceNormal, tStep ) ) {
        FloatArrayF<3> UDn( UDnormal );
        Tensor1_3d UDnormalT( UDn );
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - UDnormalT( i_3 ) * UDnormalT( j_3 );
    }

    F( i_3, j_3 ) = Fcurrent( i_3, j_3 );

    double gamma_t1, omega_t1, alpha_s_t1;
    //if ( this->gamma_ltf == 0 ) {
    //    gamma_t1 = this->gamma1;
    //    omega_t1 = this->omega1;
    //} else {
    //    gamma_t1 = this->gaminit + this->gamma1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //    omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //}

    //if ( this->fixGamma) {
    //    gamma_t1 = this->gamma1; // omega_t1 not increased gradually
    //}  
    //if ( this->fixOmega ) {
    //    omega_t1 = this->omega1; // omega_t1 not increased gradually
    //}

    gamma_t1 = this->giveGamma( tStep );
    omega_t1 = this->giveOmega( tStep );
    alpha_s_t1 = this->giveAlpha_s( tStep );
    

    // construct some useful objects
    FloatArrayF<27> vG( status->giveTempGVector() );
    FloatArrayF<27> dI( status->giveTempdIVector() );
    FloatArrayF<9> vM( status->giveTempMVector() );

    Tensor2_3d MK( vM ); // second order stress conjugate to kappa
    Tensor3_3d G( vG ), gradIt( dI ); // grad. F and grad. of surface identity

    Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
    Tensor2_3d B;
    B( j_3, m_3 )    = Finv( j_3, k_3 ) * Finv( m_3, k_3 );
    Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
    double J         = this->compute_surface_determinant( F );

    // Compute curvature measures
    Tensor2_3d K, K0, b, b0; // gradient of deformed normal nad bending strain
    K( j_3, k_3 ) = -n_def( i_3 ) * G( i_3, j_3, k_3 ); // Deformed curvature in referential frame
    b( j_3, l_3 ) = -Finv( k_3, j_3 ) * K( k_3, m_3 ) * Finv( m_3, l_3 );
    double H      = 0.5 * b( j_3, j_3 );

    double H0;
    if ( this->use_H0 == 1 ) {
        K0( j_3, k_3 ) = -Nt( i_3 ) * gradIt( i_3, j_3, k_3 ); // gradient of undeformed normal
        b0( j_3, l_3 ) = -K0( j_3, l_3 );
        H0             = 0.5 * b0( j_3, j_3 );
    } else if ( this->use_H0 == 2 ) {
        H0 = this->H0_val; 
    } else {
        H0 = 0.;
    }

    // Compute dedF
    Tensor2_3d dedF1, dedF2, dedF, dedF_alpha;
    double e          = J * ( gamma_t1 + 0.5 * omega_t1 * ( H - H0 ) * ( H - H0 ) ); 
    dedF1( i_3, j_3 ) = e * Finv( j_3, i_3 );
    dedF_alpha( i_3, j_3 ) = alpha_s_t1*(J-1.) * J * Finv( j_3, i_3 );
    dedF2( i_3, l_3 ) = J * omega_t1 * ( H - H0 ) * Finv( j_3, i_3 ) * K( j_3, k_3 ) * B( k_3, l_3 );
    dedF( i_3, j_3 )       = dedF1( i_3, j_3 ) + dedF2( i_3, j_3 ) + dedF_alpha( i_3, j_3 );

    // Compute curvature part 
    Tensor2_3d P2 = this->giveFirstPKSurfaceStressBending( F, IdentityUD, Nt, G, MK );

    Pfull( i_3, j_3 ) = dedF( i_3, j_3 ) + P2( i_3, j_3 );

    // compute Cauchy stress vector
    S( i_3, j_3 ) = ( 1 / this->compute_surface_determinant( F ) ) * Pfull( i_3, k_3 ) * F( j_3, k_3 );

    auto vP = Pfull.to_voigt_form();
    auto vS = S.to_voigt_form();
    // update gp
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );
    status->letTempCVectorBe( vS );

    return vP;
}


FloatArrayF<27>
IsotropicPolyconvexHyperelasticSurfaceMaterial2::giveSecondOrderSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns 27 components of the second order stress conjugate to gradient of F
{
    // create identities
    FloatArrayF<9> indent   = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Identity3d( indent3d );
    Tensor2_3d IdentityUD( indent );

    FloatArrayF<3> UDn( normal );
    Tensor1_3d Nt( UDn ); // normal tensor
    // if it is IGA element
    auto iga_el = dynamic_cast<IGAElement *>( gp->giveElement() );
    if ( iga_el ) {
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
    }

    // Get material parameter
    double omega_t1;
    //if ( this->gamma_ltf == 0 ) {
    //    omega_t1 = this->omega1;
    //} else {
    //    omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //}

    //if ( this->fixOmega ) {
    //    omega_t1 = this->omega1; // omega_t1 not increased gradually
    //}

    omega_t1 = this->giveOmega( tStep );

    if ( abs( omega_t1 ) > 1e-8 ) { // check if the interpolation is spline and at least quadratic
        int intOrd = iga_el->giveInterpolation()->giveInterpolationOrder();
        if ( !iga_el || iga_el->giveInterpolation()->giveInterpolationOrder() <= 1 ) {
            OOFEM_ERROR( "For second order surface elements the interpolation must be C1" );
        }
    }

    // construct some useful objects
    Tensor2_3d F( vF ); // deformation gradient
    Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
    Tensor2_3d B;
    B( j_3, m_3 )    = Finv( j_3, k_3 ) * Finv( m_3, k_3 );
    Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
    Tensor3_3d A, G( vG ), gradIt( gradI );
    double J = this->compute_surface_determinant( F );

    // Compute curvature measures
    // formulation based on Dort or Tomassetti using bending strain
    Tensor2_3d K, K0, b, b0; // gradient of deformed normal nad bending strain
    K( j_3, k_3 ) = -n_def( i_3 ) * G( i_3, j_3, k_3 ); // Deformed curvature in referential frame
    b( j_3, l_3 ) = -Finv( k_3, j_3 ) * K( k_3, m_3 ) * Finv( m_3, l_3 );
    double H      = 0.5 * b( j_3, j_3 );

    double H0;
    if ( this->use_H0 == 1 ) {
        K0( j_3, k_3 ) = -Nt( i_3 ) * gradIt( i_3, j_3, k_3 ); // gradient of undeformed normal
        b0( j_3, l_3 ) = -K0( j_3, l_3 );
        H0             = 0.5 * b0( j_3, j_3 );
    } else if ( this->use_H0 == 2 ) {
        H0 = this->H0_val;
    } else {
        H0 = 0.;
    }

    // compute the second order tensor conjugate to K, MK = dedK
    Tensor2_3d MK;
    MK( i_3, j_3 ) = -0.5 * omega_t1 * J * ( H - H0) * B( i_3, j_3 );

    // obtain the third order tensor from Mk conjugarte to grad. F
    A( i_3, j_3, k_3 ) = -n_def( i_3 ) * MK( j_3, k_3 );
    auto vA            = A.to_voigt_form_27();

    // Save object to status, note that this function needs to be called before the standar surface PK stress is evaluated
    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    status->letTempFVectorBe( vF );
    status->letTempGVectorBe( vG );
    status->letTempdIVectorBe( gradI );
    status->letTempMVectorBe( MK.to_voigt_form() );
    status->letTempKVectorBe( K.to_voigt_form() );

    return vA;
}

FloatMatrixF<9, 9>
IsotropicPolyconvexHyperelasticSurfaceMaterial2::give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    FloatArrayF<9> indent   = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d IdentityUD( indent ), Identity3d( indent3d );

    FloatArrayF<3> UDn( normal );
    Tensor1_3d Nt( UDn );
    // if it is IGA element
    if ( dynamic_cast<IGAElement *>( gp->giveElement() ) ) {
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
    }

    // get element normal in the original undeformed configuration
    FloatArray UDnormal;
    if ( gp->giveElement()->giveIPValue( UDnormal, gp, IST_ElementReferenceNormal, tStep ) ) {
        FloatArrayF<3> UDn( UDnormal );
        Tensor1_3d UDnormalT( UDn );
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - UDnormalT( i_3 ) * UDnormalT( j_3 );
    }

    Tensor2_3d F( vF );

    // Tensor2_3d F( vF );
    Tensor4_3d A, Afull, Anh;
    double gamma_t1, omega_t1, alpha_s_t1;
    //
    //if ( this->gamma_ltf == 0 ) {
    //    gamma_t1 = this->gamma1;
    //    omega_t1 = this->omega1;
    //} else {
    //    gamma_t1 = this->gaminit + this->gamma1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //    omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //}
    //if ( this->fixGamma ) {
    //    gamma_t1 = this->gamma1; // omega_t1 not increased gradually
    //}
    //if ( this->fixOmega ) {
    //    omega_t1 = this->omega1; // omega_t1 not increased gradually
    //}

    gamma_t1 = this->giveGamma( tStep );
    omega_t1 = this->giveOmega( tStep );
    alpha_s_t1 = this->giveAlpha_s( tStep );


    // construct some useful objects
    Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
    double J         = this->compute_surface_determinant( F );
    Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
    Tensor2_3d B;
    B( j_3, m_3 ) = Finv( j_3, k_3 ) * Finv( m_3, k_3 );

    FloatArrayF<27> vG( status->giveTempGVector() );
    FloatArrayF<9> vM( status->giveTempMVector() );
    FloatArrayF<27> dI( status->giveTempdIVector() );
    FloatMatrixF<9, 9> vdAK( status->giveTempD2eDK2_Matrix() ); // second gradient

    Tensor3_3d G( vG ), gradIt( dI ); // Gradient of F 
    Tensor2_3d Mk( vM ); // Second order stress conjugate to bending strain
    Tensor4_3d d2e_dKdK( vdAK );

    Tensor2_3d K, K0, b, b0; // gradient of deformed normal nad bending strain
    K( j_3, k_3 ) = -n_def( i_3 ) * G( i_3, j_3, k_3 ); // Deformed curvature in referential frame
    b( j_3, l_3 ) = -Finv( k_3, j_3 ) * K( k_3, m_3 ) * Finv( m_3, l_3 );
    double H      = 0.5 * b( j_3, j_3 );

    double H0;
    if ( this->use_H0 == 1 ) {
        K0( j_3, k_3 ) = -Nt( i_3 ) * gradIt( i_3, j_3, k_3 ); // gradient of undeformed normal
        b0( j_3, l_3 ) = -K0( j_3, l_3 );
        H0             = 0.5 * b0( j_3, j_3 );
    } else if ( this->use_H0 == 2 ) {
        H0 = this->H0_val;
    } else {
        H0 = 0.;
    }

    double e = J * ( gamma_t1 + 0.5 * omega_t1 * ( H - H0 ) * ( H - H0 ) ); 

    // compuite dedF
    Tensor2_3d dHdF;
    dHdF( i_3, l_3 ) = Finv( j_3, i_3 ) * K( j_3, k_3 ) * B( k_3, l_3 );

    Tensor2_3d dedF1, dedF2, dedF;
    dedF1( i_3, j_3 ) = e * Finv( j_3, i_3 );
    dedF2( i_3, l_3 ) = J * omega_t1 * ( H - H0 ) * dHdF( i_3, l_3 );
    dedF( i_3, j_3 )  = dedF1( i_3, j_3 ) + dedF2( i_3, j_3 );

    // compute stiffnesses
    // Compute d2e_dFdF
    Tensor4_3d d2H_dFdF, d2e_dFdF;
    d2H_dFdF( a_3, b_3, c_3, d_3 ) = ( -Finv( o_3, c_3 ) * Finv( d_3, a_3 ) + B( o_3, d_3 ) * n_def( a_3 ) * n_def( c_3 ) ) * K( o_3, p_3 ) *B( p_3, b_3 )+
         Finv( o_3, a_3 )* K( o_3, p_3 )* ( -Finv( p_3, c_3 ) * B( b_3, d_3 ) - B( p_3, d_3 ) *Finv( b_3,c_3 ) ) ;

    d2e_dFdF( a_3, b_3, c_3, d_3 ) = Finv( b_3, a_3 ) * dedF( c_3, d_3 ) + 
        e*( -Finv( b_3, c_3 ) * Finv( d_3, a_3 ) + B( b_3, d_3 ) * n_def( a_3 ) * n_def( c_3 ) ) +
        J * omega_t1 * dHdF( a_3, b_3 ) * dHdF( c_3, d_3 ) + 
        J * omega_t1 * ( H - H0 ) * dHdF( a_3, b_3 ) * Finv( d_3, c_3 ) + 
        J * omega_t1 * ( H - H0 ) * d2H_dFdF( a_3, b_3, c_3, d_3 );

    // part related to term alpha_s*(J-1)^2, first get part related to surface tension
    Tensor2_3d dedF_st, dedF_alpha;
    Tensor4_3d d2e_dFdF_st, d2e_dFdF_alpha;
    dedF_st( i_3, j_3 ) = J * Finv( j_3, i_3 );
    d2e_dFdF_st( a_3, b_3, c_3, d_3 ) = Finv( b_3, a_3 ) * dedF_st( c_3, d_3 ) + J * ( -Finv( b_3, c_3 ) * Finv( d_3, a_3 ) + B( b_3, d_3 ) * n_def( a_3 ) * n_def( c_3 ) );
    // alpha part
    //dedF_alpha( i_3, j_3 ) = alpha_s_t1 * ( J - 1. ) * dedF_st( i_3, j_3 );
    d2e_dFdF_alpha( a_3, b_3, c_3, d_3 ) = alpha_s_t1 * ( J - 1. ) * d2e_dFdF_st( a_3, b_3, c_3, d_3 ) + alpha_s_t1 * dedF_st( a_3, b_3 ) * dedF_st( c_3, d_3 );


    // Compute d2e_dFdK and d2e_dKdF
    Tensor4_3d d2e_dKdF, d2e_dFdK;
    d2e_dKdF( a_3, b_3, c_3, d_3 ) = -0.5 * J * omega_t1 * ( H - H0 ) * B( a_3, b_3 ) * Finv( d_3, c_3 ) + 
        0.5 * J * omega_t1 * ( H - H0 ) * ( Finv( a_3, c_3 ) * B( b_3, d_3 ) + B( a_3, d_3 ) * Finv( b_3, c_3 ) ) - 
        0.5 * J * omega_t1 * B( a_3, b_3 ) * dHdF( c_3, d_3 ); 

    d2e_dFdK( a_3, b_3, c_3, d_3 ) = d2e_dKdF( c_3, d_3, a_3, b_3 );

    //// Compute d2e_dKdK (This could be loaded)
    //Tensor4_3d d2e_dKdK;
    //d2e_dKdK( a_3, b_3, c_3, d_3 ) = 0.25 * J * omega_t1 * B( a_3, b_3 ) * B( c_3, d_3 ); 

    // compute tensor T 
    Tensor4_3d T;
    T( i_3, j_3, p_3, q_3 ) = ( n_def( p_3 ) * Finv( q_3, k_3 ) ) * G( k_3, i_3, j_3 );

    // compute dTdF
    Tensor6_3d dTdF;
    dTdF( i_3, j_3, k_3, l_3, m_3, n_3 ) = ( -Finv( l_3, m_3 ) * Finv( n_3, a_3 ) + B( l_3, n_3 ) * n_def( a_3 ) * n_def( m_3 ) ) * (G( a_3, i_3, j_3 ) * n_def( k_3 )) -
        (Finv( l_3, a_3 )*Finv( n_3, k_3 ))*(G( a_3, i_3, j_3 )*n_def( m_3 ));

    Afull( k_3, l_3, m_3, n_3 ) = d2e_dFdF( k_3, l_3, m_3, n_3 ) + d2e_dFdF_alpha( k_3, l_3, m_3, n_3 ) + d2e_dFdK( k_3, l_3, o_3, p_3 ) * T( o_3, p_3, m_3, n_3 ) + 
        ( d2e_dKdF( i_3, j_3, m_3, n_3 ) + d2e_dKdK( i_3, j_3, o_3, p_3 ) * T( o_3, p_3, m_3, n_3 ) ) * T( i_3, j_3, k_3, l_3 )  + 
        Mk( i_3, j_3 ) * dTdF( i_3, j_3, k_3, l_3, m_3, n_3 );


    return Afull.to_voigt_form();
}


FloatMatrixF<27, 27>
IsotropicPolyconvexHyperelasticSurfaceMaterial2::give3dSurfaceMaterialStiffnessMatrix_dAddF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns the 27x27 second order tangent stiffness matrix - dA/ddF (derivative of second order stress wrt to gradient of F)
{
    // load objects from status
    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    FloatArrayF<27> vG( status->giveTempGVector() ); // second gradient

    // create identities
    FloatArrayF<9> indent   = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d IdentityUD( indent ), Identity3d( indent3d );

    // if it is IGA element
    FloatArrayF<3> UDn( normal );
    Tensor1_3d Nt( UDn );
    if ( dynamic_cast<IGAElement *>( gp->giveElement() ) ) {
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
    }

    // Get material parameter
    double omega_t1;
    //if ( this->gamma_ltf == 0 ) {
    //    omega_t1 = this->omega1;
    //} else {
    //    omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //}

    //if ( this->fixOmega ) {
    //    omega_t1 = this->omega1; // omega_t1 not increased gradually
    //}

    omega_t1 = this->giveOmega( tStep );

    Tensor2_3d F( vF );
    Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
    double J         = this->compute_surface_determinant( F );
    Tensor2_3d B;
    B( j_3, m_3 ) = Finv( j_3, k_3 ) * Finv( m_3, k_3 );


    // stiffnes conjugate to MK
    Tensor4_3d d2e_dKdK;
    d2e_dKdK( a_3, b_3, c_3, d_3 ) = 0.25 * J * omega_t1 * B( a_3, b_3 ) * B( c_3, d_3 ); 

    // compute deformed normal
    Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal

    // Compute second order stiffness
    Tensor6_3d dA;
    // created in a stupid way by contraction
    //dA( i_3, j_3, k_3, l_3, m_3, n_3 ) = ( ( n_def( i_3 ) * n_def( l_3 ) ) * Identity3d( j_3, s_3 ) ) * d2e_dKdK( s_3, k_3, m_3, n_3 );
    dA( i_3, j_3, k_3, l_3, m_3, n_3 ) = ( n_def( i_3 ) * n_def( l_3 ) ) * d2e_dKdK( j_3, k_3, m_3, n_3 );
    /////////////////////

    status->lettempD2eDK2_MatrixBe( d2e_dKdK.to_voigt_form() );
    return dA.to_voigt_form();
}


FloatMatrixF<27, 9>
IsotropicPolyconvexHyperelasticSurfaceMaterial2::give3dSurfaceMaterialStiffnessMatrix_dAdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns the 27x9 mixed tangent stiffness matrix - dA/dF (derivative of second order stress wrt to F)
{
    // load objects from status
    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    FloatArrayF<9> vM( status->giveTempMVector() ); // second order stress conjugate to curvature tensor from Silhous
    FloatArrayF<27> vG( status->giveTempGVector() ); // second gradient
    FloatMatrixF<9, 9> vdAK( status->giveTempD2eDK2_Matrix() ); // second gradient
    FloatArrayF<27> dI( status->giveTempdIVector() );

    // create identities
    FloatArrayF<9> indent   = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d IdentityUD( indent ), Identity3d( indent3d );

    // if it is IGA element
    FloatArrayF<3> UDn( normal );
    Tensor1_3d Nt( UDn );
    if ( dynamic_cast<IGAElement *>( gp->giveElement() ) ) {
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
    }

    // Get material parameter
    double gamma_t1, omega_t1;
    //
    //if ( this->gamma_ltf == 0 ) {
    //    gamma_t1 = this->gamma1;
    //    omega_t1 = this->omega1;
    //} else {
    //    gamma_t1 = this->gaminit + this->gamma1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //    omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //}
    //if ( this->fixGamma ) {
    //    gamma_t1 = this->gamma1; // omega_t1 not increased gradually
    //}
    //if ( this->fixOmega ) {
    //    omega_t1 = this->omega1; // omega_t1 not increased gradually
    //}

    gamma_t1 = this->giveGamma( tStep );
    omega_t1 = this->giveOmega( tStep );


    Tensor3_3d G( vG ), gradIt( dI ); // Gradient of F
    Tensor2_3d Mk( vM ); // Second order stress conjugate to bending strain
    Tensor2_3d F( vF ); // get second order stress tensor and F
    Tensor4_3d dAK( vdAK );

    double J         = this->compute_surface_determinant( F );
    Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
    Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
    Tensor2_3d B;
    B( j_3, m_3 ) = Finv( j_3, k_3 ) * Finv( m_3, k_3 );

    Tensor2_3d K, K0, b, b0; // gradient of deformed normal nad bending strain
    K( j_3, k_3 )  = -n_def( i_3 ) * G( i_3, j_3, k_3 ); // Deformed curvature in referential frame
    b( j_3, l_3 )  = -Finv( k_3, j_3 ) * K( k_3, m_3 ) * Finv( m_3, l_3 );
    double H       = 0.5 * b( j_3, j_3 );

    double H0;
    if ( this->use_H0 == 1 ) {
        K0( j_3, k_3 ) = -Nt( i_3 ) * gradIt( i_3, j_3, k_3 ); // gradient of undeformed normal
        b0( j_3, l_3 ) = -K0( j_3, l_3 );
        H0             = 0.5 * b0( j_3, j_3 );
    } else if ( this->use_H0 == 2 ) {
        H0 = this->H0_val;
    } else {
        H0 = 0.;
    }
    

    double e = J * ( gamma_t1 + 0.5 * omega_t1 * ( H - H0 ) * ( H - H0 ) ); 
    
    Tensor2_3d dHdF;
    dHdF( i_3, l_3 ) = Finv( j_3, i_3 ) * K( j_3, k_3 ) * B( k_3, l_3 );

    // compute tensor T
    Tensor4_3d T;
    T( i_3, j_3, p_3, q_3 ) = ( n_def( p_3 ) * Finv( q_3, k_3 ) ) * G( k_3, i_3, j_3 );

    // Compute d2e_dKdF
    Tensor4_3d d2e_dKdF;
    d2e_dKdF( a_3, b_3, c_3, d_3 ) = -0.5 * J * omega_t1 * ( H - H0 ) * B( a_3, b_3 ) * Finv( d_3, c_3 ) + 
        0.5 * J * omega_t1 * ( H - H0 ) * ( Finv( a_3, c_3 ) * B( b_3, d_3 ) + B( a_3, d_3 ) * Finv( b_3, c_3 ) ) -
        0.5 * J * omega_t1 * B( a_3, b_3 ) * dHdF( c_3, d_3 ); 


    Tensor5_3d dA;
    dA( i_3, j_3, k_3, a_3, b_3 ) = ( Mk( j_3, k_3 ) * n_def( a_3 ) ) * Finv( b_3, i_3 ) - 
        ( d2e_dKdF( j_3, k_3, a_3, b_3 ) + dAK( j_3, k_3, o_3, p_3 ) * T( o_3, p_3, a_3, b_3 ) ) * n_def( i_3 );

    return dA.to_voigt_form_27x9();
}



void IsotropicPolyconvexHyperelasticSurfaceMaterial2::initializeFrom( InputRecord &ir )
{
    IsotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( ir );

    IR_GIVE_OPTIONAL_FIELD( ir, use_H0, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_use_H0 );
    IR_GIVE_OPTIONAL_FIELD( ir, H0_val, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_H0 );


    // arclength parameter
    int temp = 0;;
    IR_GIVE_OPTIONAL_FIELD( ir, temp, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_ArcLengthParameter );
    this->arcLengthParamater = (CALM_parameter) temp;
    if ( this->arcLengthParamater != calm_par_no ) {
        this->hasALP = true;
    }
}


double IsotropicPolyconvexHyperelasticSurfaceMaterial2::giveGamma( TimeStep *tStep ) const
{
    if ( this->hasArcLengthParameter() ) { // one of the parameters is the arc length loading
        if ( arcLengthParamater == calm_par_gamma ) { // gamma is the CALM load parameter
            if ( this->returnAll ) { // Classical internal forces are assembled. 
                return this->lambdaCALM * this->gamma1; 
            } else { // only the term multiplying lambda in CALM is assembled
                return this->gamma1; // as if lambdaCALM = 0
            }
        } else { // gamma is normal parameter
            if ( this->returnAll ) { // Classical internal forces are assembled.
                return this->gamma1; // gamma is fixed classical parameter
            } else { // only the term multiplying lambda in CALM is assembled
                return 0.; // as if lambdaCALM = 0
            } 
        }
        

    } else { // no arclength, classical formulation
        if ( this->gamma_ltf == 0 || this->fixGamma ) {
            return this->gamma1;
        } else {
            return this->gaminit + this->gamma1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        }

        //if ( this->fixGamma ) {
        //    return this->gamma1; // omega_t1 not increased gradually
        //}
    
    }
    
};

double IsotropicPolyconvexHyperelasticSurfaceMaterial2::giveOmega( TimeStep *tStep ) const
{
    if ( this->hasArcLengthParameter() ) { // one of the parameters is the arc length loading
        if ( arcLengthParamater == calm_par_omega ) { // omega is the CALM load parameter
            if ( this->returnAll ) { // Classical internal forces are assembled.
                return this->lambdaCALM * this->omega1;
            } else { // only the term multiplying lambda in CALM is assembled
                return this->omega1; // as if lambdaCALM = 0
            }
        } else { // gamma is normal parameter
            if ( this->returnAll ) { // Classical internal forces are assembled.
                return this->omega1; // gamma is fixed classical parameter
            } else { // only the term multiplying lambda in CALM is assembled
                return 0.; // as if lambdaCALM = 0
            }
        }


    } else { // no arclength, classical formulation
        if ( this->gamma_ltf == 0 || this->fixOmega ) {
            return this->omega1;
        } else {
            return this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        }

        //if ( this->fixOmega ) {
        //    return this->omega1; // omega_t1 not increased gradually
        //}
    }
};


double IsotropicPolyconvexHyperelasticSurfaceMaterial2::giveAlpha_s( TimeStep *tStep ) const
{
    if ( this->hasArcLengthParameter() ) { // one of the parameters is the arc length loading
        OOFEM_ERROR( "This model is not implemented for arc-length method" )
    } else { // no arclength, classical formulation
        if ( this->gamma_ltf == 0 || this->fixAlpha_s ) {
            return this->alpha_s1;
        } else {
            return this->alpha_s1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        }
    }
};

} // end namespace oofem
