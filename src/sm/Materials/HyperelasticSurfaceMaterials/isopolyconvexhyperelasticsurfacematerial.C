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
    Tensor2_3d F, Fcurrent( vF ), P, Pfull, S, Identity3d( indent3d ), Identity( indent ), Pgamma1, Pgamma2, Palpha1, Palpha2, Pdelta;
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

    double gamma_t1, alpha_t1, omega_t1;
    if(this->gamma_ltf == 0) {
        gamma_t1 = this->gamma1;
        alpha_t1 = this->alpha1;
        omega_t1 = this->omega1;
    } else {
        gamma_t1 = this->gaminit + this->gamma1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        alpha_t1 = this->alpha1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    //alpha_t1 = this->alpha1; // Alpha not increased gradually
    //gamma_t1 = this->gamma1; // Gamma not increased gradually

    // compute the first Piola-Kirchhoff
    Pgamma1( i_3, j_3 ) = gamma_t1 * this->compute_surface_cofactor( F, IdentityUD )( i_3, j_3 ); 
    Palpha1( i_3, j_3 ) =  alpha_t1 * this->compute_surface_d_normF_dF( F )( i_3, j_3 ) ;

    P( i_3, j_3 ) = Pgamma1( i_3, j_3 ) + Palpha1( i_3, j_3 ) ;
    
    //double J0     = this->compute_surface_determinant( F0 );
    //Pfull( i_3, j_3 ) = 1. / J0 * P( i_3, k_3 ) * F0( j_3, k_3 ); 

    ///////////////////////////////////////
    // Get the second order part
    if ( omega_t1 > 1e-6 ) {
        // load objhect from status
        FloatArrayF<27> vG( status->giveTempGVector() );
        FloatArrayF<27> dI( status->giveTempdIVector() );
        FloatArrayF<9> vM( status->giveTempMVector() );
        Tensor3_3d G( vG ), gradIt( dI ); // grad. F and grad. of surface identity

        Tensor2_3d MK( vM ); // second order stress conjugate to kappa
        Tensor2_3d P2 = this->giveFirstPKSurfaceStressBending( F, IdentityUD, Nt, G, MK );
        ////////////////////////////////////////

        Pfull( i_3, j_3 ) = P( i_3, j_3 ) + P2( i_3, j_3 );
    } else {
        Pfull( i_3, j_3 ) = P( i_3, j_3 );
    }
    
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
    Tensor2_3d Fcurrent( vF ), F;
    // Tensor2_3d F0( indent );
    // if ( gp->giveElement()->giveIPValue( vF0_temp, gp, IST_PrestrainDeformationGradient, tStep ) ) {
    //     FloatArrayF<9> vF0( vF0_temp );
    //     F0  = Tensor2_3d::Tensor2_3d( vF0 );
    // }
    // F( i_3, j_3 ) = Fcurrent( i_3, k_3 ) * F0( k_3, j_3 );
    F( i_3, j_3 ) = Fcurrent( i_3, j_3 );

    // Tensor2_3d F( vF );
    Tensor4_3d A, Afull;
    double gamma_t1, alpha_t1, omega_t1;
    //
    if ( this->gamma_ltf == 0 ) {
        gamma_t1 = this->gamma1;
        alpha_t1 = this->alpha1;
        omega_t1 = this->omega1;
    } else {
        gamma_t1 = this->gaminit + this->gamma1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        alpha_t1 = this->alpha1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    // alpha_t1 = this->alpha1; // Alpha not increased gradually
    // gamma_t1 = this->gamma1; // Gamma not increased gradually

    A( i_3, j_3, k_3, l_3 ) = gamma_t1 * this->compute_surface_dCof_dF( F, IdentityUD )( i_3, j_3, k_3, l_3 ) + alpha_t1 * this->compute_surface_d2_normF_dF2( F, IdentityUD )( i_3, j_3, k_3, l_3 );

    ///////////////////////////////////////
    // Get the second order part 
    if ( omega_t1 > 1e-6 ) {
        // construct some useful objects
        Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
        Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
        FloatArrayF<27> vG( status->giveTempGVector() );
        Tensor3_3d G( vG ); // Gradient of F
        FloatArrayF<9> vM( status->giveTempMVector() );
        Tensor2_3d Mk( vM ); // Second order stress conjugate to bending strain

        // Tensor2_3d idef;
        // idef( j_3, k_3 ) = F( j_3, m_3 ) * Finv( m_3, k_3 );


        // Compute stiffnesses
        FloatMatrixF<9, 9> vdAK( status->giveTempD2eDK2_Matrix() ); // second gradient
        Tensor4_3d dAK( vdAK ), A2, A3, A4, temp, dFinvdF;

        temp( p_3, q_3, i_3, j_3 ) = ( n_def( p_3 ) * Finv( q_3, k_3 ) ) * G( k_3, i_3, j_3 );
        A2( p_3, q_3, a_3, b_3 )   = temp( p_3, q_3, i_3, j_3 ) * dAK( i_3, j_3, r_3, s_3 ) * temp( a_3, b_3, r_3, s_3 ); // classical part

        Tensor3_3d dndF;
        dndF( p_3, a_3, b_3 ) = -Finv( b_3, p_3 ) * n_def( a_3 ); // derivative of deformed normal

        double Js = this->compute_surface_determinant( F );
        // derivative of pseudoinverse
        dFinvdF( i_3, j_3, k_3, l_3 ) = 1. / Js * this->compute_surface_dCof_dF( F, IdentityUD )( j_3, i_3, k_3, l_3 ) - 1. / Js / Js * this->compute_surface_cofactor( F, IdentityUD )( j_3, i_3 ) * this->compute_surface_cofactor( F, IdentityUD )( k_3, l_3 );

        A3( p_3, q_3, a_3, b_3 ) = dFinvdF( q_3, k_3, a_3, b_3 ) * ( G( k_3, i_3, j_3 ) * Mk( i_3, j_3 ) ) * n_def( p_3 ); // part from second variation of kappa
        A4( p_3, q_3, a_3, b_3 ) = dndF( p_3, a_3, b_3 ) * ( Finv( q_3, k_3 ) * ( G( k_3, i_3, j_3 ) * Mk( i_3, j_3 ) ) ); // part from second variation of kappa

        // Tensor4_3d A34;
        // A34( i_3, j_3, k_3, l_3 ) = A3( i_3, j_3, k_3, l_3 ) + A4( i_3, j_3, k_3, l_3 );
        // auto A34m = FloatMatrix( A34.to_voigt_form() );
        ////  A34m.printYourself();

        ////////////////////////////////////////

        ////// Transform to to original reference frame, just for testing
        //// Tensor4_3d A_rotated;
        //// Tensor2_3d F_rotated;
        //// this->TransformTensor2( F, F_rotated, gp );
        //// this->TransformTensor4( A, A_rotated, gp );

        //// double J0 = this->compute_surface_determinant( F0 );
        //// Afull( i_3, j_3, k_3, l_3 ) = 1./J0* F0( j_3, m_3 ) * A( i_3, m_3, k_3, n_3 ) * F0( l_3, n_3 );

        Afull( i_3, j_3, k_3, l_3 ) = A( i_3, j_3, k_3, l_3 ) + A2( i_3, j_3, k_3, l_3 ) + A3( i_3, j_3, k_3, l_3 ) + A4( i_3, j_3, k_3, l_3 );
    } else {
        Afull( i_3, j_3, k_3, l_3 ) = A( i_3, j_3, k_3, l_3 );
    }


    return Afull.to_voigt_form();
}

/////////
// Second order stress and stiffnesses

FloatArrayF<27>
IsotropicPolyconvexHyperelasticSurfaceMaterial::giveSecondOrderSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns 27 components of the second order stress conjugate to gradient of F
{    
    // create identities
    FloatArrayF<9> indent            = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    FloatArrayF<9> indent3d          = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Identity3d( indent3d );
    Tensor2_3d IdentityUD( indent );

    FloatArrayF<3> UDn( normal );
    Tensor1_3d Nt( UDn ); // normal tensor
    // if it is IGA element
    if ( dynamic_cast<IGAElement *>( gp->giveElement() ) ) {
        IdentityUD( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
    }

    // Get material parameter
    double omega_t1;
    if ( this->gamma_ltf == 0 ) {
        omega_t1 = this->omega1;
    } else {
        omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }

    // construct some useful objects
    Tensor2_3d F( vF ); // deformation gradient
    Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
    Tensor1_3d n_def   = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
    Tensor3_3d A, G( vG ), gradIt( gradI );

    // Compute curvature measures
    // formulation based on Dort or Tomassetti using bending strain
    Tensor2_3d K, K0; // gradient of deformed normal nad bending strain
    K( j_3, k_3 ) =  - n_def( i_3 ) * G( i_3, j_3, k_3 ) ; // Deformed curvature in referential frame
    K0( j_3, k_3 ) = - Nt( i_3 ) * gradIt( i_3, j_3, k_3 ); // gradient of undeformed normal

    // compute the second order tensor conjugate to K
    Tensor2_3d MK;
    //MK( i_3, j_3 ) = omega_t1 * ( K( l_3, l_3 ) - K0( k_3, k_3 ) ) * IdentityUD( i_3, j_3 );
    MK( i_3, j_3 ) = omega_t1 * ( K( l_3, l_3 ) ) * IdentityUD( i_3, j_3 );

    // obtain the third order tensor from Mk conjugarte to grad. F
    A( i_3, j_3, k_3 ) = -n_def( i_3 ) * MK( j_3, k_3 ); 
    auto vA = A.to_voigt_form_27();

    // Save object to status, note that this function needs to be called before the standar surface PK stress is evaluated
    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    status->letTempFVectorBe( vF );
    status->letTempGVectorBe( vG );
    status->letTempdIVectorBe( gradI );
    status->letTempMVectorBe( MK.to_voigt_form() );

    return vA;
}


FloatMatrixF<27, 27>
IsotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dAddF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
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
    if ( this->gamma_ltf == 0 ) {
        omega_t1 = this->omega1;
    } else {
        omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }

    // stiffnes conjugate to MK
    Tensor4_3d dAK;
    dAK( i_3, j_3, k_3, l_3 ) = omega_t1 * IdentityUD( i_3, j_3 ) * IdentityUD( k_3, l_3 );
    //dAK( i_3, j_3, k_3, l_3 ) = omega_t1 * Identity3d( i_3, j_3 ) * Identity3d( k_3, l_3 );

    // compute deformed normal
    Tensor2_3d F( vF );
    Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
    Tensor3_3d G( vG ); // Get gradient of F

    // Compute second order stiffness
    Tensor6_3d dA;
    // created in a stupid way by contraction
    dA( i_3, j_3, k_3, l_3, m_3, n_3 ) = ( (n_def( i_3 ) * n_def( l_3 )) * Identity3d( j_3, s_3 ) ) * dAK( s_3, k_3, m_3, n_3 ); 
    /////////////////////

    status->lettempD2eDK2_MatrixBe( dAK.to_voigt_form() );
    return dA.to_voigt_form();
}



FloatMatrixF<27, 9>
IsotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dAdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns the 27x9 mixed tangent stiffness matrix - dA/dF (derivative of second order stress wrt to F)
{
    // load objects from status
    StructuralSurfaceMaterialStatus *status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    FloatArrayF<9> vM( status->giveTempMVector() ); // second order stress conjugate to curvature tensor from Silhous
    FloatArrayF<27> vG( status->giveTempGVector() ); // second gradient
    FloatMatrixF<9,9> vdAK( status->giveTempD2eDK2_Matrix() ); // second gradient

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
    if ( this->gamma_ltf == 0 ) {
        omega_t1 = this->omega1;
    } else {
        omega_t1 = this->omega1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }

    // stiffnes conjugate to MK
    Tensor4_3d dAK( vdAK );
    //dAK( i_3, j_3, k_3, l_3 ) = omega_t1 * IdentityUD( i_3, j_3 ) * IdentityUD( k_3, l_3 );

    // construct some useful objects
    Tensor2_3d F( vF ), Mk( vM ); // get second order stress tensor and F
    Tensor2_3d Finv  = this->compute_surface_pseudoinverse( F, IdentityUD ); // compute pseudoinverse
    Tensor1_3d n_def = compute_surface_deformed_normal( F, Nt ); // compute deformed normal
    Tensor3_3d G( vG ); // Get gradient of F

    // Compute the mixed stiffness
    // creted in a stupid way with the extra identity!
    Tensor5_3d dA, dA1, dA2, dA3, dA4;
    dA1( k_3, i_3, j_3, p_3, q_3 ) = 1. * ( Mk( i_3, j_3 ) * n_def( r_3 ) ) * ( Identity3d( r_3, p_3 ) * Finv( q_3, k_3 ) ); 
    dA2( k_3, i_3, j_3, a_3, b_3 ) = -1. * ( dAK( i_3, j_3, r_3, s_3 ) * G( c_3, r_3, s_3 ) ) * ( n_def( k_3 ) * n_def( a_3 ) * Finv( b_3, c_3 ) );

    dA( p_3, q_3, l_3, i_3, j_3 ) = dA1( p_3, q_3, l_3, i_3, j_3 ) + dA2( p_3, q_3, l_3, i_3, j_3 );
    return dA.to_voigt_form_27x9();
}



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
    IR_GIVE_OPTIONAL_FIELD( ir, delta, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_delta );
    IR_GIVE_OPTIONAL_FIELD( ir, gamma_ltf, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gammaLTF );
    IR_GIVE_OPTIONAL_FIELD( ir, gaminit, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gaminit );
    IR_GIVE_OPTIONAL_FIELD( ir, omega1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_omega1 );

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

} // end namespace oofem
