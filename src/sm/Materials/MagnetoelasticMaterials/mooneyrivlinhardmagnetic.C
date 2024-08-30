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

#include "mooneyrivlinhardmagnetic.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "domain.h"
#include "function.h"
#include "tensor/tensor3.h"



namespace oofem {
REGISTER_Material(MooneyRivlinHardMagnetic);

FloatArrayF< 9 >
MooneyRivlinHardMagnetic::giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

    //get stress from Mooney-Rivlin hyperelasticity
    FloatArrayF<9> vP_mr = MooneyRivlinCompressibleMaterial::giveFirstPKStressVector_3d( vF, gp, tStep );

    // calculate the magnetoelastic 1PK
    FloatArrayF<9> vP_me;
    switch ( materialMode ) {
    case 1:
        vP_me = computeFirstPKStressVector_3d_consistent( vF, gp, tStep );
        break;
    case 2:
        vP_me = computeFirstPKStressVector_3d_mit( vF, gp, tStep );
        break;
    case 3:
        vP_me = computeFirstPKStressVector_3d_consistentfix( vF, gp, tStep );
        break;
    case 4:
        vP_me = computeFirstPKStressVector_3d_ogdenpullback( vF, gp, tStep );
        break;
    case 5:
        vP_me = computeFirstPKStressVector_3d_multiphysics( vF, gp, tStep );
        break;
    default:
        OOFEM_ERROR( "Unknown hard magnetic material mode %i", materialMode );
        break;
    }

    
    //add magnetoelastic and mooney-rivlin stress together
    auto vP    = vP_mr + vP_me;

    // update gp
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );
    
    return vP;
}

FloatMatrixF< 9, 9 >
MooneyRivlinHardMagnetic::give3dMaterialStiffnessMatrix_dPdF( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    // get stiffness from Mooney-Rivlin hyperelasticity
    FloatMatrixF<9, 9> vD_mr = MooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dPdF( matResponseMode, gp, tStep );

    // calculate magnetoelastic stiffness
    FloatMatrixF<9,9> vD_me;
    switch ( materialMode ) {
    case 1:
        vD_me = compute3dMaterialStiffnessMatrix_dPdF_consistent( matResponseMode, gp, tStep );
        break;
    case 2:
        vD_me = compute3dMaterialStiffnessMatrix_dPdF_mit( matResponseMode, gp, tStep );
        break;
    case 3:
        vD_me = compute3dMaterialStiffnessMatrix_dPdF_consistentfix( matResponseMode, gp, tStep );
        break;
    case 4:
        vD_me = compute3dMaterialStiffnessMatrix_dPdF_ogdenpullback( matResponseMode, gp, tStep );
        break;
    case 5:
        vD_me = compute3dMaterialStiffnessMatrix_dPdF_multiphysics( matResponseMode, gp, tStep );
        break;
    default:
        OOFEM_ERROR( "Unknown hard magnetic material mode %i", materialMode );
        break;
    }

    //add magnetoelastic and Mooney-Rivlin stiffness together
    auto vD = vD_me + vD_mr;

    return vD;
}

FloatArrayF<3> MooneyRivlinHardMagnetic::giveMagneticInduction_3d( const FloatArrayF<3> &vH, GaussPoint *gp, TimeStep *tStep ) const
{
    if ( materialMode < 5 ) {
        OOFEM_ERROR( "Material not in multiphysics mode." );
    }

    return computeMagneticInduction_3d_multiphysics( vH, gp, tStep );
}

FloatMatrixF<3, 9> MooneyRivlinHardMagnetic::give3dMagnetoelasticCouplingTensor_dBdF( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    if ( materialMode < 5 ) {
        OOFEM_ERROR( "Material not in multiphysics mode." );
    }

    return compute3dMagnetoelasticCouplingTensor_dBdF_multiphysics( matResponseMode, gp, tStep );
}

FloatMatrixF<3, 3> MooneyRivlinHardMagnetic::give3dPermeabilityMatrix_dBdH( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    if ( materialMode < 5 ) {
        OOFEM_ERROR( "Material not in multiphysics mode." );
    }

    return compute3dPermeabilityMatrix_dBdH_multiphysics( matResponseMode, gp, tStep );
}

void MooneyRivlinHardMagnetic::initializeFrom( InputRecord &ir )
{   
    MooneyRivlinCompressibleMaterial::initializeFrom( ir );

    FloatArray B_app_temp, B_res_temp, H_temp;

    //default values
    mu_0 = BASE_VACUUM_PERMEABILITY_MU_0;
    materialMode = 1;

    IR_GIVE_FIELD( ir, B_app_temp, _IFT_MooneyRivlinHardMagnetic_B_app );
    IR_GIVE_FIELD( ir, B_res_temp, _IFT_MooneyRivlinHardMagnetic_B_res );
    IR_GIVE_FIELD( ir, H_temp, _IFT_MooneyRivlinHardMagnetic_H );
    IR_GIVE_OPTIONAL_FIELD( ir, mu_0, _IFT_MooneyRivlinHardMagnetic_mu_0 );
    IR_GIVE_FIELD( ir, ltf_index, _IFT_MooneyRivlinHardMagnetic_ltf );
    IR_GIVE_OPTIONAL_FIELD( ir, materialMode, _IFT_MooneyRivlinHardMagnetic_mode );

    referenceB = false;
    IR_GIVE_OPTIONAL_FIELD( ir, referenceB, _IFT_MooneyRivlinHardMagnetic_referenceB );

    B_app = FloatArrayF<3>( B_app_temp );
    B_res = FloatArrayF<3>( B_res_temp );
    H = FloatArrayF<3>( H_temp );
}

FloatArrayF<9> MooneyRivlinHardMagnetic::computeFirstPKStressVector_3d_consistent( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
{

    double load_level            = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    FloatArrayF<3> B_app_at_time = load_level * B_app;

    Tensor2_3d F( vF ), P_me, delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );

    Tensor1_3d Bapp( B_app_at_time ), Bres( B_res ), Bappref;
    auto [J, cofF] = F.compute_determinant_and_cofactor();

    if ( !referenceB ) {
        Bappref( j_3 ) = Bapp( k_3 ) * cofF( k_3, j_3 );
    } else {
        Bappref = Bapp;
    }

    P_me( k_3, l_3 ) = -(1/(2*J*J*mu_0))*(Bappref(i_3) - 2.0*Bres(i_3))*F(m_3,i_3)*F(m_3,j_3)*Bappref(j_3)*cofF(k_3, l_3)
        + (1/(2*J*mu_0)) * ( Bappref( i_3 ) - 2.0 * Bres( i_3 ) ) * Bappref(j_3) * (delta(i_3, l_3)*F(k_3,j_3) + delta(j_3, l_3)*F(k_3,i_3));

    auto vP_me = P_me.to_voigt_form();

    return vP_me;
}

FloatMatrixF<9, 9> MooneyRivlinHardMagnetic::compute3dMaterialStiffnessMatrix_dPdF_consistent( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );

    
    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    FloatArrayF<3> B_app_at_time = load_level * B_app;

    Tensor2_3d F( vF );
    Tensor4_3d D_me;

    Tensor2_3d delta(1., 0., 0., 0., 1., 0., 0., 0., 1. );

    Tensor1_3d Bapp ( B_app_at_time ), Bres( B_res ), Bappref;
    auto [J, cofF] = F.compute_determinant_and_cofactor();

    if ( !referenceB ) {
        Bappref( j_3 ) = Bapp( k_3 ) * cofF( k_3, j_3 );
    } else {
        Bappref = Bapp;
    }

    D_me( k_3, l_3, p_3, q_3 ) = (1/(J*J*J*mu_0)) * ( Bappref(i_3) - 2.0 * Bres(i_3))*F(m_3,i_3)*F(m_3,j_3)*Bappref(j_3)*cofF(k_3, l_3)*cofF(p_3, q_3)
        - (1/(2*J*J*mu_0)) * ( Bappref(i_3) - 2.0 * Bres(i_3)) * Bappref(j_3) * (delta(i_3, q_3)*F(p_3,j_3) + delta(j_3, q_3)*F(p_3,i_3)) * cofF(k_3, l_3) 
        - (1/(2*J*J*mu_0)) * ( Bappref(i_3) - 2.0 * Bres(i_3)) * F(m_3,i_3)*F(m_3,j_3) * Bappref(j_3) * F.compute_tensor_cross_product()(k_3, l_3, p_3, q_3)
        - (1/(2*J*J*mu_0)) * ( Bappref(i_3) - 2.0 * Bres(i_3)) * Bappref(j_3) * (delta(i_3, l_3)*F(k_3,j_3) + delta(j_3, l_3)*F(k_3,i_3)) * cofF(p_3, q_3)
        + (1/(2*J*mu_0)) * delta(k_3,p_3) * (( Bappref(l_3) - 2.0 * Bres(l_3))*Bappref(q_3) + ( Bappref(q_3) - 2.0 * Bres(q_3))*Bappref(l_3));
    
    auto vD_me = D_me.to_voigt_form();

    return vD_me;
}

FloatArrayF<9> MooneyRivlinHardMagnetic::computeFirstPKStressVector_3d_mit( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
{

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    FloatArrayF<3> B_app_at_time = load_level * B_app;

    Tensor2_3d F( vF ), P_me;
    Tensor1_3d Bapp( B_app_at_time ), Bres( B_res );

    P_me( k_3, l_3 ) = -1.0 / mu_0 * Bapp( k_3 ) * Bres( l_3 );

    auto vP_me = P_me.to_voigt_form();

    return vP_me;
}

FloatMatrixF<9, 9> MooneyRivlinHardMagnetic::compute3dMaterialStiffnessMatrix_dPdF_mit( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    return FloatMatrixF<9, 9>();
}

FloatArrayF<9> MooneyRivlinHardMagnetic::computeFirstPKStressVector_3d_consistentfix( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
{

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    FloatArrayF<3> B_app_at_time = load_level * B_app;

    Tensor2_3d F( vF ), P_me;
    
    Tensor1_3d Bapp( B_app_at_time ), Bres( B_res );
    auto [J, cofF] = F.compute_determinant_and_cofactor();


    P_me( k_3, l_3 ) = (1/(2*mu_0))*Bapp( j_3 ) * Bapp(j_3) * cofF(k_3, l_3) - (1/mu_0)*Bapp(k_3)*Bres(l_3);

    auto vP_me = P_me.to_voigt_form();

    return vP_me;
}

FloatMatrixF<9, 9> MooneyRivlinHardMagnetic::compute3dMaterialStiffnessMatrix_dPdF_consistentfix( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );

    
    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    FloatArrayF<3> B_app_at_time = load_level * B_app;

    Tensor2_3d F( vF );
    Tensor4_3d D_me;

    Tensor2_3d delta(1., 0., 0., 0., 1., 0., 0., 0., 1. );

    Tensor1_3d Bapp ( B_app_at_time ), Bres( B_res );
    auto [J, cofF] = F.compute_determinant_and_cofactor();

    D_me( k_3, l_3, p_3, q_3 ) = (1/(2*mu_0))*Bapp( j_3 ) * Bapp(j_3) * F.compute_tensor_cross_product() (k_3, l_3, p_3, q_3);
    
    auto vD_me = D_me.to_voigt_form();

    return vD_me;
}

FloatArrayF<9> MooneyRivlinHardMagnetic::computeFirstPKStressVector_3d_ogdenpullback( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
{

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    FloatArrayF<3> B_app_at_time = load_level * B_app;

    Tensor2_3d F( vF ), P_me;
    
    Tensor1_3d Bapp( B_app_at_time ), Bres( B_res );
    auto [J, cofF] = F.compute_determinant_and_cofactor();
    auto invF      = F.compute_inverse();

    P_me( k_3, l_3 ) = (1/(2*mu_0))*Bapp( j_3 ) * Bapp(j_3) * cofF(k_3, l_3)
                     + (J/mu_0)*Bres(i_3)*Bapp(j_3)*(invF(i_3,k_3)*invF(l_3, j_3) - invF(l_3,k_3)*invF(i_3, j_3));

    auto vP_me = P_me.to_voigt_form();

    return vP_me;
}

FloatMatrixF<9, 9> MooneyRivlinHardMagnetic::compute3dMaterialStiffnessMatrix_dPdF_ogdenpullback( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );

    
    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    FloatArrayF<3> B_app_at_time = load_level * B_app;

    Tensor2_3d F( vF );
    Tensor4_3d D_me;

    Tensor2_3d delta(1., 0., 0., 0., 1., 0., 0., 0., 1. );

    Tensor1_3d Bapp ( B_app_at_time ), Bres( B_res );
    auto [J, cofF] = F.compute_determinant_and_cofactor();
    auto invF      = F.compute_inverse();

    D_me( k_3, l_3, p_3, q_3 ) = (1/(2*mu_0))*Bapp( j_3 ) * Bapp(j_3) * F.compute_tensor_cross_product() (k_3, l_3, p_3, q_3)
                               + (1/mu_0)*Bres(i_3)*Bapp(j_3)*(invF(i_3,k_3)*invF(l_3, j_3) - invF(l_3,k_3)*invF(i_3, j_3))*cofF(p_3, q_3)
                               - (J/mu_0)*Bres(i_3)*Bapp(j_3)*invF(i_3, p_3)*invF(q_3, k_3)*invF(l_3, j_3)
                               - (J/mu_0)*Bres(i_3)*Bapp(j_3)*invF(i_3, k_3)*invF(l_3, p_3)*invF(q_3, j_3)
                               + (J/mu_0)*Bres(i_3)*Bapp(j_3)*invF(l_3, p_3)*invF(i_3, j_3)*invF(q_3, k_3)
                               + (J/mu_0)*Bres(i_3)*Bapp(j_3)*invF(l_3, k_3)*invF(i_3, p_3)*invF(q_3, j_3);
    
    auto vD_me = D_me.to_voigt_form();

    return vD_me;
}

//--------------------------------------------------------------------------- Multiphysics ----------------------------------------------------------

FloatArrayF<9> MooneyRivlinHardMagnetic::computeFirstPKStressVector_3d_multiphysics( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
{
    MagnetoelasticMaterialStatus *status = static_cast<MagnetoelasticMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<3> vH( status->giveTempHVector() );

    FloatArrayF<3> vQ = vH + B_res*(1/mu_0); //Q = (H+M)

    Tensor1_3d Q( vQ );
    Tensor2_3d F( vF ), P_me;
    auto [J, cofF] = F.compute_determinant_and_cofactor();
    auto Fcross    = F.compute_tensor_cross_product();

    P_me( i_3, j_3 ) = -( mu_0 / J ) * Fcross( i_3, j_3, m_3, n_3 ) * cofF( m_3, q_3 ) * Q( q_3 ) * Q( n_3 )
                       + (mu_0 / (2*J*J)) * cofF(m_3, n_3) * Q(n_3) * cofF(m_3, k_3) * Q(k_3) * cofF(i_3, j_3);

    return P_me.to_voigt_form();
}

FloatMatrixF<9, 9> MooneyRivlinHardMagnetic::compute3dMaterialStiffnessMatrix_dPdF_multiphysics( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    MagnetoelasticMaterialStatus *status = static_cast<MagnetoelasticMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<3> vH( status->giveTempHVector() );
    FloatArrayF<9> vF( status->giveTempFVector() );

    FloatArrayF<3> vQ = vH + B_res*(1/mu_0); //Q = (H+M)

    Tensor1_3d Q( vQ );
    Tensor2_3d F( vF ), P_me, GQQ;
    auto [J, cofF] = F.compute_determinant_and_cofactor();
    auto Fcross    = F.compute_tensor_cross_product();
    GQQ( i_3, j_3 ) = cofF( i_3, k_3 ) * Q( k_3 ) * Q( j_3 );
    auto GQQcross   = GQQ.compute_tensor_cross_product();

    Tensor4_3d D_me;
    
    D_me( i_3, j_3, k_3, l_3 ) = (-mu_0/J)*(GQQcross(i_3, j_3, k_3, l_3) + ( Fcross( i_3, j_3, m_3, n_3 ) * Q(n_3) ) * ( Fcross( m_3, q_3, k_3, l_3 ) * Q( q_3 ) ))
                               + ( mu_0 / ( 2 * J * J ) ) * ( cofF( m_3, n_3 ) * Q( n_3 ) * cofF( m_3, k_3 ) * Q( k_3 ) * Fcross( i_3, j_3, k_3, l_3 ) )
                               - ( mu_0 / ( J * J * J ) ) * (cofF(m_3, n_3) * Q(n_3) * cofF(m_3, k_3) * Q(k_3) * cofF(i_3, j_3) * cofF(k_3, l_3)
                               + ( mu_0 / (  J * J ) ) * cofF(i_3, j_3) * (Fcross(k_3, l_3, m_3, n_3) * cofF(m_3,q_3) * Q(q_3) * Q(n_3)));

    return D_me.to_voigt_form();

}

FloatArrayF<3> MooneyRivlinHardMagnetic::computeMagneticInduction_3d_multiphysics( const FloatArrayF<3> &vH, GaussPoint *gp, TimeStep *tStep ) const
{
    MagnetoelasticMaterialStatus *status = static_cast<MagnetoelasticMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );

    FloatArrayF<3> vQ = vH + B_res*(1/mu_0); //Q = (H+M)

    Tensor1_3d Q( vQ ), B_me;
    Tensor2_3d F( vF );
    auto [J, cofF] = F.compute_determinant_and_cofactor();
    
    B_me( i_3 ) = mu_0 / J * cofF( j_3, i_3 ) * cofF( j_3, k_3 ) * Q( k_3 );

    return B_me.to_voigt_form();
}

FloatMatrixF<3, 9> MooneyRivlinHardMagnetic::compute3dMagnetoelasticCouplingTensor_dBdF_multiphysics( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    MagnetoelasticMaterialStatus *status = static_cast<MagnetoelasticMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<3> vH( status->giveTempHVector() );
    FloatArrayF<9> vF( status->giveTempFVector() );

    FloatArrayF<3> vQ = vH + B_res*(1/mu_0); //Q = (H+M)

    Tensor1_3d Q( vQ );
    Tensor2_3d F( vF ), P_me, delta(1., 0., 0., 0., 1., 0., 0., 0., 1. );
    Tensor3_3d A_me;
    auto [J, cofF] = F.compute_determinant_and_cofactor();
    auto Fcross    = F.compute_tensor_cross_product();
    
    A_me( i_3, k_3, l_3 ) = mu_0 / J * ( delta( i_3, n_3 ) * cofF( m_3, q_3 ) * Q( q_3 ) + G( m_3, i_3 ) * Q( n_3 ) ) * Fcross( m_3, n_3, k_3, l_3 )
                          - mu_0 / ( J * J ) * ( cofF( p_3, i_3 ) * cofF( p_3, q_3 ) * Q( q_3 ) * G( k_3, l_3 ) );

    return A_me.to_voigt_form_3x9();
}

FloatMatrixF<3, 3> MooneyRivlinHardMagnetic::compute3dPermeabilityMatrix_dBdH_multiphysics( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    MagnetoelasticMaterialStatus *status = static_cast<MagnetoelasticMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    //independent of H, actually

    Tensor2_3d F( vF ), C_me;
    auto [J, cofF] = F.compute_determinant_and_cofactor();

    C_me( i_3, k_3 ) = mu_0 / J * cofF( j_3, i_3 ) * cofF( j_3, k_3 );

    return C_me.to_voigt_form();
}

// TODO remove, temporary debug
int
MooneyRivlinHardMagnetic :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_CrackVector ) {
        answer = B_app;
        return 1;
    } else if ( type == IST_2ndCrackVector ) {
        answer = H;
        return 1;
    }
    {
        return MooneyRivlinCompressibleMaterial::giveIPValue( answer, gp, type, tStep );
  }
}

} // end namespace oofem
