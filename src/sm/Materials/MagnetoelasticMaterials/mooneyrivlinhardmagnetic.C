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
REGISTER_Material( MooneyRivlinHardMagnetic );

FloatArrayF<9>
MooneyRivlinHardMagnetic::giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

    // get stress from Mooney-Rivlin hyperelasticity
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
    default:
        OOFEM_ERROR( "Unknown hard magnetic material mode %i", materialMode );
        break;
    }


    // add magnetoelastic and mooney-rivlin stress together
    auto vP = vP_mr + vP_me;

    // update gp
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );

    return vP;
}

FloatMatrixF<9, 9>
MooneyRivlinHardMagnetic::give3dMaterialStiffnessMatrix_dPdF( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    // get stiffness from Mooney-Rivlin hyperelasticity
    FloatMatrixF<9, 9> vD_mr = MooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dPdF( matResponseMode, gp, tStep );

    // calculate magnetoelastic stiffness
    FloatMatrixF<9, 9> vD_me;
    switch ( materialMode ) {
    case 1:
        vD_me = compute3dMaterialStiffnessMatrix_dPdF_consistent( matResponseMode, gp, tStep );
        break;
    case 2:
        vD_me = compute3dMaterialStiffnessMatrix_dPdF_mit( matResponseMode, gp, tStep );
        break;
    default:
        OOFEM_ERROR( "Unknown hard magnetic material mode %i", materialMode );
        break;
    }

    // add magnetoelastic and Mooney-Rivlin stiffness together
    auto vD = vD_me + vD_mr;

    return vD;
}


void MooneyRivlinHardMagnetic::initializeFrom( InputRecord &ir )
{
    MooneyRivlinCompressibleMaterial::initializeFrom( ir );

    FloatArray tempb0, tempM0;

    // default values
    mu_0         = BASE_VACUUM_PERMEABILITY_MU_0;
    materialMode = 1;

    IR_GIVE_FIELD( ir, tempb0, _IFT_MooneyRivlinHardMagnetic_b0 );
    IR_GIVE_FIELD( ir, tempM0, _IFT_MooneyRivlinHardMagnetic_M0 );
    // IR_GIVE_FIELD( ir, H_temp, _IFT_MooneyRivlinHardMagnetic_H );
    IR_GIVE_OPTIONAL_FIELD( ir, mu_0, _IFT_MooneyRivlinHardMagnetic_mu_0 );
    IR_GIVE_FIELD( ir, ltf_index, _IFT_MooneyRivlinHardMagnetic_ltf );
    IR_GIVE_FIELD( ir, mltf_index, _IFT_MooneyRivlinHardMagnetic_mltf );
    IR_GIVE_OPTIONAL_FIELD( ir, materialMode, _IFT_MooneyRivlinHardMagnetic_mode );

    referenceB = false;
    IR_GIVE_OPTIONAL_FIELD( ir, referenceB, _IFT_MooneyRivlinHardMagnetic_referenceB );

    b0 = FloatArrayF<3>( tempb0 );
    M0 = FloatArrayF<3>( tempM0 );


    auto pullBackTypeRecord = 0; // default
    IR_GIVE_OPTIONAL_FIELD( ir, pullBackTypeRecord, _IFT_MooneyRivlinHardMagnetic_PullBackType );

    if ( pullBackTypeRecord == 0 ) {
        this->pullBackType = PBT_F;
        this->exp          = 1.0;
    } else if ( pullBackTypeRecord == 1 ) {
        this->pullBackType = PBT_R;
        this->exp          = 0.5;
    } else if ( pullBackTypeRecord == 2 ) {
        this->pullBackType = PBT_iFt;
        this->exp          = 0.0;
    } else if ( pullBackTypeRecord == 3 ) {
        this->pullBackType = PBT_iFtnoJ;
        this->exp          = 0.0;
    } else {
        throw ValueInputException( ir, _IFT_MooneyRivlinHardMagnetic_PullBackType, "Unknown pull back type" );
    }
}

FloatArrayF<9> MooneyRivlinHardMagnetic::computeFirstPKStressVector_3d_consistent( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
{

    double load_level  = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    double mload_level = this->giveDomain()->giveFunction( mltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    /*    FloatArrayF<3> b = load_level  * b0;
    FloatArrayF<3> m = mload_level * m0;
    */
    Tensor2_3d F( vF ), P_me, delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. ), C;
    Tensor1_3d b( load_level * b0 ), M( mload_level * M0 ), B;
    // b(i_3) += mu_0 * F(i_3, k_3) * M(k_3);
    C( i_3, j_3 )  = F( k_3, i_3 ) * F( k_3, j_3 );
    auto [J, cofF] = F.compute_determinant_and_cofactor();
    B( j_3 )       = b( k_3 ) * cofF( k_3, j_3 );
    //
    Tensor4_3d dCm_dC;
    dCm_dC = C.compute_dCm_dC( this->exp );
    Tensor2_3d t;
    t( p_3, j_3 ) = ( M( k_3 ) * B( l_3 ) ) * dCm_dC( k_3, l_3, p_3, j_3 );
    //
    P_me( i_3, j_3 ) = 1. / J * F( i_3, p_3 ) * ( 1. / mu_0 * B( p_3 ) * B( j_3 ) - ( M( k_3 ) * B( l_3 ) ) * dCm_dC( k_3, l_3, p_3, j_3 ) );
    //- 1. / J / J * (0.5 / mu_0 * B(k_3) * F(m_3,k_3)*F(m_3,l_3) * B(l_3) - M(k_3) * C.computeTensorPower(this->exp)(k_3,l_3) * B(l_3) ) * cofF(i_3,j_3);
    //
    auto vP_me = P_me.to_voigt_form();
    //
    return vP_me;
}


FloatMatrixF<9, 9> MooneyRivlinHardMagnetic::compute3dMaterialStiffnessMatrix_dPdF_consistent( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    //
    Tensor2_3d F( vF );
    Tensor4_3d D_me;
    //
    double load_level  = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    double mload_level = this->giveDomain()->giveFunction( mltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //
    Tensor2_3d P_me, delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
    Tensor1_3d b( load_level * b0 ), M( mload_level * M0 ), B;
    auto [J, cofF] = F.compute_determinant_and_cofactor();
    B( j_3 )       = b( k_3 ) * cofF( k_3, j_3 );
    //
    D_me( i_3, j_3, k_3, l_3 ) = 0.0 * ( 1. / J * delta( i_3, k_3 ) * ( 1. / mu_0 * B( l_3 ) * B( j_3 ) - M( l_3 ) * B( j_3 ) - B( l_3 ) * M( j_3 ) ) - -1. / J / J * F( i_3, p_3 ) * ( 1. / mu_0 * B( p_3 ) * B( j_3 ) - M( p_3 ) * B( j_3 ) - B( p_3 ) * M( j_3 ) ) * cofF( k_3, l_3 ) );
    //
    //    auto vD_me = ;
    return D_me.to_voigt_form();
}

FloatArrayF<9> MooneyRivlinHardMagnetic::computeFirstPKStressVector_3d_mit( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
{

    double load_level  = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    double mload_level = this->giveDomain()->giveFunction( mltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    //
    Tensor2_3d F( vF ), P_me, delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. ), Finv, U, invU, Y, R;
    Tensor1_3d b( load_level * b0 ), M( mload_level * M0 ), B;
    Tensor4_3d RR;

    int J, dY;

    switch ( this->pullBackType ) {

    case PBT_F:
        P_me( k_3, l_3 ) = -1.0 * b( k_3 ) * M( l_3 );
        break;
    case PBT_R:
        U = this->computeU(F);
        Y( i_3, j_3 ) = U( k_3, l_3 ) * delta( k_3, l_3 ) * delta( i_3, j_3 ) - U( i_3, j_3 );
        
        dY = Y.compute_determinant();
        invU = U.compute_inverse();
        R( i_3, j_3 )            = F( i_3, k_3 ) * invU( k_3, j_3 );
        RR( k_3, l_3, m_3, n_3 ) = 1 / dY * ( R( k_3, p_3 ) * ( Y( p_3, q_3 ) * R( m_3, q_3 ) * Y( n_3, l_3 ) - Y( p_3, n_3 ) * R( m_3, q_3 ) * Y( q_3, l_3 ) ) );
        P_me( k_3, l_3 ) = -1.0 * RR( i_3, q_3, k_3, l_3 ) * M( q_3 ) * b( i_3 );
        break;
    case PBT_iFt:
        OOFEM_ERROR( "Pullback type not implemented." )
        break;
    case PBT_iFtnoJ:
        J = F.compute_determinant_and_cofactor().first;
        Finv = F.compute_inverse();

        P_me( i_3, j_3 ) = - J * Finv( j_3, i_3 ) * Finv( q_3, p_3 ) * M( q_3 ) * b( p_3 )
                           + J * Finv( q_3, i_3 ) * Finv( j_3, p_3 ) * M( q_3 ) * b( p_3 ); 
        break;
    default:
        OOFEM_ERROR( "Unknown pullback type" );
        break;
    }

    auto vP_me = P_me.to_voigt_form();
    return vP_me;
}


FloatMatrixF<9, 9> MooneyRivlinHardMagnetic::compute3dMaterialStiffnessMatrix_dPdF_mit( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
    FloatMatrixF<9, 9> vD_me;

    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveFVector() );

    double load_level  = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    double mload_level = this->giveDomain()->giveFunction( mltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );

    Tensor1_3d M( mload_level * M0 ), b( load_level * b0 );
    Tensor2_3d F( vF ), Finv;
    Tensor4_3d D;
    int J;

    switch ( this->pullBackType ) {

    case PBT_F:
        // stiffness is zero for this case (the original MIT case)
        break;
    case PBT_R:
        // stiffness not implemented for now
        break;
    case PBT_iFt:
        OOFEM_ERROR( "Pullback type not implemented." )
        break;
    case PBT_iFtnoJ:
        J = F.compute_determinant_and_cofactor().first;
        Finv = F.compute_inverse();

        D( i_3, j_3, k_3, l_3 ) = - J * M( q_3 ) * b( p_3 ) * Finv( q_3, p_3 ) * Finv( l_3, k_3 ) * Finv( j_3, i_3 )
            + J * M( q_3 ) * b( p_3 ) * Finv( q_3, p_3 ) * Finv( j_3, k_3 ) * Finv( l_3, i_3 )
            - J * M( q_3 ) * b( p_3 ) * Finv( j_3, i_3 ) * Finv( q_3, k_3 ) * Finv( l_3, p_3 )
            + J * M( q_3 ) * b( p_3 ) * Finv( l_3, k_3 ) * Finv( q_3, i_3 ) * Finv( j_3, p_3 )
            - J * M( q_3 ) * b( p_3 ) * Finv( q_3, k_3 ) * Finv( l_3, i_3 ) * Finv( j_3, p_3 )
            + J * M( q_3 ) * b( p_3 ) * Finv( q_3, i_3 ) * Finv( j_3, k_3 ) * Finv( l_3, p_3 );

        vD_me = D.to_voigt_form();
        break;
    default:

        OOFEM_ERROR( "Unknown pullback type" );
        break;
    }
    return vD_me;
}


int MooneyRivlinHardMagnetic::giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep )
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

    if ( type == IST_CauchyStressTensor || type == IST_StrainTensorTemp || type == IST_StressTensorTemp ) {
        FloatArrayF<9> vF( status->giveFVector() );
        Tensor2_3d F( vF );
        //
        FloatArrayF<9> vP( status->givePVector() );
        Tensor2_3d P( vP ), sigma, sigma_mech, sigma_mag;
        //
        auto [J, cofF]        = F.compute_determinant_and_cofactor();
        sigma( i_3, j_3 )     = P( i_3, k_3 ) * cofF.compute_inverse()( k_3, j_3 );
        FloatArrayF<9> vPmech = MooneyRivlinCompressibleMaterial::giveFirstPKStressVector_3d( vF, gp, tStep );
        Tensor2_3d Pmech( vPmech );
        sigma_mech( i_3, j_3 ) = Pmech( i_3, k_3 ) * cofF.compute_inverse()( k_3, j_3 );
        if ( type == IST_CauchyStressTensor ) {
            answer = sigma.to_voigt_form();
        } else if ( type == IST_StrainTensorTemp ) {
            answer = sigma_mech.to_voigt_form();
        } else if ( type == IST_StressTensorTemp ) {
            sigma_mag( i_3, j_3 ) = sigma( i_3, j_3 ) - sigma_mech( i_3, j_3 );
            answer                = sigma_mag.to_voigt_form();
        }
        return 1;
    } else {
        return StructuralMaterial::giveIPValue( answer, gp, type, tStep );
    }
}

Tensor2_3d MooneyRivlinHardMagnetic::computeU( const Tensor2_3d &F ) const
{
    Tensor2_3d C;
    C( i_3, j_3 ) = F( k_3, i_3 ) * F( k_3, j_3 );

    FloatMatrixF<3, 3> mC;
    mC = C.to_matrix_form();
    //
    FloatArray eVals;
    FloatMatrix eVecs, mmC( mC );
    mmC.jaco_( eVals, eVecs, 15 );
    //
    FloatMatrixF<3, 3> mU;
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            mU.at( i, j ) = sqrt( eVals.at( 1 ) ) * eVecs.at( i, 1 ) * eVecs.at( j, 1 ) + sqrt( eVals.at( 2 ) ) * eVecs.at( i, 2 ) * eVecs.at( j, 2 ) + sqrt( eVals.at( 3 ) ) * eVecs.at( i, 3 ) * eVecs.at( j, 3 );
        }
    }
    Tensor2_3d U( mU );
    return U;
}


} // end namespace oofem
