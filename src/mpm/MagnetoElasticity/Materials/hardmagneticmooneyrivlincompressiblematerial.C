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

#include "hardmagneticmooneyrivlincompressiblematerial.h"
#include "mpm/MagnetoElasticity/Materials/magnetoelasticmaterialstatus.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "tensor/tensor3.h"
#include "tensor/tensor1.h"
#include "tensor/tensor2.h"
#include "tensor/tensor4.h"
#include "tensor/tensor5.h"
#include "tensor/tensor6.h"
#include "tensor/tensor7.h"
#include "domain.h"
#include "function.h"

//debug
#include "tensor/FTensor/Tensor4/Tensor4_times_Tensor2_outer.hpp"


namespace oofem {

  
REGISTER_Material( HardMagneticMooneyRivlinCompressibleMaterial );

HardMagneticMooneyRivlinCompressibleMaterial::HardMagneticMooneyRivlinCompressibleMaterial( int n, Domain *d ) :
    MagnetoElasticMaterial( n, d ), BaseHyperElasticMaterial()
{
}


std::tuple<FloatArrayF<9>, FloatArrayF<3> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_FirstPKStressVector_MagneticInduction_3d( const FloatArrayF<9> &vF, const FloatArrayF<3> &vH, GaussPoint *gp, TimeStep *tStep ) const
{

  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  
  Tensor1_3d H( vH ), M_full;
  Tensor2_3d F( vF );

  double m_level = this->giveDomain()->giveFunction( m_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
  M_full( i_3 ) = m_level * M( i_3 );

  // compute everything w.r.t. the original reference configuration
  auto [P, B] = this->computeFirstPKStressMagneticInductionTensors_3d(F, H, M_full);

  auto vP = P.to_voigt_form();
  auto vB = B.to_voigt_form();
  // update gp
  status->letTempFVectorBe( vF );
  status->letTempPVectorBe( vP );
  status->letTempHVectorBe( vH );
  status->letTempBVectorBe( vB );
  //
  return std::make_tuple( vP, vB );
}


std::tuple<FloatArrayF<5>, FloatArrayF<2> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_FirstPKStressVector_MagneticInduction_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<2> &vH, GaussPoint *gp, TimeStep *tStep ) const
{
  auto [P, B] = give_FirstPKStressVector_MagneticInduction_3d( assemble<9>( vF, { 0, 1, 2, 5, 8 } ), assemble<3>( vH, { 0, 1 } ), gp, tStep );
  //
  return std::make_tuple( P[{ 0, 1, 2, 5, 8 }], B[{ 0, 1 }] );
}


std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 3>, FloatMatrixF<3, 9>, FloatMatrixF<3, 3> >
HardMagneticMooneyRivlinCompressibleMaterial ::giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  FloatArrayF<9> vF( status->giveTempFVector() );
  FloatArrayF<3> vH( status->giveTempHVector() );
  Tensor1_3d H( vH ), M_full;
  Tensor2_3d F( vF );

  double m_level = this->giveDomain()->giveFunction( m_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
  M_full( i_3 ) = m_level * M( i_3 );

  auto [dPdF, dPdH, dBdF, dBdH] = this->computeStiffnessTensors_dPdF_dBdH_dPdH_3d( F, H, M_full);

  /*dPdH( i_3, j_3, k_3 ) = 0. * dPdH( i_3, j_3, k_3 );
  dBdF( i_3, j_3, k_3 ) = 0. * dPdH( i_3, j_3, k_3 );*/

  //auto [A, B, C, D] = this->compute_stiff_num( vF, vH, gp, tStep );
  //
  return std::make_tuple( dPdF.to_voigt_form(), dPdH.to_voigt_form_9x3(), dBdF.to_voigt_form_3x9(), dBdH.to_matrix_form() );
}


std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 3>, FloatMatrixF<3, 9>, FloatMatrixF<3, 3> >
HardMagneticMooneyRivlinCompressibleMaterial ::compute_stiff_num( const FloatArray &vF, const FloatArray &vH, GaussPoint *gp, TimeStep *tStep )
{
  auto [P, B] = give_FirstPKStressVector_MagneticInduction_3d( vF, vH, gp, tStep );
  FloatMatrix Kuu( 9, 9 ), Kuh( 9, 3 ), Khu( 3, 9 ), Khh( 3, 3 );
  auto pert = 1.e-9;
  for ( int i = 1; i <= 9; i++ ) {
    auto vFp = vF;
    vFp.at( i ) += pert;
    auto [Pp1, Bp1] = give_FirstPKStressVector_MagneticInduction_3d( vFp, vH, gp, tStep );
    auto dP1 = Pp1 - P;
    auto dB1 = Bp1 - B;
    Kuu.setColumn( dP1, i );
    Khu.setColumn( dB1, i );
    if ( i <= 3 ) {
      auto vHp = vH;
      vHp.at( i ) += pert;
      auto [Pp2, Bp2] = give_FirstPKStressVector_MagneticInduction_3d( vF, vHp, gp, tStep );
      auto dP2 = Pp2 - P;
      auto dB2 = Bp2 - B;
      Kuh.setColumn( dP2, i );
      Khh.setColumn( dB2, i );
    }
  }
  Kuu.times( 1. / pert );
  Kuh.times( 1. / pert );
  Khu.times( 1. / pert );
  Khh.times( 1. / pert );
  return std::make_tuple( Kuu, Kuh, Khu, Khh );
}

std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 2>, FloatMatrixF<2, 5>, FloatMatrixF<2, 2> >
HardMagneticMooneyRivlinCompressibleMaterial ::giveConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  auto [dPdF, dPdH, dBdF, dBdH] = this->giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d( mode, gp, tStep );
  auto dPdFred = dPdF( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto dPdHred = dPdH( { 0, 1, 2, 5, 8 }, { 0, 1 } );
  auto dBdFred = dBdF( { 0, 1 }, { 0, 1, 2, 5, 8 } );
  auto dBdHred = dBdH( { 0, 1 }, { 0, 1 } );
  return std::make_tuple( dPdFred, dPdHred, dBdFred, dBdHred );
}

std::tuple<Tensor2_3d, Tensor1_3d> HardMagneticMooneyRivlinCompressibleMaterial::computeFirstPKStressMagneticInductionTensors_3d( const Tensor2_3d &F, const Tensor1_3d &H, const Tensor1_3d &M) const
{
  Tensor2_3d P;
  Tensor1_3d B, Q;

  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );

  auto [J, cofF] = F.compute_determinant_and_cofactor();

  if (pullBackType == PBT_iFt ) {
    Q( i_3 ) = H( i_3 ) + M( i_3 );
    P( i_3, j_3 ) =
      // mooney rivlin part \/
      C1 * this->compute_dI1_Cdev_dF( F )( i_3, j_3 ) + C2 * this->compute_dI2_Cdev_dF( F )( i_3, j_3 ) + this->compute_dVolumetricEnergy_dF( F )( i_3, j_3 )
      // magneto part \/
      - ( mu_0 / J ) * F.compute_tensor_cross_product()( i_3, j_3, m_3, n_3 ) * ( cofF( m_3, q_3 ) * Q( q_3 ) * Q( n_3 ) )
      + ( mu_0 / ( 2 * J * J ) ) * cofF( m_3, n_3 ) * Q( n_3 ) * cofF( m_3, k_3 ) * Q( k_3 ) * cofF( i_3, j_3 )
      // isole part \/
      + lambda_isole * delta( i_3, m_3 ) * ( F( k_3, k_3 ) - 3. ) * cofF( m_3, j_3 ) + mu_isole * ( F( i_3, m_3 ) + F( m_3, i_3 ) - 2. * delta( i_3, m_3 ) ) * cofF( m_3, j_3 );

    B( i_3 ) = mu_0 / J * cofF( j_3, i_3 ) * cofF( j_3, k_3 ) * Q( k_3 );

  } else if ( pullBackType == PBT_F ) {

    P( i_3, j_3 ) =
      // mooney rivlin part \/
      C1 * this->compute_dI1_Cdev_dF( F )( i_3, j_3 ) + C2 * this->compute_dI2_Cdev_dF( F )( i_3, j_3 ) + this->compute_dVolumetricEnergy_dF( F )( i_3, j_3 )
      // magneto part \/
      - ( mu_0 / J ) * F.compute_tensor_cross_product()( i_3, j_3, m_3, n_3 ) * ( cofF( m_3, q_3 ) * H( q_3 ) * H( n_3 ) )
      + ( mu_0 / ( 2 * J * J ) ) * cofF( m_3, n_3 ) * H( n_3 ) * cofF( m_3, k_3 ) * H( k_3 ) * cofF( i_3, j_3 )
      //  -- M*C*M term
      - ( mu_0 / J ) * F( i_3, k_3) * M( k_3 ) * M( j_3 )  
      + ( mu_0 / ( 2. * J * J ) ) * F( m_3, n_3 ) * M( n_3 ) * F( m_3, k_3 ) * M( k_3 ) * cofF( i_3, j_3 )
      // isole part \/
      + lambda_isole * delta( i_3, m_3 ) * ( F( k_3, k_3 ) - 3. ) * cofF( m_3, j_3 ) + mu_isole * ( F( i_3, m_3 ) + F( m_3, i_3 ) - 2. * delta( i_3, m_3 ) ) * cofF( m_3, j_3 );

    B( i_3 ) = mu_0 / J * cofF( j_3, i_3 ) * cofF( j_3, k_3 ) * H( k_3 ) + mu_0 * M(i_3);
  } else {
    OOFEM_ERROR("Wrong pull back type");

  }

  return std::make_tuple( P, B );
}

std::tuple<Tensor4_3d, Tensor3_3d, Tensor3_3d, Tensor2_3d> HardMagneticMooneyRivlinCompressibleMaterial::computeStiffnessTensors_dPdF_dBdH_dPdH_3d( const Tensor2_3d &F, const Tensor1_3d &H, const Tensor1_3d &M ) const
{
  Tensor4_3d dPdF;
  Tensor3_3d dPdH, dBdF, dPdH_in_H;
  Tensor2_3d dBdH, GQQ;
  Tensor1_3d Q;
  //  
  auto [J, cofF] = F.compute_determinant_and_cofactor();
  auto Fcross = F.compute_tensor_cross_product();
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  if(pullBackType == PBT_iFt) {  
  
    Q( i_3 ) = H( i_3 ) + M( i_3 );


    GQQ( i_3, j_3 ) = cofF( i_3, k_3 ) * Q(k_3 ) * Q( j_3 );
    auto GQQcross = GQQ.compute_tensor_cross_product();

    dPdF( i_3, j_3, k_3, l_3 ) =
      // mooney rivlin part \/
      C1 * this->compute_d2I1_Cdev_dF2( F )( i_3, j_3, k_3, l_3 )
      + C2 * this->compute_d2I2_Cdev_dF2( F )( i_3, j_3, k_3, l_3 )
      + this->compute_d2VolumetricEnergy_dF2( F )( i_3, j_3, k_3, l_3 )
      // magneto part \/
        - mu_0 / J * ( GQQcross( i_3, j_3, k_3, l_3 ) + ( Fcross( i_3, j_3, m_3, n_3 ) * Q( n_3 ) ) * ( Fcross( m_3, q_3, k_3, l_3 ) * Q( q_3 ) ) )
      + mu_0 / J / J * Fcross( i_3, j_3, m_3, n_3 ) * GQQ( m_3, n_3 ) * cofF( k_3, l_3 )
      + ( mu_0 / ( 2 * J * J ) ) * ( cofF( m_3, n_3 ) * Q( n_3 ) * cofF( m_3, k_3 ) * Q( k_3 ) * Fcross( i_3, j_3, k_3, l_3 ) )
      - ( mu_0 / ( J * J * J ) ) * ( cofF( m_3, n_3 ) * Q( n_3 ) * cofF( m_3, k_3 ) * Q( k_3 ) * cofF( i_3, j_3 ) * cofF( k_3, l_3 ) )
      + ( mu_0 / ( J * J ) ) * cofF( i_3, j_3 ) * ( Fcross( k_3, l_3, m_3, n_3 ) * cofF( m_3, q_3 ) * Q( q_3 ) * Q( n_3 ) )
      // isole part  
      + lambda_isole * delta( i_3, m_3 ) * delta( k_3, l_3 ) * cofF( m_3, j_3 ) + mu_isole * ( delta( i_3, k_3 ) * delta( m_3, l_3 ) + delta( i_3, l_3 ) * delta( m_3, k_3 ) ) * cofF(m_3, j_3)
      + lambda_isole * delta( i_3, m_3 ) * ( F( k_3, k_3 ) - 3. ) * Fcross( m_3, j_3, k_3, l_3 ) + mu_isole * ( F( i_3, m_3 ) + F( m_3, i_3 ) - 2. * delta( i_3, m_3 ) ) * Fcross( m_3, j_3, k_3, l_3 );
      

    dBdH( i_3, k_3 ) = -mu_0 / J * cofF( j_3, i_3 ) * cofF( j_3, k_3 );
    //
    // dBdF( i_3, k_3, l_3 ) = mu_0 / J * ( G( j_3, m_3 ) * Q( m_3 ) * Fcross( j_3, i_3, k_3, l_3 ) ) + G( j_3, i_3 ) * ( Q( m_3 ) * Fcross( j_3, m_3, k_3, l_3 ) ) - mu_0 / ( J * J ) * ( G( j_3, i_3 ) * G( j_3, m_3 ) * Q( m_3 ) * G( k_3, l_3 ) );
    //
    dPdH( i_3, j_3, k_3 ) = mu_0 / J * ( Fcross( i_3, j_3, m_3, n_3 ) * Q( n_3 ) * cofF( m_3, k_3 ) + Fcross( i_3, j_3, m_3, k_3 ) * cofF( m_3, q_3 ) * Q( q_3 ) ) - mu_0 / J / J * ( cofF( m_3, k_3 ) * cofF( m_3, o_3 ) * Q( o_3 ) ) * cofF( i_3, j_3 );
    //
  } else if(pullBackType == PBT_F) {
    //
    Q( i_3 ) = H( i_3 );
    //
    GQQ( i_3, j_3 ) = cofF( i_3, k_3 ) * Q(k_3 ) * Q( j_3 );
    auto GQQcross = GQQ.compute_tensor_cross_product();
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dPdF( i_3, j_3, k_3, l_3 ) = C1 * this->compute_d2I1_Cdev_dF2( F )( i_3, j_3, k_3, l_3 )
      + C2 * this->compute_d2I2_Cdev_dF2( F )( i_3, j_3, k_3, l_3 )
        + this->compute_d2VolumetricEnergy_dF2( F )( i_3, j_3, k_3, l_3 )
        // isole part
        + lambda_isole * delta( i_3, m_3 ) * delta( k_3, l_3 ) * cofF( m_3, j_3 ) + mu_isole * ( delta( i_3, k_3 ) * delta( m_3, l_3 ) + delta( i_3, l_3 ) * delta( m_3, k_3 ) ) * cofF( m_3, j_3 )
        + lambda_isole * delta( i_3, m_3 ) * ( F( k_3, k_3 ) - 3. ) * Fcross( m_3, j_3, k_3, l_3 ) + mu_isole * ( F( i_3, m_3 ) + F( m_3, i_3 ) - 2. * delta( i_3, m_3 ) ) * Fcross( m_3, j_3, k_3, l_3 );
      //
            - mu_0 / J * ( GQQcross( i_3, j_3, k_3, l_3 ) + ( Fcross( i_3, j_3, m_3, n_3 ) * Q( n_3 ) ) * ( Fcross( m_3, q_3, k_3, l_3 ) * Q( q_3 ) ) )
      + mu_0 / J / J * Fcross( i_3, j_3, m_3, n_3 ) * GQQ( m_3, n_3 ) * cofF( k_3, l_3 )
      + ( mu_0 / ( 2 * J * J ) ) * ( cofF( m_3, n_3 ) * Q( n_3 ) * cofF( m_3, k_3 ) * Q( k_3 ) * Fcross( i_3, j_3, k_3, l_3 ) )
      - ( mu_0 / ( J * J * J ) ) * ( cofF( m_3, n_3 ) * Q( n_3 ) * cofF( m_3, k_3 ) * Q( k_3 ) * cofF( i_3, j_3 ) * cofF( k_3, l_3 ) )
      + ( mu_0 / ( J * J ) ) * cofF( i_3, j_3 ) * ( Fcross( k_3, l_3, m_3, n_3 ) * cofF( m_3, q_3 ) * Q( q_3 ) * Q( n_3 ) )
      // MCM term
      - mu_0 / J * delta(i_3, k_3) * (M(j_3) * M(l_3))
      + mu_0 / J / J * F(i_3, p_3) * M(p_3) * M(j_3) * cofF(k_3, l_3)
      +  mu_0 / ( 2. * J * J )  * ( F( m_3, n_3 ) * M( n_3 ) * F( m_3, o_3 ) * M( o_3 )) * Fcross( i_3, j_3, k_3, l_3 ) 
      - mu_0 / J / J / J * F( m_3, n_3 ) * M( n_3 ) * F( m_3, o_3 ) * M( o_3 ) * cofF( i_3, j_3 )* cofF( k_3, l_3 )
      + mu_0 / J / J * cofF(i_3, j_3) * F(k_3, p_3) * (M(p_3) * M(l_3));
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
    dBdH( i_3, k_3 ) = -mu_0 / J * cofF( j_3, i_3 ) * cofF( j_3, k_3 );
    //
    // dBdF( i_3, k_3, l_3 ) = mu_0 / J * ( G( j_3, m_3 ) * Q( m_3 ) * Fcross( j_3, i_3, k_3, l_3 ) ) + G( j_3, i_3 ) * ( Q( m_3 ) * Fcross( j_3, m_3, k_3, l_3 ) ) - mu_0 / ( J * J ) * ( G( j_3, i_3 ) * G( j_3, m_3 ) * Q( m_3 ) * G( k_3, l_3 ) );
    //
    dPdH( i_3, j_3, k_3 ) = mu_0 / J * ( Fcross( i_3, j_3, m_3, n_3 ) * Q( n_3 ) * cofF( m_3, k_3 ) + Fcross( i_3, j_3, m_3, k_3 ) * cofF( m_3, q_3 ) * Q( q_3 ) ) - mu_0 / J / J * ( cofF( m_3, k_3 ) * cofF( m_3, o_3 ) * Q( o_3 ) ) * cofF( i_3, j_3 );
    //
  }
    dBdF( i_3, k_3, l_3 ) = dPdH( k_3, l_3, i_3 );
    //
  //
    

    
  return std::make_tuple( dPdF, dPdH, dBdF, dBdH);
}

Tensor1_3d
HardMagneticMooneyRivlinCompressibleMaterial:: giveLagrangianMagnetization(TimeStep *tStep)
{

    double m_level = this->giveDomain()->giveFunction( m_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    Tensor1_3d M_time;
    M_time( i_3 ) = m_level * M( i_3 );
    return M_time;
}







MaterialStatus *
HardMagneticMooneyRivlinCompressibleMaterial::CreateStatus( GaussPoint *gp ) const
{
  return new MagnetoElasticMaterialStatus( gp );
}


void HardMagneticMooneyRivlinCompressibleMaterial::initializeFrom( InputRecord &ir )
{
  Material::initializeFrom( ir );
  BaseHyperElasticMaterial::initializeFrom( ir );
  IR_GIVE_FIELD( ir, C1, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_c1 );
  IR_GIVE_FIELD( ir, C2, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_c2 );
  FloatArray m, h_app(3);
  IR_GIVE_FIELD( ir, m, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_magnetization );
  IR_GIVE_OPTIONAL_FIELD( ir, h_app, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_happ );
  IR_GIVE_FIELD( ir, m_ltf, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_m_ltf );
  if ( m.giveSize() != 3 ) {
    OOFEM_ERROR( "Magnetization has to be vector of size 3" );
  } else {
    M = Tensor1_3d( FloatArrayF<3>( m ) );
  }
  if ( h_app.giveSize() != 3 ) {
    OOFEM_ERROR( "H applied has to be vector of size 3" );
  } else {
    this->h_app = Tensor1_3d( FloatArrayF<3>( h_app ) );
    if(h_app.computeNorm() != 0) {
      IR_GIVE_FIELD( ir, hload_ltf, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_hload_ltf );
    }
  }

  //isole simulating part
  double E_isole = 0., nu_isole = 0.;
  IR_GIVE_OPTIONAL_FIELD( ir, E_isole, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_Eisole );
  IR_GIVE_OPTIONAL_FIELD( ir, nu_isole, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_nuisole );
  this->lambda_isole = E_isole*nu_isole/((1.+nu_isole)*(1.-2.*nu_isole));
  this->mu_isole = E_isole/(2.*(1.+nu_isole));

  IR_GIVE_OPTIONAL_FIELD( ir, this->kappaGradFGradF, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_kappaGradFGradF );
  IR_GIVE_OPTIONAL_FIELD( ir, this->kappaGradJGradJ, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_kappaGradJGradJ );
  IR_GIVE_OPTIONAL_FIELD( ir, this->kappaGradRGradR, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_kappaGradRGradR );
  IR_GIVE_OPTIONAL_FIELD( ir, this->kappaFbar, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_kappaFbar );
  IR_GIVE_OPTIONAL_FIELD( ir, this->kappaJbar, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_kappaJbar );

  std::string pullBackTypeString; // default
  IR_GIVE_OPTIONAL_FIELD( ir, pullBackTypeString, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_PullBackType );

  if ( pullBackTypeString == "F" ) {
    this->pullBackType = PBT_F;
  } else if ( pullBackTypeString == "R" ) {
    this->pullBackType = PBT_R;
  } else if ( pullBackTypeString == "iFt" ) {
    this->pullBackType = PBT_iFt;
  } else {
    throw ValueInputException( ir, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_PullBackType, "Unknown pull back type" );
  }

  
  
}


int HardMagneticMooneyRivlinCompressibleMaterial ::giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep )
{

  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  if ( type == IST_DeformationGradientTensor ) {
    answer = status->giveFVector();
    return 1;
  } else if ( type == IST_FirstPKStressTensor ) {
    answer = status->givePVector();
    return 1;
  } else if ( type == IST_LagrangianMagneticInductionVector ) {
    answer = status->giveBVector();
    return 1;
  } else if ( type == IST_LagrangianMagneticFieldVector ) {
    answer = status->giveHVector();
    return 1;
  } else if ( type == IST_LagrangianMagnetizationVector ) {
    double m_level = this->giveDomain()->giveFunction( m_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    Tensor1_3d M_time;
    M_time( i_3 ) = m_level * M( i_3 );
    answer = M_time.to_voigt_form();
    return 1;
  } else if ( type == IST_CauchyStressTensor ) {
    FloatArrayF<9> vF( status->giveFVector() );
    FloatArrayF<9> vP( status->givePVector() );
    Tensor2_3d F( vF ), P( vP ), sigma;
    auto [J, G] = F.compute_determinant_and_cofactor();
    sigma( i_3, j_3 ) = 1. / J * P( i_3, k_3 ) * F( j_3, k_3 );
    answer = sigma.to_voigt_form();
    return 1;
  } else if ( type == IST_EulerianMagneticInductionVector ) {
    FloatArrayF<3> vB( status->giveBVector() );
    FloatArrayF<9> vF( status->giveFVector() );
    Tensor1_3d B( vB ), b;
    Tensor2_3d F( vF );
    auto [J, G] = F.compute_determinant_and_cofactor();
    b( i_3 ) = 1. / J * F( i_3, j_3 ) * B( j_3 );
    answer = b.to_voigt_form();
    return 1;
  } else if ( type == IST_EulerianMagneticFieldVector ) {
    FloatArrayF<3> vH( status->giveHVector() );
    FloatArrayF<9> vF( status->giveFVector() );
    Tensor1_3d H( vH ), h;
    Tensor2_3d F( vF );
    auto [J, G] = F.compute_determinant_and_cofactor();
    h( i_3 ) = (1./J) * G( i_3, j_3 ) * H( j_3 );
    answer = h.to_voigt_form();
    return 1;
  } else if ( type == IST_EulerianMagnetizationVector ) {
    FloatArrayF<3> vH( status->giveHVector() );
    FloatArrayF<3> vB( status->giveBVector() );
    FloatArrayF<9> vF( status->giveFVector() );
    Tensor1_3d H( vH ), h, B(vB), b, m;
    Tensor2_3d F( vF );
    auto [J, G] = F.compute_determinant_and_cofactor();
    h( i_3 ) = ( 1. / J ) * G( i_3, j_3 ) * H( j_3 );
    b( i_3 ) = 1. / J * F( i_3, j_3 ) * B( j_3 );
    m( i_3 ) = 1. / mu_0 * b( i_3 ) - h( i_3 ); 
    answer = m.to_voigt_form();
    return 1;
  } else if ( type == IST_EulerianHFromMVector ) {
    double h_load_level = this->giveDomain()->giveFunction( hload_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    FloatArrayF<3> vH( status->giveHVector() );
    FloatArrayF<9> vF( status->giveFVector() );
    Tensor1_3d H( vH ), h, h_from_M;
    Tensor2_3d F( vF );
    auto [J, G] = F.compute_determinant_and_cofactor();
    h_from_M( i_3 ) = ( 1. / J ) * G( i_3, j_3 ) * H( j_3 ) - h_load_level * h_app( i_3 );
    answer = h_from_M.to_voigt_form();
    return 1;
  } else {
    return Material::giveIPValue( answer, gp, type, tStep );
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////// Second Gradient ///////////////////////////////////////////////////



FloatArrayF< 8 >
HardMagneticMooneyRivlinCompressibleMaterial :: give_SecondGradient_FirstPKStressVector_PlaneStrain(const FloatArrayF<8> &vG, GaussPoint *gp, TimeStep *tStep)
  {
    auto vT = kappaGradFGradF * vG;
    return vT;
  }


FloatArrayF< 27 >
HardMagneticMooneyRivlinCompressibleMaterial :: give_SecondGradient_FirstPKStressVector_3d(const FloatArrayF<27> &vG, GaussPoint *gp, TimeStep *tStep)
  {
    auto vT = kappaGradFGradF * vG;
    return vT;
  }


FloatMatrixF<1,1>
HardMagneticMooneyRivlinCompressibleMaterial :: give_SecondGradient_ConstitutiveMatrix_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  return kappaGradFGradF;
}

FloatMatrixF<1,1>
HardMagneticMooneyRivlinCompressibleMaterial :: give_SecondGradient_ConstitutiveMatrix_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  return kappaGradFGradF;
}

////////////////////////FIRST AND SECOND GRADIENTS TOGETHER////////////////////////////

std::tuple<FloatArrayF<5>, FloatArrayF<8> >
HardMagneticMooneyRivlinCompressibleMaterial::give_FirstSecondGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep )
{
  FloatArrayF<5> vP_ans;
  FloatArrayF<8> vT_ans;

  if ( kappaGradJGradJ != 0.0 ) {
    auto [vP_gradJ, vT_gradJ] = this->give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( vF, vGradF, gp, tStep );
    vP_ans += vP_gradJ;
    vT_ans += vT_gradJ;
  }
  if ( kappaGradRGradR != 0.0 ) {
    auto [vP_gradR, vT_gradR] = this->give_RotationGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( vF, vGradF, gp, tStep );
    vP_ans += vP_gradR;
    vT_ans += vT_gradR;
  }

  return std::make_tuple( vP_ans, vT_ans );
}

std::tuple<FloatArrayF<9>, FloatArrayF<27> >
HardMagneticMooneyRivlinCompressibleMaterial::give_FirstSecondGradient_FirstPKStressVector_SecondOrderStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep )
{
  FloatArrayF<9> vP_ans;
  FloatArrayF<27> vT_ans;

  if ( kappaGradJGradJ != 0.0 ) {
    auto [vP_gradJ, vT_gradJ] = this->give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_3d( vF, vGradF, gp, tStep );
    vP_ans += vP_gradJ;
    vT_ans += vT_gradJ;
  }
  if ( kappaGradRGradR != 0.0 ) {
    OOFEM_ERROR( "gradRgradR regularization unsupported outside of Plane Strain formulations." );
  }

  return std::make_tuple( vP_ans, vT_ans );
}

std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> >
HardMagneticMooneyRivlinCompressibleMaterial::give_FirstSecondGradient_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  FloatMatrixF<9, 9> vdPdF_ans;
  FloatMatrixF<9, 27> vdPdGradF_ans;
  FloatMatrixF<27, 9> vdTdF_ans;
  FloatMatrixF<27, 27> vdTdGradF_ans;

  if ( kappaGradJGradJ != 0.0 ) {
    auto [vdPdF_gradJ, vdPdGradF_gradJ, vdTdF_gradJ, vdTdGradF_gradJ] = this->give_JacobianGradient_ConstitutiveMatrices_3d( mode, gp, tStep );
    vdPdF_ans += vdPdF_gradJ;
    vdPdGradF_ans += vdPdGradF_gradJ;
    vdTdF_ans += vdTdF_gradJ;
    vdTdGradF_ans += vdTdGradF_gradJ;
  }
  if ( kappaGradRGradR != 0.0 ) {
    OOFEM_ERROR( "gradRgradR regularization unsupported outside of Plane Strain formulations." );
  }

  return std::make_tuple( vdPdF_ans, vdPdGradF_ans, vdTdF_ans, vdTdGradF_ans );
}

std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> >
HardMagneticMooneyRivlinCompressibleMaterial::give_FirstSecondGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  FloatMatrixF<5, 5> vdPdF_ans;
  FloatMatrixF<5, 8> vdPdGradF_ans;
  FloatMatrixF<8, 5> vdTdF_ans;
  FloatMatrixF<8, 8> vdTdGradF_ans;

  if ( kappaGradJGradJ != 0.0 ) {
    auto [vdPdF_gradJ, vdPdGradF_gradJ, vdTdF_gradJ, vdTdGradF_gradJ] = this->give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( mode, gp, tStep );
    vdPdF_ans += vdPdF_gradJ;
    vdPdGradF_ans += vdPdGradF_gradJ;
    vdTdF_ans += vdTdF_gradJ;
    vdTdGradF_ans += vdTdGradF_gradJ;
  }
  if ( kappaGradRGradR != 0.0 ) {
    auto [vdPdF_gradR, vdPdGradF_gradR, vdTdF_gradR, vdTdGradF_gradR] = this->give_RotationGradient_ConstitutiveMatrices_PlaneStrain( mode, gp, tStep );
    vdPdF_ans += vdPdF_gradR;
    vdPdGradF_ans += vdPdGradF_gradR;
    vdTdF_ans += vdTdF_gradR;
    vdTdGradF_ans += vdTdGradF_gradR;
  }

  return std::make_tuple( vdPdF_ans, vdPdGradF_ans, vdTdF_ans, vdTdGradF_ans );
}


////////////////////////////////////////////////////////////////////////////////////////////////////// JacobianGradient /////////////////////////////////////////////////

std::tuple<FloatArrayF<5>, FloatArrayF<8> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep )
{
  // these elements should give all Tensor3 components without the third dimension in them
  //  according to the established Voigt notation, see Tensor3_3d::to_voigt_form_27()
  auto [vP, vT] = give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_3d( assemble<9>( vF, { 0, 1, 2, 5, 8 } ), assemble<27>( vGradF, { 0, 1, 5, 8, 9, 10, 14, 17 } ), gp, tStep );
  //
  return std::make_tuple( vP[{ 0, 1, 2, 5, 8 }], vT[{ 0, 1, 5, 8, 9, 10, 14, 17 }] );
}


std::tuple<FloatArrayF<9>, FloatArrayF<27> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep )
{
  Tensor2_3d P, F(vF);
  Tensor3_3d T, gradF(vGradF);

  auto gradJ = compute_gradJ_3d(F, gradF);
  auto [dGradJdF, dGradJdGradF] = compute_gradJ_derivatives_3d(F, gradF);

  P(r_3, s_3) = kappaGradJGradJ * gradJ(m_3) * dGradJdF(m_3,r_3,s_3);
  T(r_3, s_3, t_3) = kappaGradJGradJ * gradJ(m_3) * dGradJdGradF(m_3, r_3, s_3, t_3);

  //save gradients
  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  status->letTempFVectorBe( vF );
  status->letTempGradFVectorBe( vGradF );

  return std::make_tuple(P.to_voigt_form(), T.to_voigt_form_27());
}


std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_JacobianGradient_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  //load F and GradF
  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  FloatArrayF<9> vF = status->giveTempFVector();
  FloatArrayF<27> vGradF = status->giveTempGradFVector();
  
  Tensor2_3d F(vF);
  Tensor3_3d gradF(vGradF);

  Tensor4_3d D_FF;
  Tensor5_3d C_FG;
  Tensor6_3d C_GG;

  auto gradJ = compute_gradJ_3d( F, gradF );
  auto [dGradJdF, dGradJdGradF] = compute_gradJ_derivatives_3d( F, gradF );
  auto [d2GradJdFdF, d2GradJdFdGradF, d2GradJdGradFdF, d2GradJdGradFdGradF] = compute_gradJ_secondDerivatives_3d(F, gradF);

  D_FF( r_3, s_3, u_3, v_3 ) = kappaGradJGradJ * dGradJdF( m_3, u_3, v_3 ) * dGradJdF( m_3, r_3, s_3 ) + kappaGradJGradJ * gradJ( m_3 ) * d2GradJdFdF(m_3, r_3, s_3, u_3, v_3);
  C_FG( r_3, s_3, u_3, v_3, w_3 ) = kappaGradJGradJ * dGradJdGradF( m_3, u_3, v_3, w_3 ) * dGradJdF( m_3, r_3, s_3 ) + kappaGradJGradJ * gradJ( m_3 ) * d2GradJdFdGradF(m_3, r_3, s_3, u_3, v_3, w_3);
  C_GG( r_3, s_3, t_3, u_3, v_3, w_3 ) = kappaGradJGradJ * dGradJdGradF( m_3, u_3, v_3, w_3 ) * dGradJdGradF( m_3, r_3, s_3, t_3 ) + kappaGradJGradJ * gradJ( m_3 ) * d2GradJdGradFdGradF(m_3, r_3, s_3, t_3, u_3, v_3, w_3);

  FloatMatrix vC_FG( C_FG.to_voigt_form_9x27()), vC_GF;
  vC_GF.beTranspositionOf(vC_FG);

  return std::make_tuple(D_FF.to_voigt_form(), vC_FG, vC_GF, C_GG.to_voigt_form());

}

std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  // these elements should give all Tensor3 components without the third dimension in them
  //  according to the established Voigt notation, see Tensor3_3d::to_voigt_form_27()
  auto [vD_FF_3d, vC_FG_3d, vC_GF_3d, vC_GG_3d] = this->give_JacobianGradient_ConstitutiveMatrices_3d( mode, gp, tStep );
  auto vD_FF_red = vD_FF_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto vC_FG_red = vC_FG_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 5, 8, 9, 10, 14, 17 } );
  auto vC_GF_red = vC_GF_3d( { 0, 1, 5, 8, 9, 10, 14, 17 }, { 0, 1, 2, 5, 8 } );
  auto vC_GG_red = vC_GG_3d( { 0, 1, 5, 8, 9, 10, 14, 17 }, { 0, 1, 5, 8, 9, 10, 14, 17 } );
  return std::make_tuple( vD_FF_red, vC_FG_red, vC_GF_red, vC_GG_red );
}

Tensor1_3d HardMagneticMooneyRivlinCompressibleMaterial::compute_gradJ_3d( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor1_3d gradJ;
  Tensor3_3d eps;
  eps.be_Levi_Civita();

  gradJ( m_3 ) = 0.5 * eps( i_3, k_3, p_3 ) * (eps( j_3, l_3, q_3 ) * F( k_3, l_3 )) * F( p_3, q_3 ) * gradF( i_3, j_3, m_3 );

  return gradJ;
}

std::tuple<Tensor3_3d, Tensor4_3d> HardMagneticMooneyRivlinCompressibleMaterial::compute_gradJ_derivatives_3d( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor3_3d eps, dGradJdF;
  Tensor4_3d dGradJdGradF;
  eps.be_Levi_Civita();

  dGradJdF( m_3, r_3, s_3 ) = eps( i_3, r_3, p_3 ) * (eps( j_3, s_3, q_3 ) * F(p_3,q_3)) * gradF(i_3,j_3,m_3);
  dGradJdGradF( m_3, r_3, s_3, t_3 ) = 0.5 * eps( r_3, k_3, p_3 ) * (eps( s_3, l_3, q_3 ) * F( k_3, l_3 )) * (F( p_3, q_3 ) * delta( m_3, t_3 ));

  return std::make_tuple(dGradJdF,dGradJdGradF);
}

std::tuple<Tensor5_3d, Tensor6_3d, Tensor6_3d, Tensor7_3d> HardMagneticMooneyRivlinCompressibleMaterial::compute_gradJ_secondDerivatives_3d( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor3_3d eps;
  Tensor5_3d d2GradJdFdF;
  Tensor6_3d d2GradJdFdGradF, d2GradJdGradFdF;
  Tensor7_3d d2GradJdGradFdGradF;
  eps.be_Levi_Civita();

  d2GradJdFdF( m_3, r_3, s_3, u_3, v_3 ) = eps( i_3, r_3, u_3 ) * ( eps( j_3, s_3, v_3 ) * gradF( i_3, j_3, m_3 ) );
  d2GradJdFdGradF( m_3, r_3, s_3, u_3, v_3, w_3 ) = eps( u_3, r_3, p_3 ) * ( eps( v_3, s_3, q_3 ) * ( F( p_3, q_3 ) * delta( m_3, w_3 ) ) );
  d2GradJdGradFdF( m_3, r_3, s_3, t_3, u_3, v_3 ) = eps( r_3, u_3, p_3 ) * ( eps( s_3, v_3, q_3 ) * ( F( p_3, q_3 ) * delta( m_3, t_3 ) ) );
  d2GradJdGradFdGradF( m_3, r_3, s_3, t_3, u_3, v_3, w_3 ) = 0.;

  return std::make_tuple( d2GradJdFdF, d2GradJdFdGradF, d2GradJdGradFdF, d2GradJdGradFdGradF );
}

////////////////////////////////////////////////////////////////////////////////////////////////////// RotationGradient /////////////////////////////////////////////////

std::tuple<FloatArrayF<5>, FloatArrayF<8> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_RotationGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep )
{
  //get tensors
  FloatArrayF<9> vF_3d = assemble<9>( vF, { 0, 1, 2, 5, 8 } );
  FloatArrayF<27> vGradF_3d = assemble<27>( vGradF, { 0, 1, 5, 8, 9, 10, 14, 17 } );
  Tensor2_3d F(vF_3d), P;
  Tensor3_3d gradF(vGradF_3d), T;

  auto gradPhi = compute_gradPhi_PlaneStrain(F, gradF);
  auto [dGradPhidF, dGradPhidGradF] = compute_gradPhi_derivatives_PlaneStrain(F, gradF);

  P( i_3, j_3 ) = kappaGradRGradR * 2 * gradPhi( m_3 ) * dGradPhidF( m_3, i_3, j_3 );
  T( i_3, j_3, k_3 ) = kappaGradRGradR * 2 * gradPhi( m_3 ) * dGradPhidGradF( m_3, i_3, j_3, k_3 );

  // save gradients
  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  status->letTempFVectorBe( vF_3d );
  status->letTempGradFVectorBe( vGradF_3d );

  //output
  auto vP_3d = P.to_voigt_form();
  auto vT_3d = T.to_voigt_form_27();
  return std::make_tuple( vP_3d[{ 0, 1, 2, 5, 8 }], vT_3d[{ 0, 1, 5, 8, 9, 10, 14, 17 }] );
}


std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_RotationGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  // load F and GradF
  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  FloatArrayF<9> vF_3d = status->giveTempFVector();
  FloatArrayF<27> vGradF_3d = status->giveTempGradFVector();

  // get tensors
  Tensor2_3d F( vF_3d );
  Tensor3_3d gradF( vGradF_3d );
  Tensor4_3d D_FF;
  Tensor5_3d C_FG;
  Tensor6_3d C_GG;

  auto gradPhi = compute_gradPhi_PlaneStrain( F, gradF );
  auto [dGradPhidF, dGradPhidGradF] = compute_gradPhi_derivatives_PlaneStrain( F, gradF );
  auto [d2GradPhidFdF, d2GradPhidFdGradF, d2GradPhidGradFdGradF] = compute_gradPhi_secondDerivatives_PlaneStrain( F, gradF );

  D_FF( i_3, j_3, r_3, s_3 ) = kappaGradRGradR * 2 * dGradPhidF( m_3, r_3, s_3 ) * dGradPhidF( m_3, i_3, j_3 ) + kappaGradRGradR * 2 * gradPhi( m_3 ) * d2GradPhidFdF( m_3, i_3, j_3, r_3, s_3 );
  C_FG( i_3, j_3, r_3, s_3, t_3 ) = kappaGradRGradR * 2 * dGradPhidGradF( m_3, r_3, s_3, t_3 ) * dGradPhidF( m_3, i_3, j_3 ) + kappaGradRGradR * 2 * gradPhi( m_3 ) * d2GradPhidFdGradF( m_3, i_3, j_3, r_3, s_3, t_3 );
  C_GG( i_3, j_3, k_3, r_3, s_3, t_3 ) = kappaGradRGradR * 2 * dGradPhidGradF( m_3, r_3, s_3, t_3 ) * dGradPhidGradF( m_3, i_3, j_3, k_3 ) + kappaGradRGradR * 2 * gradPhi( m_3 ) * d2GradPhidGradFdGradF( m_3, i_3, j_3, k_3, r_3, s_3, t_3 );
  
  //output
  auto vD_FF_3d = D_FF.to_voigt_form();
  FloatMatrix vC_FG_3d( C_FG.to_voigt_form_9x27() ), vC_GF_3d;
  vC_GF_3d.beTranspositionOf( vC_FG_3d );
  auto vC_GG_3d = C_GG.to_voigt_form();

  auto vD_FF = vD_FF_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto vC_FG = ((FloatMatrixF<9,27>) vC_FG_3d).operator()( { 0, 1, 2, 5, 8 }, { 0, 1, 5, 8, 9, 10, 14, 17 } );
  auto vC_GF = ((FloatMatrixF<27,9>) vC_GF_3d).operator()( { 0, 1, 5, 8, 9, 10, 14, 17 }, { 0, 1, 2, 5, 8 } );
  auto vC_GG = vC_GG_3d( { 0, 1, 5, 8, 9, 10, 14, 17 }, { 0, 1, 5, 8, 9, 10, 14, 17 } );
  return std::make_tuple( vD_FF, vC_FG, vC_GF, vC_GG );
}


  

Tensor1_3d HardMagneticMooneyRivlinCompressibleMaterial::compute_gradPhi_PlaneStrain( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor1_3d gradPhi;

  auto Z = compute_Z_fraction_PlaneStrain( F );
  auto dZdF = compute_Z_fraction_derivative_PlaneStrain( F );

  gradPhi( m_3 ) = ( 1. / ( Z * Z + 1 ) ) * dZdF( p_3, q_3 ) * gradF( p_3, q_3, m_3 );

  return gradPhi;
}

std::tuple<Tensor3_3d, Tensor4_3d> HardMagneticMooneyRivlinCompressibleMaterial::compute_gradPhi_derivatives_PlaneStrain( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor3_3d dGradPhidF;
  Tensor4_3d dGradPhidGradF;

  auto Z = compute_Z_fraction_PlaneStrain( F );
  auto dZdF = compute_Z_fraction_derivative_PlaneStrain( F );
  auto d2ZdFdF = compute_Z_fraction_secondDerivative_PlaneStrain( F );

  dGradPhidF( m_3, i_3, j_3 ) = ( -2. * Z / ( Z * Z + 1. ) * ( Z * Z + 1. ) ) * dZdF( i_3, j_3 ) * dZdF( p_3, q_3 ) * gradF( p_3, q_3, m_3 ) + ( 1. / ( Z * Z + 1 ) ) * d2ZdFdF( p_3, q_3, i_3, j_3 ) * gradF( p_3, q_3, m_3 );
  dGradPhidGradF( m_3, i_3, j_3, k_3 ) = ( 1. / ( Z * Z + 1 ) ) * dZdF( i_3, j_3 ) * delta( m_3, k_3 );
 
  return std::make_tuple(dGradPhidF, dGradPhidGradF);
}

std::tuple<Tensor5_3d, Tensor6_3d, Tensor7_3d> HardMagneticMooneyRivlinCompressibleMaterial::compute_gradPhi_secondDerivatives_PlaneStrain( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor5_3d d2GradPhidFdF;
  Tensor6_3d d2GradPhidFdGradF;
  Tensor7_3d d2GradPhidGradFdGradF;

  auto Z = compute_Z_fraction_PlaneStrain( F );
  auto dZdF = compute_Z_fraction_derivative_PlaneStrain( F );
  auto d2ZdFdF = compute_Z_fraction_secondDerivative_PlaneStrain( F );
  auto d3ZdFdFdF = compute_Z_fraction_thirdDerivative_PlaneStrain(F);

  d2GradPhidFdF( m_3, i_3, j_3, r_3, s_3 ) = ( 8. * Z * Z / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * dZdF( r_3, s_3 ) * dZdF( i_3, j_3 ) * dZdF( p_3, q_3 ) * gradF( p_3, q_3, m_3 )
      - ( 2. / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * dZdF( r_3, s_3 ) * dZdF( i_3, j_3 ) * dZdF( p_3, q_3 ) * gradF( p_3, q_3, m_3 )
      - ( 2. * Z / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * (d2ZdFdF( i_3, j_3, r_3, s_3 ) * dZdF( p_3, q_3 ) * gradF( p_3, q_3, m_3 ))
      - ( 2. * Z / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * dZdF( i_3, j_3 ) * d2ZdFdF( p_3, q_3, r_3, s_3 ) * gradF( p_3, q_3, m_3 )
      - ( 2. * Z / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * dZdF( r_3, s_3 ) * d2ZdFdF( p_3, q_3, i_3, j_3 ) * gradF( p_3, q_3, m_3 )
      + ( 1. / ( Z * Z + 1 ) ) * d3ZdFdFdF( p_3, q_3, i_3, j_3, r_3, s_3 ) * gradF( p_3, q_3, m_3 );

  d2GradPhidFdGradF( m_3, i_3, j_3, r_3, s_3, t_3 ) = -( 2. * Z / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * dZdF( i_3, j_3 ) * ( dZdF( r_3, s_3 ) * delta( m_3, t_3 ) )
      + ( 1. / ( Z * Z + 1 ) ) * d2ZdFdF( r_3, s_3, i_3, j_3 ) * delta( m_3, t_3 );


  d2GradPhidGradFdGradF (m_3, i_3, j_3, k_3, r_3, s_3, t_3) = 0.;

  return std::make_tuple(d2GradPhidFdF, d2GradPhidFdGradF, d2GradPhidGradFdGradF);
}

double HardMagneticMooneyRivlinCompressibleMaterial::compute_Z_fraction_PlaneStrain( const Tensor2_3d &F ) const
{
  return (F(0,1) - F(1,0))/(F(0,0)+F(1,1));
}

Tensor2_3d HardMagneticMooneyRivlinCompressibleMaterial::compute_Z_fraction_derivative_PlaneStrain( const Tensor2_3d &F ) const
{
  Tensor2_3d dZdF, delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );

  dZdF( p_3, q_3 ) = -( F( 0, 1 ) - F( 1, 0 ) ) / (( F( 0, 0 ) + F( 1, 1 ) )*( F( 0, 0 ) + F( 1, 1 ) )) * (delta(p_3,0)*delta(q_3,0) + delta(p_3,1)*delta(q_3,1))
      + ( delta( p_3, 0 ) * delta( q_3, 1 ) - delta( p_3, 1 ) * delta( q_3, 0 ) ) / (F(0,0)+F(1,1));

  return dZdF;
}

Tensor4_3d HardMagneticMooneyRivlinCompressibleMaterial::compute_Z_fraction_secondDerivative_PlaneStrain( const Tensor2_3d &F ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor4_3d d2ZdFdF;

  d2ZdFdF( p_3, q_3, i_3, j_3 ) = - ( delta( i_3, 0 ) * delta( j_3, 1 ) - delta( i_3, 1 ) * delta( j_3, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) )
      + 2. * ( F( 0, 1 ) - F( 1, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) ) * ( delta( i_3, 0 ) * delta( j_3, 0 ) + delta( i_3, 1 ) * delta( j_3, 1 ) )
      - ( delta( p_3, 0 ) * delta( q_3, 1 ) - delta( p_3, 1 ) * delta( q_3, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( i_3, 0 ) * delta( j_3, 0 ) + delta( i_3, 1 ) * delta( j_3, 1 ) );

  return d2ZdFdF;
}

Tensor6_3d HardMagneticMooneyRivlinCompressibleMaterial::compute_Z_fraction_thirdDerivative_PlaneStrain( const Tensor2_3d &F ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor6_3d d3ZdFdFdF;

  d3ZdFdFdF( p_3, q_3, i_3, j_3, r_3, s_3 ) = 2. * ( delta( i_3, 0 ) * delta( j_3, 1 ) - delta( i_3, 1 ) * delta( j_3, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) ) * ( delta( r_3, 0 ) * delta( s_3, 0 ) + delta( r_3, 1 ) * delta( s_3, 1 ) )
                                            + 2. * ( delta( r_3, 0 ) * delta( s_3, 1 ) - delta( r_3, 1 ) * delta( s_3, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) ) * ( delta( i_3, 0 ) * delta( j_3, 0 ) + delta( i_3, 1 ) * delta( j_3, 1 ) )
                                                                                            - 6. * ( F( 0, 1 ) - F( 1, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) ) * ( delta( i_3, 0 ) * delta( j_3, 0 ) + delta( i_3, 1 ) * delta( j_3, 1 ) ) * ( delta( r_3, 0 ) * delta( s_3, 0 ) + delta( r_3, 1 ) * delta( s_3, 1 ) )
                                            + 2. * ( delta( p_3, 0 ) * delta( q_3, 1 ) - delta( p_3, 1 ) * delta( q_3, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( r_3, 0 ) * delta( s_3, 0 ) + delta( r_3, 1 ) * delta( s_3, 1 ) ) * ( delta( i_3, 0 ) * delta( j_3, 0 ) + delta( i_3, 1 ) * delta( j_3, 1 ) );

  return d3ZdFdFdF;
}


////////////////////////FBAR////////////////////////////

std::tuple<FloatArrayF<5>, FloatArrayF<5> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_Fbar_FirstPKStressVector_FPKbarStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<5> &vFbar, GaussPoint *gp, TimeStep *tStep )
{
  // these elements should give all Tensor3 components without the third dimension in them
  //  according to the established Voigt notation, see Tensor3_3d::to_voigt_form_27()
  auto [vP, vPbar] = give_Fbar_FirstPKStressVector_FPKbarStressVector_3d( assemble<9>( vF, { 0, 1, 2, 5, 8 } ), assemble<9>( vFbar, { 0, 1, 2, 5, 8 } ), gp, tStep );
  //
  return std::make_tuple( vP[{ 0, 1, 2, 5, 8 }], vPbar[{ 0, 1, 2, 5, 8 }] );
}


std::tuple<FloatArrayF<9>, FloatArrayF<9> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_Fbar_FirstPKStressVector_FPKbarStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<9> &vFbar, GaussPoint *gp, TimeStep *tStep )
{
  Tensor2_3d F( vF ), Fbar( vFbar );

  Tensor2_3d P, Pbar;

  auto [J, cofF] = F.compute_determinant_and_cofactor();
  auto [Jbar, cofFbar] = Fbar.compute_determinant_and_cofactor();

  if ( kappaJbar != 0.0 ) {
    P( r_3, s_3 ) += kappaJbar * ( J - Jbar ) * cofF( r_3, s_3 );
    Pbar( r_3, s_3 ) += -kappaJbar * ( J - Jbar ) * cofFbar( r_3, s_3 );
  }
  if ( kappaFbar != 0.0 ) {
    P( r_3, s_3 ) += kappaFbar * ( F( r_3, s_3 ) - Fbar( r_3, s_3 ) );
    Pbar( r_3, s_3 ) += -kappaFbar * ( F( r_3, s_3 ) - Fbar( r_3, s_3 ) );
  }
  // save gradients
  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  status->letTempFVectorBe( vF );
  status->letTempFbarVectorBe( vFbar );

  return std::make_tuple( P.to_voigt_form() , Pbar.to_voigt_form() );
}

std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_Fbar_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  // load F and GradF
  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  FloatArrayF<9> vF = status->giveTempFVector();
  FloatArrayF<9> vFbar = status->giveTempFbarVector();

  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. ), F( vF ), Fbar( vFbar );

  Tensor4_3d dPdF, dPdFbar, dPbardF, dPbardFbar;

  auto [J, cofF] = F.compute_determinant_and_cofactor();
  auto [Jbar, cofFbar] = Fbar.compute_determinant_and_cofactor();
  auto Fcross = F.compute_tensor_cross_product();
  auto Fbarcross = Fbar.compute_tensor_cross_product();

  if ( kappaJbar != 0.0 ) {
    dPdF( r_3, s_3, u_3, v_3 ) += kappaJbar * ( cofF( u_3, v_3 ) * cofF( r_3, s_3 ) + ( J - Jbar ) * Fcross( r_3, s_3, u_3, v_3 ) );
    dPdFbar( r_3, s_3, u_3, v_3 ) += -kappaJbar * cofFbar( u_3, v_3 ) * cofF( r_3, s_3 );
    dPbardF( r_3, s_3, u_3, v_3 ) += -kappaJbar * cofF( u_3, v_3 ) * cofFbar( r_3, s_3 );
    dPbardFbar( r_3, s_3, u_3, v_3 ) += kappaJbar * ( cofFbar( u_3, v_3 ) * cofFbar( r_3, s_3 ) - ( J - Jbar ) * Fbarcross( r_3, s_3, u_3, v_3 ) );
  }
  if ( kappaFbar != 0.0 ) {
    dPdF( r_3, s_3, u_3, v_3 ) += kappaFbar * delta( r_3, u_3 ) * delta( s_3, v_3 );
    dPdFbar( r_3, s_3, u_3, v_3 ) += -kappaFbar * delta( r_3, u_3 ) * delta( s_3, v_3 );
    dPbardF( r_3, s_3, u_3, v_3 ) += -kappaFbar * delta( r_3, u_3 ) * delta( s_3, v_3 );
    dPbardFbar( r_3, s_3, u_3, v_3 ) += kappaFbar * delta( r_3, u_3 ) * delta( s_3, v_3 );
  }
  return std::make_tuple( dPdF.to_voigt_form(), dPdFbar.to_voigt_form(), dPbardF.to_voigt_form(), dPbardFbar.to_voigt_form() );
}

std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 5>, FloatMatrixF<5, 5>, FloatMatrixF<5, 5> >
HardMagneticMooneyRivlinCompressibleMaterial ::give_Fbar_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  auto [vdPdF_3d, vdPdFbar_3d, vdPbardF_3d, vdPbardFbar_3d] = this->give_Fbar_ConstitutiveMatrices_3d( mode, gp, tStep );
  auto vdPdF_red = vdPdF_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto vdPdFbar_red = vdPdFbar_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto vdPbardF_red = vdPbardF_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto vdPbardFbar_red = vdPbardFbar_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  return std::make_tuple( vdPdF_red, vdPdFbar_red, vdPbardF_red, vdPbardFbar_red );
}




} // end namespace oofem
