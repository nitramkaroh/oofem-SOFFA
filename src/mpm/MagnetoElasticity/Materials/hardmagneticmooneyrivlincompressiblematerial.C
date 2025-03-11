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
#include "domain.h"
#include "function.h"


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
  //
  return std::make_tuple( dPdF.to_voigt_form(), dPdH.to_voigt_form_9x3(), dBdF.to_voigt_form_3x9(), dBdH.to_matrix_form() );
}


std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 3>, FloatMatrixF<3, 9>, FloatMatrixF<3, 3> >
HardMagneticMooneyRivlinCompressibleMaterial ::compute_stiff_num( const FloatArray &vF, const FloatArray &vH, GaussPoint *gp, TimeStep *tStep )
{
  auto [P, B] = give_FirstPKStressVector_MagneticInduction_3d( vF, vH, gp, tStep );
  FloatMatrix Kuu( 9, 9 ), Kuh( 9, 3 ), Khu( 3, 9 ), Khh( 3, 3 );
  auto pert = 1.e-5;
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

std::tuple<Tensor2_3d, Tensor1_3d> HardMagneticMooneyRivlinCompressibleMaterial::computeFirstPKStressMagneticInductionTensors_3d( const Tensor2_3d F, const Tensor1_3d H, const Tensor1_3d M) const
{
  Tensor2_3d P;
  Tensor1_3d B, Q_for_P, Q_for_B;

  auto [J, cofF] = F.compute_determinant_and_cofactor();

  Q_for_P( i_3 ) = H( i_3 ) +  M( i_3 );
  Q_for_B( i_3 ) = H( i_3 ) +  M( i_3 );
  

  P( i_3, j_3 ) = C1 * this->compute_dI1_Cdev_dF( F )( i_3, j_3 ) + C2 * this->compute_dI2_Cdev_dF( F )( i_3, j_3 ) + this->compute_dVolumetricEnergy_dF( F )( i_3, j_3 )
      - ( mu_0 / J ) * F.compute_tensor_cross_product()( i_3, j_3, m_3, n_3 ) * ( cofF( m_3, q_3 ) * Q_for_P( q_3 ) * Q_for_P( n_3 ) )
      + ( mu_0 / ( 2 * J * J ) ) * cofF( m_3, n_3 ) * Q_for_P( n_3 ) * cofF( m_3, k_3 ) * Q_for_P( k_3 ) * cofF( i_3, j_3 );

  B( i_3 ) = mu_0 / J * cofF( j_3, i_3 ) * cofF( j_3, k_3 ) * Q_for_B( k_3 );

  return std::make_tuple( P, B );
}

std::tuple<Tensor4_3d, Tensor3_3d, Tensor3_3d, Tensor2_3d> HardMagneticMooneyRivlinCompressibleMaterial::computeStiffnessTensors_dPdF_dBdH_dPdH_3d( const Tensor2_3d F, const Tensor1_3d H, const Tensor1_3d M ) const
{
  Tensor4_3d dPdF;
  Tensor3_3d dPdH, dBdF, dPdH_in_H;
  Tensor2_3d dBdH, GQQ;
  Tensor1_3d Q_for_P, Q_for_B;

  auto [J, cofF] = F.compute_determinant_and_cofactor();
  auto Fcross = F.compute_tensor_cross_product();

  Q_for_P( i_3 ) = H( i_3 ) + M( i_3 );
  Q_for_B( i_3 ) = H( i_3 ) + M( i_3 );

  GQQ( i_3, j_3 ) = cofF( i_3, k_3 ) * Q_for_P (k_3 ) * Q_for_P( j_3 );
  auto GQQcross = GQQ.compute_tensor_cross_product();

  dPdF( i_3, j_3, k_3, l_3 ) = C1 * this->compute_d2I1_Cdev_dF2( F )( i_3, j_3, k_3, l_3 )
      + C2 * this->compute_d2I2_Cdev_dF2( F )( i_3, j_3, k_3, l_3 )
      + this->compute_d2VolumetricEnergy_dF2( F )( i_3, j_3, k_3, l_3 )
      - mu_0 / J * ( GQQcross( i_3, j_3, k_3, l_3 ) + ( Fcross( i_3, j_3, m_3, n_3 ) * Q_for_P( n_3 ) ) * ( Fcross( m_3, q_3, k_3, l_3 ) * Q_for_P( q_3 ) ) )
      + mu_0 / J / J * Fcross( i_3, j_3, m_3, n_3 ) * GQQ( m_3, n_3 ) * cofF( k_3, l_3 )
      + ( mu_0 / ( 2 * J * J ) ) * ( cofF( m_3, n_3 ) * Q_for_P( n_3 ) * cofF( m_3, k_3 ) * Q_for_P( k_3 ) * Fcross( i_3, j_3, k_3, l_3 ) )
      - ( mu_0 / ( J * J * J ) ) * ( cofF( m_3, n_3 ) * Q_for_P( n_3 ) * cofF( m_3, k_3 ) * Q_for_P( k_3 ) * cofF( i_3, j_3 ) * cofF( k_3, l_3 ) )
      + ( mu_0 / ( J * J ) ) * cofF( i_3, j_3 ) * ( Fcross( k_3, l_3, m_3, n_3 ) * cofF( m_3, q_3 ) * Q_for_P( q_3 ) * Q_for_P( n_3 ) );

  dBdH( i_3, k_3 ) = -mu_0 / J * cofF( j_3, i_3 ) * cofF( j_3, k_3 );
  //
  // dBdF( i_3, k_3, l_3 ) = mu_0 / J * ( G( j_3, m_3 ) * Q( m_3 ) * Fcross( j_3, i_3, k_3, l_3 ) ) + G( j_3, i_3 ) * ( Q( m_3 ) * Fcross( j_3, m_3, k_3, l_3 ) ) - mu_0 / ( J * J ) * ( G( j_3, i_3 ) * G( j_3, m_3 ) * Q( m_3 ) * G( k_3, l_3 ) );
  //
  dPdH( i_3, j_3, k_3 ) = mu_0 / J * ( Fcross( i_3, j_3, m_3, n_3 ) * Q_for_P( n_3 ) * cofF( m_3, k_3 ) + Fcross( i_3, j_3, m_3, k_3 ) * cofF( m_3, q_3 ) * Q_for_P( q_3 ) ) - mu_0 / J / J * ( cofF( m_3, k_3 ) * cofF( m_3, o_3 ) * Q_for_P( o_3 ) ) * cofF( i_3, j_3 );
  dPdH_in_H( i_3, j_3, k_3 ) = mu_0 / J * ( Fcross( i_3, j_3, m_3, n_3 ) * Q_for_B( n_3 ) * cofF( m_3, k_3 ) + Fcross( i_3, j_3, m_3, k_3 ) * cofF( m_3, q_3 ) * Q_for_P( q_3 ) ) - mu_0 / J / J * ( cofF( m_3, k_3 ) * cofF( m_3, o_3 ) * Q_for_B( o_3 ) ) * cofF( i_3, j_3 );
  //
  dBdF( i_3, k_3, l_3 ) = dPdH_in_H( k_3, l_3, i_3 );

  return std::make_tuple( dPdF, dPdH, dBdF, dBdH);
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
  FloatArray m;
  IR_GIVE_FIELD( ir, m, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_magnetization );
  IR_GIVE_FIELD( ir, m_ltf, _IFT_HardMagneticMooneyRivlinCompressibleMaterial_m_ltf );
  if ( m.giveSize() != 3 ) {
    OOFEM_ERROR( "Magnetization has to be vector of size 3" );
  } else {
    M = Tensor1_3d( FloatArrayF<3>( m ) );
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
    h( i_3 ) = J * G( i_3, j_3 ) * H( j_3 );
    answer = h.to_voigt_form();
    return 1;
  } else {
    return Material::giveIPValue( answer, gp, type, tStep );
  }
}


} // end namespace oofem
