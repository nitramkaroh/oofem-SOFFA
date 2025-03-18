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

#include "hardmagneticmooneyrivlincompressiblematerialprestrain.h"
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
REGISTER_Material( HardMagneticMooneyRivlinCompressibleMaterialPrestrain );

HardMagneticMooneyRivlinCompressibleMaterialPrestrain::HardMagneticMooneyRivlinCompressibleMaterialPrestrain( int n, Domain *d ) :
    HardMagneticMooneyRivlinCompressibleMaterial( n, d )
{
}


std::tuple<FloatArrayF<9>, FloatArrayF<3> >
HardMagneticMooneyRivlinCompressibleMaterialPrestrain ::give_FirstPKStressVector_MagneticInduction_3d( const FloatArrayF<9> &vF, const FloatArrayF<3> &vH, GaussPoint *gp, TimeStep *tStep ) const
{
  // prestrain
  FloatArray vF0_temp;
  gp->giveElement()->giveIPValue( vF0_temp, gp, IST_PrestrainDeformationGradient, tStep );
  FloatArrayF<9> vF0( vF0_temp );
  Tensor2_3d F0( vF0 );

  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  Tensor1_3d H_interm( vH ), B, M_full, H_full, H_mag;
  Tensor2_3d Fc( vF ), F_full, P;
  F_full( i_3, j_3 ) = Fc( i_3, k_3 ) * F0( k_3, j_3 );

  // Q = (H+M)
  double m_level = this->giveDomain()->giveFunction( m_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
  M_full( i_3 ) = m_level * M( i_3 );

  H_full( i_3 ) = F0( j_3, i_3 ) * H_interm( j_3 );

  // compute everything w.r.t. the original reference configuration
  auto [P_full, B_full] = this->computeFirstPKStressMagneticInductionTensors_3d( F_full, H_full, M_full );
  // prestrain correction
  double J0 = F0.compute_determinant();
  P( i_3, j_3 ) = 1. / J0 * P_full( i_3, k_3 ) * F0( j_3, k_3 );
  B( i_3 ) = 1. / J0 * B_full( j_3 ) * F0( i_3, j_3 );
  //
  auto vP = P.to_voigt_form();
  auto vB = B.to_voigt_form();
  // update gp
  status->letTempFVectorBe( F_full.to_voigt_form() );
  status->letTempPVectorBe( vP );
  status->letTempHVectorBe( H_full.to_voigt_form() );
  status->letTempBVectorBe( vB );
  status->letHmagTensorBe( H_mag );
  //
  return std::make_tuple( vP, vB );
}

std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 3>, FloatMatrixF<3, 9>, FloatMatrixF<3, 3> >
HardMagneticMooneyRivlinCompressibleMaterialPrestrain ::giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  MagnetoElasticMaterialStatus *status = static_cast<MagnetoElasticMaterialStatus *>( this->giveStatus( gp ) );
  FloatArrayF<9> vF( status->giveTempFVector() );
  FloatArrayF<3> vH( status->giveTempHVector() );
  Tensor1_3d H_full( vH ), Q, B, M_full, H_mag;
  Tensor2_3d F_full( vF ), P, GQQ;
  Tensor4_3d dPdF;
  Tensor3_3d dPdH, dBdF;
  Tensor2_3d dBdH;

  // prestrain
  FloatArray vF0_temp;
  gp->giveElement()->giveIPValue( vF0_temp, gp, IST_PrestrainDeformationGradient, tStep );
  FloatArrayF<9> vF0( vF0_temp );
  Tensor2_3d F0( vF0 );

  double m_level = this->giveDomain()->giveFunction( m_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
  M_full( i_3 ) = m_level * M( i_3 );
  //
  auto [dPdF_full, dPdH_full, dBdF_full, dBdH_full] = this->computeStiffnessTensors_dPdF_dBdH_dPdH_3d( F_full, H_full, M_full );
  // prestrain corrections
  double J0 = F0.compute_determinant();
  dPdF( i_3, j_3, k_3, l_3 ) = ( 1. / J0 ) * F0( j_3, p_3 ) * dPdF_full( i_3, p_3, k_3, q_3 ) * F0( l_3, q_3 );
  dBdF( i_3, k_3, l_3 ) = ( 1. / J0 ) * dBdF_full( p_3, k_3, q_3 ) * F0( i_3, p_3 ) * F0( l_3, q_3 );
  dPdH( i_3, j_3, k_3 ) = dBdF( k_3, i_3, j_3 );
  dBdH( i_3, k_3 ) = ( 1. / J0 ) * dBdH_full( p_3, q_3 ) * F0( i_3, p_3 ) * F0( k_3, q_3 );
  //
  return std::make_tuple( dPdF.to_voigt_form(), dPdH.to_voigt_form_9x3(), dBdF.to_voigt_form_3x9(), dBdH.to_matrix_form() );
}

} // end namespace oofem
