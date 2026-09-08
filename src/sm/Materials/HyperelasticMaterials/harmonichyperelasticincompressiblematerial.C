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

#include "harmonichyperelasticincompressiblematerial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"


namespace oofem {
REGISTER_Material( HarmonicIncompressibleHyperelasticMaterial );

HarmonicIncompressibleHyperelasticMaterial::HarmonicIncompressibleHyperelasticMaterial( int n, Domain *d ) :
    StructuralMaterial( n, d ), BaseHyperElasticMaterial()
{
}

FloatArrayF<5>
HarmonicIncompressibleHyperelasticMaterial::giveFirstPKStressVector_PlaneStrain( const FloatArrayF<5> &vF_ps, GaussPoint *gp, TimeStep *tStep ) const
// returns 5 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
  StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

  // from plane strain to 3d
  auto vF_3d = assemble<9>( vF_ps, { 0, 1, 2, 5, 8 } );

  Tensor2_3d F( vF_3d ), P;
  double dfdI = this->compute_dFunctiondI_from_F(F);
  Tensor2_3d dIdF = this->compute_dI_of_U_2d_dF_from_F(F);
  auto [ J, cofF ] = F.compute_determinant_and_cofactor();

  // compute the first Piola-Kirchhoff: harmonic core + volumetric penalty
  P( i_3, j_3 ) = 2 * mu * ( dfdI * dIdF( i_3, j_3 ) - cofF( i_3, j_3 ) )
      + lam * ( J - 1. ) * cofF( i_3, j_3 );

  auto vP_3d = P.to_voigt_form();
  // update gp
  status->letTempFVectorBe( vF_3d );
  status->letTempPVectorBe( vP_3d );
  //
  auto vP_red = vP_3d[ { 0, 1, 2, 5, 8 } ];
  return vP_red;
}


FloatMatrixF<5, 5>
HarmonicIncompressibleHyperelasticMaterial::givePlaneStrainStiffnessMatrix_dPdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const
// returns the 5x5 tangent stiffness matrix - dP/dF
{
  StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

  FloatArrayF<9> vF_3d( status->giveTempFVector() );

  Tensor2_3d F( vF_3d );
  Tensor4_3d A;
  auto [ J, cofF ] = F.compute_determinant_and_cofactor();
  auto Fcross = F.compute_tensor_cross_product();
  double dfdI = this->compute_dFunctiondI_from_F(F);
  double d2fdI2 = this->compute_d2FunctiondIdI_from_F(F);
  Tensor2_3d dIdF = this->compute_dI_of_U_2d_dF_from_F(F);
  Tensor4_3d d2IdFdF = this->compute_d2I_of_U_2d_dFdF_from_F(F);

  // harmonic core + volumetric penalty
  // d/dF [ lam (J-1) cofF ] = lam ( cofF x cofF + (J-1) Fcross )
  A( i_3, j_3, k_3, l_3 ) = 2. * mu * ( d2fdI2 * dIdF( k_3, l_3 ) * dIdF( i_3, j_3 ) + dfdI * d2IdFdF( i_3, j_3, k_3, l_3 ) - Fcross( i_3, j_3, k_3, l_3 ) )
      + lam * ( cofF( i_3, j_3 ) * cofF( k_3, l_3 ) + ( J - 1. ) * Fcross( i_3, j_3, k_3, l_3 ) );

  //back to plane strain
  auto vA_3d = A.to_voigt_form();
  auto vA_red = vA_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  return vA_red;
}


MaterialStatus *
HarmonicIncompressibleHyperelasticMaterial::CreateStatus( GaussPoint *gp ) const
{
  return new StructuralMaterialStatus( gp );
}


void HarmonicIncompressibleHyperelasticMaterial::initializeFrom( InputRecord &ir )
{
  StructuralMaterial::initializeFrom( ir );
  BaseHyperElasticMaterial::initializeFrom( ir );
  IR_GIVE_FIELD( ir, alpha, _IFT_HarmonicIncompressibleHyperelasticMaterial_alpha );
  IR_GIVE_FIELD( ir, beta, _IFT_HarmonicIncompressibleHyperelasticMaterial_beta );
  IR_GIVE_FIELD( ir, mu, _IFT_HarmonicIncompressibleHyperelasticMaterial_mu );
  IR_GIVE_FIELD( ir, lam, _IFT_HarmonicIncompressibleHyperelasticMaterial_lam );
}

double HarmonicIncompressibleHyperelasticMaterial::compute_dFunctiondI_from_F( const Tensor2_3d &F ) const
{
  double I = this->compute_I_of_U_2d_from_F( F );
  return 1. / ( 4. * alpha ) * ( I + std::sqrt( I * I - 16 * alpha * beta ) );
}

double HarmonicIncompressibleHyperelasticMaterial::compute_d2FunctiondIdI_from_F( const Tensor2_3d &F ) const
{
  double I = this->compute_I_of_U_2d_from_F( F );
  return 1. / ( 4. * alpha ) * ( 1. + I / std::sqrt( I * I - 16 * alpha * beta ) );
}


double HarmonicIncompressibleHyperelasticMaterial::compute_I_of_U_2d_from_F( const Tensor2_3d &F ) const
{
  return std::sqrt( F( p_3, q_3 ) * F( p_3, q_3 ) - 1. + 2. * F.compute_determinant() );
}
Tensor2_3d HarmonicIncompressibleHyperelasticMaterial::compute_dI_of_U_2d_dF_from_F( const Tensor2_3d &F ) const
{
  Tensor2_3d answer;
  auto [J, cofF] = F.compute_determinant_and_cofactor();

  answer( i_3, j_3 ) = ( F( i_3, j_3 ) + cofF( i_3, j_3 ) ) / std::sqrt( F( p_3, q_3 ) * F( p_3, q_3 ) - 1. + 2. * J );
  return answer;
}

Tensor4_3d HarmonicIncompressibleHyperelasticMaterial::compute_d2I_of_U_2d_dFdF_from_F( const Tensor2_3d &F ) const
{
  Tensor4_3d answer;
  auto cofF = F.compute_cofactor();
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  double I = this->compute_I_of_U_2d_from_F( F ); // the square root which repeats again and again
  Tensor4_3d Fcross = F.compute_tensor_cross_product();

  answer( i_3, j_3, k_3, l_3 ) = -( F( i_3, j_3 ) + cofF( i_3, j_3 ) ) * ( F( k_3, l_3 ) + cofF( k_3, l_3 ) ) * ( 1. / ( I * I * I ) )
      + ( delta( i_3, k_3 ) * delta( j_3, l_3 ) + Fcross( i_3, j_3, k_3, l_3 ) ) * ( 1. / I );

  return answer;
}


} // end namespace oofem
