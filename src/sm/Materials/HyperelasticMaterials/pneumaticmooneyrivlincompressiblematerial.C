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

#include "pneumaticmooneyrivlincompressiblematerial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "domain.h"
#include "function.h"


namespace oofem {
REGISTER_Material( PneumaticMooneyRivlinCompressibleMaterial );

PneumaticMooneyRivlinCompressibleMaterial::PneumaticMooneyRivlinCompressibleMaterial( int n, Domain *d ) :
    MooneyRivlinCompressibleMaterial( n, d )
{
}

FloatArrayF<9>
PneumaticMooneyRivlinCompressibleMaterial::giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
  // mooney rivlin solution
  auto vP = MooneyRivlinCompressibleMaterial::giveFirstPKStressVector_3d( vF, gp, tStep );

  // necessary data
  StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
  Tensor2_3d F( vF ), P_pneum;
  auto cofF = F.compute_cofactor();
  double pressure_level = this->giveDomain()->giveFunction( pressure_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
  double dp = pressure * pressure_level;

  // compute pneumatic term
  P_pneum( i_3, j_3 ) = dp * cofF( i_3, j_3 );

  // add together
  vP += P_pneum.to_voigt_form();

  // update gp
  status->letTempFVectorBe( vF );
  status->letTempPVectorBe( vP );

  //
  return vP;
}


FloatMatrixF<9, 9>
PneumaticMooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dPdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
  // mooney rivlin solution
  auto vA = MooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dPdF( mode, gp, tStep );

  // get data
  StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
  FloatArrayF<9> vF = status->giveTempFVector();
  Tensor2_3d F( vF );
  Tensor4_3d A_pneum;
  auto Fcross = F.compute_tensor_cross_product();
  double pressure_level = this->giveDomain()->giveFunction( pressure_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
  double dp = pressure * pressure_level;

  // compute pneumatic contribution
  A_pneum(i_3, j_3, k_3, l_3) = dp * Fcross(i_3, j_3, k_3, l_3);

  //add together and return
  vA += A_pneum.to_voigt_form();
  return vA;
}

void PneumaticMooneyRivlinCompressibleMaterial::initializeFrom( InputRecord &ir )
{
  MooneyRivlinCompressibleMaterial::initializeFrom( ir );
  IR_GIVE_FIELD( ir, pressure, _IFT_PneumaticMooneyRivlinCompressibleMaterial_pressure );
  IR_GIVE_FIELD( ir, pressure_ltf, _IFT_PneumaticMooneyRivlinCompressibleMaterial_pressure_ltf );
}
} // end namespace oofem
