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

#include "mooneyrivlincompressiblematerial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"



namespace oofem {
REGISTER_Material(MooneyRivlinCompressibleMaterial);

  MooneyRivlinCompressibleMaterial::MooneyRivlinCompressibleMaterial(int n, Domain *d) : StructuralMaterial(n, d), BaseHyperElasticMaterial(), LargeStrainMixedPressureMaterialExtensionInterface(d)
{ }

FloatArrayF< 9 >
MooneyRivlinCompressibleMaterial::giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    Tensor2_3d F(vF), P;
    // compute the first Piola-Kirchhoff
    P(i_3, j_3) =  C1 * this->compute_dI1_Cdev_dF(F)(i_3, j_3) + C2 * this->compute_dI2_Cdev_dF(F)(i_3, j_3) + this->compute_dVolumetricEnergy_dF(F)(i_3, j_3);
    auto vP = P.to_voigt_form();
    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(vP);
    //
    return vP;
}



FloatMatrixF< 9, 9 >
MooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dPdF(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    FloatArrayF< 9 >vF(status->giveTempFVector() );
    Tensor2_3d F(vF);
    Tensor4_3d A;
    A(i_3, j_3, k_3, l_3) = C1 * this->compute_d2I1_Cdev_dF2(F)(i_3, j_3, k_3, l_3) + C2 * this->compute_d2I2_Cdev_dF2(F)(i_3, j_3, k_3, l_3) + this->compute_d2VolumetricEnergy_dF2(F)(i_3, j_3, k_3, l_3);
    return A.to_voigt_form();
}




MaterialStatus *
MooneyRivlinCompressibleMaterial::CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(gp);
}


void
MooneyRivlinCompressibleMaterial::initializeFrom(InputRecord &ir)
{
    StructuralMaterial::initializeFrom(ir);
    BaseHyperElasticMaterial::initializeFrom(ir);
    IR_GIVE_FIELD(ir, C1, _IFT_MooneyRivlinCompressibleMaterial_c1);
    IR_GIVE_FIELD(ir, C2, _IFT_MooneyRivlinCompressibleMaterial_c2);
}

  /////////////////////////// Functions for Mixed Pressure formulation

std::tuple<FloatArrayF<9>, double>
MooneyRivlinCompressibleMaterial ::  giveFirstPKStressVector_3d(const FloatArrayF<9> &vF, double pressure, GaussPoint *gp, TimeStep *tStep)
{

  StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
  //
  Tensor2_3d F(vF), P;
  //
  auto [ detF, cofF ] = F.compute_determinant_and_cofactor();
  //
  P(i_3, j_3) =  C1 * this->compute_dI1_Cdev_dF(F)(i_3, j_3) + C2 * this->compute_dI2_Cdev_dF(F)(i_3, j_3) + pressure * cofF(i_3, j_3);
  auto vP = P.to_voigt_form();
  // update gp
  status->letTempFVectorBe(vF);
  status->letTempPVectorBe(vP);
  //
  auto dJ = (detF - 1.) - pressure/this->K;
  return std::make_tuple(vP, dJ);
}

std::tuple<FloatArrayF<5>, double>
MooneyRivlinCompressibleMaterial ::  giveFirstPKStressVector_PlaneStrain(const FloatArrayF<5> &vF, double pressure, GaussPoint *gp, TimeStep *tStep)
{
  auto [vP, dJ] = this->giveFirstPKStressVector_3d(assemble< 9 >(vF, { 0, 1, 2, 5, 8 }), pressure, gp, tStep);
  //
  return std::make_tuple(vP[ { 0, 1, 2, 5, 8 } ], dJ);
}

  

  

std::tuple<FloatMatrixF<9,9>, FloatArrayF<9>, double>
MooneyRivlinCompressibleMaterial :: giveLargeStrainMixedPressureConstitutiveMatrices_3d(double pressure, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
  FloatArrayF< 9 >vF(status->giveTempFVector() );
  //  double pressure = status->giveTempPressure();
  Tensor2_3d F(vF);
  Tensor4_3d A;
  A(i_3, j_3, k_3, l_3) = C1 * this->compute_d2I1_Cdev_dF2(F)(i_3, j_3, k_3, l_3) + C2 * this->compute_d2I2_Cdev_dF2(F)(i_3, j_3, k_3, l_3) + pressure * this->compute_d2J_dF2(F)(i_3, j_3, k_3, l_3);
;
  //
  auto [ detF, cofF ] = F.compute_determinant_and_cofactor();
  //  
  return std::make_tuple(A.to_voigt_form(), cofF.to_voigt_form(), -this->K);
  
}


std::tuple<FloatMatrixF<5,5>, FloatArrayF<5>, double>
MooneyRivlinCompressibleMaterial :: giveLargeStrainMixedPressureConstitutiveMatrices_PlaneStrain(double pressure, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  auto [A, cofF, minusK] = this->giveLargeStrainMixedPressureConstitutiveMatrices_3d(pressure, rMode, gp, tStep);
  //
  return std::make_tuple(A({ 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 }), cofF[ { 0, 1, 2, 5, 8 } ], minusK);
  
}

  

  
} // end namespace oofem
