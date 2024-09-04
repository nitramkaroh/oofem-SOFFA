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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#include "mpm/CrossSections/magnetoelasticcrosssection.h"
#include "mpm/Materials/magnetoelasticmaterial.h"
#include "mpm/Materials/magnetoelasticms.h"
#include "gausspoint.h"
#include "element.h"
#include "floatarray.h"

namespace oofem {


std::tuple<FloatArray, FloatArray>
MagnetoElasticCrossSection::giveFluxVectors(const FloatArray &vGrad, GaussPoint *gp, TimeStep *tStep) const
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
      auto vF = vGrad.at(0);
      auto vH = vGrad.at(1);
      return this->giveFluxVectors_P_B_3d(vF,vH, gp, tStep);
    } else if ( mode == _PlaneStrain ) {
      auto vF = vGrad.at(0);
      auto vH = vGrad.at(1);
      return this->giveFluxVectors_P_B_PlaneStrain(vF, vH, gp, tStep);
    } else {
        OOFEM_ERROR("unsupported mode");
    }
}


std::tuple<FloatArrayF< 9 >,FloatArrayF< 3 >>
MagnetoElasticCrossSection :: giveFluxVectors_P_B_3d(const FloatArrayF< 9 > &vF, const FloatArrayF< 3 > &vH, GaussPoint *gp, TimeStep *tStep) const
{
  auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
  return mat->give_1PKStressVector_MagneticInduction_3d(vF, vH, gp, tStep);
}

  
std::tuple<FloatArrayF< 5 >,FloatArrayF< 2 >>
SimpleCrossSection :: giveFluxVectors_P_B_3d(const FloatArrayF< 5 > &vF, const FloatArrayF< 2 > &vH, GaussPoint *gp, TimeStep *tStep) const
{
  auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
  return mat->give_1PKStressVector_MagneticInduction_PlameStrain(vF, vH, gp, tStep);
}
 
  
void
MagnetoElasticCrossSection :: giveCharMaterialMatrix_dPdF_dBdH_dPdH(std::vector<FloatMatrix> &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode mode = gp->giveMaterialMode();
  if ( mode == _3dMat ) {
    auto [dPdF, dBdH, dPdH] = give_dPdF_dBdH_dPdH_3d(mode, gp, tStep);
    answer.emplace_back(dPdF);
    answer.emplace_back(dBdH);
    answer.emplace_back(dPdH);
  } else if ( mode == _PlaneStrain ) {
    auto [dPdF, dBdH, dPdH] = give_dPdF_dBdH_dPdH_PlaneStrain(mode, gp, tStep);
    answer.emplace_back(dPdF);
    answer.emplace_back(dBdH);
    answer.emplace_back(dPdH);
  } else {
    OOFEM_ERROR("unsupported mode");
  }
}

  st::tuple<FloatMatrixF< 9, 9 >,FloatMatrixF< 3, 3 >,FloatMatrixF< 9, 3 >>
MagnetoElasticCrossSection :: giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
    return mat->giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d(rMode, gp, tStep);
}


FloatMatrixF< 5, 5 >
MagnetoElasticCrossSection :: giveConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
    return mat->giveConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain(rMode, gp, tStep);
}



  



} // end namespace oofem
