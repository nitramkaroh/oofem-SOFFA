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

#include "mpm/MagnetoElasticity/CrossSections/magnetoelasticcrosssection.h"
#include "mpm/MagnetoElasticity/Materials/magnetoelasticmaterial.h"
#include "mpm/MagnetoElasticity/Materials/magnetoelasticmaterialstatus.h"
#include "gausspoint.h"
#include "element.h"
#include "floatarray.h"
#include "classfactory.h"

namespace oofem {
REGISTER_CrossSection(MagnetoElasticCrossSection);

std::vector<FloatArray>
MagnetoElasticCrossSection::give_FluxVectors(const std::vector<FloatArray> &vGrad, GaussPoint *gp, TimeStep *tStep) const
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode == _3dMat ) {
      auto vF = vGrad.at(0);
      auto vH = vGrad.at(1);
      auto[P, B] = this->give_FluxVectors_P_B_3d(vF,vH, gp, tStep);
      std::vector<FloatArray> answer = {P, B};
      return answer;
      
    } else if ( mode == _PlaneStrain ) {
      auto vF = vGrad.at(0);
      auto vH = vGrad.at(1);
      auto[P, B] = this->give_FluxVectors_P_B_PlaneStrain(vF,vH, gp, tStep);
      std::vector<FloatArray> answer = {P, B};
      return answer;
    } else {
        OOFEM_ERROR("unsupported mode");
    }
}


std::tuple<FloatArrayF< 9 >,FloatArrayF< 3 >>
MagnetoElasticCrossSection :: give_FluxVectors_P_B_3d(const FloatArrayF< 9 > &vF, const FloatArrayF< 3 > &vH, GaussPoint *gp, TimeStep *tStep) const
{
  auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
  return mat->give_FirstPKStressVector_MagneticInduction_3d(vF, vH, gp, tStep);
}

  
std::tuple<FloatArrayF< 5 >,FloatArrayF< 2 >>
MagnetoElasticCrossSection :: give_FluxVectors_P_B_PlaneStrain(const FloatArrayF< 5 > &vF, const FloatArrayF< 2 > &vH, GaussPoint *gp, TimeStep *tStep) const
{
  auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
  return mat->give_FirstPKStressVector_MagneticInduction_PlaneStrain(vF, vH, gp, tStep);
}
 
  
void
MagnetoElasticCrossSection :: give_dFlux_dGrad_Matrices(std::vector<FloatMatrix> &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode mode = gp->giveMaterialMode();
  if ( mode == _3dMat ) {
    auto [dPdF, dPdH, dBdF, dBdH] = this->give_ConstitutiveMatrices_dPdF_dBdH_dPdH_3d(rmode, gp, tStep);
    answer.emplace_back(dPdF);
    answer.emplace_back(dPdH);
    answer.emplace_back(dBdF);
    answer.emplace_back(dBdH);
  } else if ( mode == _PlaneStrain ) {
    auto [dPdF, dPdH, dBdF, dBdH] = this->give_ConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain(rmode, gp, tStep);
    answer.emplace_back(dPdF);
    answer.emplace_back(dPdH);
    answer.emplace_back(dBdF);
    answer.emplace_back(dBdH);
  } else {
    OOFEM_ERROR("unsupported mode");
  }
}

  std::tuple<FloatMatrixF< 9, 9 >,FloatMatrixF< 9, 3 >, FloatMatrixF< 3, 9 >, FloatMatrixF< 3, 3 >>
MagnetoElasticCrossSection :: give_ConstitutiveMatrices_dPdF_dBdH_dPdH_3d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
    return mat->giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d(rMode, gp, tStep);
}


  std::tuple<FloatMatrixF< 5, 5 >, FloatMatrixF< 5, 2 >, FloatMatrixF< 2, 5 >, FloatMatrixF< 2, 2 >>
MagnetoElasticCrossSection :: give_ConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
    return mat->giveConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain(rMode, gp, tStep);
}


void
MagnetoElasticCrossSection :: initializeFrom(InputRecord &ir)
{
    CrossSection::initializeFrom(ir);
    this->materialNumber = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->materialNumber, _IFT_MagnetoElasticCrossSection_MaterialNumber);


}

void
MagnetoElasticCrossSection :: createMaterialStatus(GaussPoint &iGP)
{
    Material *mat = domain->giveMaterial(materialNumber);
    MaterialStatus *matStat = mat->CreateStatus(& iGP);
    iGP.setMaterialStatus(matStat);
}



int
MagnetoElasticCrossSection :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveMaterial(gp)->packUnknowns(buff, tStep, gp);
}

int
MagnetoElasticCrossSection :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveMaterial(gp)->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
MagnetoElasticCrossSection :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->giveMaterial(gp)->estimatePackSize(buff, gp);
}

Material *
MagnetoElasticCrossSection :: giveMaterial(IntegrationPoint *ip) const
{
    if ( this->giveMaterialNumber() ) {
        return this->giveDomain()->giveMaterial(this->giveMaterialNumber() );
    } else {
        return ip->giveElement()->giveMaterial();
    }
}

  ////////////////////////////////////////////////////////////////////////////////////////////////
FloatArray
MagnetoElasticCrossSection::give_SecondGradient_FluxVector(const FloatArray &vGrad, GaussPoint *gp, TimeStep *tStep) const
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode == _3dMat ) {
      return this->give_SecondGradient_FluxVector_3d(vGrad, gp, tStep);
    } else if ( mode == _PlaneStrain ) {
      return this->give_SecondGradient_FluxVector_PlaneStrain(vGrad, gp, tStep);
    } else {
        OOFEM_ERROR("unsupported mode");
    }
}


FloatArrayF< 27 >
MagnetoElasticCrossSection :: give_SecondGradient_FluxVector_3d(const FloatArrayF< 27 > &vG, GaussPoint *gp, TimeStep *tStep) const
{
  auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
  return mat->give_SecondGradient_FirstPKStressVector_3d(vG, gp, tStep);
}

  
FloatArrayF< 8 >
MagnetoElasticCrossSection :: give_SecondGradient_FluxVector_PlaneStrain(const FloatArrayF< 8 > &vG, GaussPoint *gp, TimeStep *tStep) const
{
  auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
  return mat->give_SecondGradient_FirstPKStressVector_PlaneStrain(vG, gp, tStep);
}
 
  
void
MagnetoElasticCrossSection :: give_SecondGradient_dFlux_dGrad(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode mode = gp->giveMaterialMode();
  if ( mode == _3dMat ) {
    answer = this->give_SecondGradient_ConstitutiveMatrix_3d(rmode, gp, tStep);
  } else if ( mode == _PlaneStrain ) {
    answer = this->give_SecondGradient_ConstitutiveMatrix_PlaneStrain(rmode, gp, tStep);
  } else {
    OOFEM_ERROR("unsupported mode");
  }
}


FloatMatrixF< 1, 1 >
MagnetoElasticCrossSection :: give_SecondGradient_ConstitutiveMatrix_3d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
    return mat->give_SecondGradient_ConstitutiveMatrix_3d(rMode, gp, tStep);
}


FloatMatrixF< 1, 1 >
MagnetoElasticCrossSection :: give_SecondGradient_ConstitutiveMatrix_PlaneStrain(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< MagnetoElasticMaterial * >( this->giveMaterial(gp) );
    return mat->give_SecondGradient_ConstitutiveMatrix_PlaneStrain(rMode, gp, tStep);
}










  

} // end namespace oofem
