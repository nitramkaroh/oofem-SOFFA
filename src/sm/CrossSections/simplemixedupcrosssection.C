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

#include "sm/CrossSections/simplemixedupcrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Elements/structuralelement.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "datastream.h"
#include "contextioerr.h"
#include "engngm.h"
#include "../src/sm/Materials/MixedPressure/mixedpressurematerialextensioninterface.h"


namespace oofem {
REGISTER_CrossSection(SimpleMixedUPCrossSection);



std::tuple<FloatArray, double>
SimpleMixedUPCrossSection::giveRealStresses(const FloatArray &strain, double p, GaussPoint *gp, TimeStep *tStep) const
{
  auto mode = gp->giveMaterialMode();
  if ( mode == _3dMat ) {
      return this->giveRealStress_3d(strain, p, gp, tStep);
    } else if ( mode == _PlaneStrain ) {
      return this->giveRealStress_PlaneStrain(strain, p, gp, tStep);
    } else {
      OOFEM_ERROR("unsupported mode");
    }
}
  
std::tuple<FloatArrayF< 6 >, double>
SimpleMixedUPCrossSection::giveRealStress_3d(const FloatArrayF< 6 > &strain, double p, GaussPoint *gp, TimeStep *tStep) const
{
    SmallStrainMixedPressureMaterialExtensionInterface *mixedPressureMat = static_cast< SmallStrainMixedPressureMaterialExtensionInterface * >( this->giveMaterialInterface(SmallStrainMixedPressureMaterialExtensionInterfaceType, gp) );
    if ( !mixedPressureMat ) {
        OOFEM_ERROR("Material doesn't implement the required Mixed Pressure interface!");
    }

    return mixedPressureMat->giveRealStressVector_3d(strain, p, gp, tStep);
}


  std::tuple<FloatArrayF< 4 >, double>
SimpleMixedUPCrossSection::giveRealStress_PlaneStrain(const FloatArrayF< 4 > &strain, double p, GaussPoint *gp, TimeStep *tStep) const
{
      SmallStrainMixedPressureMaterialExtensionInterface *mixedPressureMat = static_cast<  SmallStrainMixedPressureMaterialExtensionInterface * >( this->giveMaterialInterface(SmallStrainMixedPressureMaterialExtensionInterfaceType, gp) );
    if ( !mixedPressureMat ) {
        OOFEM_ERROR("Material doesn't implement the required Mixed Pressure interface!");
    }

    return mixedPressureMat->giveRealStressVector_PlaneStrain(strain, p, gp, tStep);
}


std::tuple<FloatMatrixF<6,6>, FloatArrayF<6>, double>
SimpleMixedUPCrossSection::giveSmallStrainMixedPressureConstitutiveMatrices_3d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
      SmallStrainMixedPressureMaterialExtensionInterface *mixedPressureMat = static_cast<  SmallStrainMixedPressureMaterialExtensionInterface * >( this->giveMaterialInterface(SmallStrainMixedPressureMaterialExtensionInterfaceType, gp) );
    if ( !mixedPressureMat ) {
        OOFEM_ERROR("Material doesn't implement the required Mixed Pressure interface!");
    }
    return mixedPressureMat->giveSmallStrainMixedPressureConstitutiveMatrices_3d(rMode, gp, tStep);
}


std::tuple<FloatMatrixF<4,4>, FloatArrayF<4>, double>
SimpleMixedUPCrossSection::giveSmallStrainMixedPressureConstitutiveMatrices_PlaneStrain(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    SmallStrainMixedPressureMaterialExtensionInterface *mixedPressureMat = static_cast<  SmallStrainMixedPressureMaterialExtensionInterface * >( this->giveMaterialInterface(SmallStrainMixedPressureMaterialExtensionInterfaceType, gp) );
    if ( !mixedPressureMat ) {
        OOFEM_ERROR("Material doesn't implement the required Mixed Pressure interface!");
    }

    return mixedPressureMat->giveSmallStrainMixedPressureConstitutiveMatrices_PlaneStrain(rMode, gp, tStep);
}

/// Large strain
std::tuple<FloatArray, double>
SimpleMixedUPCrossSection::giveFirstPKStresses(const FloatArray &reducedF, double p, GaussPoint *gp, TimeStep *tStep) const
{
    auto mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
      return this->giveFirstPKStress_3d(reducedF, p, gp, tStep);
    } else if ( mode == _PlaneStrain ) {
      return this->giveFirstPKStress_PlaneStrain(reducedF, p, gp, tStep);
    } else {
        OOFEM_ERROR("unsupported mode");
    }
}


std::tuple<FloatArrayF< 9 >, double >
SimpleMixedUPCrossSection::giveFirstPKStress_3d(const FloatArrayF< 9 > &vF, double p, GaussPoint *gp, TimeStep *tStep) const
{
    LargeStrainMixedPressureMaterialExtensionInterface *mixedPressureMat = static_cast<  LargeStrainMixedPressureMaterialExtensionInterface * >( this->giveMaterialInterface(LargeStrainMixedPressureMaterialExtensionInterfaceType, gp) );
    if ( !mixedPressureMat ) {
        OOFEM_ERROR("Material doesn't implement the required Mixed Pressure interface!");
    }

    return mixedPressureMat->giveFirstPKStressVector_3d(vF, p, gp, tStep);
}


std::tuple< FloatArrayF< 5 >, double >
SimpleMixedUPCrossSection::giveFirstPKStress_PlaneStrain(const FloatArrayF< 5 > &vF, double p, GaussPoint *gp, TimeStep *tStep) const
{
  
      LargeStrainMixedPressureMaterialExtensionInterface *mixedPressureMat = static_cast<  LargeStrainMixedPressureMaterialExtensionInterface * >( this->giveMaterialInterface(LargeStrainMixedPressureMaterialExtensionInterfaceType, gp) );
    if ( !mixedPressureMat ) {
        OOFEM_ERROR("Material doesn't implement the required Mixed Pressure interface!");
    }
    return mixedPressureMat->giveFirstPKStressVector_PlaneStrain(vF, p, gp, tStep);
}



void
SimpleMixedUPCrossSection :: giveSmallStrainStiffnessMatrix(FloatMatrix &answer_uu, FloatArray &answer_up, double& K, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
      std::tie(answer_uu, answer_up, K) = this->giveSmallStrainMixedPressureConstitutiveMatrices_3d(rMode, gp, tStep);
    }  else if ( mode == _PlaneStrain ) {
       std::tie(answer_uu, answer_up, K) = this->giveSmallStrainMixedPressureConstitutiveMatrices_PlaneStrain(rMode, gp, tStep);
    }
}




void
SimpleMixedUPCrossSection::giveLargeStrainStiffnessMatrix(FloatMatrix &answer_uu, FloatArray &answer_up, double& K, double pressure, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
      std::tie(answer_uu, answer_up, K) = this->giveLargeStrainMixedPressureConstitutiveMatrices_3d(pressure, rMode, gp, tStep);
    }  else if ( mode == _PlaneStrain ) {
      std::tie(answer_uu, answer_up, K) = this->giveLargeStrainMixedPressureConstitutiveMatrices_PlaneStrain(pressure, rMode, gp, tStep);
    }
}
  



std::tuple<FloatMatrixF<9,9>, FloatArrayF<9>, double>
SimpleMixedUPCrossSection::giveLargeStrainMixedPressureConstitutiveMatrices_3d(double pressure, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
      LargeStrainMixedPressureMaterialExtensionInterface *mixedPressureMat = static_cast<  LargeStrainMixedPressureMaterialExtensionInterface * >( this->giveMaterialInterface(LargeStrainMixedPressureMaterialExtensionInterfaceType, gp) );
    if ( !mixedPressureMat ) {
        OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
    }


    return mixedPressureMat->giveLargeStrainMixedPressureConstitutiveMatrices_3d(pressure, rMode, gp, tStep);
}


std::tuple<FloatMatrixF<5,5>, FloatArrayF<5>, double>
SimpleMixedUPCrossSection::giveLargeStrainMixedPressureConstitutiveMatrices_PlaneStrain(double pressure, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  
      LargeStrainMixedPressureMaterialExtensionInterface *mixedPressureMat = static_cast<  LargeStrainMixedPressureMaterialExtensionInterface * >( this->giveMaterialInterface(LargeStrainMixedPressureMaterialExtensionInterfaceType, gp) );
    if ( !mixedPressureMat ) {
        OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
    }

    return mixedPressureMat->giveLargeStrainMixedPressureConstitutiveMatrices_PlaneStrain(pressure, rMode, gp, tStep);
}




} // end namespace oofem
