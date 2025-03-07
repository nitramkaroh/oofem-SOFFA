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



#include "domain.h"
#include "mixedpressurematerialextensioninterface.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "../sm/Materials/structuralmaterial.h"

#include <list>


namespace oofem {
// constructor
MixedPressureMaterialExtensionInterface :: MixedPressureMaterialExtensionInterface(Domain *d)  : Interface()
{
    dom = d;
}


SmallStrainMixedPressureMaterialExtensionInterface :: SmallStrainMixedPressureMaterialExtensionInterface(Domain *d)  : MixedPressureMaterialExtensionInterface(d)
{
} 

  
std::tuple<FloatMatrixF<4,4>, FloatArrayF<4>, double>
SmallStrainMixedPressureMaterialExtensionInterface :: giveSmallStrainMixedPressureConstitutiveMatrices_PlaneStrain(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  OOFEM_ERROR("Not implemented");
}

std::tuple<FloatArrayF<4>, double>
SmallStrainMixedPressureMaterialExtensionInterface :: giveRealStressVector_PlaneStrain(const FloatArrayF<4> &reducedStrain, double pressure, GaussPoint *gp, TimeStep *tStep)
{
  OOFEM_ERROR("Not implemented");

}

LargeStrainMixedPressureMaterialExtensionInterface :: LargeStrainMixedPressureMaterialExtensionInterface(Domain *d)  : MixedPressureMaterialExtensionInterface(d)
{
} 

  

std::tuple<FloatArrayF<5>, double>
LargeStrainMixedPressureMaterialExtensionInterface :: giveFirstPKStressVector_PlaneStrain(const FloatArrayF<5> &F, double pressure, GaussPoint *gp, TimeStep *tStep)
{
  OOFEM_ERROR("Not implemented");
}


std::tuple<FloatMatrixF<5,5>, FloatArrayF<5>, double>
LargeStrainMixedPressureMaterialExtensionInterface :: giveLargeStrainMixedPressureConstitutiveMatrices_PlaneStrain(double pressure, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  OOFEM_ERROR("Not implemented");
}



  
  

  

  
  
}
