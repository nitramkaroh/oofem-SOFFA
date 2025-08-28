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

#pragma once

#include "HyperelasticMaterials/constantstiffnesshyperelasticmaterial.h"
#include "thirdmediummaterial.h"
#include "thirdmediummaterialstatus.h"
#include "classfactory.h"


///@name Input fields for MooneyRivlinMaterial
//@{
#define _IFT_ThirdMediumConstantStiffnessMaterial_Name "thirdmediumconstantstiffnessmat"
//@}

namespace oofem {
/**
 * This class implements a constant stiffness material in conjunction with 3M functionality
 *
 * @author Ondrej Faltus
 */
class ThirdMediumConstantStiffnessMaterial : public ConstantStiffnessHyperElasticMaterial, public ThirdMediumMaterial
{

  public:
  ThirdMediumConstantStiffnessMaterial( int n, Domain *d ) :
          ConstantStiffnessHyperElasticMaterial( n, d ) {};

  void initializeFrom( InputRecord &ir ) override
  {
    ThirdMediumMaterial::initializeFrom( ir );
    ConstantStiffnessHyperElasticMaterial::initializeFrom( ir );
  }

  MaterialStatus *CreateStatus( GaussPoint *gp ) const override {
    return new ThirdMediumMaterialStatus(gp);
  }


  const char *giveInputRecordName() const override { return _IFT_ThirdMediumConstantStiffnessMaterial_Name; }
  const char *giveClassName() const override { return "ThirdMediumConstantStiffnessMaterial"; }
};

REGISTER_Material( ThirdMediumConstantStiffnessMaterial );

} // end namespace oofem
