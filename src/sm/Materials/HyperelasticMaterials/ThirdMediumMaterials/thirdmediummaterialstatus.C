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

#include "thirdmediummaterialstatus.h"
#include "gausspoint.h"

namespace oofem {
ThirdMediumMaterialStatus ::ThirdMediumMaterialStatus( GaussPoint *g ) :
  StructuralMaterialStatus(g), gradFVector(), tempGradFVector(), FbarVector(), tempFbarVector()
{
  gradFVector.resize( 27 );
  tempGradFVector = gradFVector;
  FbarVector.resize( 9 );
  FbarVector.at(1) = FbarVector.at(2) = FbarVector.at(3) = 1.;
  tempFbarVector = FbarVector;    

}

void ThirdMediumMaterialStatus ::updateYourself( TimeStep *tStep )
// Performs end-of-step updates.
{
    StructuralMaterialStatus :: updateYourself(tStep);

    gradFVector = tempGradFVector;
    FbarVector = tempFbarVector;
}


void ThirdMediumMaterialStatus ::initTempStatus()
//
// initialize record at the begining of new load step
//
{
    StructuralMaterialStatus :: initTempStatus();

    tempGradFVector = gradFVector;
    tempFbarVector = FbarVector;

}

void ThirdMediumMaterialStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
  StructuralMaterialStatus::copyStateVariables( iStatus );
  
  const ThirdMediumMaterialStatus &structStatus = static_cast<const ThirdMediumMaterialStatus &>( iStatus );
  gradFVector = structStatus.giveGradFVector();
  tempGradFVector = structStatus.giveTempGradFVector();
  FbarVector = structStatus.giveFbarVector();
  tempFbarVector = structStatus.giveTempFbarVector();
}

} // end namespace oofem
